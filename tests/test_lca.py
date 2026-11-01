import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import lca  # noqa: E402


class TestLCA(unittest.TestCase):
    def test_ar6_characterisation(self):
        inv = {st: np.zeros(5) for st in lca.STAGES}
        inv["Filling"] = np.array([1.0, 1.0, 1.0, 0, 0])
        self.assertAlmostEqual(float(lca.characterise(inv, "AR6")["Total"]["GWP"]), 1 + 29.8 + 273)
        self.assertAlmostEqual(float(lca.characterise(inv, "AR5")["Total"]["GWP"]), 1 + 30 + 265)

    def test_mass_balance_closes(self):
        for s in lca.SYSTEMS:
            m = lca.inventory(s)["_mass"]
            self.assertAlmostEqual(m["new"], m["recycled"] + m["waste"], places=9)

    def test_results_scale_linearly_with_functional_unit(self):
        for s in lca.SYSTEMS:
            a = lca.characterise(lca.inventory(s, fu_litres=1000))["Total"]["GWP"]
            b = lca.characterise(lca.inventory(s, fu_litres=2000))["Total"]["GWP"]
            self.assertAlmostEqual(float(b), 2 * float(a), places=6)

    def test_more_trips_lower_impact(self):
        self.assertLess(lca.total("Glass refillable", {"trips": 30}), lca.total("Glass refillable", {"trips": 10}))

    def test_substitution_gives_credit_when_recycling_exceeds_content(self):
        for s in lca.SYSTEMS:
            if lca.SYSTEMS[s]["rr"] > lca.SYSTEMS[s]["rc"]:
                self.assertLess(lca.total(s, method="substitution"), lca.total(s, method="cutoff"))

    def test_monte_carlo_centred_on_deterministic(self):
        mc = lca.monte_carlo(4000, seed=1)
        for s in lca.SYSTEMS:
            ratio = np.median(mc[s]["GWP"]) / lca.total(s)
            self.assertTrue(0.9 < ratio < 1.15, (s, ratio))

    def test_monte_carlo_restores_unit_processes(self):
        before = dict(lca.PROC); lca.monte_carlo(100)
        self.assertEqual(before, lca.PROC)


if __name__ == "__main__":
    unittest.main()
