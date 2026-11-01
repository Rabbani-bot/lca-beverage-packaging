"""Comparative, ISO 14040/44-structured LCA of beverage packaging.

Functional unit: delivery of 1,000 L of carbonated soft drink to retail.
Systems: PET bottle 0.5 L, single-use glass 0.33 L, refillable glass 0.33 L, aluminium can 0.33 L.
Impact categories: climate change (GWP100, IPCC AR6; AR5 as sensitivity), cumulative energy demand (MJ),
water use (L).

INVENTORY NOTE. Process data are rounded values within ranges typical of public LCA literature.
They are NOT drawn from ecoinvent or any licensed database and must not be used for real claims.
"""
from __future__ import annotations
import numpy as np

GWP = {"AR6": {"co2": 1.0, "ch4": 29.8, "n2o": 273.0},   # fossil CH4, AR6 WG1 Table 7.15
       "AR5": {"co2": 1.0, "ch4": 30.0, "n2o": 265.0}}   # fossil CH4 incl. carbon feedback, AR5
FLOWS = ["co2", "ch4", "n2o", "mj", "water"]

# unit process inventory: kg CO2, kg CH4, kg N2O, MJ primary energy, L water per reference unit
PROC = {
    "pet_virgin":    (2.05, 0.0055, 5e-5, 75.0, 30.0),    # per kg
    "pet_recycled":  (0.55, 8e-4, 2e-5, 15.0, 5.0),
    "glass_virgin":  (0.85, 1.2e-3, 2e-5, 12.0, 3.0),
    "glass_cullet":  (0.55, 8e-4, 1e-5, 8.0, 2.0),
    "al_primary":    (11.5, 0.012, 2e-4, 160.0, 60.0),
    "al_recycled":   (0.55, 1e-3, 2e-5, 10.0, 3.0),
    "electricity":   (0.48, 1.2e-3, 1e-5, 9.5, 2.0),      # per kWh, generic grid mix
    "natgas_heat":   (0.056, 9e-5, 1e-6, 1.1, 0.02),      # per MJ, incl. upstream methane
    "truck":         (0.095, 1e-5, 3e-6, 1.4, 0.0),       # per t*km
    "incin_pet":     (2.29, 0.0, 4e-5, 0.0, 0.0),         # per kg, fossil carbon combusted
    "landfill_inert": (0.01, 0.0, 0.0, 0.1, 0.0),         # per kg
    "tap_water":     (0.0003, 0.0, 0.0, 0.004, 1.0),      # per L
}

SYSTEMS = {
    #                vol L, mass kg, material, recycled content, forming kWh/kg, recycling rate, trips
    "PET 0.5 L":            dict(vol=0.5, mass=0.022, mat="pet", rc=0.25, form_kwh=1.2, rr=0.45, trips=1),
    "Glass single-use":     dict(vol=0.33, mass=0.200, mat="glass", rc=0.50, form_kwh=0.5, rr=0.70, trips=1),
    "Glass refillable":     dict(vol=0.33, mass=0.350, mat="glass", rc=0.50, form_kwh=0.5, rr=0.90, trips=25),
    "Aluminium can":        dict(vol=0.33, mass=0.013, mat="al", rc=0.70, form_kwh=2.0, rr=0.75, trips=1),
}
VIRGIN = {"pet": "pet_virgin", "glass": "glass_virgin", "al": "al_primary"}
RECYC = {"pet": "pet_recycled", "glass": "glass_cullet", "al": "al_recycled"}
COMMON = dict(dist_km=300.0, fill_kwh=0.01, wash_mj=0.15, wash_kwh=0.01, wash_water=0.5,
              incin_share=0.4, quality=0.9)
STAGES = ["Material production", "Forming", "Filling", "Distribution", "Washing & return", "End of life", "Recycling credit"]


def _vec(name, scale):
    """Scale a unit process; works for scalar or array parameters and scalar or array (Monte Carlo) flows."""
    return np.stack([np.asarray(scale, float) * np.asarray(f, float) for f in PROC[name]], axis=-1)


def inventory(sys_name: str, p: dict | None = None, method: str = "cutoff", fu_litres: float = 1000.0) -> dict:
    """Return {stage: array(..., 5 flows)} for one system. Parameters may be numpy arrays (Monte Carlo)."""
    s = {**SYSTEMS[sys_name], **COMMON, **(p or {})}
    m = s["mat"]
    n_fill = fu_litres / s["vol"]                     # fillings needed
    n_new = n_fill / s["trips"]                       # new containers produced
    mass_new = n_new * s["mass"]                      # kg of new packaging
    rc = s["rc"]
    inv = {st: 0.0 for st in STAGES}
    inv["Material production"] = _vec(VIRGIN[m], mass_new * (1 - rc)) + _vec(RECYC[m], mass_new * rc)
    inv["Forming"] = _vec("electricity", mass_new * s["form_kwh"])
    inv["Filling"] = _vec("electricity", n_fill * s["fill_kwh"])
    product_t = fu_litres * 1.04 / 1000                # beverage density 1.04 kg/L
    pkg_t_full = n_fill * s["mass"] / 1000
    inv["Distribution"] = _vec("truck", (product_t + pkg_t_full) * s["dist_km"])
    if SYSTEMS[sys_name]["trips"] > 1:
        returns = n_fill - n_new
        inv["Washing & return"] = (_vec("natgas_heat", returns * s["wash_mj"]) + _vec("electricity", returns * s["wash_kwh"])
                                   + _vec("tap_water", returns * s["wash_water"]) + _vec("truck", returns * s["mass"] / 1000 * s["dist_km"]))
    else:
        inv["Washing & return"] = _vec("tap_water", 0.0 * n_fill)
    waste = mass_new * (1 - s["rr"])
    if m == "pet":
        inv["End of life"] = _vec("incin_pet", waste * s["incin_share"]) + _vec("landfill_inert", waste * (1 - s["incin_share"]))
    else:
        inv["End of life"] = _vec("landfill_inert", waste)
    inv["Recycling credit"] = _vec("landfill_inert", 0.0 * mass_new)
    if method == "substitution":
        # avoided burden: recycled output displaces virgin material, net of recycled content already credited
        net_out = mass_new * (s["rr"] - rc) * s["quality"]
        inv["Recycling credit"] = -np.stack([np.asarray(net_out, float) * (np.asarray(a) - np.asarray(b))
                                             for a, b in zip(PROC[VIRGIN[m]], PROC[RECYC[m]])], axis=-1)
    inv["_mass"] = dict(new=mass_new, recycled=mass_new * s["rr"], waste=waste)
    return inv


def characterise(inv: dict, gwp: str = "AR6") -> dict:
    cf = GWP[gwp]
    out = {}
    for st in STAGES:
        a = np.asarray(inv[st])
        out[st] = {"GWP": a[..., 0] * cf["co2"] + a[..., 1] * cf["ch4"] + a[..., 2] * cf["n2o"],
                   "CED": a[..., 3], "Water": a[..., 4]}
    out["Total"] = {k: sum(out[st][k] for st in STAGES) for k in ["GWP", "CED", "Water"]}
    return out


def total(sys_name, p=None, method="cutoff", gwp="AR6", indicator="GWP"):
    return characterise(inventory(sys_name, p, method), gwp)["Total"][indicator]


# ---------- uncertainty ----------
GSD = {"mass": 1.05, "form_kwh": 1.2, "dist_km": 1.4, "trips": 1.35, "wash_mj": 1.3}


def sample_params(sys_name, n, rng):
    s = SYSTEMS[sys_name]
    p = {k: s[k] * rng.lognormal(0, np.log(g), n) for k, g in GSD.items() if k in s}
    p["dist_km"] = COMMON["dist_km"] * rng.lognormal(0, np.log(GSD["dist_km"]), n)
    p["wash_mj"] = COMMON["wash_mj"] * rng.lognormal(0, np.log(GSD["wash_mj"]), n)
    if s["trips"] > 1:
        p["trips"] = np.clip(p["trips"], 2, None)
    else:
        p["trips"] = np.ones(n)
    p["rr"] = np.clip(s["rr"] + rng.normal(0, 0.07, n), 0, 0.98)
    p["rc"] = np.clip(s["rc"] + rng.normal(0, 0.05, n), 0, 0.95)
    return p


def monte_carlo(n=10000, seed=7, method="cutoff"):
    """Background unit processes share one random draw across systems (common random numbers),
    so pairwise comparisons are correlated correctly; foreground parameters are sampled per system."""
    rng = np.random.default_rng(seed)
    base = dict(PROC)
    res = {}
    try:
        for k in base:
            f = rng.lognormal(0, np.log(1.15), n)               # GSD 1.15 on each unit process
            PROC[k] = tuple(f * v for v in base[k])
        for sname in SYSTEMS:
            c = characterise(inventory(sname, sample_params(sname, n, rng), method))
            res[sname] = {k: c["Total"][k] for k in ["GWP", "CED", "Water"]}
    finally:
        PROC.clear(); PROC.update(base)
    return res
