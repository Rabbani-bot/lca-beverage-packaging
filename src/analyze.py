"""Run the comparative LCA: tables and figures."""
from pathlib import Path
import json
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(ROOT.parent / "common"))
import lca  # noqa: E402
import style  # noqa: E402

FIG, OUT = ROOT / "outputs" / "figures", ROOT / "outputs"
SYS = list(lca.SYSTEMS)
COL = dict(zip(SYS, [style.PALETTE[1], style.PALETTE[3], style.PALETTE[5], style.PALETTE[0]]))
STAGE_COL = dict(zip(lca.STAGES, ["#1D4E89", "#4F7CAC", "#9DB4C0", "#F6AE2D", "#00A6A6", "#5C677D", "#8AC926"]))


def main():
    style.apply(); FIG.mkdir(parents=True, exist_ok=True)
    res = {}
    # ---------- deterministic ----------
    rows, stage_rows = [], []
    for s in SYS:
        for method in ["cutoff", "substitution"]:
            for g in ["AR6", "AR5"]:
                c = lca.characterise(lca.inventory(s, method=method), g)
                rows.append({"system": s, "method": method, "gwp_set": g, **{k: float(v) for k, v in c["Total"].items()}})
                if method == "cutoff" and g == "AR6":
                    for st in lca.STAGES:
                        stage_rows.append({"system": s, "stage": st, "GWP": float(c[st]["GWP"])})
    det = pd.DataFrame(rows); stages = pd.DataFrame(stage_rows)
    det.round(2).to_csv(OUT / "table2_deterministic_results.csv", index=False)
    piv = stages.pivot(index="system", columns="stage", values="GWP")[lca.STAGES].loc[SYS]
    (piv.div(piv.sum(axis=1), axis=0) * 100).round(1).to_csv(OUT / "table3_stage_contribution_pct.csv")
    pd.DataFrame([{"system": k, **v} for k, v in lca.SYSTEMS.items()]).to_csv(OUT / "table1_foreground_parameters.csv", index=False)
    # mass balance check
    mb = []
    for s in SYS:
        m = lca.inventory(s)["_mass"]
        mb.append({"system": s, "new_kg": m["new"], "recycled_kg": m["recycled"], "waste_kg": m["waste"],
                   "balance_error_kg": m["new"] - m["recycled"] - m["waste"]})
    pd.DataFrame(mb).round(4).to_csv(OUT / "table6_mass_balance.csv", index=False)

    # ---------- Monte Carlo ----------
    mc = lca.monte_carlo(10000)
    summ = []
    for s in SYS:
        for k in ["GWP", "CED", "Water"]:
            x = mc[s][k]
            summ.append({"system": s, "indicator": k, "mean": x.mean(), "sd": x.std(), "P5": np.percentile(x, 5),
                         "median": np.median(x), "P95": np.percentile(x, 95), "CV_%": 100 * x.std() / x.mean()})
    pd.DataFrame(summ).round(2).to_csv(OUT / "table4_monte_carlo_summary.csv", index=False)
    prob = pd.DataFrame([[np.mean(mc[a]["GWP"] < mc[b]["GWP"]) * 100 if a != b else np.nan for b in SYS] for a in SYS], index=SYS, columns=SYS)
    prob.round(1).to_csv(OUT / "table5_probability_lower_gwp.csv")
    res["deterministic_gwp_cutoff_AR6"] = {s: float(det.query("system==@s and method=='cutoff' and gwp_set=='AR6'").GWP.iloc[0]) for s in SYS}
    res["mc_gwp_P5_median_P95"] = {s: [float(np.percentile(mc[s]["GWP"], q)) for q in (5, 50, 95)] for s in SYS}
    res["prob_refillable_below_pet_pct"] = float(prob.loc["Glass refillable", "PET 0.5 L"])
    res["prob_refillable_below_can_pct"] = float(prob.loc["Glass refillable", "Aluminium can"])

    fig_stages(piv)
    fig_mc(mc)
    fig_prob(prob)
    res["tornado"] = fig_tornado()
    res["breakeven_trips_vs_pet"] = fig_breakeven()
    fig_recycling()
    fig_heatmap(det)
    fig_methods(det)
    (OUT / "results.json").write_text(json.dumps(res, indent=2))
    return res


def fig_stages(piv):
    fig, ax = plt.subplots(figsize=(10, 5))
    left = np.zeros(len(piv))
    for st in lca.STAGES:
        v = piv[st].to_numpy()
        if np.allclose(v, 0):
            continue
        ax.barh(piv.index, v, left=left, color=STAGE_COL[st], label=st, height=0.6, edgecolor="white")
        left += v
    for i, t in enumerate(piv.sum(axis=1)):
        ax.text(t + 8, i, f"{t:.0f}", va="center", fontweight="semibold")
    ax.invert_yaxis(); ax.set_xlabel("kg CO2e per 1,000 L delivered (IPCC AR6 GWP100, cut-off)")
    ax.set_title("Climate impact by life-cycle stage")
    ax.legend(ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.14), fontsize=9)
    style.footnote(fig, style.CREDIT); fig.savefig(FIG / "01_stage_contribution.png"); plt.close(fig)


def fig_mc(mc):
    fig, ax = plt.subplots(figsize=(10, 5))
    data = [mc[s]["GWP"] for s in SYS]
    parts = ax.violinplot(data, orientation="horizontal", showextrema=False, widths=0.8)
    for pc, s in zip(parts["bodies"], SYS):
        pc.set_facecolor(COL[s]); pc.set_alpha(0.55); pc.set_edgecolor("none")
    for i, x in enumerate(data, 1):
        p5, p25, p50, p75, p95 = np.percentile(x, [5, 25, 50, 75, 95])
        ax.plot([p5, p95], [i, i], color="#222", lw=1); ax.plot([p25, p75], [i, i], color="#222", lw=5)
        ax.scatter(p50, i, color="white", edgecolor="#222", zorder=4, s=40)
        ax.text(p95 + 15, i, f"median {p50:.0f}  [P5 {p5:.0f}, P95 {p95:.0f}]", va="center", fontsize=9)
    ax.set_yticks(range(1, len(SYS) + 1)); ax.set_yticklabels(SYS); ax.invert_yaxis()
    ax.set_xlim(0, 1150); ax.set_xlabel("kg CO2e per 1,000 L")
    ax.set_title("Monte Carlo uncertainty (10,000 runs, lognormal parameters, shared background draws)")
    style.footnote(fig, style.CREDIT); fig.savefig(FIG / "02_monte_carlo.png"); plt.close(fig)


def fig_prob(prob):
    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    im = ax.imshow(prob.to_numpy(), cmap="RdYlGn", vmin=0, vmax=100)
    for i in range(len(SYS)):
        for j in range(len(SYS)):
            if i != j:
                ax.text(j, i, f"{prob.iloc[i, j]:.0f}%", ha="center", va="center", fontweight="semibold")
    ax.set_xticks(range(len(SYS))); ax.set_xticklabels(SYS, rotation=25, ha="right"); ax.set_yticks(range(len(SYS))); ax.set_yticklabels(SYS)
    ax.grid(False); fig.colorbar(im, ax=ax, label="P(row has lower GWP than column)")
    ax.set_title("Probability of lower climate impact\n(paired Monte Carlo comparison)")
    style.footnote(fig, style.CREDIT); fig.savefig(FIG / "03_probability_matrix.png"); plt.close(fig)


def fig_tornado():
    s = "Glass refillable"
    base = float(lca.total(s))
    tests = {"Trips per bottle (15 / 40)": ("trips", 15, 40), "Bottle mass (300 / 400 g)": ("mass", 0.30, 0.40),
             "Distribution distance (150 / 600 km)": ("dist_km", 150, 600), "Washing heat (0.08 / 0.30 MJ)": ("wash_mj", 0.08, 0.30),
"Recycled content (30 / 70%)": ("rc", 0.30, 0.70),
             "Filling electricity (0.005 / 0.02 kWh)": ("fill_kwh", 0.005, 0.02)}
    rows = []
    for lab, (k, lo, hi) in tests.items():
        rows.append((lab, float(lca.total(s, {k: lo})) - base, float(lca.total(s, {k: hi})) - base))
    rows.sort(key=lambda r: abs(r[2] - r[1]))
    fig, ax = plt.subplots(figsize=(10, 4.8))
    for i, (lab, a, b) in enumerate(rows):
        ax.barh(i, a, color=style.PALETTE[0], height=0.6); ax.barh(i, b, color=style.PALETTE[3], height=0.6)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows]); ax.axvline(0, color="#222", lw=1)
    ax.set_xlabel(f"Change in kg CO2e per 1,000 L vs base case ({base:.0f})")
    ax.set_title("Sensitivity (tornado): refillable glass\nblue = first value, red = second value")
    style.footnote(fig, style.CREDIT); fig.savefig(FIG / "04_tornado_refillable.png"); plt.close(fig)
    return {r[0]: [round(r[1], 1), round(r[2], 1)] for r in rows}


def fig_breakeven():
    trips = np.arange(1, 41); dists = np.array([50, 100, 200, 300, 500, 800, 1200])
    pet = np.array([float(lca.total("PET 0.5 L", {"dist_km": d})) for d in dists])
    grid = np.array([[float(lca.total("Glass refillable", {"trips": t, "dist_km": d})) for t in trips] for d in dists])
    ratio = grid / pet[:, None]
    fig, ax = plt.subplots(figsize=(11, 4.8))
    im = ax.pcolormesh(trips, np.arange(len(dists)), ratio, cmap="RdBu_r", vmin=0.4, vmax=1.6, shading="nearest")
    cs = ax.contour(trips, np.arange(len(dists)), ratio, levels=[1.0], colors="#222", linewidths=2)
    ax.clabel(cs, fmt={1.0: "break-even"}, fontsize=9)
    ax.set_yticks(range(len(dists))); ax.set_yticklabels([f"{d} km" for d in dists]); ax.grid(False)
    ax.set_xlabel("Trips per refillable bottle"); ax.set_ylabel("Distribution distance")
    fig.colorbar(im, ax=ax, label="Refillable GWP / PET GWP")
    ax.set_title("When does refillable glass beat PET? Break-even depends on trips and distance")
    style.footnote(fig, style.CREDIT); fig.savefig(FIG / "05_breakeven_heatmap.png"); plt.close(fig)
    out = {}
    for d, row in zip(dists, ratio):
        below = np.where(row < 1)[0]
        out[f"{d} km"] = int(trips[below[0]]) if len(below) else "not reached within 40 trips"
    return out


def fig_recycling():
    rr = np.linspace(0.2, 0.95, 16)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), sharey=True)
    for ax, method in zip(axes, ["cutoff", "substitution"]):
        for s in SYS:
            y = [float(lca.total(s, {"rr": r}, method=method)) for r in rr]
            ax.plot(rr * 100, y, color=COL[s], lw=2.2, label=s)
            ax.scatter(lca.SYSTEMS[s]["rr"] * 100, float(lca.total(s, method=method)), color=COL[s], s=40, zorder=4, edgecolor="white")
        ax.set_title(f"{'Cut-off (recycled content)' if method == 'cutoff' else 'Substitution (avoided burden)'}")
        ax.set_xlabel("End-of-life recycling rate (%)")
    axes[0].set_ylabel("kg CO2e per 1,000 L"); axes[0].legend(fontsize=9)
    fig.suptitle("Recycling-rate scenarios: the end-of-life method changes how much recycling 'counts' (dots = base case)",
                 x=0.01, ha="left", fontweight="semibold")
    style.footnote(fig, style.CREDIT); fig.tight_layout(); fig.savefig(FIG / "06_recycling_scenarios.png"); plt.close(fig)


def fig_heatmap(det):
    d = det.query("method=='cutoff' and gwp_set=='AR6'").set_index("system")[["GWP", "CED", "Water"]].loc[SYS]
    norm = d / d.min()
    fig, ax = plt.subplots(figsize=(8, 4.6))
    im = ax.imshow(norm.to_numpy(), cmap="YlOrRd", vmin=1, vmax=4)
    units = {"GWP": "kg CO2e", "CED": "MJ", "Water": "L"}
    for i in range(len(SYS)):
        for j, k in enumerate(d.columns):
            ax.text(j, i, f"{d.iloc[i, j]:,.0f} {units[k]}\n({norm.iloc[i, j]:.1f}x best)", ha="center", va="center", fontsize=9)
    ax.set_xticks(range(3)); ax.set_xticklabels(["Climate change", "Cumulative energy demand", "Water use"])
    ax.set_yticks(range(len(SYS))); ax.set_yticklabels(SYS); ax.grid(False)
    ax.set_title("Multi-indicator results, normalised to the best system\nNo single winner: water ranks differently from GWP")
    style.footnote(fig, style.CREDIT); fig.savefig(FIG / "07_multi_indicator_heatmap.png"); plt.close(fig)


def fig_methods(det):
    fig, ax = plt.subplots(figsize=(10, 4.6))
    variants = [("cutoff", "AR6", "o", "Cut-off, AR6 (base)"), ("substitution", "AR6", "s", "Substitution, AR6"), ("cutoff", "AR5", "^", "Cut-off, AR5")]
    for i, s in enumerate(SYS):
        vals = [float(det.query("system==@s and method==@m and gwp_set==@g").GWP.iloc[0]) for m, g, _, _ in variants]
        ax.plot([min(vals), max(vals)], [i, i], color="#BBB", lw=3, zorder=1)
        for v, (_, _, mk, lab) in zip(vals, variants):
            ax.scatter(v, i, marker=mk, s=70, color=COL[s], edgecolor="#222", zorder=3, label=lab if i == 0 else None)
    ax.set_yticks(range(len(SYS))); ax.set_yticklabels(SYS); ax.invert_yaxis()
    ax.set_xlabel("kg CO2e per 1,000 L"); ax.legend(fontsize=9, loc="lower right")
    ax.set_title("Methodological choices: end-of-life allocation and GWP set (AR5 vs AR6)")
    style.footnote(fig, style.CREDIT); fig.savefig(FIG / "08_method_sensitivity.png"); plt.close(fig)


if __name__ == "__main__":
    print(json.dumps(main(), indent=1))
