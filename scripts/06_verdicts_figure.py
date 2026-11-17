"""Verdict table and the summary figure."""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from tssc import RES  # noqa: E402

FIG = ROOT / "figures"
INK, MUTED, GRID = "#1f2328", "#5b616b", "#e3e5e8"
EVID, CLIN, BASE = "#2f6db3", "#c4622d", "#8a5fb0"


def verdicts():
    p1 = pd.read_csv(RES / "p1_concordance.csv").iloc[0]
    auc = pd.read_csv(RES / "auroc.csv").set_index("score")
    dd = pd.read_csv(RES / "auroc_differences.csv").set_index("comparison")
    p7 = pd.read_csv(RES / "p7_organ_match.csv").iloc[0]
    p8 = pd.read_csv(RES / "p8_audit.csv").iloc[0]
    d4 = dd.iloc[0]
    d5 = dd.iloc[1]
    d6 = dd.iloc[3]
    rows = [
        ("P1", "warnings shared within targets",
         f"{p1.concordant_share:.3f} vs null {p1.null_share:.3f}, p = {p1.p_one_sided:.4f}",
         p1.p_one_sided < 0.01),
        ("P2", "FAERS AUROC >= 0.65", f"{auc.loc['FAERS terms (raw)', 'auroc']:.3f}",
         auc.loc["FAERS terms (raw)", "auroc"] >= 0.65),
        ("P3", "drug count AUROC >= 0.60", f"{auc.loc['drug count (raw)', 'auroc']:.3f}",
         auc.loc["drug count (raw)", "auroc"] >= 0.60),
        ("P4", "M1 - M0 > 0, interval excludes 0",
         f"{d4['diff']:.3f} ({d4['low']:.3f} to {d4['high']:.3f})", d4["low"] > 0),
        ("P5", "M2 - M1 <= 0.03", f"{d5['diff']:.3f} ({d5['low']:.3f} to {d5['high']:.3f})",
         d5["diff"] <= 0.03),
        ("P6", "clinical minus experimental > 0.03",
         f"{d6['diff']:.3f} ({d6['low']:.3f} to {d6['high']:.3f})", d6["diff"] > 0.03),
        ("P7", "organ score higher, p < 0.05",
         f"{p7.weighted_mean_difference:.3f} SD, p = {p7.p_one_sided:.4f}", p7.p_one_sided < 0.05),
        ("P8", "sentences and numbers >= 95% traceable",
         f"{p8.share_sentences_pass:.1%} of {int(p8.sentences)}; {p8.share_numbers_supported:.1%} of {int(p8.numbers_in_facts)}",
         p8.share_sentences_pass >= 0.95 and p8.share_numbers_supported >= 0.95),
    ]
    v = pd.DataFrame(rows, columns=["prediction", "criterion", "result", "met"])
    v.to_csv(RES / "verdicts.csv", index=False)
    print(v.to_string(index=False))
    return auc


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelcolor=INK, labelsize=8.5)
    ax.xaxis.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)


def main():
    auc = verdicts()
    order = [
        ("drug count (raw)", "Drug count", BASE),
        ("FAERS terms (raw)", "FAERS signal terms", CLIN),
        ("curated event, clinical sources", "Curated events, clinical sources", CLIN),
        ("curated event, experimental sources", "Curated events, experimental", EVID),
        ("equal-weight evidence (raw)", "Equal-weight pre-drug score", EVID),
        ("risk_mouseKOScore", "  mouse knockout score", EVID),
        ("risk_tissueDistribution", "  tissue distribution", EVID),
        ("risk_tissueSpecificity", "  tissue specificity", EVID),
        ("risk_geneticConstraint", "  genetic constraint", EVID),
        ("risk_isCancerDriverGene", "  cancer driver gene", EVID),
        ("risk_paralogMaxIdentityPercentage", "  paralog identity", EVID),
    ]
    p7 = pd.read_csv(RES / "p7_organ_by_class.csv").sort_values("difference")

    fig, (a, b) = plt.subplots(1, 2, figsize=(11.5, 4.9), gridspec_kw={"width_ratios": [1.15, 1]})
    for k, (key, lab, col) in enumerate(order):
        r = auc.loc[key]
        yk = len(order) - 1 - k
        a.plot([r["low"], r["high"]], [yk, yk], color=col, lw=2, solid_capstyle="round")
        a.plot(r["auroc"], yk, "o", ms=6.5, color=col, mec="white", mew=1.2)
    a.axvline(0.5, color=MUTED, lw=1, ls=(0, (3, 3)))
    a.set_yticks(range(len(order)))
    a.set_yticklabels([o[1] for o in order][::-1])
    a.set_xlim(0.44, 0.78)
    a.set_xlabel("AUROC for a boxed warning or withdrawal (95% CI)", color=INK, fontsize=9)
    a.set_title("A  Which evidence ranks 438 targets", loc="left", fontsize=10.5, color=INK)
    style(a)
    a.yaxis.grid(False)
    from matplotlib.lines import Line2D
    a.legend(handles=[Line2D([], [], color=c, marker="o", lw=2, label=l) for c, l in
                      ((BASE, "baseline"), (CLIN, "built from drug experience"), (EVID, "available before a drug"))],
             loc="lower right", fontsize=8, frameon=False, labelcolor=INK)

    ys = range(len(p7))
    cols = [EVID if d > 0 else MUTED for d in p7["difference"]]
    b.hlines(list(ys), 0, p7["difference"], color=cols, lw=2)
    b.scatter(p7["difference"], list(ys), s=46, c=cols, edgecolors="white", linewidths=1.2, zorder=3)
    for yk, (d, n) in enumerate(zip(p7["difference"], p7["n_class"])):
        b.text(d + (0.05 if d >= 0 else -0.05), yk, f"n = {int(n)}", va="center",
               ha="left" if d >= 0 else "right", fontsize=7.5, color=MUTED)
    b.axvline(0, color=MUTED, lw=1)
    b.set_yticks(list(ys))
    b.set_yticklabels([c.replace(" toxicity", "") for c in p7["toxicityClass"]])
    b.set_xlim(-0.9, 1.5)
    b.set_xlabel("Expression z-score in the organ, difference", color=INK, fontsize=9)
    b.set_title("B  Do warnings fall where the target is expressed?", loc="left", fontsize=10.5, color=INK)
    style(b)
    b.yaxis.grid(False)
    fig.tight_layout()
    FIG.mkdir(exist_ok=True)
    fig.savefig(FIG / "evidence_and_organs.png", dpi=160, facecolor="white", metadata={"Software": None})
    print("figure written")


if __name__ == "__main__":
    main()
