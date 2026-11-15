# -*- coding: utf-8 -*-
"""Figure 2, rebuilt from the actual analysis outputs.
Panel A: Spearman correlation (miRNA vs target enzyme) in TCGA-BRCA.
Panel B: Cox hazard ratios for overall survival, with 95% CIs derived
exactly from HR and the Wald p-value (se = |ln HR| / z, z = Phi^-1(1-p/2))."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap
from matplotlib.patches import Rectangle
import numpy as np
from statistics import NormalDist

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "DejaVu Sans"]
plt.rcParams["mathtext.fontset"] = "dejavusans"

import csv, os, re, sys
RES = sys.argv[1] if len(sys.argv) > 1 else "../04_Analysis_results"
FIGOUT = sys.argv[2] if len(sys.argv) > 2 else "."

CORR = []
with open(os.path.join(RES, "02_tcga_correlation.csv"), newline="") as fh:
    for row in csv.DictReader(fh):
        m = row["mirna"].replace("hsa-", "")
        e = row["target"]; e = "OGA" if e in ("MGEA5", "OGA") else e
        CORR.append((m, e, float(row["rho"]), float(row["fdr"])))

SURV = []
with open(os.path.join(RES, "03_survival.csv"), newline="") as fh:
    for row in csv.DictReader(fh):
        feat = re.sub(r"\(.*\)$", "", row["feature"]).replace("hsa-", "")
        grp = "mir" if feat.startswith("miR") else "enzyme"
        SURV.append((feat, float(row["HR"]), float(row["HR_p"]), grp))
N_PAIRS = 1076
PRIORITY = ("miR-200c-3p", "EOGT")

def stars(p):
    return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 5e-2 else "ns"

def ci_from_p(hr, p):
    if p >= 1 or hr <= 0:
        return (np.nan, np.nan)
    b = np.log(hr)
    z = NormalDist().inv_cdf(1 - p / 2.0)
    if z == 0:
        return (np.nan, np.nan)
    se = abs(b) / z
    return (np.exp(b - 1.959964 * se), np.exp(b + 1.959964 * se))

ENZ = ["OGT", "OGA", "EOGT"]
mirs = {}
for m, e, r, f in CORR:
    mirs.setdefault(m, {})[e] = (r, f)
# order rows by most-negative rho achieved
row_order = sorted(mirs, key=lambda m: min(v[0] for v in mirs[m].values()))
_pri_y = [None]

# ---------------- figure ----------------
fig = plt.figure(figsize=(7.2, 11.4))
gs = fig.add_gridspec(2, 1, height_ratios=[len(row_order) + 2.2, len(SURV) + 2.0],
                      hspace=0.28, left=0.24, right=0.9, top=0.955, bottom=0.06)

# custom blue-white-red so negative rho is unambiguously blue
cmap = LinearSegmentedColormap.from_list("bwr2", ["#2166AC", "#4393C3", "#F7F7F7", "#D6604D", "#B2182B"])
norm = TwoSlopeNorm(vmin=-0.30, vcenter=0.0, vmax=0.30)

# ===== Panel A: correlation heatmap =====
axA = fig.add_subplot(gs[0])
nR, nC = len(row_order), len(ENZ)
for i, m in enumerate(row_order):
    for j, e in enumerate(ENZ):
        y = nR - 1 - i
        if e in mirs[m]:
            r, f = mirs[m][e]
            axA.add_patch(Rectangle((j, y), 1, 1, facecolor=cmap(norm(r)), edgecolor="white", lw=1.5))
            tcol = "white" if (m, e) == PRIORITY else "#111111"
            axA.text(j + 0.5, y + 0.60, f"{r:+.2f}".replace("-", "\u2212"), ha="center", va="center",
                     fontsize=10.5, fontweight="bold", color=tcol)
            axA.text(j + 0.5, y + 0.28, stars(f), ha="center", va="center",
                     fontsize=9, color=tcol)
            if (m, e) == PRIORITY:
                axA.add_patch(Rectangle((j, y), 1, 1, fill=False, edgecolor="#C9960A", lw=3.6, zorder=5))
        else:
            axA.add_patch(Rectangle((j, y), 1, 1, facecolor="#F0F0F0", edgecolor="white", lw=1.5))
            axA.text(j + 0.5, y + 0.5, "n/t", ha="center", va="center", fontsize=8.5, color="#B8B8B8")

axA.set_xlim(0, nC); axA.set_ylim(0, nR)
axA.set_xticks([j + 0.5 for j in range(nC)]); axA.set_xticklabels(ENZ, fontsize=12, fontweight="bold")
axA.xaxis.set_ticks_position("top"); axA.xaxis.set_label_position("top")
axA.set_yticks([nR - 1 - i + 0.5 for i in range(nR)]); axA.set_yticklabels(row_order, fontsize=11)
for s in axA.spines.values():
    s.set_visible(False)
axA.tick_params(length=0)
axA.set_title("A   miRNA\u2013enzyme correlation in TCGA-BRCA (n = 1076)",
              fontsize=12.5, fontweight="bold", loc="left", pad=22)

cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axA,
                    orientation="vertical", fraction=0.045, pad=0.05, ticks=[-0.3, -0.15, 0, 0.15, 0.3])
cbar.set_label("Spearman \u03c1", fontsize=10.5)
cbar.ax.tick_params(labelsize=9)
axA.text(0.0, -0.055, "* FDR<0.05    ** FDR<0.01    *** FDR<0.001     n/t, not tested",
         transform=axA.transAxes, fontsize=8.6, color="#444444")
axA.add_patch(Rectangle((0.0, -0.128), 0.032, 0.026, transform=axA.transAxes,
              fill=False, edgecolor="#C9960A", lw=2.4, clip_on=False, zorder=6))
axA.text(0.048, -0.115, "priority candidate (miR-200c-3p to EOGT)",
         transform=axA.transAxes, fontsize=8.6, color="#444444")

# ===== Panel B: HR forest =====
axB = fig.add_subplot(gs[1])
rows = [s for s in SURV if s[3] == "enzyme"] + [None] + sorted([s for s in SURV if s[3] == "mir"], key=lambda r: -r[1])
ypos, labels, ycur = [], [], len(rows)
for s in rows:
    ycur -= 1
    if s is None:
        continue
    feat, hr, p, grp = s
    lo, hi = ci_from_p(hr, p)
    sig = p < 0.05
    col = "#B2182B" if sig else "#5A5A5A"
    axB.plot([lo, hi], [ycur, ycur], color=col, lw=1.8, zorder=2, solid_capstyle="round")
    axB.plot([lo, lo], [ycur - 0.14, ycur + 0.14], color=col, lw=1.6)
    axB.plot([hi, hi], [ycur - 0.14, ycur + 0.14], color=col, lw=1.6)
    axB.scatter([hr], [ycur], s=58, color=col, zorder=3, edgecolors="white", linewidths=0.8)
    axB.text(1.86, ycur, f"{hr:.2f} ({lo:.2f}\u2013{hi:.2f}){'  '+stars(p) if sig else ''}",
             va="center", ha="left", fontsize=9, color=col)
    labels.append((ycur, feat, grp))

axB.axvline(1.0, ls="--", color="#999999", lw=1.2, zorder=1)
axB.set_xscale("log")
axB.set_xlim(0.6, 1.65)
axB.set_xticks([0.7, 0.85, 1.0, 1.2, 1.5])
axB.set_xticklabels(["0.7", "0.85", "1.0", "1.2", "1.5"], fontsize=9.5)
from matplotlib.ticker import NullFormatter
axB.xaxis.set_minor_formatter(NullFormatter())
axB.set_ylim(-0.7, len(rows) + 0.5)
axB.set_yticks([y for y, _, _ in labels])
axB.set_yticklabels([f for _, f, _ in labels], fontsize=11)
# color enzyme labels distinctly
for tick, (_, feat, grp) in zip(axB.get_yticklabels(), labels):
    tick.set_fontweight("bold" if grp == "enzyme" else "normal")
    tick.set_color("#1F4E79" if grp == "enzyme" else "#000000")
axB.set_xlabel("Hazard ratio for overall survival (95% CI, log scale)", fontsize=11)
for s in ["top", "right"]:
    axB.spines[s].set_visible(False)
axB.tick_params(length=3)
axB.set_title("B   Prognostic association (Cox, nominal)", fontsize=12.5, fontweight="bold", loc="left", pad=16)
# group separation without colliding with the feature labels
enz_ys = [y for y, _, g in labels if g == "enzyme"]
mir_ys = [y for y, _, g in labels if g == "mir"]
gap_y = (min(enz_ys) + max(mir_ys)) / 2.0
axB.axhline(gap_y, color="#DADADA", lw=1.0, zorder=0)
# vertical group labels placed in the empty strip just inside the left spine
axB.text(0.632, np.mean(enz_ys), "Enzymes", rotation=90, va="center", ha="center",
         fontsize=8.5, color="#1F4E79", fontweight="bold")
axB.text(0.632, np.mean(mir_ys), "microRNAs", rotation=90, va="center", ha="center",
         fontsize=8.5, color="#555555", fontweight="bold")

os.makedirs(FIGOUT, exist_ok=True)
fig.savefig(f"{FIGOUT}/Figure_2.pdf", bbox_inches="tight", facecolor="white")
fig.savefig(f"{FIGOUT}/Figure_2.png", dpi=400, bbox_inches="tight", facecolor="white")
fig.savefig(f"{FIGOUT}/_f2_tmp.png", dpi=1200, bbox_inches="tight", facecolor="white")
from PIL import Image
Image.MAX_IMAGE_PIXELS = None
Image.open(f"{FIGOUT}/_f2_tmp.png").convert("RGB").save(f"{FIGOUT}/Figure_2.tiff", dpi=(1200, 1200), compression="tiff_lzw")
os.remove(f"{FIGOUT}/_f2_tmp.png")
print("Figure 2 written to", FIGOUT)
