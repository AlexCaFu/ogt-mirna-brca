# -*- coding: utf-8 -*-
"""Graphical Abstract. Numbers are read from the result tables.
Usage: python3 ga.py <results_dir> <out_dir>
"""
import csv, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

RES = sys.argv[1] if len(sys.argv) > 1 else "../04_Analysis_results"
OUT = sys.argv[2] if len(sys.argv) > 2 else "."
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "DejaVu Sans"]
plt.rcParams["mathtext.fontset"] = "dejavusans"

def rd(name): return list(csv.DictReader(open(os.path.join(RES, name))))
pri = [r for r in rd("02_tcga_correlation.csv") if r["mirna"].endswith("miR-200c-3p") and r["target"] == "EOGT"][0]
rho, fdr, n = float(pri["rho"]), float(pri["fdr"]), int(pri["n"])
part = [r for r in rd("06A_partial_correlation.csv") if r["model"].startswith("Partial")][0]
G8 = rd("08_geo_named_arm_foldchange.csv")
geo = {r["feature"]: r for r in G8 if r["series"] == "GSE26459"}
mA = [r for r in G8 if r["platform"] == "GPL7722" and r["feature"] == "miR-200c-3p"][0]
mB = [r for r in G8 if r["platform"] == "GPL7723" and r["feature"] == "miR-200c-3p"][0]
tn = rd("06D_EOGT_tumor_vs_normal.csv")[0]
M = lambda v, d=2: f"{float(v):+.{d}f}".replace("-", "−")
mant, expo = f"{fdr:.1e}".split("e")
FDR_TXT = f"{mant} × 10$^{{−{abs(int(expo))}}}$"

C_PRI, C_BLUE, C_EDGE = "#B01245", "#1F3B5C", "#2E6DA4"
fig = plt.figure(figsize=(13.2, 8.6))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 13.2); ax.set_ylim(0, 8.6); ax.axis("off")

def box(x, y, w, h, text, fc="#DCEBF7", ec=C_EDGE, fs=10.5, fw="bold", tc=C_BLUE, lw=1.6, ls="-", rs=0.12, z=3):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={rs}", fc=fc, ec=ec, lw=lw, ls=ls, zorder=z))
    if text: ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, fontweight=fw, color=tc, zorder=4, linespacing=1.25)
def arrow(p, q, col="#222222", lw=1.8, style="-|>", ls="-"):
    ax.annotate("", xy=q, xytext=p, arrowprops=dict(arrowstyle=style, color=col, lw=lw, mutation_scale=15, ls=ls, shrinkA=0, shrinkB=0), zorder=2)

# ---------- left: metabolic input ----------
box(2.55, 7.85, 3.5, 0.55, "Nutrient excess (glucose, glutamine)")
box(2.55, 6.95, 3.5, 0.55, "Hexosamine biosynthetic pathway flux")
box(3.3, 6.05, 2.0, 0.55, "UDP-GlcNAc")
arrow((4.3, 7.85), (4.3, 7.52)); arrow((4.3, 6.95), (4.3, 6.62))

# ---------- enzymes in two compartments ----------
box(0.35, 4.15, 5.2, 1.5, "", fc="#F4F7FB", ec="#7F8FA6", lw=1.1, ls=(0, (5, 3)), rs=0.2, z=1)
ax.text(2.95, 4.23, "Nucleus and cytoplasm: O-GlcNAc cycling", ha="center", va="bottom", fontsize=8.4, style="italic", color="#44546A", zorder=4)
box(0.75, 4.65, 1.7, 0.75, "OGT", fc="#D9C9EC", ec="#444444", fs=13, tc="#111111")
box(3.45, 4.65, 1.7, 0.75, "OGA", fc="#C5E3D6", ec="#444444", fs=13, tc="#111111")
arrow((2.5, 5.02), (3.4, 5.02), style="<|-|>", lw=1.8)
box(6.0, 4.15, 2.5, 1.5, "", fc="#FBF6EF", ec="#A6927F", lw=1.1, ls=(0, (5, 3)), rs=0.2, z=1)
ax.text(7.25, 4.23, "ER lumen: EGF-like domains", ha="center", va="bottom", fontsize=8.4, style="italic", color="#6B5A45", zorder=4)
box(6.4, 4.65, 1.7, 0.75, "EOGT", fc="#F8D2AE", ec="#444444", fs=13, tc="#111111")
arrow((3.9, 6.05), (2.95, 5.67)); arrow((4.7, 6.05), (7.0, 5.67))

# ---------- microRNA inputs ----------
box(0.35, 6.2, 2.0, 0.95, "Candidate microRNAs\n(12 pairs tested, Table 2)", fc="white", ec="#555555", fs=8.6, tc="#222222", fw="normal")
arrow((1.35, 6.2), (1.6, 5.67), col="#8A8A8A", lw=1.8, style="-", ls=(0, (1, 2)))
ax.text(0.3, 5.86, "node-specific,\nmostly unvalidated", fontsize=7.6, color="#555555", va="center", ha="left")
box(6.25, 6.55, 2.0, 0.6, "miR-200c-3p", fc="white", ec=C_PRI, fs=11.5, tc=C_PRI, lw=2.4)
ax.plot([7.25, 7.25], [6.55, 5.52], color=C_PRI, lw=3.2, zorder=5, solid_capstyle="butt")
ax.plot([6.98, 7.52], [5.5, 5.5], color=C_PRI, lw=3.6, zorder=5, solid_capstyle="butt")
ax.text(7.42, 6.08, "inverse association\n(candidate, not validated)", fontsize=8.0, color=C_PRI, va="center", ha="left")
ax.scatter([8.5], [6.85], marker="*", s=260, color=C_PRI, zorder=6)

# ---------- downstream ----------
box(0.35, 3.0, 5.2, 0.62, "Protein O-GlcNAcylation of resistance effectors")
arrow((2.95, 4.15), (2.95, 3.64))
box(6.0, 3.0, 2.5, 0.62, "Role in breast cancer\nresistance: hypothesis", fc="white", ec="#A6927F", fs=8.6, tc="#6B5A45", fw="normal", ls=(0, (4, 2)))
arrow((7.25, 4.15), (7.25, 3.64), col="#A6927F", ls=(0, (4, 2)))
eff = [("ERα", "↓ ERα expression\n↓ tamoxifen response"),
       ("TOP2A", "↑ TOP2A activity\n↓ doxorubicin response"),
       ("MITF", "↑ MITF activation\n↓ palbociclib response"),
       ("RUNX1", "↑ RUNX1 stability\n↓ palbociclib response")]
for i, (nm, tx) in enumerate(eff):
    x = 0.35 + i * 2.06
    box(x, 1.25, 1.92, 1.3, "", fc="#EAF2FB", ec=C_EDGE, lw=1.3)
    ax.text(x + 0.96, 2.27, nm, ha="center", va="center", fontsize=11, fontweight="bold", color=C_BLUE, zorder=4)
    ax.text(x + 0.96, 1.72, tx, ha="center", va="center", fontsize=8.1, color="#222222", zorder=4, linespacing=1.4)
    arrow((2.95, 3.0), (x + 0.96, 2.57), lw=1.3)

# ---------- banner ----------
box(0.35, 0.25, 8.15, 0.72, "", fc="#F6B21B", ec="#B07A00", lw=1.6)
ax.text(4.42, 0.74, "Node-specific candidate relationships, not a global rheostat", ha="center", va="center", fontsize=11, fontweight="bold", color="#3B2A00", zorder=4)
ax.text(4.42, 0.44, "miR-200c-3p/EOGT is prioritized for functional validation", ha="center", va="center", fontsize=9.4, color="#3B2A00", zorder=4)

# ---------- right: evidence panel ----------
box(8.95, 0.25, 4.0, 8.15, "", fc="white", ec=C_PRI, lw=2.0, rs=0.15)
ax.text(10.95, 8.02, "Evidence for the\nmiR-200c-3p/EOGT candidate node", ha="center", va="center", fontsize=11.2, fontweight="bold", color=C_PRI, zorder=4, linespacing=1.3)
ax.plot([9.15, 12.75], [7.6, 7.6], color=C_PRI, lw=1.0, zorder=4)
e = geo["EOGT"]
sig = lambda r: "FDR < 0.05" if float(r["adj_p"]) < 0.05 else "ns"
items = [
 ("Supports prioritization", None),
 ("Strongest inverse correlation of 12 tested\npairs in TCGA-BRCA",
  f"Spearman ρ = {M(rho)}, FDR = {FDR_TXT}, n = {n}"),
 ("Persists after adjustment for tumor purity\nand PAM50 subtype", f"partial ρ = {M(part['rho'])}, n = {part['n']}"),
 ("EOGT higher in tamoxifen-resistant cells", f"log2 fold change {M(e['log2FC'])}, FDR {float(e['adj_p']):.3f}"),
 ("miR-200c-3p lower in resistant LY2 cells;\nsize differs between two array batches", f"log2 fold change {M(mA['log2FC'])} ({sig(mA)}) and {M(mB['log2FC'])} ({sig(mB)})"),
 ("Limits interpretation", None),
 ("EOGT and miR-200c-3p measured in different\ncell models", None),
 ("Not a seed-predicted target of miR-200c-3p;\ndirect targeting not tested", None),
 ("EOGT lower in tumor than in normal tissue", f"log2 fold change {M(tn['log2FC'])}"),
 ("No association with overall survival", None),
]
y = 7.3
for head, sub in items:
    if sub is None and head in ("Supports prioritization", "Limits interpretation"):
        ax.text(9.2, y, head, fontsize=9.6, fontweight="bold", color=C_BLUE, va="center", zorder=4); y -= 0.42; continue
    nl = head.count("\n")
    ax.scatter([9.32], [y - 0.01], s=26, color=C_PRI, zorder=4)
    ax.text(9.52, y + 0.11, head, fontsize=8.9, fontweight="bold", color="#111111", va="top", zorder=4, linespacing=1.25)
    yy = y + 0.11 - 0.24 * (nl + 1) - 0.06
    if sub:
        ax.text(9.52, yy, sub, fontsize=8.5, color="#333333", va="top", zorder=4); yy -= 0.24
    y = yy - 0.3

os.makedirs(OUT, exist_ok=True)
fig.savefig(os.path.join(OUT, "Graphical_Abstract.png"), dpi=300, facecolor="white")
fig.savefig(os.path.join(OUT, "Graphical_Abstract.pdf"), facecolor="white")
fig.savefig(os.path.join(OUT, "Graphical_Abstract.tiff"), dpi=600, facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
print("Graphical Abstract written to", OUT, "| last y:", round(y, 2))
