# -*- coding: utf-8 -*-
"""Figure 1. Evidence map of candidate microRNA relationships with OGT, OGA and EOGT.
All numbers are read from the result tables:
  02_tcga_correlation.csv  (Spearman rho and BH-FDR, 12 tested pairs, n = 1076)
  08_geo_named_arm_foldchange.csv  (resistant vs sensitive log2FC of the named arm; GSE28267 batches GPL7722 and GPL7723)
Usage: python3 fig1.py <results_dir> <out_dir>
"""
import csv, os, sys, math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
from matplotlib.lines import Line2D

RES = sys.argv[1] if len(sys.argv) > 1 else "../04_Analysis_results"
OUT = sys.argv[2] if len(sys.argv) > 2 else "."
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "DejaVu Sans"]
plt.rcParams["mathtext.fontset"] = "dejavusans"

# ---------------- data ----------------
corr = {}
for r in csv.DictReader(open(os.path.join(RES, "02_tcga_correlation.csv"))):
    m = r["mirna"].replace("hsa-", ""); e = "OGA" if r["target"] in ("MGEA5", "OGA") else r["target"]
    corr[(m, e)] = (float(r["rho"]), float(r["fdr"]))
geo, geoA, geoB = {}, {}, {}
for r in csv.DictReader(open(os.path.join(RES, "08_geo_named_arm_foldchange.csv"))):
    v = (float(r["log2FC"]), float(r["adj_p"])) if r["adj_p"] else None
    if r["series"] == "GSE26459": geo[r["feature"]] = v
    elif r["platform"] == "GPL7722": geoA[r["feature"]] = v
    elif r["platform"] == "GPL7723": geoB[r["feature"]] = v

# literature facts carried from Table 2 of the manuscript (validation level, drug context)
VALID = {("miR-24-3p", "OGT"): "breast cancer", ("miR-7-5p", "OGT"): "colorectal cancer",
         ("miR-181d-5p", "OGT"): "ovarian cancer"}
DRUG = {"miR-24-3p": "T", "miR-148a-3p": "T", "miR-200c-3p": "T", "miR-152-3p": "T",
        "miR-7-5p": "D", "miR-222-3p": "D", "miR-181b-5p": "D", "miR-181d-5p": "D", "miR-195-5p": "D"}
PRIORITY = ("miR-200c-3p", "EOGT")

C_INV, C_POS, C_NS, C_PRI = "#1F5FA8", "#C0392B", "#8A8A8A", "#B01245"
C_T, C_D = "#11A7C7", "#6B3FA0"

def fmt_fc(name):
    if name in geo:
        fc, q = geo[name]
        s = f"{fc:+.2f}".replace("-", "−")
        return "resistant model: " + f"{s}, " + ("ns" if q >= 0.05 else f"FDR {q:.3f}")
    one = lambda v: f"{v[0]:+.2f}".replace("-", "−") + ("*" if v[1] < 0.05 else "")
    return "LY2 vs MCF-7: " + one(geoA[name]) + " | " + one(geoB[name])

# ---------------- layout ----------------
ENZ = {"OGT": (4.6, 5.6), "OGA": (9.6, 5.6), "EOGT": (12.6, 1.55)}
MIR = {"miR-24-3p": (1.55, 8.75), "miR-181d-5p": (4.6, 9.45), "miR-7-5p": (7.75, 8.75),
       "miR-181b-5p": (11.0, 9.45), "miR-222-3p": (14.6, 6.3), "miR-148a-3p": (1.2, 5.6),
       "miR-195-5p": (1.55, 2.5), "miR-200c-3p": (7.2, 1.55), "miR-152-3p": (15.2, 3.7)}
EW, EH = 1.42, 0.72      # enzyme half-size
MW, MH = 1.28, 0.47      # microRNA half-size

fig = plt.figure(figsize=(11.6, 9.4))
ax = fig.add_axes([0.01, 0.235, 0.98, 0.755]); ax.set_xlim(-0.35, 16.75); ax.set_ylim(0.05, 10.25); ax.axis("off")

def edge_point(c, hw, hh, toward, pad=0.0):
    """Point where the segment from box center c toward `toward` leaves the box (plus pad)."""
    dx, dy = toward[0] - c[0], toward[1] - c[1]
    t = min((hw + pad) / abs(dx) if dx else 1e9, (hh + pad) / abs(dy) if dy else 1e9)
    return (c[0] + dx * t, c[1] + dy * t)

# compartments
ax.add_patch(FancyBboxPatch((2.75, 4.35), 8.7, 2.55, boxstyle="round,pad=0.02,rounding_size=0.25",
             fc="#F4F7FB", ec="#7F8FA6", lw=1.1, ls=(0, (5, 3)), zorder=0))
ax.text(7.0, 4.5, "Nucleus and cytoplasm: O-GlcNAc cycling", ha="center", va="bottom",
        fontsize=9.3, color="#44546A", style="italic")
ax.add_patch(FancyBboxPatch((10.75, 0.42), 3.7, 2.38, boxstyle="round,pad=0.02,rounding_size=0.25",
             fc="#FBF6EF", ec="#A6927F", lw=1.1, ls=(0, (5, 3)), zorder=0))
ax.text(12.6, 0.52, "ER lumen: O-GlcNAc on EGF-like domains", ha="center", va="bottom",
        fontsize=9.0, color="#6B5A45", style="italic")

# OGT <-> OGA cycling arrow
ax.annotate("", xy=(ENZ["OGA"][0] - EW - 0.08, 5.6), xytext=(ENZ["OGT"][0] + EW + 0.08, 5.6),
            arrowprops=dict(arrowstyle="<|-|>", color="#222222", lw=2.0, mutation_scale=18), zorder=2)

# edges
def draw_edge(m, e):
    rho, fdr = corr[(m, e)]
    a = edge_point(MIR[m], MW, MH, ENZ[e], pad=0.05)
    b = edge_point(ENZ[e], EW, EH, MIR[m], pad=0.12)
    pri = (m, e) == PRIORITY
    if fdr < 0.05 and rho < 0:
        col, ls, lw = (C_PRI if pri else C_INV), "-", (3.6 if pri else 2.3)
    elif fdr < 0.05 and rho > 0:
        col, ls, lw = C_POS, (0, (5, 2.5)), 2.0
    else:
        col, ls, lw = C_NS, (0, (1, 2.2)), 2.0
    ax.plot([a[0], b[0]], [a[1], b[1]], color=col, lw=lw, ls=ls, zorder=1, solid_capstyle="butt")
    if fdr < 0.05 and rho < 0:      # T-bar at the enzyme end
        dx, dy = b[0] - a[0], b[1] - a[1]; L = math.hypot(dx, dy); px, py = -dy / L, dx / L; h = 0.26 if pri else 0.2
        ax.plot([b[0] - px * h, b[0] + px * h], [b[1] - py * h, b[1] + py * h], color=col, lw=lw + 0.4, zorder=3, solid_capstyle="butt")
    # rho label
    f = 0.5
    if (m, e) == ("miR-7-5p", "OGA"): f = 0.42
    if (m, e) == ("miR-152-3p", "EOGT"): f = 0.3
    mx, my = a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f
    if (m, e) == ("miR-148a-3p", "OGT"): mx, my = 2.05, 6.36
    if (m, e) == ("miR-222-3p", "EOGT"): mx, my = a[0] + (b[0] - a[0]) * 0.28, a[1] + (b[1] - a[1]) * 0.28
    lab = ("ρ = " + f"{rho:+.2f}".replace("-", "−")) + ("" if fdr < 0.05 else ", ns")
    if not pri:
        ax.text(mx, my, lab, ha="center", va="center", fontsize=8.6, color=col, fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.16", fc="white", ec="none"), zorder=4)
    if (m, e) in VALID:             # luciferase-validated badge near the microRNA end
        vx, vy = a[0] + (b[0] - a[0]) * 0.2, a[1] + (b[1] - a[1]) * 0.2
        ax.scatter([vx], [vy], s=230, marker="o", fc="#111111", ec="white", lw=1.2, zorder=5)
        ax.text(vx, vy - 0.01, "L", ha="center", va="center", fontsize=8.6, color="white", fontweight="bold", zorder=6)

for (m, e) in corr:
    draw_edge(m, e)

# enzyme boxes
ECOL = {"OGT": "#D9C9EC", "OGA": "#C5E3D6", "EOGT": "#F8D2AE"}
ENAME = {"OGT": "OGT", "OGA": "OGA (MGEA5)", "EOGT": "EOGT"}
for e, (x, y) in ENZ.items():
    ax.add_patch(FancyBboxPatch((x - EW, y - EH), 2 * EW, 2 * EH, boxstyle="round,pad=0.0,rounding_size=0.16",
                 fc=ECOL[e], ec="#444444", lw=1.6, zorder=7))
    ax.text(x, y + 0.2, ENAME[e], ha="center", va="center", fontsize=14.5 if e != "OGA" else 13, fontweight="bold", zorder=8)
    ax.text(x, y - 0.33, fmt_fc(e), ha="center", va="center", fontsize=7.9, color="#222222", zorder=8)

# microRNA boxes
for m, (x, y) in MIR.items():
    oc = C_T if DRUG[m] == "T" else C_D
    ax.add_patch(FancyBboxPatch((x - MW, y - MH), 2 * MW, 2 * MH, boxstyle="round,pad=0.0,rounding_size=0.12",
                 fc="white", ec=oc, lw=2.6, zorder=7))
    ax.text(x, y + 0.14, m, ha="center", va="center", fontsize=11, fontweight="bold", zorder=8)
    ax.text(x, y - 0.2, fmt_fc(m), ha="center", va="center", fontsize=7.4, color="#333333", zorder=8)
for (m, e), where in VALID.items():
    x, y = MIR[m]
    ax.text(x, y + MH + 0.1, f"luciferase-validated in {where}", ha="center", va="bottom", fontsize=7.8, color="#111111", style="italic")

# priority annotation
px = (MIR["miR-200c-3p"][0] + MW + ENZ["EOGT"][0] - EW) / 2
ax.scatter([MIR["miR-200c-3p"][0] - MW - 0.33], [1.55], marker="*", s=330, color=C_PRI, zorder=9)
rho, fdr = corr[PRIORITY]
mant, expo = f"{fdr:.1e}".split("e")
ax.text(px, 1.88, "ρ = " + f"{rho:.2f}".replace("-", "−"), ha="center", va="bottom", fontsize=9.6, color=C_PRI, fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none"), zorder=4)
ax.text(px, 1.2, f"FDR = {mant} × 10$^{{−{abs(int(expo))}}}$", ha="center", va="top", fontsize=8.8, color=C_PRI, fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none"), zorder=4)
ax.text(MIR["miR-200c-3p"][0], 0.98, "Priority candidate. No predicted site for the 3p arm;\ndirect targeting not tested.", ha="center", va="top",
        fontsize=8.4, color=C_PRI)

# ---------------- legend ----------------
lg = fig.add_axes([0.01, 0.005, 0.98, 0.225]); lg.set_xlim(0, 16.8); lg.set_ylim(0, 3.1); lg.axis("off")
lg.add_patch(Rectangle((0.12, 0.08), 16.56, 2.94, fc="#FAFAFA", ec="#555555", lw=1.0))
def head(x, t): lg.text(x, 2.72, t, fontsize=9.6, fontweight="bold", color="#1F3B5C", va="center")
head(0.4, "Connector: TCGA-BRCA co-expression (n = 1076, BH across 12 pairs)")
rows = [(C_INV, "-", 2.3, True, "Inverse correlation, FDR < 0.05 (compatible with repression; not proof of targeting)"),
        (C_PRI, "-", 3.6, True, "Priority candidate: strongest inverse correlation of the 12 tested pairs"),
        (C_POS, (0, (5, 2.5)), 2.0, False, "Positive correlation, FDR < 0.05 (opposite to a repression model)"),
        (C_NS, (0, (1, 2.2)), 2.0, False, "Not significant (not corroborated)")]
for i, (c, ls, lw, bar, txt) in enumerate(rows):
    y = 2.24 - i * 0.47
    lg.plot([0.45, 1.5], [y, y], color=c, lw=lw, ls=ls, solid_capstyle="butt")
    if bar: lg.plot([1.5, 1.5], [y - 0.13, y + 0.13], color=c, lw=lw + 0.4, solid_capstyle="butt")
    lg.text(1.75, y, txt, fontsize=8.7, va="center")
head(9.0, "Symbols and boxes")
lg.scatter([9.3], [2.24], s=210, marker="o", fc="#111111", ec="white", lw=1.2); lg.text(9.3, 2.23, "L", ha="center", va="center", fontsize=8.2, color="white", fontweight="bold")
lg.text(9.65, 2.24, "Direct interaction validated by luciferase reporter (cancer type stated)", fontsize=8.7, va="center")
lg.scatter([9.3], [1.77], marker="*", s=250, color=C_PRI); lg.text(9.65, 1.77, "Priority candidate for functional validation", fontsize=8.7, va="center")
for i, (c, t) in enumerate([(C_T, "tamoxifen context"), (C_D, "doxorubicin context")]):
    x0 = 9.05 + i * 3.7
    lg.add_patch(FancyBboxPatch((x0, 1.17), 0.55, 0.26, boxstyle="round,pad=0.0,rounding_size=0.06", fc="white", ec=c, lw=2.4))
    lg.text(x0 + 0.7, 1.3, "Box outline: " + t, fontsize=8.7, va="center")
lg.text(9.05, 0.83, "log2 fold change, resistant vs sensitive cells. Enzymes: GSE26459. microRNAs:", fontsize=8.7, va="center")
lg.text(9.05, 0.5, "GSE28267, array batch GPL7722 | batch GPL7723 (*, FDR < 0.05; ns, FDR ≥ 0.05)", fontsize=8.7, va="center")

os.makedirs(OUT, exist_ok=True)
fig.savefig(os.path.join(OUT, "Figure_1.png"), dpi=300, facecolor="white")
fig.savefig(os.path.join(OUT, "Figure_1.pdf"), facecolor="white")
fig.savefig(os.path.join(OUT, "Figure_1.tiff"), dpi=600, facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
print("Figure 1 written to", OUT)
