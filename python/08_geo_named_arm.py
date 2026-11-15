# -*- coding: utf-8 -*-
"""Resistant-versus-sensitive log2 fold change of the NAMED ARM of each candidate microRNA.

GSE28267 is analyzed separately for its two array batches (GPL7722 and GPL7723).
Each microRNA is read from the probe whose platform annotation (miRBase 12.0)
carries the miRBase accession of the named mature arm. GSE66607 (one array per
state) is read from the mature or star probe that corresponds to the same arm.
Enzyme values (GSE26459) are copied from results/04_geo_foldchange.csv (limma in R).

Requires: python >= 3.10, numpy, pandas, inmoose (limma port). Run from the repository root:
    python3 python/08_geo_named_arm.py
Downloads the GEO series matrices and platform tables into data/geo/ (not stored in the repository)
and writes results/08_geo_named_arm_foldchange.csv.
"""
import gzip, math, os, urllib.request, warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
from inmoose.limma import lmFit, eBayes, topTable

GEO = "https://ftp.ncbi.nlm.nih.gov/geo/"
FILES = {
    "GSE28267-GPL7722.txt.gz": "series/GSE28nnn/GSE28267/matrix/GSE28267-GPL7722_series_matrix.txt.gz",
    "GSE28267-GPL7723.txt.gz": "series/GSE28nnn/GSE28267/matrix/GSE28267-GPL7723_series_matrix.txt.gz",
    "GSE66607.txt.gz": "series/GSE66nnn/GSE66607/matrix/GSE66607_series_matrix.txt.gz",
    "GPL7722.soft.gz": "platforms/GPL7nnn/GPL7722/soft/GPL7722_family.soft.gz",
    "GPL7723.soft.gz": "platforms/GPL7nnn/GPL7723/soft/GPL7723_family.soft.gz",
}
CACHE = os.path.join("data", "geo")
ARM = {"miR-200c-3p": "MIMAT0000617", "miR-200c-5p": "MIMAT0004657", "miR-148a-3p": "MIMAT0000243",
       "miR-152-3p": "MIMAT0000438", "miR-24-3p": "MIMAT0000080", "miR-222-3p": "MIMAT0000279",
       "miR-7-5p": "MIMAT0000252", "miR-181b-5p": "MIMAT0000257", "miR-181d-5p": "MIMAT0002821",
       "miR-195-5p": "MIMAT0000461"}
AFFY = {"miR-200c-3p": "hsa-miR-200c_st", "miR-200c-5p": "hsa-miR-200c-star_st", "miR-148a-3p": "hsa-miR-148a_st",
        "miR-152-3p": "hsa-miR-152_st", "miR-24-3p": "hsa-miR-24_st", "miR-222-3p": "hsa-miR-222_st",
        "miR-7-5p": "hsa-miR-7_st", "miR-181b-5p": "hsa-miR-181b_st", "miR-181d-5p": "hsa-miR-181d_st",
        "miR-195-5p": "hsa-miR-195_st"}

def fetch():
    os.makedirs(CACHE, exist_ok=True)
    for name, path in FILES.items():
        dst = os.path.join(CACHE, name)
        if not os.path.exists(dst):
            print("downloading", name); urllib.request.urlretrieve(GEO + path, dst)

def read_matrix(f):
    meta, rows, ids, cols, intab = {}, [], [], None, False
    for l in gzip.open(f, "rt", errors="replace"):
        l = l.rstrip("\n")
        if l.startswith("!Sample_"):
            k, *v = l.split("\t"); meta.setdefault(k, []).append([x.strip('"') for x in v])
        elif l.startswith("!series_matrix_table_begin"): intab = True
        elif l.startswith("!series_matrix_table_end"): intab = False
        elif intab:
            p = l.split("\t")
            if cols is None: cols = [x.strip('"') for x in p[1:]]
            else:
                ids.append(p[0].strip('"'))
                rows.append([float(x) if x not in ("", "null", "NA", "NaN") else np.nan for x in p[1:]])
    return meta, cols, ids, np.array(rows)

def read_platform(f):
    tab, hdr, intab = {}, None, False
    for l in gzip.open(f, "rt", errors="replace"):
        l = l.rstrip("\n")
        if l.startswith("!platform_table_begin"): intab = True; continue
        if l.startswith("!platform_table_end"): break
        if intab:
            p = l.split("\t")
            if hdr is None: hdr = p
            else: tab[p[0]] = dict(zip(hdr, p))
    return tab

def moderated_t(ex, ids, res):
    df = pd.DataFrame(ex, index=ids)
    ok = ((~np.isnan(ex[:, res])).sum(1) >= 1) & ((~np.isnan(ex[:, ~res])).sum(1) >= 1)
    df = df[ok]
    design = pd.DataFrame({"Intercept": 1.0, "resistant": res.astype(float)}, index=df.columns)
    fit = eBayes(lmFit(df, design))
    return topTable(fit, coef="column1", number=len(df), sort_by="p")

def main():
    fetch(); out = []
    for pf in ("GPL7722", "GPL7723"):
        meta, cols, ids, ex = read_matrix(os.path.join(CACHE, f"GSE28267-{pf}.txt.gz"))
        ch = meta["!Sample_characteristics_ch1"]
        res = np.array(["resistant" in s for s in ch[2]])
        agent = np.array([a.replace("agent: ", "") for a in ch[3]])
        keep = agent != "E2"                               # estradiol arrays excluded
        gpl = read_platform(os.path.join(CACHE, f"{pf}.soft.gz"))
        tt = moderated_t(ex[:, keep], ids, res[keep])
        for nm, acc in ARM.items():
            for pid in tt.index:
                if acc in gpl.get(pid, {}).get("Accession", "").split("/"):
                    v = ex[ids.index(pid), keep]; r = res[keep]
                    out.append(dict(series="GSE28267", platform=pf,
                        comparison="LY2 (resistant) minus MCF-7 (sensitive); vehicle and 4-OHT arrays; estradiol arrays excluded",
                        n_resistant=int(r.sum()), n_sensitive=int((~r).sum()), feature=nm, accession=acc, probe=pid,
                        platform_name=[t for t in gpl[pid]["miRNA_ID_LIST"].split("/") if t.startswith("hsa")][0],
                        log2FC=round(float(np.nanmean(v[r]) - np.nanmean(v[~r])), 4),
                        p=float(tt.loc[pid, "pvalue"]), adj_p=float(tt.loc[pid, "adj_pvalue"]),
                        test="moderated t (limma method; InMoose 0.9.1 port), BH across all probes of the batch"))
    meta, cols, ids, ex = read_matrix(os.path.join(CACHE, "GSE66607.txt.gz"))
    assert meta["!Sample_title"][0] == ["MCF7 parent cell", "MCF7 tamoxifen resistant cell"]
    for nm, pid in AFFY.items():
        par, rst = ex[ids.index(pid)]
        out.append(dict(series="GSE66607", platform="GPL16384",
            comparison="MCF7 tamoxifen-resistant minus MCF7 parent; log2(signal+1)", n_resistant=1, n_sensitive=1,
            feature=nm, accession="", probe=pid, platform_name=pid,
            log2FC=round(math.log2(rst + 1) - math.log2(par + 1), 4), p=None, adj_p=None,
            test="none (single array per state)"))
    old = pd.read_csv(os.path.join("results", "04_geo_foldchange.csv"))
    probe = {"OGT": "207564_x_at", "OGA": "235868_at", "EOGT": "221935_s_at"}
    for _, r in old[old.series == "GSE26459"].iterrows():
        out.append(dict(series="GSE26459", platform="GPL570",
            comparison="resistant subclone G11 minus sensitive subclone B7; 12 vs 12 pooled over four 4 h treatments; log2(signal+1)",
            n_resistant=12, n_sensitive=12, feature=r.target, accession="", probe=probe[r.target], platform_name=r.target,
            log2FC=r.log2FC, p=None, adj_p=r.adj_p, test="limma moderated t (authors' R run), BH across all probes"))
    pd.DataFrame(out).to_csv(os.path.join("results", "08_geo_named_arm_foldchange.csv"), index=False)
    print("written results/08_geo_named_arm_foldchange.csv:", len(out), "rows")

if __name__ == "__main__":
    main()
