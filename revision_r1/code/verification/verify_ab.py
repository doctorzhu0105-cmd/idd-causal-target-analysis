"""Independent numerical verifier, adapted for isolated inputs and outputs."""
import numpy as np, pandas as pd, json, os

def run(paired, reference, output_dir):
    PAIRED = str(paired)
    OUT = str(output_dir)
    df = pd.read_csv(PAIRED, sep="\t")
    res = {"paired_shape": list(df.shape), "columns": list(df.columns)}
    res["scenarios"] = sorted(df["scenario"].dropna().unique().tolist())
    res["endpoints"] = sorted(df["endpoint"].dropna().unique().tolist())

    # --- Exposure-side MATCHED distributions per scenario ---
    exp = df[(df["endpoint"]=="EXPOSURE") & (df["status"]=="MATCHED")].copy()
    res["exposure_matched_N_per_scenario"] = exp.groupby("scenario").size().to_dict()
    tab=[]
    for sc, g in exp.groupby("scenario"):
        v = g["log_ratio_information"].to_numpy(float)
        d = g["diff_information"].to_numpy(float)
        tab.append({"scenario":sc,"N":len(v),"median":float(np.median(v)),
            "Q1":float(np.percentile(v,25)),"Q3":float(np.percentile(v,75)),
            "P5":float(np.percentile(v,5)),"P90":float(np.percentile(v,90)),
            "P95":float(np.percentile(v,95)),"P99":float(np.percentile(v,99)),
            "min":float(v.min()),"max":float(v.max()),
            "pos":int((d>1e-12).sum()),"zero":int((np.abs(d)<=1e-12).sum()),"neg":int((d<-1e-12).sum())})
    res["exposure_distributions"]=tab

    # --- Compare to sealed AB_FULL_DISTRIBUTIONS.tsv ---
    ab = pd.read_csv(str(reference),sep="\t",on_bad_lines="skip")
    ab.columns = ab.columns.str.strip()
    ab_exp = ab[ab["endpoint"]=="EXPOSURE"]
    cmp=[]
    for _,r in ab_exp.iterrows():
        sc=r["scenario"]; mine=[t for t in tab if t["scenario"]==sc][0]
        cmp.append({"scenario":sc,"sealed_median":float(r["median_log_I"]),"recomputed_median":mine["median"],
                    "abs_diff":abs(float(r["median_log_I"])-mine["median"]),
                    "sealed_pos":int(r["positive_diff"]),"recomputed_pos":mine["pos"],
                    "sealed_N":int(r["N"]),"recomputed_N":mine["N"]})
    res["median_reconciliation"]=cmp
    res["max_median_abs_diff"]=float(max(c["abs_diff"] for c in cmp))

    # --- Outcome endpoints identical I_X claim ---
    out_med={}
    for ep in ["M13_INTERVERTEB","M13_LOWBACKPAIN","M13_SCIATICA"]:
        g=df[(df["endpoint"]==ep)&(df["status"]=="MATCHED")&(df["scenario"]=="A0")]
        out_med[ep]=float(np.median(g["log_ratio_information"].to_numpy(float)))
    res["outcome_A0_medians"]=out_med
    res["outcome_identical"]=bool(max(out_med.values())-min(out_med.values())<1e-9)

    # --- Direction transitions vs A0 baseline (exposure + 3 endpoints) ---
    sub=df[df["status"]=="MATCHED"].copy()
    def sign(x): return "POS" if x>1e-12 else ("NEG" if x<-1e-12 else "ZERO")
    sub["sgn"]=sub["diff_information"].map(sign)
    base=sub[sub["scenario"]=="A0"][["ensembl_gene_id","endpoint","sgn","log_ratio_information"]].rename(columns={"sgn":"base","log_ratio_information":"bval"})
    trans={}
    for sc in ["A5","A10","B1","B4","B9"]:
        cur=sub[sub["scenario"]==sc][["ensembl_gene_id","endpoint","sgn","log_ratio_information"]].rename(columns={"sgn":"cur","log_ratio_information":"cval"})
        m=cur.merge(base,on=["ensembl_gene_id","endpoint"],how="inner")
        ct=pd.crosstab(m["base"],m["cur"])
        def cell(a,b): return int(ct.loc[a,b]) if (a in ct.index and b in ct.columns) else 0
        gt=int((m["cval"]>m["bval"]+1e-12).sum()); lt=int((m["cval"]<m["bval"]-1e-12).sum()); eq=int((np.abs(m["cval"]-m["bval"])<=1e-12).sum())
        trans[sc]={"N":len(m),"POS_POS":cell("POS","POS"),"POS_NEG":cell("POS","NEG"),
                   "ZERO_ZERO":cell("ZERO","ZERO"),"NEG_POS":cell("NEG","POS"),"NEG_NEG":cell("NEG","NEG"),
                   "scenario_gt_baseline":gt,"scenario_lt_baseline":lt,"scenario_eq_baseline":eq}
    res["direction_transitions"]=trans

    # --- A/B formula properties (analytic numeric check) ---
    F=np.array([5.,10.,20.,50.,100.,200.])
    def aA(F,lam): return np.maximum(0.5,1-lam/F)
    def aB(F,lam): return F/(F+lam)
    res["formula_check"]={
      "A0_equals_1":bool(np.allclose(aA(F,0),1.0)),
      "B0_equals_1":bool(np.allclose(aB(F,0),1.0)),
      "A_range_0.5_to_1":bool((aA(F,10)>=0.5-1e-12).all() and (aA(F,10)<=1+1e-12).all()),
      "B_range_0_to_1":bool((aB(F,9)>0).all() and (aB(F,9)<=1+1e-12).all()),
      "A_monotone_decreasing_in_lambda":bool(np.all(aA(F,10)<=aA(F,5)+1e-12) and np.all(aA(F,5)<=aA(F,0)+1e-12)),
      "B_monotone_decreasing_in_lambda":bool(np.all(aB(F,9)<=aB(F,4)+1e-12) and np.all(aB(F,4)<=aB(F,1)+1e-12)),
      "weaker_F_more_attenuated_A":bool(np.all(np.diff(aA(np.sort(F),5))>=-1e-12)),
      "weaker_F_more_attenuated_B":bool(np.all(np.diff(aB(np.sort(F),4))>=-1e-12)),
      "A_floor_active_when_F_lt_2lambda":bool(np.allclose(aA(np.array([10.,15.]),10),0.5)),
    }
    print(json.dumps(res,indent=2,ensure_ascii=False))
    with open(os.path.join(OUT,"ab_verification_results.json"),"w",encoding="utf-8") as f:
        json.dump(res,f,indent=2,ensure_ascii=False)