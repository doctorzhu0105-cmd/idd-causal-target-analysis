"""QWEN independent numerical verifier, adapted for isolated inputs and outputs."""
import numpy as np, pandas as pd, json, os

def run(source, output_dir):
    SRC = str(source)
    OUT = str(output_dir)
    df = pd.read_csv(SRC, sep="\t")
    out = {}
    out["n_rows"] = len(df)
    out["n_genes"] = int(df["ensembl_gene_id"].nunique())
    out["endpoints"] = sorted(df["endpoint"].unique().tolist())
    out["K_from_values"] = sorted(df["K_from"].unique().tolist())
    out["K_to_min"] = int(df["K_to"].min()); out["K_to_max"] = int(df["K_to"].max())
    out["all_K1_to_multi"] = bool((df["K_transition"]=="K1_to_multi").all())

    # --- Identity residual (independent recompute) ---
    phi_from = np.ones(len(df), dtype=float)
    phi_to = np.ones(len(df), dtype=float)
    from_multi = df["K_from"].to_numpy() >= 2
    to_multi = df["K_to"].to_numpy() >= 2
    phi_from[from_multi] = np.maximum(1.0, df["Q_from"].to_numpy()[from_multi] / (df["K_from"].to_numpy()[from_multi]-1))
    phi_to[to_multi] = np.maximum(1.0, df["Q_to"].to_numpy()[to_multi] / (df["K_to"].to_numpy()[to_multi]-1))
    out["phi_from_all_one"] = bool(np.allclose(phi_from,1.0))
    out["phi_to_matches_stored"] = float(np.max(np.abs(phi_to - df["phi_to_indep"])))
    se_fe_from = 1.0/np.sqrt(df["D_from"]); se_fe_to = 1.0/np.sqrt(df["D_to"])
    out["se_fe_from_maxerr"] = float(np.max(np.abs(se_fe_from - df["SE_FE_from_indep"])))
    out["se_fe_to_maxerr"]   = float(np.max(np.abs(se_fe_to   - df["SE_FE_to_indep"])))
    se_re_from = se_fe_from*np.sqrt(phi_from); se_re_to = se_fe_to*np.sqrt(phi_to)
    out["se_re_from_maxerr"] = float(np.max(np.abs(se_re_from - df["SE_RE_from_indep"])))
    out["se_re_to_maxerr"]   = float(np.max(np.abs(se_re_to   - df["SE_RE_to_indep"])))
    lhs = np.log(se_re_to/se_re_from)
    rhs = np.log(se_fe_to/se_fe_from) + 0.5*np.log(phi_to/phi_from)
    resid = lhs - rhs
    out["identity_max_abs_resid"] = float(np.max(np.abs(resid)))
    out["stored_resid_max"] = float(np.max(np.abs(df["log_decomp_residual_indep"])))
    # denominator check: phi uses K-1
    out["phi_to_ge1_all"] = bool((phi_to>=1.0-1e-12).all())
    out["phi_to_gt1_count"] = int((phi_to>1.0+1e-12).sum())
    out["phi_to_eq1_count"] = int((np.abs(phi_to-1.0)<=1e-12).sum())

    # --- D component sign ---
    Dcomp = -0.5*np.log(df["D_to"]/df["D_from"])
    out["D_component_negative_count"] = int((Dcomp<0).sum())
    out["D_component_negative_pct"] = float((Dcomp<0).mean()*100)

    # --- direction counts ---
    abs_re = se_re_to - se_re_from
    re_dec = int((abs_re < -1e-15).sum()); re_inc = int((abs_re > 1e-15).sum()); re_eq=int((np.abs(abs_re)<=1e-15).sum())
    out["re_dec"]=re_dec; out["re_inc"]=re_inc; out["re_eq"]=re_eq
    out["re_dec_pct"]=re_dec/len(df)*100; out["re_inc_pct"]=re_inc/len(df)*100
    rel_re = abs_re/se_fe_from/np.sqrt(phi_from)*100  # = abs_re/se_re_from*100
    rel_re = (se_re_to-se_re_from)/se_re_from*100
    out["med_rel_SE_RE_pct"]=float(np.median(rel_re))
    out["q25_rel_SE_RE_pct"]=float(np.percentile(rel_re,25)); out["q75_rel_SE_RE_pct"]=float(np.percentile(rel_re,75))
    rel_fe=(se_fe_to-se_fe_from)/se_fe_from*100
    out["med_rel_SE_FE_pct"]=float(np.median(rel_fe))
    out["med_rel_VAR_FE_pct"]=float(np.median(((df["VY_to"]-df["VY_from"])/df["VY_from"])*100)) if "VY_to" in df else None

    # --- endpoint-specific ---
    ep_stats={}
    for ep in out["endpoints"]:
        s=df[df["endpoint"]==ep]
        a_re=(s["SE_RE_to_indep"]-s["SE_RE_from_indep"])
        ep_stats[ep]={"n":len(s),"re_dec":int((a_re<-1e-15).sum()),"re_inc":int((a_re>1e-15).sum()),
                      "re_dec_pct":float((a_re<-1e-15).mean()*100),
                      "med_rel_SE_RE_pct":float(np.median((s["SE_RE_to_indep"]-s["SE_RE_from_indep"])/s["SE_RE_from_indep"]*100)),
                      "med_rel_SE_FE_pct":float(np.median((s["SE_FE_to_indep"]-s["SE_FE_from_indep"])/s["SE_FE_from_indep"]*100)),
                      "med_D_comp":float(np.median(-0.5*np.log(s["D_to"]/s["D_from"]))),
                      "med_phi_comp":float(np.median(0.5*np.log(s["phi_to_indep"]/s["phi_from_indep"])))}
    out["endpoint_stats"]=ep_stats

    # --- three regimes (manuscript claim verification) ---
    DbDa = df["D_to"]/df["D_from"]
    reg1 = (np.abs(phi_to-1.0)<=1e-12)
    reg3 = (phi_to > DbDa + 1e-12)
    reg2 = (~reg1) & (~reg3)
    def med_rel(mask):
        s=df[mask]; return float(np.median((s["SE_RE_to_indep"]-s["SE_RE_from_indep"])/s["SE_RE_from_indep"]*100))
    out["regime1_n"]=int(reg1.sum()); out["regime1_med_rel_SE_RE_pct"]=med_rel(reg1)
    out["regime2_n"]=int(reg2.sum()); out["regime2_med_rel_SE_RE_pct"]=med_rel(reg2)
    out["regime3_n"]=int(reg3.sum()); out["regime3_med_rel_SE_RE_pct"]=med_rel(reg3)
    out["regime1_pct"]=float(reg1.mean()*100); out["regime2_pct"]=float(reg2.mean()*100); out["regime3_pct"]=float(reg3.mean()*100)

    # --- gene concordance ---
    g=df.assign(is_dec=(df["SE_RE_direction_indep"]=="decrease").astype(int)).groupby("ensembl_gene_id")["is_dec"].sum()
    out["gene_re_dec_distribution"]={str(k):int(v) for k,v in g.value_counts().sort_index().items()}
    gp=df.groupby("ensembl_gene_id")["phi_to_gt_1_flag"].sum()
    out["gene_phi_gt1_distribution"]={str(k):int(v) for k,v in gp.value_counts().sort_index().items()}

    # --- heterogeneity overwhelm flag ---
    phicomp=0.5*np.log(phi_to/phi_from)
    out["phi_overwhelm_count"]=int((phicomp > -Dcomp + 1e-15).sum())

    print(json.dumps(out, indent=2, ensure_ascii=False))
    with open(os.path.join(OUT, "wp1_verification_results.json"), "w", encoding="utf-8") as f:
        json.dump(out,f,indent=2,ensure_ascii=False)