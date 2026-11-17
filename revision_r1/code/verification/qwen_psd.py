"""QWEN independent numerical verifier, adapted for isolated inputs and outputs."""
import numpy as np, pandas as pd, os, json, hashlib

def run(manifest, output_dir):
    P0 = os.path.realpath(os.path.dirname(str(manifest)))
    OUT = str(output_dir)
    man = pd.read_csv(os.path.join(P0,"P0_MATRIX_MANIFEST.tsv"), sep="\t")
    man = man[man["status"]=="estimable"]
    # pick a spread of K sizes
    paths = sorted(set(os.path.join(P0, p.split("/",1)[1]) for p in man["matrix_file"]))
    expected = {}
    for _, source_row in man.iterrows():
        file_path = os.path.realpath(os.path.join(P0, source_row["matrix_file"].split("/",1)[1]))
        if os.path.commonpath([P0, file_path]) != P0:
            raise ValueError("Matrix path escapes input directory")
        old_digest = expected.setdefault(file_path, source_row["matrix_file_sha256"])
        if old_digest != source_row["matrix_file_sha256"]:
            raise ValueError("Conflicting matrix file hashes")
    results=[]
    tested=0
    for p in paths:
        if tested>=40: break
        try:
            real = os.path.realpath(p)
            with open(real, "rb") as fh:
                digest = hashlib.sha256(fh.read()).hexdigest()
            if digest != expected.get(real):
                raise ValueError("Matrix hash mismatch")
            d=np.load(real, allow_pickle=False)
            R=d["R"].astype(float); K=R.shape[0]
            if K<2 or K>40: continue
        except Exception as exc:
            raise ValueError("Could not load an eligible matrix file") from exc
        tested+=1
        # use unit D for structural test (D cancels in congruence for PSD sign)
        D=np.eye(K)
        R_reg=0.95*R+0.05*np.eye(K)
        def joint(Rx,Ry,Rov,rho):
            Sx=D@Rx@D; Oy=D@Ry@D; C=rho*(D@Rov@D)
            return np.block([[Sx,C],[C.T,Oy]])
        def mineig(M): return float(np.min(np.linalg.eigvalsh(M)))
        row={"path":os.path.basename(p),"K":K}
        for rho in [0.0,0.0177,0.10,1.0]:
            row[f"sharedR_mineig_rho{rho}"]=mineig(joint(R,R,R,rho))
            row[f"R5mismatch_mineig_rho{rho}"]=mineig(joint(R_reg,R_reg,R,rho))
        # Schur complement test for R5 mismatch at rho=1: Oy - C' Sx^+ C
        rho=1.0; Sx=D@R_reg@D; Oy=D@R_reg@D; C=rho*(D@R@D)
        Schur=Oy - C.T @ np.linalg.pinv(Sx) @ C
        row["R5_schur_mineig_rho1"]=mineig(Schur)
        results.append(row)

    df=pd.DataFrame(results)
    out={
     "n_matrices_tested":len(df),
     "sharedR_min_eig_over_all_rho1":float(df["sharedR_mineig_rho1.0"].min()),
     "sharedR_any_negative_rho1":bool((df["sharedR_mineig_rho1.0"]<-1e-9).any()),
     "sharedR_min_eig_rho0.10":float(df["sharedR_mineig_rho0.1"].min()),
     "R5mismatch_min_eig_rho1":float(df["R5mismatch_mineig_rho1.0"].min()),
     "R5mismatch_any_negative_rho1":bool((df["R5mismatch_mineig_rho1.0"]<-1e-9).any()),
     "R5mismatch_n_negative_rho1":int((df["R5mismatch_mineig_rho1.0"]<-1e-9).sum()),
     "R5mismatch_min_eig_rho0.10":float(df["R5mismatch_mineig_rho0.1"].min()),
     "R5mismatch_any_negative_rho0.10":bool((df["R5mismatch_mineig_rho0.1"]<-1e-9).any()),
     "R5mismatch_n_negative_rho0.10":int((df["R5mismatch_mineig_rho0.1"]<-1e-9).sum()),
     "R5_schur_min_eig_rho1":float(df["R5_schur_mineig_rho1"].min()),
     "R5_schur_any_negative_rho1":bool((df["R5_schur_mineig_rho1"]<-1e-9).any()),
    }
    # show a few concrete violating examples
    viol=df[df["R5mismatch_mineig_rho1.0"]<-1e-9].head(5)
    out["example_violations_rho1"]=viol[["path","K","R5mismatch_mineig_rho1.0","R5mismatch_mineig_rho0.1"]].to_dict("records")
    print(json.dumps(out,indent=2))
    df.to_csv(os.path.join(OUT, "wp3_psd_actual_matrix_results.tsv"),sep="	",index=False)
    with open(os.path.join(OUT, "wp3_psd_summary.json"), "w", encoding="utf-8") as f:
        json.dump(out,f,indent=2)