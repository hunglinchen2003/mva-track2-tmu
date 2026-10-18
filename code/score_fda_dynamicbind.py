"""Score FDA DynamicBind poses for Lys1002 contact and 1002-1012 window rescue."""
import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, r"F:\BUB1B篩藥\analysis")
import score_docks as sd

HITS = Path(r"F:\DynamicBind\input\bub1b_fda_tcm.csv")
RUN_MUT = Path(r"F:\DynamicBind\runs\bub1b_fda_mut")
RUN_WT = Path(r"F:\DynamicBind\runs\bub1b_fda_wt")
OUT = Path(r"F:\2026_MVA_tract2_submission\analysis\fda_tcm_k1002\dynamicbind_scores.csv")
VINA = Path(r"F:\2026_MVA_tract2_submission\analysis\fda_tcm_k1002\dynamicbind_input_table.csv")


def names():
    with HITS.open(encoding="utf-8") as handle:
        return [row["name"] for row in csv.DictReader(handle)]


def vina_by_name():
    out = {}
    with VINA.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            key = row["name"].replace(" ", "_").replace("(", "").replace(")", "")
            out[key] = row
    return out


def score_run(run, wt_ca, label, drugs):
    rows = []
    for i, name in enumerate(drugs):
        folder = run / f"index{i}_idx_{i}"
        pose = sd.best_pose(folder) if folder.exists() else None
        if pose is None:
            rows.append({"name": name, "protein": label, "status": "missing"})
            continue
        rec, lig, lddt, aff = pose
        _, ca = sd.ca_coords(rec)
        n = min(len(ca), len(wt_ca))
        aligned = sd.transform(wt_ca[:n], ca[:n])
        rmsd = sd.window_rmsd(wt_ca[:n], aligned)
        d1002 = sd.min_dist(sd.ligand_xyz(lig), sd.residue_heavy(rec, 302))
        rows.append({
            "name": name,
            "protein": label,
            "status": "ok",
            "lddt": round(lddt, 3),
            "affinity": round(aff, 3),
            "window_rmsd": round(rmsd, 3),
            "dist_k1002": None if d1002 is None else round(d1002, 2),
        })
    return rows


def main():
    drugs = names()
    meta = vina_by_name()
    _, wt = sd.ca_coords(sd.WT_PDB)
    _, apo = sd.ca_coords(sd.APO_MUT)
    n = min(len(wt), len(apo))
    baseline = sd.window_rmsd(wt[:n], sd.transform(wt[:n], apo[:n]))
    mut_rows = score_run(RUN_MUT, wt, "N1002K", drugs)
    wt_rows = score_run(RUN_WT, wt, "WT", drugs)
    wt_by = {r["name"]: r for r in wt_rows}
    merged = []
    for r in mut_rows:
        twin = wt_by.get(r["name"], {})
        info = meta.get(r["name"], {})
        rescue = None
        if r.get("status") == "ok":
            rescue = round(baseline - r["window_rmsd"], 3)
        merged.append({
            "name": r["name"],
            "display_name": info.get("name", r["name"]),
            "symptom": info.get("symptom", ""),
            "vina_mut": info.get("mut", ""),
            "vina_nz": info.get("nz", ""),
            "db_affinity_mut": r.get("affinity"),
            "db_lddt_mut": r.get("lddt"),
            "window_rmsd_mut": r.get("window_rmsd"),
            "rescue_A": rescue,
            "dist_k1002_mut": r.get("dist_k1002"),
            "db_affinity_wt": twin.get("affinity"),
            "window_rmsd_wt": twin.get("window_rmsd"),
            "dist_k1002_wt": twin.get("dist_k1002"),
            "closer_than_8A": r.get("dist_k1002") is not None and r["dist_k1002"] <= 8,
            "moved_toward_wt": rescue is not None and rescue > 0.2,
        })
    fields = list(merged[0].keys())
    with OUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(merged)
    print(f"apo window RMSD {baseline:.3f} A")
    for r in sorted(merged, key=lambda x: -(x["rescue_A"] or -99)):
        print(
            f"{r['display_name'][:28]:28} aff {r['db_affinity_mut']} "
            f"rmsd {r['window_rmsd_mut']} rescue {r['rescue_A']} "
            f"K1002 {r['dist_k1002_mut']} A"
        )


if __name__ == "__main__":
    main()
