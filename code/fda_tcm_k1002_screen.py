"""FDA + scientific TCM constituents, Vina at Lys1002, then MVA-symptom filter.

HY-L066 mutant/wild-type scores are the completed screen
(Vina 1.2.7, 22 A box on Lys1002 NZ, exhaustiveness 8, seed 17).
Constituents that are not already in that library are prepared and docked
with the same receptor, box, and seed.
"""
from __future__ import annotations

import csv
import json
import math
import os
import subprocess
import time
import urllib.request
from pathlib import Path

import openpyxl
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
from meeko import MoleculePreparation, PDBQTWriterLegacy

ROOT = Path(r"F:\2026_MVA_tract2_submission")
OUT = ROOT / "analysis" / "fda_tcm_k1002"
OUT.mkdir(parents=True, exist_ok=True)
XLSX = ROOT / "HY-L066_FDA Approved & Pharmacopeial Drug Library_2025-02-13.xlsx"
VINA_ROOT = Path(r"F:\BUB1B篩藥\analysis\vina")
VINA = VINA_ROOT / "vina_1.2.7_win.exe"
MAPS = VINA_ROOT / "maps" / "mut"
MAPS_WT = VINA_ROOT / "maps" / "wt"
NZ = (-26.781, -7.289, 7.438)

# Characterized constituents. PubChem CID is the structure source.
# link is the MVA-related phenotype cluster used later, not a treatment claim.
TCM = [
    ("Curcumin", 969516, "Curcuma longa", "cancer"),
    ("Berberine", 2353, "Coptis chinensis", "cancer; metabolic/growth"),
    ("Artemisinin", 68827, "Artemisia annua", "cancer"),
    ("Artesunate", 6917864, "Artemisia annua", "cancer"),
    ("Baicalein", 5281605, "Scutellaria baicalensis", "cancer"),
    ("Baicalin", 64982, "Scutellaria baicalensis", "cancer"),
    ("Wogonin", 5281703, "Scutellaria baicalensis", "cancer"),
    ("Tanshinone IIA", 164676, "Salvia miltiorrhiza", "cancer"),
    ("Cryptotanshinone", 160254, "Salvia miltiorrhiza", "cancer"),
    ("Salvianolic acid B", 6451084, "Salvia miltiorrhiza", "kidney; cardiovascular"),
    ("Oridonin", 34378, "Isodon rubescens", "cancer"),
    ("Celastrol", 122724, "Tripterygium wilfordii", "cancer"),
    ("Triptolide", 107985, "Tripterygium wilfordii", "cancer"),
    ("Shikonin", 479503, "Lithospermum erythrorhizon", "cancer"),
    ("Emodin", 3220, "Rheum palmatum", "cancer; kidney"),
    ("Rhein", 10168, "Rheum palmatum", "kidney"),
    ("Honokiol", 72303, "Magnolia officinalis", "cancer"),
    ("Magnolol", 72300, "Magnolia officinalis", "cancer"),
    ("Andrographolide", 5318517, "Andrographis paniculata", "cancer; inflammation"),
    ("Matrine", 91466, "Sophora flavescens", "cancer"),
    ("Oxymatrine", 114850, "Sophora flavescens", "cancer"),
    ("Quercetin", 5280343, "widely distributed flavonoid", "cancer"),
    ("Resveratrol", 445154, "Polygonum cuspidatum", "cancer"),
    ("Epigallocatechin gallate", 65064, "Camellia sinensis", "cancer"),
    ("Luteolin", 5280445, "widely distributed flavonoid", "cancer"),
    ("Apigenin", 5280443, "widely distributed flavonoid", "cancer"),
    ("Kaempferol", 5280863, "widely distributed flavonoid", "cancer"),
    ("Genistein", 5280961, "Glycine max", "cancer"),
    ("Daidzein", 5281708, "Glycine max", "cancer"),
    ("Formononetin", 5280378, "Astragalus membranaceus", "cancer"),
    ("Icariin", 5318997, "Epimedium", "growth/bone"),
    ("Astragaloside IV", 13943297, "Astragalus membranaceus", "growth/immune"),
    ("Ginsenoside Rg3", 9918693, "Panax ginseng", "cancer"),
    ("Ginsenoside Rh2", 119307, "Panax ginseng", "cancer"),
    ("Notoginsenoside R1", 441934, "Panax notoginseng", "circulation"),
    ("Glycyrrhizin", 14982, "Glycyrrhiza uralensis", "inflammation"),
    ("Glycyrrhetinic acid", 10114, "Glycyrrhiza uralensis", "inflammation"),
    ("Isoliquiritigenin", 638278, "Glycyrrhiza uralensis", "cancer"),
    ("Paeoniflorin", 442534, "Paeonia lactiflora", "inflammation"),
    ("Puerarin", 5281807, "Pueraria lobata", "metabolic/growth"),
    ("Ferulic acid", 445858, "Ligusticum / grains", "antioxidant"),
    ("Chlorogenic acid", 1794427, "Lonicera / many herbs", "metabolic/growth"),
    ("Oleanolic acid", 10494, "widely distributed triterpene", "metabolic/liver"),
    ("Ursolic acid", 64945, "widely distributed triterpene", "cancer; muscle"),
    ("Tetramethylpyrazine", 14296, "Ligusticum chuanxiong", "kidney; circulation"),
    ("Gastrodin", 115067, "Gastrodia elata", "neurodevelopment"),
    ("Salidroside", 159278, "Rhodiola rosea", "fatigue/growth"),
    ("Cordycepin", 6303, "Cordyceps militaris", "cancer"),
    ("Osthole", 10228, "Cnidium monnieri", "cancer"),
    ("Hesperetin", 72281, "Citrus", "cancer"),
    ("Naringenin", 439246, "Citrus", "metabolic/growth"),
    ("Polydatin", 5281718, "Polygonum cuspidatum", "metabolic/growth"),
    ("Schisandrin B", 158103, "Schisandra chinensis", "liver/cytoprotection"),
]

CANCER = (
    "cancer", "tumor", "tumour", "carcinoma", "sarcoma", "leukem", "oncolog",
    "neoplas", "myeloma", "lymphoma", "rhabdomyo", "antineoplast", "anticancer",
    "anti-cancer", "antitumor", "anti-tumor", "myelofibrosis",
)
ONCO_TARGET = (
    "bcr-abl", "egfr", "erbb", "her2", "anaplastic lymphoma", "ros kinase",
    "trk", "c-kit", "c-met", "hgfr", "fgfr", "pdgfr", "vegfr", "flt3", "btk",
    "parp", "hdac", "cdk", "aromatase", "proteasome", "pi3k", "ret", "idh", "braf",
)
GROWTH = (
    "growth hormone", "somatropin", "igf", "failure to thrive", "short stature",
    "cachexia", "ghsr", "gh secretagogue",
)
KIDNEY = ("renal", "kidney", "nephro", "urolith", "urinary", "bladder")
MUSCLE = ("muscle", "myopath", "dystroph", "myasthenia", "sarcopen")
POISON = ("aurora", "polo-like", "plk1", "plk-1", "vincristine", "vinblastine", "vinorelbine", "paclitaxel", "docetaxel", "cabazitaxel", "colchicine", "barasertib", "volasertib", "alisertib")
TOXIC_TCM = {"triptolide", "celastrol"}


def norm(name: str) -> str:
    return "".join(ch for ch in name.lower() if ch.isalnum())


def load_excel():
    wb = openpyxl.load_workbook(XLSX, read_only=True, data_only=True)
    ws = wb["Compound Information"]
    by_cat = {}
    by_name = {}
    for row in ws.iter_rows(min_row=2, values_only=True):
        if not row or not row[3]:
            continue
        item = {
            "catalog": str(row[2] or ""),
            "name": str(row[3] or "").strip(),
            "synonyms": str(row[4] or ""),
            "cas": str(row[5] or ""),
            "mw": row[6],
            "nda": str(row[7] or ""),
            "target": str(row[8] or ""),
            "pathway": str(row[9] or ""),
            "approved": str(row[10] or ""),
            "pharmacopoeia": str(row[11] or ""),
            "activity": str(row[12] or ""),
        }
        by_cat[item["catalog"]] = item
        by_name.setdefault(norm(item["name"]), item)
    return by_cat, by_name


def load_vina():
    meta = json.loads((VINA_ROOT / "library_meta.json").read_text(encoding="utf-8"))
    mut = {r["ligand"]: r.get("score") for r in json.loads((VINA_ROOT / "scores_mut.json").read_text(encoding="utf-8"))}
    wt = {r["ligand"]: r.get("score") for r in json.loads((VINA_ROOT / "scores_wt.json").read_text(encoding="utf-8"))}
    rows = []
    for item in meta["ligands"]:
        key = item.get("pdbqt")
        if not key or mut.get(key) is None:
            continue
        rows.append({
            "name": item["name"],
            "catalog": item.get("catalog", ""),
            "cas": item.get("cas", ""),
            "smiles": item.get("smiles", ""),
            "target": item.get("target", ""),
            "pdbqt": key,
            "mut": mut[key],
            "wt": wt.get(key),
            "source": "HY-L066",
        })
    return rows


def pose_nz(path: Path):
    if not path.exists():
        return None
    best = None
    for line in path.read_text(errors="ignore").splitlines():
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        try:
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
        except ValueError:
            continue
        d = math.dist((x, y, z), NZ)
        best = d if best is None else min(best, d)
    return None if best is None else round(best, 2)


def blob(rec):
    return " ".join([
        rec.get("name", ""), rec.get("synonyms", ""), rec.get("target", ""),
        rec.get("pathway", ""), rec.get("activity", ""), rec.get("tcm_link", ""),
        rec.get("herb", ""),
    ]).lower()


HEME_ONCOLOGY = ("momelotinib", "ruxolitinib", "fedratinib", "pacritinib")


def tags(rec):
    text = blob(rec)
    target = (rec.get("target") or "").lower()
    name = (rec.get("name") or "").lower()
    hits = []
    onco_target = any(k in target for k in ONCO_TARGET) or any(k in name for k in HEME_ONCOLOGY)
    if onco_target:
        hits.append("cancer")
    if any(k in text for k in GROWTH):
        hits.append("growth")
    if any(k in text for k in KIDNEY):
        hits.append("kidney")
    if any(k in text for k in MUSCLE):
        hits.append("muscle")
    link = (rec.get("tcm_link") or "").lower()
    for label, keys in (
        ("cancer", ("cancer",)),
        ("growth", ("growth", "metabolic")),
        ("kidney", ("kidney",)),
        ("muscle", ("muscle",)),
    ):
        if any(k in link for k in keys) and label not in hits:
            hits.append(label)
    return hits


def incidental_anticancer(rec):
    """Library text mentions tumor activity, but the annotated target is not an oncology target."""
    if "cancer" in tags(rec):
        return False
    return any(k in blob(rec) for k in CANCER)


def is_poison(rec):
    text = blob(rec)
    return any(k in text for k in POISON)


def is_toxic_tcm(name):
    return norm(name) in {norm(x) for x in TOXIC_TCM}


def fetch_smiles(cids):
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
        + ",".join(str(c) for c in cids)
        + "/property/CanonicalSMILES,Title/JSON"
    )
    req = urllib.request.Request(url, headers={"User-Agent": "mva-fda-tcm/1.0"})
    with urllib.request.urlopen(req, timeout=60) as handle:
        data = json.loads(handle.read().decode("utf-8"))
    out = {}
    for prop in data["PropertyTable"]["Properties"]:
        smi = prop.get("CanonicalSMILES") or prop.get("ConnectivitySMILES")
        out[prop["CID"]] = smi
    return out


def prep_ligand(smiles, dest: Path):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "parse"
    frags = Chem.GetMolFrags(mol, asMols=True)
    mol = max(frags, key=lambda m: m.GetNumHeavyAtoms())
    heavy = mol.GetNumHeavyAtoms()
    if heavy < 6 or heavy > 80:
        return f"size {heavy}"
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 17
    if AllChem.EmbedMolecule(mol, params) != 0:
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, params) != 0:
            return "embed"
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=80)
    except Exception:
        pass
    setups = MoleculePreparation().prepare(mol)
    pdbqt, success, err = PDBQTWriterLegacy.write_string(setups[0])
    if not success:
        return str(err)
    dest.write_text(pdbqt if pdbqt.endswith("\n") else pdbqt + "\n", encoding="utf-8")
    return None


def dock(maps, ligand, out):
    if out.exists() and out.stat().st_size > 0:
        for line in out.read_text(errors="ignore").splitlines():
            if line.startswith("REMARK VINA RESULT:"):
                return float(line.split()[3])
    cmd = [
        str(VINA), "--maps", str(maps), "--ligand", str(ligand), "--out", str(out),
        "--exhaustiveness", "8", "--num_modes", "1", "--cpu", "1", "--seed", "17", "--verbosity", "0",
    ]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=100)
    except subprocess.TimeoutExpired:
        if out.exists():
            out.unlink()
        return None
    text = out.read_text(errors="ignore") if out.exists() else ""
    for src in (text, proc.stdout or ""):
        for line in src.splitlines():
            if line.startswith("REMARK VINA RESULT:"):
                return float(line.split()[3])
    return None


def main():
    by_cat, by_name = load_excel()
    base = load_vina()
    print(f"scored HY-L066 ligands {len(base)}", flush=True)
    known_smiles = set()
    known_names = set()
    for rec in base:
        known_names.add(norm(rec["name"]))
        mol = Chem.MolFromSmiles(rec["smiles"]) if rec["smiles"] else None
        if mol is not None:
            known_smiles.add(Chem.MolToSmiles(mol))
        info = by_cat.get(rec["catalog"]) or by_name.get(norm(rec["name"]))
        if info:
            rec.update({k: info.get(k, rec.get(k)) for k in ("synonyms", "nda", "target", "pathway", "approved", "pharmacopoeia", "activity", "cas")})
        approved = rec.get("approved") or ""
        nda = rec.get("nda") or ""
        rec["fda"] = ("FDA" in approved.upper()) or bool(nda.strip())
        rec["nz"] = pose_nz(VINA_ROOT / "out_mut" / rec["pdbqt"])

    print("fetching PubChem SMILES", flush=True)
    smiles = fetch_smiles([cid for _, cid, _, _ in TCM])
    tcm_dir = OUT / "tcm_ligands"
    tcm_dir.mkdir(exist_ok=True)
    failed_path = OUT / "dock_failed.txt"
    failed_names = set()
    if failed_path.exists():
        failed_names = {line.strip() for line in failed_path.read_text(encoding="utf-8").splitlines() if line.strip()}
    tcm_rows = []
    for name, cid, herb, link in TCM:
        if name in failed_names:
            print(f"skip previously failed {name}", flush=True)
            continue
        smi = smiles.get(cid)
        if not smi:
            print("no smiles", name, flush=True)
            continue
        mol = Chem.MolFromSmiles(smi)
        can = Chem.MolToSmiles(mol) if mol is not None else ""
        existing = None
        if norm(name) in known_names or can in known_smiles:
            for rec in base:
                if norm(rec["name"]) == norm(name) or rec.get("smiles") == can:
                    existing = rec
                    break
        if existing:
            existing["source"] = "HY-L066; TCM constituent"
            existing["herb"] = herb
            existing["tcm_link"] = link
            existing["pubchem"] = cid
            print(f"already scored {name} mut={existing['mut']}", flush=True)
            continue
        if mol is not None:
            heavy = mol.GetNumHeavyAtoms()
            nrot = rdMolDescriptors.CalcNumRotatableBonds(mol)
            if heavy > 45 or nrot > 10:
                print(f"skip flexible {name} heavy={heavy} rot={nrot}", flush=True)
                continue
        lig = tcm_dir / f"{norm(name)}.pdbqt"
        err = prep_ligand(smi, lig)
        if err:
            print(f"prep fail {name}: {err}", flush=True)
            continue
        t0 = time.time()
        mut = dock(MAPS, lig, tcm_dir / f"{norm(name)}_mut.pdbqt")
        wt = dock(MAPS_WT, lig, tcm_dir / f"{norm(name)}_wt.pdbqt")
        if mut is None or wt is None:
            print(f"dock timeout or no score {name} {time.time()-t0:.0f}s", flush=True)
            failed = OUT / "dock_failed.txt"
            with failed.open("a", encoding="utf-8") as handle:
                handle.write(name + "\n")
            continue
        row = {
            "name": name, "catalog": "", "cas": "", "smiles": can, "target": "",
            "pdbqt": lig.name, "mut": mut, "wt": wt, "source": "TCM constituent",
            "herb": herb, "tcm_link": link, "pubchem": cid, "fda": False,
            "approved": "scientific TCM constituent, not an FDA drug",
            "pathway": "", "activity": link, "synonyms": herb, "nda": "",
            "nz": pose_nz(tcm_dir / f"{norm(name)}_mut.pdbqt"),
        }
        tcm_rows.append(row)
        print(f"docked {name} mut={mut:.3f} wt={wt:.3f} nz={row['nz']} {time.time()-t0:.0f}s", flush=True)

    pool = [r for r in base if r.get("fda")] + tcm_rows
    # TCM constituents found inside HY-L066 are already in `base` and may not be FDA.
    for rec in base:
        if rec.get("tcm_link") and rec not in pool:
            pool.append(rec)
    pool.sort(key=lambda r: r["mut"])
    top = pool[:100]
    for i, rec in enumerate(top, 1):
        rec["rank"] = i
        rec["selectivity"] = None if rec.get("wt") is None else round(rec["wt"] - rec["mut"], 3)
        rec["symptom"] = "; ".join(tags(rec))
        rec["mitotic_poison"] = is_poison(rec)
        rec["severe_tcm_toxicity"] = is_toxic_tcm(rec["name"])

    fields = [
        "rank", "name", "source", "fda", "approved", "catalog", "cas", "pubchem", "herb",
        "mut", "wt", "selectivity", "nz", "target", "pathway", "symptom",
        "mitotic_poison", "severe_tcm_toxicity", "tcm_link", "activity",
    ]
    def write(path, rows):
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)

    write(OUT / "top100_vina.csv", top)
    kept = [r for r in top if r["symptom"] and not r["mitotic_poison"] and not r["severe_tcm_toxicity"]]
    write(OUT / "mva_symptom_kept.csv", kept)
    dropped_poison = [r for r in top if r["symptom"] and (r["mitotic_poison"] or r["severe_tcm_toxicity"])]
    write(OUT / "symptom_related_but_not_protective.csv", dropped_poison)
    # DynamicBind input: tighter poses first, cap at 12 so the GPU job can finish.
    db = sorted(kept, key=lambda r: (r["nz"] if r["nz"] is not None else 99, r["mut"]))[:12]
    write(OUT / "dynamicbind_input_table.csv", db)
    with (OUT / "dynamicbind_ligands.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["ligand", "name"])
        for rec in db:
            writer.writerow([rec["smiles"], rec["name"].replace(" ", "_")])
    summary = {
        "hy_l066_scored": len(base),
        "fda_in_scored_set": sum(1 for r in base if r.get("fda")),
        "tcm_newly_docked": len(tcm_rows),
        "pool": len(pool),
        "top100": len(top),
        "symptom_kept": len(kept),
        "excluded_as_mitotic_or_toxic": len(dropped_poison),
        "dynamicbind_n": len(db),
        "box": "22 A cube centered on Lys1002 NZ (-26.781, -7.289, 7.438)",
        "vina": "1.2.7 exhaustiveness 8 seed 17",
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print("--- top 15 ---")
    for rec in top[:15]:
        print(f"{rec['rank']:3} {rec['mut']:7.3f} nz={rec['nz']} {rec['name'][:40]:40} {rec['symptom'] or '-'}")
    print("--- dynamicbind ---")
    for rec in db:
        print(f"{rec['name']:30} mut={rec['mut']:.3f} nz={rec['nz']} {rec['symptom']}")


if __name__ == "__main__":
    main()
