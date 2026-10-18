"""Re-tag the saved top 100 with oncology targets and MVA symptom links."""
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fda_tcm_k1002_screen as screen

OUT = screen.OUT
META = json.loads((screen.VINA_ROOT / "library_meta.json").read_text(encoding="utf-8"))
by_cat = {item.get("catalog"): item for item in META["ligands"]}
by_name = {screen.norm(item.get("name", "")): item for item in META["ligands"]}

rows = list(csv.DictReader((OUT / "top100_vina.csv").open(encoding="utf-8")))
for rec in rows:
    rec["mitotic_poison"] = screen.is_poison(rec)
    rec["severe_tcm_toxicity"] = screen.is_toxic_tcm(rec["name"])
    rec["symptom"] = "; ".join(screen.tags(rec))
    rec["incidental_anticancer"] = screen.incidental_anticancer(rec)
    src = by_cat.get(rec["catalog"]) or by_name.get(screen.norm(rec["name"]))
    rec["smiles"] = (src or {}).get("smiles", "")
    rec["nz"] = float(rec["nz"]) if rec["nz"] else None
    rec["mut"] = float(rec["mut"])

fields = list(rows[0].keys())
kept = [r for r in rows if r["symptom"] and not r["mitotic_poison"] and not r["severe_tcm_toxicity"]]
caution = [r for r in rows if (r["symptom"] or r["incidental_anticancer"]) and (r["mitotic_poison"] or r["severe_tcm_toxicity"])]
incidental = [r for r in rows if r["incidental_anticancer"] and not r["symptom"]]
# Library "FDA" flag on a research code, not a market name.
not_market_name = {screen.norm("c-Kit-IN-1")}
near = [r for r in kept if r["nz"] is not None and r["nz"] <= 4.2 and screen.norm(r["name"]) not in not_market_name]
near.sort(key=lambda r: r["mut"])
db = near[:12]
for extra in ("Macimorelin (acetate)", "Infigratinib"):
    if any(r["name"] == extra for r in db):
        continue
    match = next((r for r in kept if r["name"] == extra), None)
    if match:
        db.append(match)

def write(path, data):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(data)

write(OUT / "top100_vina.csv", rows)
write(OUT / "mva_symptom_kept.csv", kept)
write(OUT / "symptom_related_but_not_protective.csv", caution)
write(OUT / "incidental_anticancer_mention.csv", incidental)
write(OUT / "dynamicbind_input_table.csv", db)
with (OUT / "dynamicbind_ligands.csv").open("w", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["ligand", "name"])
    for rec in db:
        writer.writerow([rec["smiles"], rec["name"].replace(" ", "_").replace("(", "").replace(")", "")])

summary = json.loads((OUT / "summary.json").read_text(encoding="utf-8"))
summary["symptom_kept"] = len(kept)
summary["excluded_as_mitotic_or_toxic"] = len(caution)
summary["dynamicbind_n"] = len(db)
summary["incidental_anticancer_only"] = len(incidental)
summary["dynamicbind_rule"] = "primary oncology target, growth, kidney, or muscle link; not a mitotic poison; heavy-atom distance to Lys1002 NZ <= 4.2 A; best 12 mutant Vina scores"
(OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
print("kept", len(kept), "caution", len(caution), "incidental", len(incidental), "db", len(db))
print("--- incidental ---")
for rec in incidental:
    print(rec["rank"], rec["name"])
for rec in kept:
    print(f"{rec['rank']:>3} {rec['mut']:7.3f} nz={rec['nz']} {rec['symptom']:22} {rec['name']}")
print("--- db ---")
for rec in db:
    print(rec["name"], rec["mut"], rec["nz"], rec["symptom"])
