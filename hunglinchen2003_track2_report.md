# TMU Track 2 report: repositioned FDA drugs for BUB1B p.Asn1002Lys

Hung-Lin Chen, on behalf of the faculty–student team  
College of Interdisciplinary Studies and College of Pharmacy, Taipei Medical University

Hugging Face account: `hunglinchen2003`  
GitHub: https://github.com/hunglinchen2003/mva-track2-tmu  
Bilingual page: `index.html` on that repository

This is a computational and cellular hypothesis. It is not a dose and not a recommendation to treat a child.

## Track 2 form

Team name: TMU faculty–student team, Taipei Medical University, College of Interdisciplinary Studies and College of Pharmacy.

Approach: One child’s GRCh38 VCF and eight HPO terms identified compound-heterozygous *BUB1B* (p.Leu737Ter and p.Asn1002Lys) as the only biallelic functional result among the MVA genes scanned. The mechanism is loss of function of the spindle-assembly checkpoint. Full-length chaperone docking and kinase-fold rescue failed. The Track 2 list is fourteen FDA-approved drugs from the HY-L066 Lys1002 Vina screen that match cancer or growth failure, ordered by DynamicBind predicted affinity for the kinase. The lead six are pralsetinib, lapatinib, tucatinib, entrectinib, capmatinib and ponatinib. The proposed experiment is reversal of the HeLa *BUB1B*-knockdown nucleus.

Generative AI (required): Cursor session, model Grok 4.7 (xAI). The subscription tier was not written into the analysis log. The language model wrote scripts and this report from local outputs. AutoDock Vina, DynamicBind, PLIP and the nuclear counts were computed locally and were not assigned by the model. Docking ran on the local workstation. This report does not claim that variant-level data were used to train the model.

Automated versus manual: Variant scanning, Vina, DynamicBind and PLIP were automated. Manual decisions were to stop full-length docking, exclude Aurora/PLK/microtubule poisons, drop the research code c-Kit-IN-1, not treat bladder drugs as nephrocalcinosis drugs, and rank approved drugs by DynamicBind affinity rather than by distance to Lys1002.

Public data: hackathon VCF and phenotypes, Ensembl, VEP, gnomAD annotation, ClinVar, HPO, OMIM, MalaCards, PubMed, ChEMBL, ZINC20, PubChem, UniProt.  
Proprietary or vendor: MedChemExpress HY-L066 table and local Vina affinity maps, which are not republished in full. Score tables are in `results/`.

Mechanism: loss of function, not gain of function. Pathway: BubR1 and the spindle-assembly checkpoint, including the kinase fold’s role in holding PP2A-B56 at the kinetochore. Downstream: chromosome mis-segregation, mosaic aneuploidy, embryonal tumour risk (rhabdomyosarcoma in this child), growth failure, and parental pregnancy loss.

Time: the campaign ran across the hackathon window. Person-hours were not kept as one timesheet.

Pitch video: the form requires a YouTube or Vimeo URL of at most about three minutes. This file does not replace that video.

## Mechanism

The null allele is chr15:40209701 T>G, p.Leu737Ter, stop_gained, HIGH, depth 46, allelic depths 21 and 25, VAF about 54%, QUAL 708.77, PASS. The stop removes the C-terminal kinase and is a plausible nonsense-mediated-decay allele. The second allele is chr15:40220612 T>G, p.Asn1002Lys, missense, MODERATE, depth 28, allelic depths 15 and 13, VAF about 46%, QUAL 344.77, PASS. SIFT 0.00, PolyPhen-2 0.99, CADD 24.5. Residue 1002 is the N-cap of helix αG, ten residues from p.Leu1012Pro. There is no parental DNA, so trans configuration is inferred. Track 1 confidence for the pair was 0.95. No cDNA HGVS was invented.

On the kinase model, aligned on residues 956–995, Lys795 moves 1.03 Å and the Cα at 1002 moves 4.07 Å. Side-chain SASA rises from 27 to 82 Å². That fits instability of BubR1 better than a dead ATP site. Ubiquitin on the new lysine was not measured.

## Methods

Vina box: 22 Å cube centred on mutant Lys1002 NZ (−26.781, −7.289, 7.438), AutoDock Vina 1.2.7, exhaustiveness 8, seed 17. HY-L066 contributed 2,432 scored unique ligands, 1,339 flagged FDA. Fifty-three PubChem constituents of Chinese materia medica were joined by name or canonical SMILES, or docked if new. The top 100 by mutant Vina energy were kept for an oncology target or a growth-hormone or urinary annotation. Celastrol was removed for *Tripterygium* toxicity. DynamicBind used published weights, 10 inference steps, two samples, no relaxation, on UniProt 701–1050. Table 6 uses the affinity on the higher-lDDT pose. A larger number is a higher predicted affinity. It is not a kcal/mol energy and distance to Lys1002 was not a filter. PLIP 3.0.1 was run afterwards on the six lead rigid poses only.

HeLa cells, pLKO.1-puro *BUB1B* shRNA versus control, Hoechst 33342, 40×, 14 October 2026. Control N = 88, knockdown N = 43. This is wild-type knockdown, not N1002K.

## Results

The two-million-compound screen (ChEMBL + ZINC, Wager-type CNS filter, Uni-Dock 1.2 / Vina, 2,427,793 docks) left ten research poses. They are Table 5 below. None is an FDA drug, so they are not the proposal.

The approved-drug rebuild is Table 6. Nilotinib had the tightest mutant Vina score in that FDA set (−8.282 kcal/mol). The DynamicBind order is different: pralsetinib 7.99, lapatinib 7.48, tucatinib 7.07. Of 53 herbal constituents, 45 received a score. Only celastrol (−6.794 kcal/mol) entered the Vina top 100, and it was excluded for toxicity.

### Table 5. Research poses, not the Track 2 list

| Rank | Compound | Mutant | Wild type | Selectivity | NZ (Å) |
| --- | --- | ---: | ---: | ---: | ---: |
| 1 | CHEMBL1445111 | −8.64 | −5.75 | +2.90 | 2.96 |
| 2 | CHEMBL229651 | −8.09 | −6.29 | +1.80 | 3.30 |
| 3 | CHEMBL5816698 | −7.96 | −6.11 | +1.85 | 3.22 |
| 4 | CHEMBL4633401 | −7.95 | −6.02 | +1.93 | 3.39 |
| 5 | CHEMBL1255082 | −7.85 | −6.06 | +1.79 | 2.98 |
| 6 | CHEMBL3781825 | −7.83 | −6.35 | +1.48 | 3.04 |
| 7 | CHEMBL1365318 | −7.81 | −5.94 | +1.87 | 3.04 |
| 8 | CHEMBL1537224 | −7.75 | −6.31 | +1.44 | 3.00 |
| 9 | CHEMBL6015444 | −7.74 | −6.69 | +1.05 | 3.05 |
| 10 | CHEMBL2042022, notoamide B | −7.74 | −7.23 | +0.51 | 3.16 |

Energies are kcal/mol. Selectivity is wild-type score minus mutant score.

### Table 6. FDA drugs ordered by DynamicBind affinity

| Rank | Drug | Link | Vina mut | Vina WT | DB mut | DB WT | Role |
| --- | --- | --- | ---: | ---: | ---: | ---: | --- |
| 1 | Pralsetinib | RET; cancer | −7.127 | −6.646 | 7.99 | 7.85 | Lead |
| 2 | Lapatinib | EGFR/HER2; cancer | −7.181 | −6.516 | 7.48 | 7.56 | Lead |
| 3 | Tucatinib | HER2; cancer | −7.718 | −6.847 | 7.07 | 7.03 | Lead |
| 4 | Entrectinib | TRK/ROS1/ALK; cancer | −7.560 | −7.151 | 6.89 | 6.63 | Lead |
| 5 | Capmatinib | MET; cancer | −7.207 | −7.043 | 6.71 | 6.76 | Lead |
| 6 | Ponatinib | BCR-ABL and others; cancer | −7.631 | −7.287 | 6.70 | 6.67 | Lead |
| 7 | Infigratinib | FGFR; cancer | −6.790 | −5.489 | 6.60 | 6.54 | Earlier cover seed |
| 8 | Nilotinib | BCR-ABL; cancer | −8.282 | −7.104 | 6.58 | 6.68 | Tightest Vina |
| 9 | Pexidartinib | CSF1R/c-Kit; cancer | −6.984 | −6.321 | 6.43 | 6.39 | Extended |
| 10 | Regorafenib | Multikinase; cancer | −7.098 | −6.607 | 6.26 | 6.61 | Extended |
| 11 | Macimorelin | GH secretagogue; growth/cachexia | −6.507 | −6.456 | 6.20 | 6.27 | Growth arm |
| 12 | Cabozantinib | MET/VEGFR; cancer | −7.012 | −6.286 | 6.07 | 6.24 | Extended |
| 13 | Momelotinib | JAK1/2; myelofibrosis | −7.393 | −6.527 | 5.99 | 6.06 | Extended |
| 14 | Sorafenib | RAF/VEGFR; cancer | −7.063 | −6.712 | 5.91 | 6.48 | Extended |

Rigid-pose PLIP, which does not set the rank: pralsetinib, lapatinib, tucatinib and capmatinib show a hydrogen bond to K1002 (nearest 2.60, 3.09, 2.06 and 2.32 Å). Entrectinib and ponatinib do not. Their contacts are on neighbouring residues. The figure is `assets/k1002_vina_top6.png`.

HeLa knockdown: nuclei 88 to 43 (−51%), mean area 288.8 to 233.7 (−19%), micronucleus-like objects 1.1% to 4.7%, pyknosis 3.4% to 9.3%. A hit moves those counts toward control. Binding BubR1 is not required.

## Proposal

Open a knockdown nuclear screen with the fourteen approved drugs, lead six first. The cancer link is the rhabdomyosarcoma phenotype, not a claim that the tumour is driven by RET, HER2 or MET. Macimorelin is the non-kinase comparator for short stature, failure to thrive and cancer cachexia. Aurora, PLK and microtubule poisons stay off the list because they can worsen a broken checkpoint. Kinase inhibitors may still be cytotoxic. The readout is the nucleus, after a dose–response, at a screening core. Nothing here is a paediatric dose.

## Limits

No parental phase, no RNA, no patient cell, no N1002K knock-in. Predictor scores are predictors. Vina is not a Kd. DynamicBind affinity is not a Kd, and the flexible poses did not restore the wild-type 1002–1012 window; affinity is used as a binding reference only. The ubiquitin hypothesis is untested. The imaging is one field.

## 中文摘要

台北醫學大學跨領域學院與藥學院師生團隊。這個 MVA 孩子唯一的雙等位功能變異是 *BUB1B* p.Leu737Ter 與 p.Asn1002Lys，機制是紡錘體檢查點的功能失去。Track 2 提案是 14 個 FDA 核准藥，依 DynamicBind 對激酶的預測親和力排序。優先六名是 pralsetinib、lapatinib、tucatinib、entrectinib、capmatinib、ponatinib。實驗讀值是 HeLa *BUB1B* knockdown 的細胞核是否回到對照。這不是治療建議。完整中英對照在同一儲存庫的 `index.html`。
