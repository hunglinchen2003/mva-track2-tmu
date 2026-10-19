# TMU Track 2 report: repositioned FDA drugs for BUB1B p.Asn1002Lys

Authors: wang3yuan; Hung-Lin Chen (corresponding author), TMU-IIH PharmOmics  
College of Interdisciplinary Studies and College of Pharmacy, Taipei Medical University

Hugging Face accounts: `wang3yuan`; `hunglinchen2003` (corresponding author)  
GitHub: https://github.com/hunglinchen2003/mva-track2-tmu  
Bilingual page: `index.html` on that repository

This is a computational and cellular hypothesis. It is not a dose and not a recommendation to treat a child.

## Track 2 form

Team name: TMU-IIH PharmOmics, Taipei Medical University, College of Interdisciplinary Studies and College of Pharmacy.

Approach: One child’s GRCh38 VCF and eight HPO terms identified compound-heterozygous *BUB1B* (p.Leu737Ter and p.Asn1002Lys) as the only biallelic functional result among the MVA genes scanned. The mechanism is loss of function of the spindle-assembly checkpoint. Full-length chaperone docking and kinase-fold rescue failed. The Track 2 list is fourteen FDA-approved drugs from the HY-L066 Lys1002 Vina screen that match cancer or growth failure, ordered by DynamicBind predicted affinity for the kinase. The lead six are pralsetinib, lapatinib, tucatinib, entrectinib, capmatinib and ponatinib. The proposed experiment is reversal of the HeLa *BUB1B*-knockdown nucleus. Knockdown, rather than an N1002K knock-in, is deliberate: a hit should be readable in other *BUB1B* alleles and in other spindle-assembly or spindle-metabolism defects. In parallel, Connectivity Map (LINCS L1000) listed 32 signatures that raise *BUB1B* mRNA. The cell hypothesis restores the defect from both sides, message and protein. The only pair staged after each single works is doxycycline with entrectinib or pralsetinib. Aurora, PLK and CHK1 inducers stay excluded. The lead six are single-agent assay candidates, not a combination and not a paediatric regimen. A tumour somatic genotype would be required before any labelled use.

Generative AI (required): Cursor session, model Grok 4.7 (xAI). The subscription tier was not written into the analysis log. The language model wrote scripts and this report from local outputs. AutoDock Vina, DynamicBind, PLIP and the nuclear counts were computed locally and were not assigned by the model. Docking ran on the local workstation. This report does not claim that variant-level data were used to train the model.

Automated versus manual: Variant scanning, Vina, DynamicBind and PLIP were automated. Manual decisions were to stop full-length docking, exclude Aurora/PLK/microtubule poisons, drop the research code c-Kit-IN-1, not treat bladder drugs as nephrocalcinosis drugs, and rank approved drugs by DynamicBind affinity rather than by distance to Lys1002.

Public data: hackathon VCF and phenotypes, Ensembl, VEP, gnomAD annotation, ClinVar, HPO, OMIM, MalaCards, PubMed, ChEMBL, ZINC20, PubChem, UniProt.  
Proprietary or vendor: MedChemExpress HY-L066 table and local Vina affinity maps, which are not republished in full. Score tables are in `results/`.

Mechanism: loss of function, not gain of function. Pathway: BubR1 and the spindle-assembly checkpoint, including the kinase fold’s role in holding PP2A-B56 at the kinetochore. Downstream: chromosome mis-segregation, mosaic aneuploidy, embryonal tumour risk (rhabdomyosarcoma in this child), growth failure, and parental pregnancy loss.

Time: the campaign ran across the hackathon window. Person-hours were not kept as one timesheet.

Pitch video: https://youtu.be/415cWfMVpw8

## Mechanism

The null allele is chr15:40209701 T>G, p.Leu737Ter, stop_gained, HIGH, depth 46, allelic depths 21 and 25, VAF about 54%, QUAL 708.77, PASS. The stop removes the C-terminal kinase and is a plausible nonsense-mediated-decay allele. The second allele is chr15:40220612 T>G, p.Asn1002Lys, missense, MODERATE, depth 28, allelic depths 15 and 13, VAF about 46%, QUAL 344.77, PASS. SIFT 0.00, PolyPhen-2 0.99, CADD 24.5. Residue 1002 is the N-cap of helix αG, ten residues from p.Leu1012Pro. There is no parental DNA, so trans configuration is inferred. Track 1 confidence for the pair was 0.95. No cDNA HGVS was invented.

On the kinase model, aligned on residues 956–995, Lys795 moves 1.03 Å and the Cα at 1002 moves 4.07 Å. Side-chain SASA rises from 27 to 82 Å². That fits instability of BubR1 better than a dead ATP site. Ubiquitin on the new lysine was not measured.

## Methods

The Track 2 method is the FDA-approved screen in Figure 1 of the manuscript and on the bilingual page. Start from compound-heterozygous *BUB1B*, dock FDA-flagged HY-L066 drugs on Lys1002, keep cancer or growth-failure uses, order the fourteen marketed names by DynamicBind affinity, and read a *BUB1B*-knockdown nucleus. Full-length docking, fold rescue, the generator, and the ChEMBL/ZINC library are not this method.

Vina box: 22 Å cube centred on mutant Lys1002 NZ (−26.781, −7.289, 7.438), AutoDock Vina 1.2.7, exhaustiveness 8, seed 17. HY-L066 contributed 2,432 scored unique ligands, 1,339 flagged FDA. Fifty-three PubChem constituents of Chinese materia medica were joined by name or canonical SMILES, or docked if new. The top 100 by mutant Vina energy were kept for an oncology target or a growth-hormone or urinary annotation. Celastrol was removed for *Tripterygium* toxicity. DynamicBind used published weights, 10 inference steps, two samples, no relaxation, on UniProt 701–1050. Table 6 uses the affinity on the higher-lDDT pose. A larger number is a higher predicted affinity. It is not a kcal/mol energy and distance to Lys1002 was not a filter. PLIP 3.0.1 was run afterwards on the six lead rigid poses only.

HeLa cells, pLKO.1-puro *BUB1B* shRNA versus control, Hoechst 33342, 40×, 14 October 2026. Control N = 88, knockdown N = 43. This is wild-type knockdown, not N1002K.

## Results

The FDA screen is the result that matters for Track 2. Of 2,432 scored HY-L066 ligands, 1,339 were flagged FDA. The top-100 cutoff was pazopanib at −6.491 kcal/mol. Thirty-nine rows matched cancer or growth failure. Fourteen marketed drugs were ordered by DynamicBind affinity (Table 6). Nilotinib had the tightest mutant Vina score in that set (−8.282 kcal/mol). The DynamicBind order is different: pralsetinib 7.99, lapatinib 7.48, tucatinib 7.07. Of 53 herbal constituents, 45 received a score. Only celastrol (−6.794 kcal/mol) entered the Vina top 100, and it was excluded for toxicity. No herbal constituent is in Table 6.

A separate ChEMBL and ZINC screen (Wager-type CNS filter, Uni-Dock 1.2 / Vina, 2,427,793 docks) left ten research poses (Table 5). None is an FDA drug. That table is not the method and not the proposal.

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

Rigid-pose PLIP, which does not set the rank: pralsetinib, lapatinib, tucatinib and capmatinib show a hydrogen bond to K1002 (nearest 2.60, 3.09, 2.06 and 2.32 Å). Entrectinib and ponatinib do not. Their contacts are on neighbouring residues. The figure is `assets/k1002_vina_top6.png`.

HeLa knockdown: nuclei 88 to 43 (−51%), mean area 288.8 to 233.7 (−19%), micronucleus-like objects 1.1% to 4.7%, pyknosis 3.4% to 9.3%. A hit moves those counts toward control. Binding BubR1 is not required.

## Proposal

Open a knockdown nuclear screen with the fourteen approved drugs, lead six first, each drug alone. The cancer link is the rhabdomyosarcoma phenotype, not a claim that the tumour is driven by RET, HER2 or MET. Macimorelin is the non-kinase comparator for short stature, failure to thrive and cancer cachexia. Aurora, PLK and microtubule poisons stay off the list because they can worsen a broken checkpoint. Kinase inhibitors may still be cytotoxic. The readout is the nucleus, after a dose–response, at a screening core. Nothing here is a paediatric dose.

The screen uses *BUB1B* knockdown rather than an N1002K knock-in on purpose. Knock-in would test only this missense allele. Knockdown tests loss of BubR1, the state shared by biallelic *BUB1B* disease. A compound that restores nuclear count, area, micronucleus-like objects or pyknosis is then a candidate for other *BUB1B* genotypes and for other defects of spindle assembly or spindle metabolism, so the same assay can be taken to other MVA patients. Binding Lys1002 is not required. The HeLa field is still not this child’s cells, and a second genotype would be needed before a hit was treated as shared.

## Use of the lead six

The six are a screening order, not a regimen. United States labels consulted in October 2026 (Gavreto, December 2025; Rozlytrek, January 2024; Tukysa, lapatinib, Tabrecta and Iclusig through 2025) give a paediatric age to only two of them, and only with a somatic driver this germline VCF does not contain.

Pralsetinib: adults, and paediatric patients 12 years and older, with advanced or metastatic RET fusion-positive thyroid cancer that needs systemic therapy and is radioactive-iodine refractory when iodine is appropriate. Not established under 12 years, and not established for paediatric RET fusion-positive lung cancer. The label says to monitor open growth plates; animal studies showed physeal dysplasia near the clinical exposure. Entrectinib: patients older than 1 month with an NTRK fusion-positive solid tumour that is metastatic or whose resection would cause severe morbidity, that has progressed or has no satisfactory alternative, and that lacks a known acquired resistance mutation. ROS1 lung cancer is an adult indication. The label reports fractures in children, including fractures without trauma. Lapatinib, tucatinib, capmatinib and ponatinib: safety and effectiveness in paediatric patients have not been established. Ponatinib’s juvenile-rat study recorded deaths at 3 mg/kg/day, about 0.32 times the 45 mg human dose on a mg/m² basis, and the adult label has a boxed warning for arterial occlusion, venous thromboembolism, heart failure and hepatotoxicity. None of the six is labelled for *BUB1B* or for mosaic variegated aneuploidy.

A tumour would need its own somatic mutation analysis before any of these labels could apply: RET fusion, NTRK fusion, HER2, MET exon 14 skipping, or *BCR-ABL1*, matched to the drug. Germline *BUB1B* does not replace that test. This study did not sequence the tumour.

The six are not a combination. Labels do not name one as a contraindication to the others, and they do not approve the six together. Approved combinations pair one kinase inhibitor with an antibody or with chemotherapy. Lapatinib and tucatinib are sequential HER2 regimens in the clinical reports, not a simultaneous pair. Tucatinib is a strong CYP3A inhibitor and its label says to avoid CYP3A substrates when a small concentration rise can be serious. The other five are CYP3A substrates. Pralsetinib’s label says to avoid moderate or strong CYP3A inhibitors and P-gp inhibitors; lapatinib inhibits CYP3A4 and P-gp. Liver, lung and QT warnings would stack. A pair is allowed in the cell assay only after both singles rescue the nucleus. That is not a fixed-dose product and not a paediatric dose.

## Clinical compatibility and toxicity exclusion

Table 6 is the affinity order. The tier below is a separate filter for this child’s recorded phenotype: skeletal muscle atrophy, nephrocalcinosis, rhabdomyosarcoma, birth at 32 weeks, birth weight about 1 kg, short stature and failure to thrive. The eight HPO terms do not include developmental delay. The tier is not a dose and does not remove the need for tumour somatic genotyping.

First tier, carried forward. Entrectinib was developed for central-nervous-system exposure. The label gives paediatric doses, by age and body surface area, for patients older than 1 month with an NTRK-fusion solid tumour. A sarcoma is inside that indication only when an NTRK fusion is present. This child’s rhabdomyosarcoma has not been sequenced, so the match is to the labelled tumour class, not a claim that the tumour is NTRK-positive. Fractures and cognitive effects are labelled risks. Pralsetinib is the more selective name in this tier. Its advantage is single-agent use without a moderate or strong CYP3A or P-gp inhibitor. Azole antifungals are in that group. Safety under 12 years is not established, and open growth plates must be monitored.

Second tier, limited. Tucatinib has central-nervous-system activity in HER2-positive brain metastases and is more HER2-selective than lapatinib, so rash and diarrhoea are less dominant than with an EGFR/HER2 inhibitor. It is a strong CYP3A inhibitor, which makes CYP3A-substrate supportive drugs, including some antiemetics, sedatives and anti-infectives, difficult to combine with it. Capmatinib is kept because MET signalling participates in skeletal-muscle biology and has been studied in rhabdomyosarcoma. That is a biological association, not a labelled use. The interstitial-lung-disease warning, in a child born at 32 weeks, is why it stays limited.

Excluded. Ponatinib’s boxed warning covers arterial occlusion, venous thromboembolism, heart failure and hepatotoxicity. Prematurity, growth failure, muscle atrophy and nephrocalcinosis do not make that warning acceptable. Lapatinib is excluded for boxed hepatotoxicity. Its label requires each dose at least one hour before or one hour after a meal, because food, including a high-fat meal, increases exposure. The rule is separation from meals, not dosing with a high-fat meal, and that separation fits poorly with failure to thrive.

If the only question is the knockdown nucleus, any of the six may still be tested singly. The names this filter would carry toward a later clinical discussion are entrectinib and pralsetinib. Tucatinib and capmatinib are backups. Ponatinib and lapatinib are not.

## CMap drugs that raise BUB1B mRNA

A LINCS L1000-to-RNAseq file, `results/bub1b_l1000_rnaseq_up.csv`, lists 32 signatures with BUB1B fold change above 1. BUB1B is an L1000 landmark gene. Each row is one cell line, one dose and one time, with connectivity about 0.03. This is not a protein stability measurement, and it does not reorder Table 6.

Mitotic poisons in that file raise the transcript because cells collect in G2/M. PF-477736, a CHK1 inhibitor, is the top row at fold change 1.94. HMN-214 (1.27) and tozasertib (1.04) are in the file, and the accompanying CMap report adds barasertib. Ponatinib is also in the file, at 1.14. These stay off the protective list. They are not combination partners for the six.

The marketed, non-mitotic names are the repurposing set: risperidone 1.43, doxycycline 1.38, duloxetine 1.36, quetiapine 1.32, haloperidol 1.28, isradipine 1.26, bezafibrate 1.12 and vorinostat 1.07. The CMap report also cites valproic acid. That drug is not a row in the 32-line file, so no fold change is given here.

A pair would mean more mRNA from the p.Asn1002Lys allele plus one first-tier kinase drug. p.Leu737Ter is a truncation. Table 6 compounds were not shown to stabilise the protein. The pair is a cell hypothesis. Test each drug alone first. The only pair worth staging next is doxycycline with entrectinib, or doxycycline with pralsetinib. Doxycycline is marketed, its fold change is among the higher approved-drug values, and it is not a strong CYP3A inhibitor. In a young child the tetracycline warning on teeth and bone still applies, so this is not a dose.

Antipsychotics and duloxetine are not partners. Quetiapine, risperidone and haloperidol are CYP3A substrates and prolong QT. Entrectinib inhibits CYP3A and has a QT warning. Tucatinib is a strong CYP3A inhibitor. Duloxetine adds a liver warning. Isradipine is a CYP3A substrate and lowers blood pressure. Bezafibrate needs renal caution beside nephrocalcinosis. Vorinostat and valproic acid are HDAC inhibitors. Valproic acid has a boxed hepatotoxicity warning in young children. They are not add-on doses for the first-tier drugs.

## Limits

No parental phase, no RNA, no patient cell. The absent N1002K knock-in is the screening design above, not a missing allele-specific control. Predictor scores are predictors. Vina is not a Kd. DynamicBind affinity is not a Kd, and the flexible poses did not restore the wild-type 1002–1012 window; affinity is used as a binding reference only. The ubiquitin hypothesis is untested. The imaging is one field.

## 中文摘要

TMU-IIH PharmOmics，台北醫學大學跨領域學院與藥學院。作者是 wang3yuan 與 Hung-Lin Chen，通訊作者是 Hung-Lin Chen（hunglinchen2003）。這個 MVA 孩子唯一的雙等位功能變異是 *BUB1B* p.Leu737Ter 與 p.Asn1002Lys，機制是紡錘體檢查點的功能失去。Track 2 提案是 14 個 FDA 核准藥，依 DynamicBind 對激酶的預測親和力排序，並用 Connectivity Map（LINCS L1000）找會提高 BUB1B mRNA 的老藥，讓 mRNA 與蛋白並行回補缺陷。優先六名是 pralsetinib、lapatinib、tucatinib、entrectinib、capmatinib、ponatinib。實驗讀值是 HeLa *BUB1B* knockdown 的細胞核是否回到對照。不用 N1002K knock-in，是為了讓讀值對其他 *BUB1B* 基因型，以及其他紡錘體組成或代謝異常的 MVA 病人仍然適用。前六名不是複方。兒童標示只涵蓋 pralsetinib 的 12 歲以上 RET 融合甲狀腺癌，以及 entrectinib 大於 1 個月且腫瘤有 NTRK 融合的情形；另外四個藥的兒童安全性未建立。臨床上還需要腫瘤自己的體細胞突變分析。臨床相容性另排一層，不取代表 6：第一梯隊是 entrectinib 與 pralsetinib；第二梯隊是 tucatinib 與 capmatinib；ponatinib 與 lapatinib 依黑框警告排除。八個 HPO 沒有發育遲緩。Entrectinib 的兒童適應症是 NTRK 融合實體瘤，不是未定序的橫紋肌肉瘤。Lapatinib 要與餐食隔開至少一小時，不是配高脂食物。CMap L1000 另有 32 個會提高 BUB1B mRNA 的簽名。CHK1、PLK、Aurora 那一類，以及清單裡的 ponatinib，不拿來合併。已上市、非有絲分裂的老藥裡，只有 doxycycline 值得在單藥有效之後，於細胞實驗跟 entrectinib 或 pralsetinib 做一對。抗精神病藥、duloxetine、vorinostat 與 valproic acid 不做成併用。這不是治療建議。完整中英對照在同一儲存庫的 `index.html`。
