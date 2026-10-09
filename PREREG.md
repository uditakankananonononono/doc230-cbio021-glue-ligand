# PREREG: CBIO021 GlueForge (Molecular Glue Discovery Pipeline), ligand-only partial slice
Written and committed before any model is scored. Date 2026-10-09.
Scope: PARTIAL. The parent is a molecular-glue discovery pipeline that needs protein structure modeling (ternary complexes). This repo has NO structure model and makes no claim about the structural core. It tests only what public ligand data can: can chemistry alone predict glue properties in MolGlueDB, under scaffold-cold validation?
Data: MolGlueDB full CSV (1840 rows; SHA256 in data_manifest). Dedupe on StdInChIKey (keep first). Features: ECFP4 2048-bit (RDKit, radius 2) + RDKit descriptors. Groups: Bemis-Murcko scaffold (generic-scaffold-free exact SMILES; empty scaffold gets its own group per molecule).
Split: 5-fold GroupKFold on scaffold ("cold scaffold"), seed-free deterministic. Model: RandomForest 300 trees (seed 0, class_weight balanced); baseline models listed per direction.
## Directions and gates (fixed now)
D1 IsActive (1 vs 0): cold-scaffold OOF AUROC. Gate: AUROC >= 0.70 AND beats a MolWeight+tPSA+clogP-only RF by >= 0.03.
D2 Recruiting protein, CRBN vs non-CRBN (rows with RecruitingProtein not Unknown): OOF AUPRC for non-CRBN. Gate: AUPRC >= prevalence + 0.20.
D3 Primary target, top 5 target classes (rows in them): OOF macro-F1. Gate: macro-F1 >= majority-class macro-F1 + 0.15.
D4 Label leakage control: D1 with labels shuffled within fold (seed 0..4). Gate: mean shuffled AUROC in [0.45, 0.55]. Also report max train-test Tanimoto (ECFP4) per fold and the share of test molecules with max sim >= 0.85 to train.
D5 Metadata-only baseline: D1 and D2 from Pharmacophore + Core + ResearchStage one-hot only. Gate (informational): fingerprint model adds >= 0.03 over metadata-only for D1.
D6 Seed stability: D1 repeated with RF seeds 0-4. Gate: AUROC std <= 0.02.
Honest reporting: every gate reported pass/fail; failures stay failures. No gate edits after scoring. Amendments only before the run and noted here.
Limits: database-curated labels (publication bias, IsActive definition per source), activity labels are heterogeneous, no structure, no experimental validation.
