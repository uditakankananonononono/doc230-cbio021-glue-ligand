# RESULTS: CBIO021 GlueForge, ligand-only partial slice
Data: MolGlueDB full CSV (1840 rows, 1629 after InChIKey dedupe and parse, 872 Bemis-Murcko scaffolds). 5-fold cold-scaffold CV, RF 300 trees on ECFP4 + RDKit descriptors. Protocol PREREG.md (Amendment 1 disclosed: first run had a 'nan' target-class bug and was kept in results/first_run; D7 robustness checks were added after seeing a suspiciously strong D2).
| Dir | Result | Gate |
|---|---|---|
| D1 IsActive | AUROC 0.847 (physchem-only 0.703, metadata-only 0.478, prevalence 0.81) | met (>=0.70 and +0.03 over physchem) |
| D2 CRBN vs non-CRBN recruiter | AUPRC 0.998 (prevalence 0.279; physchem-only 0.767; metadata-only 0.866) | met, but see caveat |
| D3 top-5 primary targets (n 511) | macro-F1 0.920 vs majority 0.095 | met, but see caveat |
| D4 shuffled-label control | AUROC 0.496 mean (0.468-0.547); mean max Tanimoto test-to-train 0.707; 14.8% of test molecules >= 0.85 | met |
| D5 fingerprint minus metadata-only on D1 | +0.369 | met (informational) |
| D6 seed stability | AUROC std 0.0009 | met |
| D7 near-duplicate-excluded (test sim < 0.85) | D1 AUROC 0.851 (n 1388); D2 AUPRC 0.9975 vs prevalence 0.298 (n 1367) | met |
All 6 gates plus D7 are met on paper. What that means: 
- Strong scores here mostly reflect that glues from one paper or chemotype series share labels (target, recruiter, activity). Cold-scaffold splitting removes identical scaffolds but not close analogs: 14.8% of test molecules have a train neighbour at Tanimoto >= 0.85, and D7 shows the result survives their removal, yet series-level structure still exists below 0.85.
- D2 is easy: CRBN glutarimide glues versus VHL/DCAF/other recruiters are different chemical classes, and even metadata-only reaches AUPRC 0.866. This is classification of known chemotypes, not discovery of new glues.
- D1 IsActive labels are heterogeneous across sources; prevalence 0.81.
- No structure, no ternary-complex modelling, no experimental validation. This does not test the parent's core claim (undruggable-target glue discovery) at all.
Caveats on counts: gates are my own prereg thresholds, not published benchmarks. D3 majority macro-F1 is the macro-F1 of always predicting the largest class (0.095), so the +0.15 bar was easy.

## Extension D8-D10 (PREREG_EXT.md and src/ext_run.py committed together at 3ca147e before any score; results/ext_results.json, ext_run.log)
Same pipeline scored all splits in one script. Scaffold-cold reference: D1 AUROC 0.847 (share of test molecules with a train neighbour >=0.85: see ext_results.json), D2 AUPRC 0.998 (prevalence 0.279).
| Dir | Result | Gate |
|---|---|---|
| D8 cluster-cold D1 (Butina 0.60) | AUROC 0.784; share test with train neighbour >=0.85: 1.0% | MET (>=0.70) |
| D9 source-cold D1 (142 publications) | AUROC 0.630 (physchem-only 0.588; prevalence 0.807); share >=0.85: 2.5% | NOT MET (>=0.65) |
| D10 source-cold D2 | AUPRC 0.974 vs prevalence 0.279, AUROC 0.992 | MET (>= prev+0.10) |
Reading: removing close analogs (D8) cuts D1 from 0.847 to 0.784, and holding out whole publications (D9) cuts it to 0.630, missing my gate. So most of the D1 signal in the original result came from shared chemotype series and papers, as the earlier caveat suspected. The recruiter label (D10) stays easy even source-cold, since it separates chemical classes. Limits: one deterministic split per direction (no variance), uneven source group sizes (1 to 97 molecules), uneven fold sizes, labels are database-curated. Ligand-only, no structure, the parent's core remains unbuilt. With this batch the repo has D1-D10.
