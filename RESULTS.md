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
