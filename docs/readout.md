# Readout

The plan is [design.md](design.md), committed before any outcome was linked to a feature. Checks added afterwards are marked post hoc.

## Short version

Across 438 human targets of 1,332 approved drugs, boxed warnings and withdrawals cluster by target: two single-target drugs on the same target share their warning status 72.8% of the time against 60.0% by chance, so a target-level safety score has something to find. But the evidence available before a drug exists finds little of it. Genetic constraint, knockout mice, tissue expression, cancer-driver status and paralogs each reached AUROC 0.51 to 0.56, and their equal-weight combination added 0.015 (interval -0.001 to 0.031) to a model that only counts how many approved drugs act on the target, which alone reached 0.70. Learned weights did worse than equal weights. Evidence built from experience with drugs did better, as expected and partly for circular reasons: curated safety events from clinical reviews beat those from experimental sources by 0.043, and FAERS signal counts reached 0.67 but added nothing beyond drug count. One mechanistic pattern held: warnings land more often in organs where the target is expressed (0.15 SD, p = 0.033), most clearly for psychiatric, neurological, skin, vascular, kidney and heart warnings, and not for liver, consistent with drug liver injury coming mostly from metabolism. Claude wrote safety rationales for 24 targets from an evidence pack; all 222 sentences cited valid evidence and all 264 numbers in factual sentences traced to the cited rows.

## Data

Open Targets Platform release 26.09: drug molecules, mechanisms of action and ChEMBL drug warnings; target annotations, target prioritisation factors, curated target safety events, mouse phenotypes and Reactome pathways; significant FAERS adverse reactions (openFDA, likelihood ratio test) as processed by Open Targets. GTEx median expression came from the Open Targets API for the 148 targets in the organ analysis.

## Cohort

1,332 approved drugs with a single-protein mechanism (1,051 small molecules, 149 proteins, 121 antibodies, 10 oligosaccharides, 1 oligonucleotide), of which 397 have a boxed warning or a withdrawal (603 boxed-warning records and 318 withdrawal records). They act on 438 targets in 14 ChEMBL target classes; 211 targets (48%) have at least one warned drug. A spot check against known cases came out as expected: rofecoxib to PTGS2 with its 2004 withdrawal, infliximab and adalimumab to TNF with infection warnings, cerivastatin withdrawn while atorvastatin carries no warning.

## P1: are warnings shared within targets?

Among 913 drugs that act on exactly one target, on 157 targets with at least two such drugs, 6,326 same-target pairs: 72.8% agree on warning status, against 60.0% when warning labels are shuffled across drugs (p = 0.0001, the smallest value 10,000 shuffles can give). Met.

## P2 to P6: which evidence ranks targets

| Score | AUROC (95% CI) |
|---|---|
| Drug count (log number of approved drugs) | 0.699 (0.655 to 0.748) |
| FAERS signal terms, most among the target's drugs | 0.673 (0.622 to 0.721) |
| Curated safety event, clinical sources | 0.598 (0.555 to 0.642) |
| Curated safety event, experimental sources | 0.554 (0.514 to 0.594) |
| Equal-weight pre-drug score | 0.567 (0.513 to 0.619) |
| mouse knockout score | 0.558 (0.502 to 0.613) |
| tissue distribution | 0.551 (0.498 to 0.600) |
| tissue specificity | 0.550 (0.494 to 0.604) |
| genetic constraint | 0.532 (0.478 to 0.584) |
| cancer driver gene | 0.530 (0.498 to 0.563) |
| paralog identity | 0.511 (0.473 to 0.548) |

Cross-validated models, folds grouped by target class:

| Model | AUROC |
|---|---|
| M0 drug count | 0.673 |
| M1 drug count + equal-weight pre-drug score | 0.688 |
| M2 drug count + learned weights on the six factors | 0.670 |
| M3 M1 + curated-event indicators | 0.688 |

- P2 met: FAERS 0.673.
- P3 met: drug count 0.699.
- P4 not met: M1 minus M0 was 0.015 (-0.001 to 0.031). In the single-target cohort (306 targets, 135 warned) the gain was 0.027 (-0.006 to 0.066).
- P5 met: learned weights scored 0.018 below equal weights (-0.047 to 0.010).
- P6 met on the point estimate: clinical minus experimental 0.043 (-0.008 to 0.093); the interval includes 0.

M0 scores 0.673 against 0.699 for the raw drug count because pooled out-of-fold predictions come from models fitted without the held-out target class; all model comparisons share that footing.

## P7: do warnings fall where the target is expressed?

148 targets have a warning in a toxicity class with a GTEx organ. For each class, targets with that warning were compared with targets whose warnings are in other classes, on the target's expression z-score in the matching organ.

| Toxicity class | Targets with it | z, with | z, others | Difference |
|---|---|---|---|---|
| psychiatric | 18 | 1.87 | 0.72 | 1.15 |
| neurotoxicity | 37 | 1.18 | 0.75 | 0.43 |
| dermatological | 17 | 0.58 | 0.19 | 0.39 |
| vascular | 34 | 0.76 | 0.37 | 0.39 |
| nephrotoxicity | 14 | 0.43 | 0.19 | 0.24 |
| cardiotoxicity | 64 | 0.52 | 0.31 | 0.21 |
| respiratory | 26 | 0.62 | 0.60 | 0.02 |
| hepatotoxicity | 74 | -0.27 | -0.22 | -0.05 |
| gastrointestinal | 22 | 1.01 | 1.10 | -0.09 |
| hematological | 29 | -0.95 | -0.53 | -0.43 |
| musculoskeletal | 7 | -1.10 | -0.59 | -0.51 |

Weighted mean difference 0.153 SD, one-sided p = 0.033. Met. Liver is the largest class and shows nothing, consistent with drug liver injury coming mostly from reactive metabolites rather than from the target. GTEx has whole blood but no bone marrow, so the hematological comparison measures the wrong tissue.

The API names tissues differently from GTEx (skeletal muscle is "gastrocnemius medialis" under "muscle tissue", and "blood vessel" contains "blood"), so the matching rule was fixed as: a keyword equal to the parent organ or inside the tissue name. This was set after seeing the tissue list for one gene and before the organ comparison was run.

## P8: LLM safety rationales

Claude Opus 5.5 (effort medium, structured JSON output, server-side refusal fallback enabled) wrote one assessment per target for the 12 targets with the highest M3 risk and 12 drawn at random. 24 reports, 222 sentences: 131 fact, 50 hypothesis, 41 gap. Every sentence cited at least one ID and every cited ID existed. All 264 numbers in fact sentences appeared in a cited row at the same rounding, and all 50 toxicity classes named in fact sentences appeared in a cited row. Met. Cost: $0.66.

The model gave "high" to 8 targets, all with warned drugs, "medium" to 15 (11 warned) and "low" to one (not warned). This is not a prediction: the pack lists each drug's warnings.

## Post hoc

- FAERS beyond drug count: adding the FAERS score to M0 changed AUROC by 0.006 (-0.008 to 0.018). The number of significant FAERS terms mostly measures how many and how widely used the drugs are.
- FAERS source data: 55 adverse-reaction profiles are attached to more than one parent molecule (125 of 3,028), 21 of them cohort drugs. Many are synonyms or prodrugs (baloxavir and baloxavir marboxil), but some are mismatches: fluvastatin and pravastatin have identical rows, as do bentiromide and pirfenidone. Name-based mapping of FAERS reports needs checking before its counts are used for a single drug.
- Numbers outside fact sentences: the same check on hypothesis and gap sentences found 54 numbers, all supported.
- A full read of three reports with their packs (SOST, NGF, HRH1) found no misstatement. One SOST sentence restating FAERS counts was labeled hypothesis, which the audit would not have checked had it been wrong; the post hoc check above covers that gap.
- The number check first accepted a match after dropping a minus sign. A unit test caught it, the check was made strict, and the audit was rerun on the stored reports with the same result.

## Verdicts

| Prediction | Result | Met |
|---|---|---|
| P1 warnings shared within targets | 72.8% vs 60.0%, p = 0.0001 | yes |
| P2 FAERS AUROC at least 0.65 | 0.673 | yes |
| P3 drug count AUROC at least 0.60 | 0.699 | yes |
| P4 pre-drug evidence adds to drug count | +0.015 (-0.001 to 0.031) | no |
| P5 learned weights not better by more than 0.03 | -0.018 | yes |
| P6 clinical minus experimental curated events above 0.03 | +0.043 (-0.008 to 0.093) | yes |
| P7 organ score higher where warned | 0.153 SD, p = 0.033 | yes |
| P8 LLM rationales traceable | 100% of 222 sentences, 100% of 264 numbers | yes |

## What this says for target safety scoring

- The question is real: warnings follow targets more than chance, so scoring at the target level is not misplaced.
- Generic pre-drug factors are weak on their own. A score built from them should be read as a prompt for specific questions, not as a ranking, and any score has to be compared with the trivial baseline of how much drug exposure a target has had.
- Organ-level reasoning (where is the target expressed, and does that organ match the liability) carried more signal than any single generic factor, which argues for scoring liabilities organ by organ rather than as one number.
- Curated events and FAERS are useful but partly restate what is already known about existing drugs; for a new target with no drugs they will be empty.
- An LLM can turn an evidence pack into a readable rationale whose every statement is checkable by code, which is the property that makes such reports reviewable.

## Limits

A boxed warning is a regulatory decision, not a measure of harm, and oncology drugs carry warnings for toxicities accepted for benefit. Many harms come from off-target activity or reactive metabolites that no target-level factor sees. The pre-drug factors are Open Targets' own scaled values and were used as given. Expression comes from healthy adult tissue, without bone marrow. The LLM sample is small, its audit checks traceability rather than judgment, and the packs include the outcome.
