# Design

Written after reading the table schemas and the category labels of `warningType`, `toxicityClass`, `maximumClinicalStage`, `drugType`, `actionType`, `targetType` and the safety-event `datasource` in Open Targets Platform 26.09, and the value ranges of the target prioritisation factors over all 78,733 genes. No outcome was linked to any feature before this plan was committed.

## Question

Target safety teams combine evidence about a drug target (genetics, knockout mice, tissue expression, curated safety events, adverse-event reports) into a risk score. With public data only:

1. Are serious drug safety problems shared by drugs that act on the same target, so that a target-level score can work at all?
2. Does evidence available before any drug exists rank targets with boxed warnings or withdrawals above other targets, beyond the number of drugs on each target?
3. Do learned weights beat equal weights?
4. How much of the signal in curated safety events is circular, because the curation itself drew on drugs already known to cause harm?
5. Do warnings fall in organs where the target is expressed?
6. Can an LLM write a target safety rationale in which every statement is traceable to the evidence?

## Cohort and outcome

- Approved drugs: `drug_molecule` rows with `maximumClinicalStage` = APPROVAL, salts mapped to their parent molecule. Drug types kept: small molecule, antibody, protein, oligonucleotide, oligosaccharide. Antibody-drug conjugates, cells, genes, enzymes, vaccine components and unknown types are dropped, since their harms often come from a payload or a replacement product rather than the target.
- Qualifying mechanisms: `targetType` = single protein, excluding action types that are not pharmacology on a target (exogenous gene or protein, vaccine antigen, hydrolytic or proteolytic enzyme, cross-linking or disrupting agent).
- Targets: human protein-coding genes named in a qualifying mechanism of at least one approved drug.
- Outcome Y: 1 if any approved drug with a qualifying mechanism on the target has a ChEMBL warning of type Black Box Warning or Withdrawn (any country, any year).
- Drug count n: the number of such drugs per target. More drugs give more chances of a warning, so n is the baseline every score has to beat.
- Sensitivity: the same analysis restricted to drugs with exactly one qualifying target, where attribution to the target is clean.

## Evidence

Pre-drug evidence, from `target_prioritisation`, where negative values are unfavourable by the Open Targets convention; risk is the negated value:

- genetic constraint (gnomAD loss-of-function)
- mouse knockout score
- cancer driver gene (missing means not a driver, risk 0)
- maximum paralog identity
- tissue specificity
- tissue distribution

Curated evidence: target safety events, split by source.

- Clinical-derived (reviews and panels built from harms seen with drugs in people): Bowes et al. (2012), Brennan et al. (2024), Urban et al. (2012), Lynch et al. (2017), Force et al. (2011), ClinPGx.
- Experimental: ToxCast, AOP-Wiki, Lamore et al. (2017).

Post-marketing evidence, the positive control: for each target, the largest number of significant FAERS adverse-reaction terms (openFDA, likelihood ratio test, from Open Targets) among its approved drugs, on the log1p scale.

## Scores and models

- Equal-weight pre-drug score: mean percentile rank (among cohort targets) of the six risk-oriented pre-drug features; a missing continuous feature gets the middle rank.
- M0: log n alone.
- M1: logistic regression on log n and the equal-weight score.
- M2: logistic regression on log n, the six features and missing-value indicators (learned weights).
- M3: M1 plus the two curated-source indicators.
- FAERS score alone.

Fitted models are scored on out-of-fold predictions from 5-fold cross-validation grouped by ChEMBL target class (level 1), so related targets are not split across training and test folds. AUROC with 1,000 bootstrap resamples of targets; differences use paired resamples.

## Mechanism check (organ match)

Toxicity classes with an organ in GTEx: hepatotoxicity (liver), cardiotoxicity (heart), nephrotoxicity (kidney), neurotoxicity and psychiatric toxicity (brain), dermatological (skin), gastrointestinal (colon, stomach, esophagus, small intestine), respiratory (lung), hematological (blood), musculoskeletal (skeletal muscle), vascular (artery). For each target with at least one warning in these classes, GTEx median expression (log1p TPM) is fetched from the Open Targets API, and the organ score is the z-score of the highest matching tissue against all GTEx tissues of that target. Within each class, targets with a warning of that class are compared with targets whose warnings are only in other mapped classes.

## LLM safety rationales

For 24 targets (the 12 with the highest M3 out-of-fold risk and 12 drawn at random from the rest, seed 7), an evidence pack lists every fact with an ID (constraint, knockout phenotypes, expression factors, curated events, Reactome top-level pathways, approved drugs with their warnings and top FAERS terms). Claude writes a short assessment as structured output: each sentence carries the IDs it rests on and is marked as fact or hypothesis. Audits are automatic:

- every sentence cites at least one ID;
- every cited ID exists in that target's pack;
- every number in a fact sentence appears in a cited evidence row (same rounding);
- every toxicity class named in a fact sentence appears in a cited row.

## Predictions

| | Prediction | Met if |
|---|---|---|
| P1 | Warnings are shared within targets | Among single-target drugs on targets with at least two such drugs, the share of same-target drug pairs with the same warning status exceeds a permutation null (10,000 shuffles of drug labels), one-sided p < 0.01 |
| P2 | FAERS works as a positive control | AUROC of the FAERS score at least 0.65 |
| P3 | Drug count is a real baseline | AUROC of log n at least 0.60 |
| P4 | Pre-drug evidence adds to drug count | AUROC of M1 minus M0 above 0, 95% interval excluding 0 |
| P5 | Learned weights do not beat equal weights | AUROC of M2 minus M1 at most 0.03 |
| P6 | Curated events are partly circular | AUROC of the clinical-derived indicator minus the experimental indicator above 0.03 |
| P7 | Warnings fall where the target is expressed | Class-stratified mean difference in organ score above 0, one-sided permutation p < 0.05 |
| P8 | LLM rationales are traceable | At least 95% of sentences pass the citation checks and at least 95% of numbers in fact sentences are supported |

P1, P3 and P4 are also reported for the single-target sensitivity cohort; verdicts use the main cohort.

## Limits known in advance

A boxed warning is a regulatory decision, not a measure of harm, and many harms come from off-target activity or reactive metabolites that no target-level feature can see. Oncology drugs carry warnings for toxicities accepted in exchange for benefit. The pre-drug features are Open Targets' own scaled factors, so their orientation is taken as given. Expression comes from healthy adult tissue.
