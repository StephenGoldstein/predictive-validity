# Target origin: full run over the 2008–2019 window (#146 step 5), 2026-11-20 to 2026-11-21

The full run started after both pilot gates passed under the aligned definitions (`target_origin_pilot.md`, last
part) and after its announcement on #146 (2026-11-20). It ran from 2026-11-20 15:04 to 2026-11-21 09:32 Pacific
time, unattended, as three invocations of one wrapper outside the repository: search until no pair had failed, then
`label`, `check` and `summary`, each over `--pairs window`, each re-invoked after a pause when it exited non-zero.
Every step is idempotent (search skips searched pairs, label reuses validated decisions, check skips checked rows),
so the second and third invocations only asked again for the pairs still unlabelled. Numbers below are from the run
log, `python -m pipeline.target_disease summary --pairs window` and `pv_evidence.llm_decision`.

**Result: 5,101 of 5,174 pairs carry a labelled origin; 73 (1.4%) are unlabelled because the API refused them on
every pass.** After the check, 1,582 of the 5,101 are unclear (31%).

## Scope

The window is every target × disease pair whose first Phase 2 or 3 entry falls in 2008–2019: 5,174 pairs
(`pair_frame("window")`). The 100 pilot pairs are inside it; their label decisions (100) and check decisions (87)
were reused unchanged, so the pilot's gate numbers stand for the same rows.

## Retrieval (`search --pairs window`)

| pass | pairs asked | searched | records stored | failed |
|---|---|---|---|---|
| 1 (15:04 to 03:54) | 5,074 | 5,063 | 235,969 | 11 |
| 2 (03:54 to 03:55) | 11 | 11 | 73 | 0 |

236,042 records in all. The 11 failures were Europe PMC 500 and 503 responses in one half hour around 17:00; a failed
pair is rolled back and left unmarked, so the second pass found exactly those 11. Single requests that failed and
were retried by the client (up to five times) appear 99 times in the log and cost nothing. 3,016 pairs (58%) have at
least one record; the median over all pairs is 1 and the largest bundle is 31,193.

## Labelling (`label --pairs window`, Opus 5.5, prompt v2, web search and fetch)

Three invocations of `label`, three passes each (`LABEL_PASSES`); every pass asks again only for the pairs whose last
answer failed validation.

| pass | asked | validated | still undecided |
|---|---|---|---|
| 1 | 5,074 | 4,791 | 283 |
| 2 | 283 | 135 | 148 |
| 3 | 148 | 25 | 123 |
| 4 | 123 | 18 | 105 |
| 5 | 105 | 13 | 92 |
| 6 | 92 | 7 | 85 |
| 7 | 85 | 6 | 79 |
| 8 | 79 | 2 | 77 |
| 9 | 77 | 4 | 73 |

6,066 model calls: 5,001 validated, 950 refused by the API, 115 broke a rule (a cited record the model was not
shown, a year that is not a cited record's, a route quote without a route). Tokens: 34.7 million in (uncached),
331.9 million read from cache, 10.5 million out; 16,087 web searches and 4,789 fetches.

**Refusals.** A refusal is the API's `stop_reason: refusal`, the streaming safety classifier, not a rule break. On
pass 1 it hit 187 of 5,074 calls (3.7%); 176 of the 187 came after the model had begun its answer (one cut off
mid-sentence while reading CD19 expression on AML blasts). The refused inputs are bigger (median 4,251 input tokens
against 959 for validated calls) and carry more pathogen, toxin, vaccine, gene-therapy, radiation and infection
language in their abstracts (each two to four times as often as in validated inputs). They are not one disease or one
target: the 118 pairs still refused after pass 3 spanned more than a hundred targets, the largest groups CSF2 (7
pairs) and FKBP1A (5), and the diseases ran from acute myeloid leukemia (7) to hemophilia, ovarian cancer and
malaria. A refusal is partly random: of the 187 refused on pass 1, 48 passed on pass 2 and 22 of the rest on pass 3;
later passes recovered 18, 13, 7, 6, 2 and 4. The run was stopped at 73 because the last three passes had converged.

**The 73 unlabelled pairs** (`target_origin_unlabelled.csv`: pair, symbol, disease, records stored, reason) have no
row in `preclin.target_disease_origin`. They are missing data, not `unclear`: an unclear origin is a finding about the
record, an unlabelled pair is a pair the method did not reach. Their record counts run from 0 to 558 (median 11).
The in-method remedy is more passes, now at a few pairs a pass; the out-of-method remedy is another model for these
pairs, recorded as such. Neither was taken; both are the owner's call.

## Check (`check --pairs window`, Sonnet 5.5, prompt v3, fetch)

4,419 graded labels (87 reused from the pilot check); 4,335 calls, 4,332 validated, 3 asked again. The first
invocation died on an NCBI efetch outage before any decision (`RuntimeError: eutils efetch failed`) and the wrapper
re-invoked it after two minutes. Tokens: 4.2 million in, 65.5 million read from cache, 1.5 million out; 2,115 fetches.

| outcome | labels |
|---|---|
| supported | 3,519 |
| unsupported | 730 |
| unverified | 170 |

A label that is unsupported or unverified is stored as `unclear` with grade `none`.

## Stored result (`preclin.target_disease_origin`, basis `labelled`, 5,101 rows)

| origin | label as read | stored after the check |
|---|---|---|
| animal_model_pharmacology | 944 | 797 |
| diseased_tissue | 925 | 850 |
| earlier_drug_other_disease | 727 | 361 |
| unclear | 682 | 1,582 |
| human_cell_perturbation | 482 | 425 |
| clinical_observation | 392 | 302 |
| observational_biomarker | 280 | 245 |
| mouse_genetics | 209 | 172 |
| human_genetics | 198 | 170 |
| tumour_genetics | 132 | 122 |
| phenotypic_screen | 113 | 66 |
| curated_pathway | 12 | 6 |
| cell_screen | 5 | 3 |

The 1,582 unclear rows are 682 the labeller left ungraded, 730 the checker rejected and 170 it could not verify.
Grades stored: stated 1,730, pattern 1,789, none 1,582. The check falls hardest on earlier_drug_other_disease (366 of
727 labels as read lost, 50%) and phenotypic_screen (47 of 113, 42%); diseased_tissue loses 75 of 925 (8%).

Discovery route of the first drug: derived_from_active 1,380, repurposed_drug 1,158, target_screen 669,
de_novo_design 228, phenotypic_screen 86, unclear 1,580.

The `summary` command's "tokens per decision" table counts reused decisions again on every run (they are copied as
rows with `prior_decision_id`), so its `n` overstates the calls; the call counts above are rows without a prior.

## Timeline

| step | from | to |
|---|---|---|
| search, pass 1 | 11-20 15:04 | 11-21 03:54 |
| search, pass 2 | 03:54 | 03:55 |
| label, passes 1 to 3 | 03:55 | 05:57 |
| check (first attempt died at 05:57; completed run) | 06:08 | 06:48 |
| summary | 06:48 | |
| label, passes 4 to 6; check; summary | 07:34 | 07:57 |
| label, passes 7 to 9; check; summary | 09:13 | 09:32 |

## What this leaves for step 6

Step 6 (the O1 test, descriptives and the figure) is the owner's. Two things it should know: the coverage split by
group over the whole window is not computed here (the coverage gate scored it on the pilot's settled-stop frame), and
the 73 unlabelled pairs are absent from the origin table rather than unclear.

## Checker sensitivity: Opus 5.5 as a second checker (2026-11-21, 09:50 to 10:51)

Asked on the day the run finished: the checker (Sonnet 5.5) is a smaller model than the labeller (Opus 5.5), so
would an Opus checker verify differently? Measured rather than argued: every one of the 4,419 graded labels was
checked again by Opus 5.5 under the v3 prompt text, recorded as `target_origin_check` prompt version v4 (runs
ddc03aac, ac2d3948, 8e9f70c6, 2a997db7, 878feeb8), in the decision log only; `preclin.target_disease_origin` was not
touched. 4,402 labels were decided by both (Opus refused 17). `target_origin_check_opus.csv` carries both verdicts
per label; `target_origin_check_disagreements.md` shows twenty disagreements with the quote, the definition and both
checkers' notes, for a human read.

| | Opus supported | Opus unsupported | Opus unverified | all |
|---|---|---|---|---|
| Sonnet supported | 3,439 | 47 | 24 | 3,510 |
| Sonnet unsupported | 590 | 125 | 9 | 724 |
| Sonnet unverified | 150 | 12 | 6 | 168 |

κ between the checkers is 0.20 (0.21 on supported against the rest). Opus rejects 184 labels (4.2%), Sonnet 724
(16.4%); they agree on 125 rejections. Under Opus the stored table would hold 223 unclear graded rows instead of 892.

**Where the disagreement sits.** Of the 590 labels Sonnet rejected and Opus supported, Opus scored the category
"partly" on 528 and "yes" on 61; 455 of the 590 are labels Sonnet rejected on the category alone. Read side by side,
the two checkers describe the same facts (the quoted pre-cutoff source shows the drug in use in the same disease and
does not show an earlier approval for an unrelated disease) and set the threshold differently: Sonnet says the
category does not follow, Opus says "partly" and supplies the missing step from its own knowledge ("already approved
for renal cancer", "an existing approved drug repurposed"). The prompt defines "partly" as support only together with
the other sources, or support for a neighbouring category as well; a link the record does not contain is not in that
definition, so on the rule as written Sonnet applied it and Opus was lenient. The leniency is largest where the
labeller over-calls most: earlier_drug_other_disease keeps 663 of 726 labels under Opus against 361 under Sonnet,
phenotypic_screen 105 of 111 against 65. An Opus checker shares the Opus labeller's reading; that is the case for a
second model, in numbers.

**Where Sonnet is weak: the year field.** Opus supported 71 of the 87 labels Sonnet rejected on the year alone, and
38 of Sonnet's 116 year-based rejections carry a note that says the year is fine while the field says no. Those are
checker errors, about 2% of the stored labels. A later prompt version can take the year out of the model's hands:
ask for the year it read for each web source and compute `year_ok` in code.

**The other side.** Opus rejects 47 labels Sonnet supported: 29 on the category, 10 on the year, 7 on the quote
(13 animal_model_pharmacology, 7 diseased_tissue, 7 earlier_drug_other_disease). Opus marks "partly" on 61% of the
labels it supports, Sonnet on 53%.

Stored origin under each checker, over the 4,402 labels both decided (the 682 the labeller left ungraded are unclear
under both):

| origin | label as read | Sonnet (stored) | Opus |
|---|---|---|---|
| animal_model_pharmacology | 935 | 793 | 884 |
| diseased_tissue | 925 | 850 | 901 |
| earlier_drug_other_disease | 726 | 361 | 663 |
| human_cell_perturbation | 481 | 424 | 464 |
| clinical_observation | 392 | 302 | 370 |
| observational_biomarker | 279 | 244 | 267 |
| mouse_genetics | 206 | 170 | 189 |
| human_genetics | 198 | 170 | 196 |
| tumour_genetics | 132 | 122 | 126 |
| phenotypic_screen | 111 | 65 | 105 |
| curated_pathway | 12 | 6 | 10 |
| cell_screen | 5 | 3 | 4 |
| unclear (graded rows) | 0 | 892 | 223 |

**Status.** The stored table keeps the Sonnet verdicts, which is the registered design. The threshold moves the
owner's variable by up to 669 rows, so the choice is the owner's, with the twenty-case sheet as the material. If the
Opus verdicts are chosen, `label --pairs window` then `check --pairs window` with `CHECK_TASK` at prompt version v4
and model `claude-opus-5-5` rewrite the rows from the stored v4 decisions without new check calls.
