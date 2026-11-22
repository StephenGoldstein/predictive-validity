# Target origin: full run over the 2008–2019 window (#146 step 5), 2026-11-20 to 2026-11-21

## State (2026-11-21 16:08 Pacific, after check method v2)

5,174 pairs. 241,118 pre-entry records; 3,016 pairs (58%) have at least one. 5,113 pairs labelled; 61 unlabelled, the
API's refusal on every pass (`target_origin_unlabelled.csv`). 4,431 graded labels checked: supported 3,654, unsupported
392 (301 of them earlier-drug labels with an approval for a related disease, 18 of those the hub's
`validated_at_entry`), disputed 380 (`target_origin_disputes.csv`, 426 rows awaiting a verdict),
unverified 5. Stored: 1,459 of the 5,113 rows are unclear (29%); the table by origin, as read and stored, is in the
last section. Gates, scored on the label as read: set b κ 0.65, set a 0.69, coverage gap 8 points, both pass.

The sections below are the dated record in the order things happened: the run as it ran (search, label, and the v1
check by Sonnet 5.5 alone), the checker sensitivity study of the next morning, and check method v2, which replaced the
v1 outcomes in the stored table that afternoon. The v1 check and the stored result under it are kept and marked
superseded. Numbers in the run sections are from the run log, `python -m pipeline.target_disease summary --pairs
window` and `pv_evidence.llm_decision`; the state above is from the stored table at the last commit.

## How it ran

The full run started after both pilot gates passed under the aligned definitions (`target_origin_pilot.md`, last
part) and after its announcement on #146 (2026-11-20). It ran from 2026-11-20 15:04 to 2026-11-21 09:32 Pacific
time, unattended, as three invocations of one wrapper outside the repository: search until no pair had failed, then
`label`, `check` and `summary`, each over `--pairs window`, each re-invoked after a pause when it exited non-zero.
Every step is idempotent (search skips searched pairs, label reuses validated decisions, check skips checked rows),
so the second and third invocations only asked again for the pairs still unlabelled. The wrapper's output is
`target_origin_run_log.txt` beside this file (personal paths removed) and, for the afternoon's check v2 passes,
`target_origin_run_v2_log.txt`; the pass tables, refusal counts and token counts below are read from them.

## Scope

The window is every target × disease pair whose first Phase 2 or 3 entry falls in 2008–2019: 5,174 pairs
(`pair_frame("window")`). The 100 pilot pairs are inside it; their label decisions (100) and check decisions (87)
were reused unchanged, so the pilot's gate numbers stand for the same rows.

## Retrieval (`search --pairs window`)

| pass | pairs asked | searched | records stored | failed |
|---|---|---|---|---|
| 1 (15:04 to 03:54) | 5,074 | 5,063 | 235,969 | 11 |
| 2 (03:54 to 03:55) | 11 | 11 | 73 | 0 |

236,042 records stored by this run; with the pilot's 100 pairs, searched on 2026-11-20 (5,076 records), the window
holds 241,118. The 11 failures were Europe PMC 500 and 503 responses in one half hour around 17:00; a failed
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
pairs) and FKBP1A (5), and the diseases ran from acute myeloid leukemia (7) to hemophilia, ovarian cancer and malaria.
A refusal is partly random: of the 187 refused on pass 1, 48 passed on pass 2 and 22 of the rest on pass 3; later
passes recovered 18, 13, 7, 6, 2 and 4. The run was stopped at 73 after nine passes; three more passes on 2026-11-21
(the check v2 section) recovered 12 at the same rate, so 61 remain.

**The unlabelled pairs** (73 at this point, 61 after the passes of 2026-11-21; `target_origin_unlabelled.csv` lists
the 61: pair, symbol, disease, records stored, reason) have no row in `preclin.target_disease_origin`. They are
missing data, not `unclear`: an unclear origin is a finding about the record, an unlabelled pair is a pair the method
did not reach. Their record counts run from 0 to 558 (median 11). The in-method remedy is more passes, now at a few
pairs a pass (taken once more on 2026-11-21); the out-of-method remedy is another model for these pairs, recorded as
such (not taken; the owner's call).

## Check v1 (`check --pairs window`, Sonnet 5.5, prompt v3, search and fetch): superseded by check method v2 below

4,419 graded labels (87 reused from the pilot check); 4,335 calls, 4,332 validated, 3 asked again. The first
invocation died on an NCBI efetch outage before any decision (`RuntimeError: eutils efetch failed`) and the wrapper
re-invoked it after two minutes. Tokens: 4.2 million in, 65.5 million read from cache, 1.5 million out; 2,115 fetches.

| outcome | labels |
|---|---|
| supported | 3,519 |
| unsupported | 730 |
| unverified | 170 |

A label that is unsupported or unverified is stored as `unclear` with grade `none`.

## Stored result after check v1 (5,101 rows; superseded, the v2 table is in the last section)

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
"partly" on 529 and "yes" on 61; 455 of the 590 are labels Sonnet rejected on the category alone. Read side by side,
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

**The other side.** Opus rejects 47 labels Sonnet supported: 29 on the category, 10 on the year, 7 on the quote, 1 on
both quote and category (13 animal_model_pharmacology, 7 diseased_tissue, 7 earlier_drug_other_disease). Opus marks
"partly" on 61% of the labels it supports, Sonnet on 53%.

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

**Status.** Written at 11:00 with the stored table on the Sonnet verdicts. Superseded the same afternoon by check
method v2 (next section): the disagreement was scored against the dated approval index, which showed each checker
wrong in its own direction on the earlier-drug category, and the stored table now follows neither checker alone.

## Check method v2: facts in code, judgment by agreement (2026-11-21, 14:11 to 16:08)

Posted on #146 at 14:08, before any row was rewritten; migration `db/148_target_disease_origin_check_v2.sql` adds
the outcome `disputed` and the column `check_fields`.

**Why.** The earlier-drug category's defining fact, a drug on the target approved before entry, is a dated fact in the
database: `report.section3.Evidence` indexes Drugs@FDA approvals by target gene for the chart's precedent lines.
Scored against it, the 727 earlier-drug labels as read split 408 (approved for an unrelated disease, which is the
definition), 301 (approved for a related disease, Lin at least 0.5, the chart's rule) and 18 (no approval on the
target before entry in the index). The v1 checker kept 219 of the 408 confirmed labels, rejected 168 (154 on the
category alone, for the cited paper not restating the approval) and left 21 unverified; it kept 134 of the 301 related
ones and 8 of the 18 unconfirmed (219 + 134 + 8 is its stored 361). Opus rejected 27 of the 408 and 25 of the 301 and
passed all 18 it cannot confirm, among them roxadustat for myelodysplastic syndrome, approved nowhere before entry,
which its own note says. The v1 check asked a literature record for a database fact, and each model failed it in its
own direction. Separately, Sonnet's year field contradicts its own note on 38 of its 116 year-based rejections.
Neither defect is model strength; both are the design asking a model for what the data or the code already hold.

**The rule** (`check_fields`, `check_outcome_v2`; `CHECK_DESCRIPTION` has it in full). Three fields per graded label.
quote: the two checkers' answers where they agree; could_not_open carries no information. year: no when a cited record
is dated in or after the cutoff year; yes when every source is a record; otherwise the two answers where they agree,
unknown carrying no information. category: for earlier_drug_other_disease the approval fact (unrelated yes; related
no, since the chart's approved_unrelated_only tier, which the owner's one-list decision makes the definition, excludes
a pair with an approval for a related disease at Lin 0.5 or more; none disputed), for every other origin the two
answers where they agree, partly counting as yes. A field the two disagree on, or that one of them did not answer, is
disputed until a human verdict from `target_origin_disputes.csv` settles it: `check` reads the sheet (`--disputes
PATH`, default that file) and rewrites it with the open rows. A verdict is `yes` (the field holds) or `no`; any other
value stops the run; a verdict applies only where the checkers left the field disputed. Outcome: unsupported if a
field is no, else disputed if one is disputed, else unverified if one is unverified, else supported. A row that is not
supported holds unclear with grade none; `check_fields` records the fields, the fact, both decision ids and the
verdicts applied. The sheet lists only labels whose outcome a verdict can change.

**The run.** `label --pairs window` first, to restore every row to the label as read: 5,101 decisions reused, the 73
refused pairs asked three more times (201 calls, 12 validated, 189 refused), so 5,113 pairs are labelled and 61 are
not (`target_origin_unlabelled.csv` is now the 61). Then `check --pairs window`: Sonnet v3 reused 4,419 decisions
and decided the 12 new labels (12 calls); Opus v4 reused 4,402 and was asked again for the labels it had refused and
the 12 new ones (78 calls over the day's three check passes, 24 validated, 54 refused; 5 labels still have one
checker only and are disputed on that account). No other model calls; the fields are computed from the stored
decisions. Three passes because the second exposed a gap: `check` read only rows still graded, so after the first
pass it saw the 3,649 supported rows and wrote the disputes sheet empty (the table was untouched). It now reads every
label graded as read from its decision and restores a supported row, so it reruns on its own to apply verdicts; the
third pass (15:49 to 16:08) reused all 4,431 stored decisions and wrote the sheet in full.

| outcome | labels (of 4,431 graded) |
|---|---|
| supported | 3,654 |
| unsupported | 392: the 301 validated at entry, and 91 on a field both checkers, or the code, answered no |
| disputed | 380 |
| unverified | 5 |

Disputes sheet: 426 rows over the 380 disputed labels (category 277, year 86, quote 63;
18 of the category rows are the approval fact). On the category rows Sonnet says no and Opus partly on
228: the same threshold the previous section describes, now on a person's desk rather than either model's.

Outcome by origin as read (the 4,431 graded labels; the share is disputed over as read):

| origin as read | labels | supported | unsupported | disputed | unverified | disputed share |
|---|---|---|---|---|---|---|
| animal_model_pharmacology | 949 | 833 | 30 | 85 | 1 | 9% |
| diseased_tissue | 927 | 870 | 9 | 46 | 2 | 5% |
| earlier_drug_other_disease | 727 | 388 | 303 | 35 | 1 | 5% |
| human_cell_perturbation | 486 | 439 | 12 | 35 | 0 | 7% |
| clinical_observation | 393 | 322 | 14 | 56 | 1 | 14% |
| observational_biomarker | 280 | 252 | 6 | 22 | 0 | 8% |
| mouse_genetics | 209 | 176 | 7 | 26 | 0 | 12% |
| human_genetics | 198 | 175 | 1 | 22 | 0 | 11% |
| tumour_genetics | 132 | 120 | 4 | 8 | 0 | 6% |
| phenotypic_screen | 113 | 69 | 3 | 41 | 0 | 36% |
| curated_pathway | 12 | 6 | 2 | 4 | 0 | 33% |
| cell_screen | 5 | 4 | 1 | 0 | 0 | 0% |

The dispute rate is 9% overall and falls unevenly: 36% of phenotypic_screen labels as read and 14% of
clinical_observation, against 5 to 12% in the other categories with more than a dozen labels. Every dispute holds
unclear, so the stored table under-counts those categories until the sheet is read; a verdict can only move a row out
of unclear, and 345 of the 380 disputed labels have one disputed field with the other two yes. The 91 unsupported
labels outside the 301 failed on the category (64), the year (24) or the quote (4); one on two fields. The 5
unverified are labels whose cited source neither checker could open.

Stored origin after v2 (5,113 labelled rows), beside v1:

| origin | v1 (Sonnet alone) | v2 |
|---|---|---|
| unclear | 1,582 | 1,459 (682 ungraded, 392 unsupported, 380 disputed, 5 unverified) |
| diseased_tissue | 850 | 870 |
| animal_model_pharmacology | 797 | 833 |
| human_cell_perturbation | 425 | 439 |
| earlier_drug_other_disease | 361 | 388 |
| clinical_observation | 302 | 322 |
| observational_biomarker | 245 | 252 |
| human_genetics | 170 | 175 |
| mouse_genetics | 172 | 176 |
| tumour_genetics | 122 | 120 |
| phenotypic_screen | 66 | 69 |
| curated_pathway | 6 | 6 |
| cell_screen | 3 | 4 |

The 301 rows with `check_fields.fact = related` hold unclear, and "related" has two thresholds in this repository.
The fact uses the chart's rule, Lin 0.5 or more (`report.section3.RELATED`, the b1 line), because that is the
threshold inside the approved_unrelated_only tier the one-list decision defines the category by. The hub's
`validated_at_entry` flag (migration 142) is stricter, Lin 0.8 or more, and marks 83 pairs in the window. The owner's
one-list comment and the header comment of migration 145 name the hub flag for the related-disease case; the stored
fact follows the chart's tier. So 18 of the 301 are validated at entry by the hub's definition and 283 carry a
related-disease approval below it. All 83 hub-flagged pairs are in the stored table (18 among the 301, 63 with
another origin as read, 2 ungraded); by the hub's own column comment their origin is not tested, and step 6 can set
them aside by the flag. Whether the 283 should count as earlier_drug_other_disease instead is put to the owner on #146
with the closing note; either answer is a rewrite from the stored decisions.

Discovery route is unchanged but for the 12 new labels (unclear 1,583, derived_from_active 1,385, repurposed_drug
1,159, target_screen 670, de_novo_design 230, phenotypic_screen 86).

**Gates, rerun on the v2 table.** Set a: kappa 0.69 on the label as read (stored 0.52), coverage gap 8.0 points,
passed. Set b: 0.65 (stored 0.55), gap 8.0, passed. The gate is scored on the label as read, which v2 does not
touch; the stored kappas are reported beside it as before.

**What this leaves.** The disputes read, the owner's by default (the #146 note offers it back); every verdict moves
one row out of unclear, and unread rows stay unclear with the count reported. Also for step 6: the coverage split by
group over the whole window is not computed (the gate scored it on the pilot's settled-stop frame), and the 61
unlabelled pairs are absent from the table. The residual the check cannot remove: two capable models disagree on at
least one field for 357 of 4,431 graded labels (on the category for 254); the other disputed labels are the 18 the
approval index is silent on and the 5 one checker refused. The audit reports that rate rather than present the check
as exact.
