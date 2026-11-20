# Target origin: pilot of 100 pairs (#146 step 4), 2026-11-19

The pilot ran once, on prompt v1, with the gate reference the owner decided on 2026-11-19: three independent blind
model readings of hand-check set a, majority taken. **The agreement gate failed and the coverage gate passed.** Nothing
beyond the pilot ran. Numbers below are from `report/audits/target_origin_gate.json`, `python -m pipeline.target_disease
summary --pairs pilot` and the decision log.

## Sample

Seed 20261119: 50 cohort pairs with a settled "not causal" stop and 50 without (pairs with any stop filed as
unresolved left out of both groups, 127 of 5,091, per #158; frame 703 and 4,261). Set a: 25 per group.

## Retrieval (`search --pairs pilot`)

| | not causal | other | all |
|---|---|---|---|
| pairs | 50 | 50 | 100 |
| pairs with at least one pre-entry record | 20 (0.40) | 24 (0.48) | 44 |
| records stored | | | 4,583 (median per pair 0, max 1,446) |

Coverage gap 8 points: **passes** (at most 10). Fifty-six pairs have no record that names both the gene and the
disease before the cutoff year. The query is correct; the names are too few. The hub's HGNC table carries symbols only,
so a receptor of a well-known ligand is searched under its HGNC name alone (PTGIR: "prostaglandin I2 receptor" and the
two-letter alias "IP", while papers say "prostacyclin receptor" or "IP receptor"; the gene alone finds 38 pre-2017
abstracts, the disease alone 1,141, both none).

## Labelling (`label --pairs pilot`, Opus 5.5, web search and fetch)

99 of 100 pairs decided over three passes (94, 3, 2). Rejected answers: 5 refusals (two pairs, one never decided) and
5 rule breaks (a PMID cited from memory rather than retrieved, 3; a year count mismatch, 1; an origin year not among
the cited sources, 1). Tokens per decision: 4,318 in, 2,144 out.

| origin (labeller, before the check) | n |
|---|---|
| animal_model_pharmacology | 11 |
| earlier_drug_other_disease | 7 |
| observational_biomarker | 7 |
| diseased_tissue | 7 |
| human_cell_perturbation | 4 |
| clinical_observation | 4 |
| human_genetics | 3 |
| tumour_genetics | 2 |
| mouse_genetics | 2 |
| curated_pathway | 1 |
| unclear | 10 |

Grades: stated 22, pattern 26, none 10 (the 41 rows the check rewrote are counted below). Discovery route: derived
from an existing active 25, repurposed drug 25, target screen 12, de novo design 5, unclear 32.

## Check (`check --pairs pilot`, Sonnet 5.5, fetch only)

89 graded labels checked: supported 48, unverified 37, unsupported 4. Of the 37 unverified, 33 are "could not open":
the checker's fetch failed on PubMed ("too many requests"), ResearchGate, Semantic Scholar and blocked hosts, and 34 of
the 37 labels cite URLs only, so no shown record could carry the quote. Unverified is a fetch failure, not a reading.
Under the design an unverified or unsupported row is rewritten as unclear, so after the check 51 of 99 labelled rows
are unclear.

## Reference reading (`blind --pairs set_a`, Opus 5.5, records only, three replicates)

150 readings, all valid. Consistency: 46 of 50 pairs unanimous, 3 by majority, 1 with no majority; pairwise
agreement 0.94. Reference origins: unclear 39, human_genetics 3, earlier_drug_other_disease 2, and one each of
human_cell_perturbation, diseased_tissue, observational_biomarker, tumour_genetics, curated_pathway.

## Gate (`gate --set a`)

| | pairs | κ labeller raw | κ stored after check | reference unclear | labeller unclear |
|---|---|---|---|---|---|
| all compared | 49 | 0.21 | **0.42** | 39 | 6 |
| pairs with a record | 16 | 0.58 | 0.77 | 6 | 1 |
| pairs with no record | 33 | 0.00 | 0.00 | 33 | 5 |

Agreement gate **fails** (0.42 against 0.6). Collapsed-layer κ 0.41 (stored) and 0.20 (raw); the reference's origin
is the labeller's primary or secondary origin for 0.31 of pairs.

## Reading of the result

- The failure is structural, not a reading disagreement. The reference reads the stored records and nothing else; 33
  of the 49 compared pairs have no record, so the reference says unclear for all 33 while the labeller, with web
  search, labels 28 of them. On the 16 pairs where both saw records, κ is 0.58 raw and 0.77 after the check, above the
  gate.
- The check's "unverified" outcome is dominated by fetch failures (33 of 37), and it removes 37 labels from the tests.
- Two defects of this pilot are therefore the reference's access and the checker's reach, not the labeller's prompt.
  The registered design allows one prompt revision; the revision worth making is elsewhere.

## Options for the owner (none taken)

1. Give the reference reader the labeller's access (web search and fetch, still no provisional origin and no program,
   three replicates) and re-read set a; the labeller's prompt is unchanged, so set a is still unshaped by it. Set b stays
   reserved for a labeller revision.
2. Add the UniProt and HGNC protein names to the gene names and re-run `search`; fewer pairs without records.
3. Let the checker read a PubMed URL's record through E-utilities instead of fetching the page, throttle its fetches,
   and count a label with a cited but unopenable source as `unverified` only when every source failed; keep the
   registered "unverified as labelled" sensitivity.

Posted to #146 as a dated comment for decision; the full run waits for both gates.
