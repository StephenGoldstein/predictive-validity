# Polygenicity and failure (#147): announcements and the dated definitions

The text of each #147 comment this work posts, with its id and time once posted. Registration: Melissa's comment of
2026-11-19 (IC_kwDOTg_ees8AAAABou6xvw); decisions of 2026-11-23: no-score rule and secondary measures
(IC_kwDOTg_ees8AAAABqDb0FA), infectious diseases out (IC_kwDOTg_ees8AAAABqFfZeg), heritability source, prevalence
and the four sensitivities (issuecomment-7123250021).

## Announcement (before any database write) — posted 2026-11-23 23:10 PST, issuecomment-7124701902

**Database writes and tools for #147, starting now.**

- Migration 160: four new tables, nothing existing changed. `preclin.gwas_catalog_study` (one row per GWAS Catalog study: mapped MONDO ids, cases, controls, European share, summary-statistics path), `preclin.ldsc_run` (one row per LDSC run: observed-scale h2, intercept, clumped loci, heritability concentration), `preclin.gbd_prevalence` (GBD 2023 prevalence by cause, 2021), `preclin.disease_polygenicity` (one row per candidate disease: the profile, K and its source, the score, the secondary measures).
- Loads, in order: the Catalog study table (release 2026-11-19), the GBD file Melissa sent, then one LDSC run per chosen study (about 208, hours, resumable), then the disease table. Each load registers its source in `preclin.evidence_source_release`.
- Tools: `ldsc` 2.0.1 (the Python 3 port) in its own Python 3.11 venv; PLINK 1.90b7.7 from bioconda; EUR LD scores (Zenodo 8182036); HapMap3 SNP list (Zenodo 7773502); 1000 Genomes Phase 3 European genotypes (Zenodo 6614170) cut to HapMap3 SNPs for clumping at r² < 0.1 within 1 Mb.
- No outcome is joined to any of this until the dated definitions comment is posted.

### Downloads (for Stephen's approval, not part of the comment)

| File | Source | Size |
|---|---|---|
| eur_w_ld_chr.tar.gz | zenodo.org/records/8182036 | 33 MB |
| w_hm3.snplist.gz | zenodo.org/records/7773502 | 5 MB |
| 1000G_EUR.bed, .bim, .fam, .afreq | zenodo.org/records/6614170 | 231 + 52 + 0 + 56 MB, cut to HapMap3 SNPs then deleted |
| gwas-catalog-download-ancestries-v1.0.3.1.tsv | ftp.ebi.ac.uk/pub/databases/gwas/releases/latest | about 20 MB |
| O'Connor 2019 Supplementary Data S1 (mmc2.xlsx) | cell.com, PMC | under 1 MB |
| ldsc 2.0.1 wheel and its dependencies; plink 1.90b7.7 conda package | PyPI; bioconda | about 60 MB together |
| harmonised summary statistics, about 208 files | ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics | 50 MB to 2 GB each, one at a time, deleted after use |

## Scored count (Melissa's request, before the tertiles) — posted 2026-11-25 14:48 PST, issuecomment-7176802793

**Scored count for #147 (build 5945a7e4, 2026-11-25).** Of the 531 diseases (2,321 US Phase 2 entries 2008–2019 on the rebuilt programmes, 2,275 before PR #188), 112 have a score (1,092 entries): a GWAS Catalog study of the disease or a subtype with public harmonised summary statistics, ≥ 90% European ancestry, case-control, a genome-wide platform, an LDSC heritability estimate at z ≥ 4 and a prevalence from GBD (48 diseases) or Pan-UKB (64). The rest by reason: `no_study` 195 (347 entries), `low_z` 94 (345 entries), `no_sumstats` 81 (293 entries), `no_k` 16 (51 entries), `targeted_array` 12 (18 entries), `zero_loci` 9 (62 entries), `not_european` 7 (29 entries), `no_effect` 3 (76 entries), `no_case_control` 2 (8 entries). Tertiles follow in the dated definitions comment.

## Two questions before the gate — posted 2026-11-25 17:09 PST, issuecomment-7179483929

**Two questions before the gate (Stephen, 2026-11-25).** Everything else in the definitions comment is registered or ruled. These two change what gets fixed, so we hold the gate for your answer. No outcome has been joined.

1. **Unusable files.** The largest eligible study for Alzheimer disease (48 Phase 2 entries) is GCST90444373, Belloy 2024's X-chromosome-only study. For Crohn disease (27 entries) it is GCST000879, a 2010 file with p-values only. Both run as no-score (`no_effect`). We propose: a chosen study whose file yields no usable autosomal effect is excluded and the next eligible study is chosen (Alzheimer → GCST90027158, Bellenguez 2022; Crohn → GCST003044, Liu 2015). Juvenile idiopathic arthritis (1 entry) has the same problem and no fallback. Yes or no?

2. **P3's cell minimum.** On this unit the bottom-tertile cell with genetic support holds 9 target sets (12 Phase 2 entries) against the registered 15. P3 is therefore not estimable in any power replicate (other cells: 27 / 60, 216 / 380, 138 / 179). As registered we report it not estimable and keep Holm m = 3. If you want a lower minimum, now is the time. We have read no outcome by tertile.

Our pick: yes on 1, keep 15 on 2.

Melissa's answer, 2026-11-25 17:18 PST, issuecomment-7179634211: (1) yes, the fallback is part of the selection rule for every disease alike; (2) keep 15, report P3 not estimable with the reason (the bottom × support cell is small by construction: fewest loci, fewest supported targets), and post the Mendelian-form split's four P3 cell counts with the gate.

## Dated definitions (the gate), full method — the short form below was posted 2026-11-25 19:14 PST, issuecomment-7181532783

**Dated definitions for #147 (2026-11-25, Stephen; before any outcome is joined).** Everything below is fixed as of this comment. Tertile membership: `report/audits/polygenicity_tertiles.csv` (commit 5fb4a30).

**Build and unit**
- Programme build `5945a7e4-1b7d-4f1c-9767-2c3424b76512` (2026-11-23 19:42 PST, PR #188's rebuild).
- Unit: non-oncology, non-infectious diseases with a US-developed Phase 2 entry in 2008–2019. 531 diseases, 2,321 Phase 2 entries, all profiled; 474 engaged-without-benefit stops, each with its history row.

**Primary measure (as registered, source A)**
- Study: the GWAS Catalog study (release 2026-11-19) of the disease, else of a MONDO subtype, with public harmonised summary statistics, ≥ 90% European ancestry, cases and controls reported, platform > 200,000 SNPs, never exome-wide. Largest n_eff = 4 / (1/cases + 1/controls); ties break on accession. A chosen study whose run merges fewer than 200,000 HapMap3 SNPs or yields no usable autosomal effect is excluded and the next eligible study is chosen (your ruling of 2026-11-25): Alzheimer disease → GCST90027158 (Bellenguez 2022); Crohn disease → GCST90473822 (your named GCST003044 is Immunochip data, 48,002 merged SNPs, so the merged-SNP rule excludes it); juvenile idiopathic arthritis stays no-score.
- h2: LDSC 2.0.1, EUR LD scores (Zenodo 8182036), HapMap3 SNPs. Floor z = h2 / SE ≥ 4 (your default); below it, no-score (`low_z`).
- Liability scale: h2_obs × K²(1−K)² / (P(1−P) φ(t)²), P = the study's case fraction, t = Φ⁻¹(1−K).
- K: GBD 2021 global prevalence of the matching cause (your export of 2026-11-23; map `pipeline/data/gbd_cause_mondo.csv`), else the Pan-UKB European case fraction of the ICD-10 code, four-character codes through their three-character parent (45 scored diseases). Source and code recorded per disease.
- Loci: PLINK 1.90b7.7 clumping at p < 5 × 10⁻⁸, r² < 0.1 within 1,000 kb, 1000 Genomes EUR (Zenodo 6614170) cut to HapMap3; the extended MHC (chr6:25–34 Mb) counts once; significant SNPs the panel lacks count when > 1 Mb from every index SNP and each other; p < 1 × 10⁻³⁰⁰ floored for munge.
- Score = loci / h2_liab, kept when 0 < h2_liab ≤ 1. Zero loci above the floor is no-score (`zero_loci`), as you ruled.
- Scored: 114 diseases (1,167 entries). No-score, reported as a group: 417 diseases (1,154 entries): `no_study` 195, `low_z` 94, `no_sumstats` 81, `no_k` 16, `targeted_array` 13, `zero_loci` 9, `not_european` 7, `no_case_control` 2.

**Tertiles**
- Cutpoints at the 33.3rd and 66.7th percentiles of the score: 164.37 and 796.29; a disease at a cutpoint goes down (none is).
- Bottom 38 diseases (250 entries), middle 38 (350), top 38 (567).

**Tests (items 4–7 of our list, unanswered, so stated as used)**
- P1: Phase 2 and 3 entries; outcome the "target not causal" mode of the current `engaged_modes` (PR #183); COVARS_23 as in `target_origin`.
- P2, P3: Phase 2 entries with COVARS. P3 = the top × genetics (g3 > 0) interaction on S4; estimable only when each cell holds ≥ 15 target sets and an event.
- All: top against bottom, middle reported; entries with a resolved primary target only, the rest reported per tertile; one-sided p in the registered direction before Holm (as the #140 confirmation), two-sided beside; errors clustered by disease; Holm keeps m = 3 when a test is not estimable.

**Sensitivities (your 2026-11-24 amendment, as we read it)**
- (1) GWAS size: + log10 n_eff; P1–P3 with Holm within.
- (2) Mendelian form over all 531: Orphanet, ClinGen, Genomics England or Gene2Phenotype evidence in the full Open Targets 26.09 direct-association file for the disease, a subtype or Lin ≥ 0.8 (OMIM is not a datasource there); "no form" in place of the top tertile; 333 of 531 have a form; P1–P3 with Holm within.
- (3) Negative control: P1's model, outcome a stop traced in `preclin.stop` to `sponsor_decision` or `operational`, safety stops left out, PK/dose/formulation stops on the biology side. Denominator (confirmed by you today): entries with a resolved primary target in the top and bottom tertiles, minus safety stops. One two-sided test outside any family.
- Dropped as you ruled: continuous score, concentration.

**Descriptives (item 8, no Holm family):** as amended; O'Connor's Me covers 5 scored diseases.

**Checks before the join (your two)**
- Within-area variance: 79% of the log10 score's variance over the 114 scored diseases lies within disease area (13 areas).
- Power: 200 replicates on the real entries, tertiles, genetic support, covariates and clusters; each outcome drawn at its rate pooled over the top and bottom tertiles together, never by tertile (not-causal stop 9.4%, Phase 2 advance 17.7%), disease random effect SD 0.5 on the log scale; P3 assumes a twofold S4 advantage of genetic support times the ratio in the top tertile. One-sided 0.05, Holm's 0.05/3 in parentheses:

| Test | RR | Power (Holm) |
|---|---|---|
| P1 not-causal stop, top vs bottom | 1.25 / 1.5 / 2.0 | 0.30 (0.20) / 0.53 (0.36) / 0.86 (0.71) |
| P2 Phase 2 advance | 0.8 / 0.67 / 0.5 | 0.28 (0.18) / 0.50 (0.38) / 0.81 (0.64) |
| P3 genetics × top on advance | any | not estimable in any replicate |

- P3's bottom-tertile cell with genetic support holds 12 target sets (15 entries) against the 15 required (top with support 33 / 74, top without 233 / 409, bottom without 143 / 193; counts in `report/audits/polygenicity_p3_primary_cells.json`, outcome columns dropped before counting). No event was read. P3 is reported as not estimable and the family keeps Holm m = 3 (your ruling of 2026-11-25): the cell is small by construction, since the score's numerator is the locus count, so the least polygenic diseases have the fewest loci and the fewest supported targets.
- P3 under the Mendelian-form sensitivity, same minimum of 15, outcome columns dropped before counting: no form × support 5 target sets (6 Phase 2 entries), no form × none 314 (539), form × support 153 (278), form × none 501 (1,125). Not estimable either: diseases without a Mendelian form carry almost no genetic support. Counts in `report/audits/polygenicity_p3_mendelian_cells.json`.
- No outcome was read by tertile for either check: `report/audits/polygenicity_check.json`, `report/audits/polygenicity_power.json`.

**Provenance:** genetic support from evidence-line build `37f24080-9f7a-47a7-8656-f287befabc6e` (2026-11-24 23:03 UTC, commit ee87412); tools and references per run in `preclin.ldsc_run`.

**Caveats fixed now:** K is a point prevalence, not lifetime risk, so h2_liab runs low for episodic diseases; the Pan-UKB fallback is a cohort case fraction.

## Gate comment as posted — 2026-11-25 19:14 PST, issuecomment-7181532783

**Dated definitions for #147 (2026-11-25, Stephen; before any outcome is joined).** Fixed as of this comment. Full method and tertile membership: `report/audits/polygenicity_definitions.md` and `report/audits/polygenicity_tertiles.csv` (commit a23538f).

**Unit**
- Build `5945a7e4-1b7d-4f1c-9767-2c3424b76512`.
- 531 non-oncology, non-infectious diseases with a US Phase 2 entry in 2008–2019; 2,321 Phase 2 entries.

**Score (source A, as registered)**
- Study: GWAS Catalog (2026-11-19) harmonised summary statistics, ≥ 90% European, cases and controls reported, platform > 200,000 SNPs, never exome-wide; largest n_eff wins.
- Fallback (your ruling today): a run that merges < 200,000 HapMap3 SNPs or has no usable autosomal effect is skipped for the next eligible study. Alzheimer → GCST90027158. Crohn → GCST90473822 (GCST003044 is Immunochip, 48,002 merged SNPs). Juvenile idiopathic arthritis stays no-score.
- h2: LDSC 2.0.1, EUR LD scores, HapMap3 SNPs; floor z ≥ 4; liability scale with K from GBD 2021, else the Pan-UKB case fraction (45 scored diseases through the three-character ICD-10 parent).
- Loci: PLINK clumping at p < 5 × 10⁻⁸, r² < 0.1 within 1,000 kb, 1000 Genomes EUR; the extended MHC counts once.
- Score = loci / h2_liab, kept when 0 < h2_liab ≤ 1; zero loci is no-score.
- Scored: 114 diseases, 1,167 entries. No-score 417: no_study 195, low_z 94, no_sumstats 81, no_k 16, targeted_array 13, zero_loci 9, not_european 7, no_case_control 2.

**Tertiles**
- Cutpoints 164.37 and 796.29 (33.3rd and 66.7th percentiles; a tie goes down, none occurs).
- Bottom 38 diseases (250 entries), middle 38 (350), top 38 (567).

**Tests**
- P1: not-causal stop, Phase 2 and 3, COVARS_23. P2: Phase 2 advance, COVARS. P3: top × genetic support on S4, estimable only when every cell holds ≥ 15 target sets and an event.
- Top vs bottom, middle reported; resolved primary target only; one-sided in the registered direction, two-sided beside; errors clustered by disease; Holm m = 3, kept when a test is not estimable.

**Sensitivities (your 2026-11-24 amendment)**
- GWAS size: + log10 n_eff.
- Mendelian form (Orphanet, ClinGen, Genomics England, Gene2Phenotype in Open Targets 26.09): 333 of 531 have one; "no form" replaces the top tertile.
- Negative control: P1's model, outcome a sponsor_decision or operational stop, safety stops out of the denominator; one two-sided test outside the family.
- Dropped: continuous score, concentration. Descriptives (item 8) as amended, no Holm family.

**Pre-join checks (no outcome read by tertile)**
- Within-area variance: 79% of the log10 score's variance lies within disease area.
- Power, 200 replicates, outcomes drawn at the rate pooled over the top and bottom tertiles together:

| Test | RR | Power (Holm) |
|---|---|---|
| P1 not-causal stop | 1.25 / 1.5 / 2.0 | 0.30 (0.20) / 0.53 (0.36) / 0.86 (0.71) |
| P2 Phase 2 advance | 0.8 / 0.67 / 0.5 | 0.28 (0.18) / 0.50 (0.38) / 0.81 (0.64) |
| P3 genetics × top | any | not estimable |

- P3 cells (outcome columns dropped): bottom × support holds 12 target sets against the 15 required, so P3 is not estimable and Holm keeps m = 3 (your ruling today; the cell is small by construction, since the score's numerator is the locus count). Mendelian split: no form × support holds 5. Counts in `report/audits/polygenicity_p3_*_cells.json`.

**Provenance:** genetic support from evidence-line build `37f24080-9f7a-47a7-8656-f287befabc6e`; tools per run in `preclin.ldsc_run`. K is a point prevalence, so h2_liab runs low for episodic diseases.

## Readings — posted 2026-11-25 19:45 PST, issuecomment-7181992251

**Readings for #147 (2026-11-25, Stephen; outcomes joined after the gate, issuecomment-7181532783).** Tables and figures come with the PR.

**Registered tests (top against bottom tertile, one-sided in the registered direction, Holm m = 3)**

| Test | Outcome | Top | Bottom | RR (95% CI) | One-sided p | Holm | Reading |
|---|---|---|---|---|---|---|---|
| P1 | not-causal stop, Phase 2 and 3 | 59 / 673 | 31 / 286 | 1.00 (0.64–1.55) | 0.51 | 1.00 | not replicated |
| P2 | Phase 2 advance | 89 / 483 | 33 / 208 | 1.33 (0.74–2.39) | 0.83 | 1.00 | not replicated |
| P3 | genetics × top on advance | 18 / 74 with support, 71 / 409 without | 1 / 15, 32 / 193 | not estimable | | | 12 target sets in bottom × support against 15 |

- P2's registered direction was fewer advances in the top tertile; the point estimate runs the other way.
- P3 as you ruled: not estimable, Holm keeps m = 3.

**Sensitivities**
- GWAS size (+ log10 n_eff): P1 RR 1.04, same direction, not confirmed; P2 RR 1.45, not replicated; P3 not estimable.
- Mendelian form (no form against form, all 531): P1 RR 1.18 (0.89–1.56), one-sided p 0.13, Holm 0.38, same direction, not confirmed; P2 RR 0.99, same direction, not confirmed; P3 not estimable.
- Negative control (sponsor or operational stops, two-sided): RR 0.82 (0.55–1.22), p 0.33, no association.

**Rates by tertile** (not-causal stop over Phase 2 and 3 entries; advance over Phase 2 entries): bottom 9.1% and 14.0%; middle 11.1% and 15.1%; top 7.6% and 16.0%; no-score 10.0% and 18.0%.

**Burden descriptive (item 8)**: of 2021 noncommunicable DALYs, 51% fall in causes mapped to a disease with a Mendelian form, 53% to one with a drugged large-effect gene, 59% to either (333 and 188 of the 531 diseases).

Neither registered hypothesis replicates. The power table put the mid effects at 0.3 to 0.5, and both point estimates sit at or across the null. The PR follows with tables, figures and narrative.
