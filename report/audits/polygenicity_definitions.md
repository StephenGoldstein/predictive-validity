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

## Scored count (Melissa's request, before the tertiles) — to post once the LDSC runs finish

**Scored count for #147 (build 5945a7e4, 2026-11-24).** Of the 531 diseases (2,321 US Phase 2 entries 2008–2019 on the rebuilt programmes, 2,275 before PR #188), {{SCORED_DISEASES}} have a score ({{SCORED_ENTRIES}} entries): a GWAS Catalog study of the disease or a subtype with public harmonised summary statistics, ≥ 90% European ancestry, case-control, a genome-wide platform, an LDSC heritability estimate at z ≥ 4 and a prevalence from GBD ({{K_GBD}} diseases) or Pan-UKB ({{K_UKB}}). The rest by reason: {{NO_SCORE_REASONS}}. Tertiles follow in the dated definitions comment.

## Dated definitions (the gate) — to post after the scored count

**Dated definitions for #147 (2026-11-24, Stephen; before any outcome is joined).** Everything below is fixed as of this comment; the tertile membership is committed at `report/audits/polygenicity_tertiles.csv` (commit {{COMMIT}}).

**Build and unit.** Programme build `5945a7e4-1b7d-4f1c-9767-2c3424b76512` (2026-11-23 19:42 PST, PR #188's rebuild, 50,228 programmes). Unit: non-oncology, non-infectious diseases with a US-developed Phase 2 entry in 2008–2019: 531 diseases, 2,321 Phase 2 entries. All 531 are profiled; 474 engaged-without-benefit stops in the unit, every one with its history row.

**Primary measure (as registered, source A).** Per disease: the GWAS Catalog study (release 2026-11-19) of the disease, else of a MONDO subtype, with public harmonised summary statistics, ≥ 90% European ancestry (the Catalog's ancestry table), cases and controls both reported, platform > 200,000 SNPs, taking the largest by effective sample size (4 / (1/cases + 1/controls)). LDSC 2.0.1 (Python 3 port) with the EUR LD scores (Zenodo 8182036) on HapMap3 SNPs, observed-scale h2 and its SE from the same run. **Heritability floor:** the estimate counts at z = h2 / SE ≥ 4 (your default; we had no reasoned value); below it the disease is no-score (`low_z`). Liability scale: h2_obs × K²(1−K)² / (P(1−P) φ(t)²), P = the study's case fraction, t = Φ⁻¹(1−K), K = GBD 2021 global prevalence of the matching cause (all ages, both sexes; GBD 2023 round, Melissa's export of 2026-11-23; hand map `pipeline/data/gbd_cause_mondo.csv`), else the Pan-UKB European case fraction of the disease's ICD-10 code, a four-character code matching through its three-character parent ({{N_K_PARENT}} diseases take K that way); K's source and code are recorded per disease. Loci: PLINK 1.90b7.7 clumping of the same summary statistics at p < 5 × 10⁻⁸, r² < 0.1 within 1,000 kb, against 1000 Genomes EUR (Zenodo 6614170) cut to HapMap3 SNPs; the extended MHC (chr6:25–34 Mb) counts once. Score = loci / h2_liab, kept when 0 < h2_liab ≤ 1. A disease at or above the floor with **zero loci** is no-score (`zero_loci`), not score 0, as you ruled. Scored: {{SCORED_DISEASES}} diseases ({{SCORED_ENTRIES}} entries). No-score group (outside the tertiles and P1–P3, reported as a group): {{NO_SCORE_DISEASES}} diseases ({{NO_SCORE_ENTRIES}} entries), by reason {{NO_SCORE_REASONS}}.

**Tertiles.** Cutpoints at the 33.3rd and 66.7th percentiles of the score over the scored diseases: {{CUT_LO}} and {{CUT_HI}}; a disease at a cutpoint goes down. Bottom {{N_BOTTOM}} diseases ({{E_BOTTOM}} entries), middle {{N_MIDDLE}} ({{E_MIDDLE}}), top {{N_TOP}} ({{E_TOP}}). The full membership with each disease's study, h2, K and its source and code, loci, score and Mendelian form: `report/audits/polygenicity_tertiles.csv`.

**Tests (items 4–7 of our list, unanswered, so stated here as used).** P1 on Phase 2 and 3 entries, the "target not causal" failure mode of the current `engaged_modes` (PR #183's three groups; `not_causal_basis` reported as a descriptive), COVARS_23 as in `target_origin`; P2 and P3 on Phase 2 entries with COVARS; top tertile against bottom, middle reported; entries with a resolved primary target only, entries without one reported per tertile; one-sided p in the registered direction before Holm (as the #140 confirmation), two-sided beside; errors clustered by disease. P3 = the top × genetics (g3 > 0) interaction on S4, with the genetics ratio within each tertile beside it; estimable only when each cell holds ≥ 15 target sets and an event.

**Sensitivities (as amended on 2026-11-24).** Each of the first two repeats P1–P3 with Holm within itself; the third is one two-sided test outside any family. (1) GWAS size: + log10 n_eff, n_eff = 4 / (1/cases + 1/controls). (2) Mendelian form over all 531 diseases, the no-score group included: a gene with Orphanet, ClinGen, Genomics England or Gene2Phenotype evidence for the disease, a subtype or a disease at Lin ≥ 0.8 in the full Open Targets 26.09 direct-association file (every gene; OMIM is not a datasource there, so it enters only through those four); "no Mendelian form" in place of the top tertile. {{N_MENDELIAN}} of 531 have a form. (3) Negative control: P1's model (Phase 2 and 3 entries, COVARS_23, clustered by disease) with the outcome a stop whose traced cause in `preclin.stop` is `sponsor_decision` or `operational`; entries whose traced stop is a safety stop are left out; PK, dose and formulation stops are efficacy in that table and so count on the biology side. Denominator as we read it: every entry with a resolved primary target in the top and bottom tertiles, minus the safety stops. Reported two-sided; no association expected. The continuous-score and concentration sensitivities are dropped, as you ruled.

**Descriptives, none in a Holm family (item 8).** Gene counts from the full Open Targets association file, every gene (GWAS-linked genes per 10,000 cases, rare/Mendelian genes), the Catalog's all-study loci, inclusion at z ≥ 2, GBD-K diseases alone, O'Connor's Me for the {{N_OCONNOR}} scored diseases it covers, tertile by disease area, the no-score group's rates, and the DALY share (WHO GHE 2021) in diseases with a Mendelian form or a drugged large-effect gene.

**Rules you ruled on (2026-11-24, before any run was trusted).** Exome-wide studies (exome sequencing or exome arrays, by the Catalog's genotyping technology) are never chosen: their files hold exome variants only, and 27 diseases would otherwise have lost a score that 16 of them have in a genome-wide study. Genome-wide significant SNPs the clumping panel lacks (about a tenth of HapMap3) count as loci when they lie more than 1 Mb from every index SNP and from each other, taken by p; the count is in each run's log. The extended MHC counts once among the loci. p below 1 × 10⁻³⁰⁰ is floored so munge keeps the SNP. Pan-UKB ICD-10 phenotypes are three-character codes, so a disease's four-character codes match through their parent (count above). A disease with zero loci is no-score. Holm keeps the family of three when a test is not estimable. Ties on effective sample size break on the accession.

**Checks before the join (your two).** (1) Within-area variance: {{WITHIN_AREA_SHARE}} of the log10 score's variance over the {{SCORED_DISEASES}} scored diseases lies within disease area ({{N_AREAS}} areas). (2) Power: {{POWER_REPS}} replicates on the real entries, tertiles, genetic support, covariates and disease clusters; each outcome drawn from its rate pooled over the whole unit (not-causal stop {{BASE_NC}}, Phase 2 advance {{BASE_S4}}), a disease random effect of SD 0.5 on the log scale, and the posited ratio; P3 assumes a twofold S4 advantage of genetic support, multiplied by the posited ratio in the top tertile. Power at one-sided 0.05 (and at 0.05/3, Holm's threshold for the smallest p): {{POWER_TABLE}}. No outcome was read by tertile for either check. `report/output/polygenicity/check.json` and `power.json`.

**Provenance.** Genetic support in P3 (g3 > 0) comes from evidence-line build `{{ELB_ID}}` ({{ELB_AT}}, commit {{ELB_COMMIT}}); the unit from programme build `5945a7e4-1b7d-4f1c-9767-2c3424b76512`; tools and references per run in `preclin.ldsc_run`.

**Caveats fixed now.** K is a point prevalence (GBD "prevalence", 2021), not lifetime risk, so h2_liab is lower than a lifetime-risk conversion would give for episodic diseases; the Pan-UKB fallback is a cohort case fraction. {{N_RUN_FAILED}} chosen studies have no usable estimate (no effect column in the file, or the run failed) and sit in the no-score group. Tools and references are recorded per run in `preclin.ldsc_run`.
