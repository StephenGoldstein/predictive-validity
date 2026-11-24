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
