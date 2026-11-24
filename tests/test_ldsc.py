"""Polygenicity (#147): LDSC heritability, clumping and heritability concentration per study (pipeline.ldsc)."""
import numpy as np
import pandas as pd

from pipeline.ldsc import concentration, loci, munge_input, parse_clumped, parse_h2_log


def test_munge_input_uses_beta_then_log_odds_ratio_and_hapmap_columns():
    h = pd.DataFrame({"hm_rsid": ["rs1", "rs2", "rs3"], "hm_effect_allele": ["A", "C", "G"],
                      "hm_other_allele": ["G", "T", "A"], "hm_beta": [0.1, np.nan, np.nan],
                      "hm_odds_ratio": [np.nan, 1.5, np.nan], "p_value": [1e-3, 1e-9, 0.5]})
    m = munge_input(h)
    assert list(m.columns) == ["SNP", "A1", "A2", "BETA", "P"] and len(m) == 2
    assert abs(m.BETA.iloc[1] - np.log(1.5)) < 1e-9


def test_munge_input_is_none_without_any_effect():   # Review Focus 2
    h = pd.DataFrame({"hm_rsid": ["rs1"], "hm_effect_allele": ["A"], "hm_other_allele": ["G"], "hm_beta": [np.nan],
                      "hm_odds_ratio": [np.nan], "p_value": [0.1]})
    assert munge_input(h) is None


LOG = """After merging with regression SNP LD, 1150000 SNPs remain.
Total Observed scale h2: 0.1234 (0.0123)
Lambda GC: 1.0912
Mean Chi^2: 1.1402
Intercept: 1.0123 (0.0089)
Ratio: 0.0877 (0.0634)
"""


def test_parse_h2_log():
    assert parse_h2_log(LOG) == {"h2_obs": 0.1234, "h2_obs_se": 0.0123, "intercept": 1.0123, "ratio": 0.0877,
                                 "lambda_gc": 1.0912, "mean_chi2": 1.1402, "n_snps": 1150000}


def test_parse_h2_log_with_a_ratio_note():
    text = LOG.replace("Ratio: 0.0877 (0.0634)", "Ratio < 0 (usually indicates GC correction).")
    assert parse_h2_log(text)["ratio"] is None and parse_h2_log(text)["h2_obs"] == 0.1234


CLUMPED = """ CHR    F         SNP         BP        P    TOTAL   NSIG    S05    S01   S001  S0001    SP2
   1    1      rs10    1000000  1.0e-12       10      0      0      0      0     10 rs11(1),rs12(1)
   6    1      rs60   26000000  1.0e-20       50      0      0      0      0     50 NONE
   6    1      rs61   33000000  1.0e-09        3      0      0      0      0      3 NONE

"""


def test_parse_clumped_and_mhc_counts_once():
    c = parse_clumped(CLUMPED)
    assert list(c.SNP) == ["rs10", "rs60", "rs61"] and list(c.CHR) == [1, 6, 6] and loci(c) == 2
    assert loci(parse_clumped("")) == 0


def test_concentration_is_top10_per_snp_h2_over_observed_h2():
    z2 = pd.Series([101.0, 51.0, 26.0])   # (z^2 - 1)/N = 0.1, 0.05, 0.025 -> 0.175 over 0.5
    share, top = concentration(z2, 1000, 0.5)
    assert abs(share - 0.35) < 1e-9 and len(top) == 3
    assert concentration(z2, 1000, 0.0) == (None, [])
    share12, top12 = concentration(pd.Series([11.0] * 12), 1000, 1.0)
    assert len(top12) == 10 and abs(share12 - 0.1) < 1e-9


def test_munge_input_accepts_the_plain_harmonised_columns():
    """Newer harmonised files (GWAS-SSF) have rsid, effect_allele, beta... without the hm_ prefix."""
    h = pd.DataFrame({"rsid": ["rs1", "rs2"], "effect_allele": ["A", "C"], "other_allele": ["G", "T"],
                      "beta": [0.1, np.nan], "odds_ratio": [np.nan, 2.0], "p_value": [1e-3, 1e-9]})
    m = munge_input(h)
    assert list(m.SNP) == ["rs1", "rs2"] and abs(m.BETA.iloc[1] - np.log(2.0)) < 1e-9


def test_munge_input_without_an_effect_column_is_none():
    h = pd.DataFrame({"rsid": ["rs1"], "effect_allele": ["A"], "other_allele": ["G"], "p_value": [1e-3]})
    assert munge_input(h) is None
