"""Polygenicity (#147): the disease-profile pipeline (pipeline.polygenicity) and its tables."""
import pytest

from pipeline import llm


@pytest.mark.db
def test_migration_156_tables_exist():
    with llm.connect() as conn, conn.cursor() as cur:
        cur.execute("""SELECT table_name FROM information_schema.tables WHERE table_schema = 'preclin'
                       AND table_name IN ('gwas_catalog_study', 'ldsc_run', 'gbd_prevalence', 'disease_polygenicity')""")
        assert {r[0] for r in cur.fetchall()} == {"gwas_catalog_study", "ldsc_run", "gbd_prevalence",
                                                   "disease_polygenicity"}
