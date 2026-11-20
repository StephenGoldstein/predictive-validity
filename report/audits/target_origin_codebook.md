# Target origin hand check: codebook

One row per target × disease pair. The question: what first proposed this target for this disease, before any drug on
it was tested in patients for it? Read the records in the package (the same ones the model saw: the earliest 25 and
the latest 15 dated before the cutoff year, abstracts cut at 1,500 characters; the link opens the full record). You
may also follow the two search links, PubMed and Google Patents, which are already restricted to dates before the
cutoff year. Anything you find that way counts only if its own date is before the cutoff year. You never see the
model's answer.

## Origin (one primary; the definitions are the database's own)

- human_genetics: germline variants link the target to the disease (GWAS, gene burden, Mendelian)
- tumour_genetics: the gene is recurrently mutated in the cancer
- diseased_tissue: the gene's RNA level differs in diseased tissue
- mouse_genetics: knocked-out or mutant mice show the disease
- curated_pathway: a curated pathway shows a disease mutation disrupting the target
- cell_screen: a gene-knockout screen in cells
- earlier_drug_other_disease: a drug on the target was approved for another disease
- observational_biomarker: the target or its product tracks the disease in patients
- human_cell_perturbation: perturbing the target in human cells or tissue changes a disease phenotype
- animal_model_pharmacology: modulating the target in an animal disease model changes it
- clinical_observation: an effect seen in people, such as a side effect
- phenotypic_screen: a molecule found without a target hypothesis came first
- unclear: no source settles it

The origin is the kind of evidence that first proposed the pair, not the strongest evidence and not the latest. Where
the sources show two kinds, the earlier is the origin and the other goes in secondary_origin. The provisional origin
shown on each card is the earliest dated item in the evidence databases; the model saw it too. Confirm or correct it
from the sources; do not copy it without a source.

## Grade

- stated: a source says why the target was chosen for the disease (a rationale in an introduction or discussion, a
  review written before the cutoff, a patent's background).
- pattern: no source says so, but a source shows the defining observation of the category, published before any
  other kind of evidence for the pair.
- none: the sources would leave the origin to inference. Then the origin is unclear, and quote and source stay empty.

A label graded stated or pattern carries a verbatim quote from its source, the source's id as shown on the card
(PMID:… or the Europe PMC id) or its URL, and the source's own publication year.

## Discovery route (optional, descriptive)

Where did the first program's drug come from? phenotypic_screen (a screen with no target hypothesis), target_screen
(a high-throughput, fragment or virtual screen against the target; immunization or display for a biologic),
de_novo_design (antisense, siRNA, structure-based design), derived_from_active (a known drug, an endogenous ligand, a
natural product), repurposed_drug (an existing drug tried in this disease), unclear. Fill it in when a record settles
it; leave it blank otherwise.

## Filling the form

Every pair needs an origin and a grade. Quote, source and source year are required when the grade is stated or
pattern. Note is free text for anything the next reader should know. Answers save in the browser as you type; when
done, use "Copy answers as CSV" and send the text, or "Download CSV".
