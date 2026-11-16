"""Small simulated teaching data sets, written on demand by the web interface
(and usable from Python). Standard library only; every writer produces the
same bytes on every run.

- `association.write_demo(folder)`: a candidate-gene cohort for
  `genotype_association` and `association_power` (screen 30).
- `contaminants.write_demo(folder)`: references, a "raw material" and a
  reagent blank for `build_screen_database` and `contaminant_screen`
  (screen 31).
- `rnaseq_fastq.write_demo(folder)`: a transcriptome and six paired-end
  samples for `build_transcript_index`, `quantify_samples` and then
  `differential_expression` (screen 22).
"""
