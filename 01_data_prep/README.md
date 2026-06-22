# Reference and Control Genome Preparation

This directory contains scripts used to prepare the reference and control genomes required for the EM-seq analysis workflow.

## Reference genome

The reference genome is used by the main Snakemake workflow.

1. Navigate to the directory that will contain the reference genome for the main Snakemake workflow.
2. Run:

```bash
bash prep_ref_genome.sh
```

This script downloads the reference genome, renames the FASTA file, and prepares the Bismark genome indices.

The completed reference genome directory should remain in the location used by the main Snakemake workflow.

---

## Control genomes

The control genomes are used by the control Snakemake workflow.

Save the control genome sequences as:

```text
pUC19.fa
lambda.fa
```

Place both FASTA files in the directory used by the control Snakemake workflow.

From that directory, run:

```bash
bash prep_control_genomes.sh
```

This script prepares the Bismark genome indices required for alignment to the control genomes.

---

## Notes

* Run each preparation script only once for each genome.
* Ensure that Bismark and its dependencies are available before running the scripts.
* The generated Bismark index files should remain with their corresponding reference genome directories for use by the Snakemake workflows.
