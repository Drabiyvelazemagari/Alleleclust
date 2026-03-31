AlleleClust: Deterministic Exact-Sequence Allele Clustering for CDS Datasets

AlleleClust is a standalone Python tool designed for the reproducible grouping of coding sequences (CDSs) into allele sets based on exact nucleotide sequence identity. The tool integrates sequence data with externally provided genomic metadata to enable structured analysis of allele diversity, host distribution, and replicon-associated dissemination in bacterial genomic datasets.

Overview and Design Principles

AlleleClust operates under a strictly deterministic framework in which allele definitions are based exclusively on exact nucleotide sequence identity. Unlike similarity-based clustering approaches, the tool does not rely on alignment, distance thresholds, or heuristic grouping strategies. Instead, each unique CDS sequence is treated as a distinct allele, ensuring maximal transparency and reproducibility.

The software is designed to function downstream of CDS extraction pipelines, particularly in conjunction with BAKO, from which it consumes both sequence data and associated metadata. All annotation is externally defined, and AlleleClust does not perform taxonomic or functional inference independently.

The core functionality consists of three stages: (i) parsing and validation of sequence and metadata inputs, (ii) deterministic clustering of sequences into allele groups, and (iii) generation of structured summaries and representative sequence sets. All outputs are explicitly defined and reproducible given identical inputs.

Input Parsing and Metadata Integration

AlleleClust requires two primary inputs: a CDS FASTA file and a corresponding metadata table (e.g., BAKO SUMMARY TSV). FASTA headers are parsed to extract accession identifiers, which are used to map sequences to metadata entries.

The metadata table is treated as the sole source of contextual information, including genus, species, and replicon type. The tool enforces strict consistency between FASTA and metadata inputs and reports discrepancies in accession coverage to prevent silent annotation errors.

Exact-Sequence Clustering

CDS sequences are grouped into allele clusters based on exact nucleotide sequence identity. Each unique sequence defines a single allele, and sequences are assigned to clusters using direct string comparison without alignment or similarity scoring.

Allele identifiers are assigned deterministically in the order of sequence processing, with optional zero-padding for consistent formatting. This ensures reproducibility within a dataset, although identifiers are not globally stable across independent datasets.

Allele Annotation and Characterization

For each allele cluster, AlleleClust computes a range of descriptive statistics and contextual annotations, including:

Number of member sequences
Sequence length (nucleotide and inferred amino acid length)
Set of contributing accessions
Distribution of genus and species
Replicon composition (chromosome, plasmid, phage, or unknown)

These annotations enable downstream analyses of allele distribution, host specificity, and potential dissemination patterns across genomic contexts.

Representative Sequence Selection

AlleleClust implements a structured and deterministic approach to representative sequence selection. Representatives are defined at multiple levels:

Species-level representatives, selecting one sequence per unique (genus, species) combination
Genus-level representatives, used when species-level annotation is absent or ambiguous
Replicon-specific representatives, capturing chromosome-, plasmid-, and phage-associated subsets
Single representative per allele, selected deterministically (minimum accession fallback) for use in downstream analyses

This multi-level representation facilitates flexible integration with phylogenetic, comparative, and structural workflows.

Output Structure and Reproducibility

AlleleClust produces a standardized set of outputs describing allele composition and membership:

Tabular summaries of allele membership and cluster characteristics
Representative mapping tables capturing taxonomic and replicon context
FASTA outputs containing allele-tagged sequences and representative subsets

All outputs are deterministically generated and consistently formatted to support reproducible downstream analysis. FASTA outputs include allele identifiers appended to sequence headers, enabling direct traceability between sequence data and cluster assignments.

GUI 

<img width="853" height="717" alt="image" src="https://github.com/user-attachments/assets/252fa097-27f3-4748-97c0-3e95e281cc3a" />


Implementation

AlleleClust is implemented as a single-file Python application (Python ≥3.9) with both command-line and optional graphical user interface (GUI) support via PyQt6. The tool requires minimal dependencies and is designed for local execution without reliance on external services or databases.

Scope and Limitations

AlleleClust is designed for exact-sequence grouping of CDS datasets and does not perform sequence alignment, similarity-based clustering, or evolutionary inference. As a result, it does not capture relationships between closely related but non-identical sequences (e.g., SNP variants) and should be interpreted as a strict definition of allele identity.

The accuracy of allele annotation depends on the quality and completeness of externally provided metadata. Inconsistent or incomplete metadata may affect downstream interpretation but does not influence clustering itself.
