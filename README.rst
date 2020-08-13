=======================================================
``qctk``: Quality Control Toolkit for Clinical NGS Data
=======================================================

This software package allows for the easy generation of reports for NGS data with a focus on clinical applications.
It allows to process variant call files (VCF/BCF), alignment files (BAM/CRAM), and raw read data (FASTQ).
From these, it generates quality control reports aimed clinical bioinformaticians and data analysts.

=========
Use Cases
=========

- Compare the sex indicated the molecular data to expected sex from meta data (pedigree file).
- Check the actual relationship indicated by the sample variants with the expected relationship from the pedigree file.
- Find cross-contaminations from a all-to-all relationship heatmap and variance in alternate allele frequency.

Right now, only germline data is supported.

=======
Reports
=======

``qctk`` can generate the following reports:

- Relationship scatter plot for checking pedigree structure.
- Relationship heatmap for for checking pedigree structure and detecting sample cross-contamination.
- Alternate allele balance boxplot for detecting sample cross-contamination.
- Variant coverage information.

The following are only available when starting from BAM/CRAM:

- Target coverage information for targeted sequencing ("how many exons have coverage above 10x/20x/...")?
- Overall coverage statistics for whole genome sequencing.

===========
Input Files
===========

- Family structure is provided in **PLINK PED** format.
- For targeted sequencing, the targeted exons are provided as **BED** files.
- Variants can be extracted form **VCF/BCF** files and the variant coverage information can be used.
- Alternatively, **BAM/CRAM** files can be queried for their variant status at (predefined) polymorphic sites.
  This is more time (and I/O) consuming but yields higher quality results as this yields comprehensive results of wild-type sites.
- Further, raw read data in **FASTQ** files can be analyzed.
  For this, an "in silico SNP-array" approach is taken (similar to the one by NGSCheckMate).

============
How It Works
============

============
Related Work
============

While functionality in ``qctk`` can also be found in other tools it focuses on the use cases most relevant to DNA NGS data analysis.
``qctk``

- is limited to targeted (exome/panel) and whole genome sequencing data.
- assumes that the sequencing and alignment process is well established and only provide summary statistics that allow to detect low quality samples.
  For further insight, tools such as ``samtools stats`` can provide more insight.
- focuses on checking variants for their quality (coverage) and to to detect sample swaps and contaminations.

Some References:

- Pedersen, Brent S., et al.
  ["Somalier: rapid relatedness estimation for cancer and germline studies using efficient genome sketches."](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-020-00761-2)
  Genome medicine 12.1 (2020): 1-9.
- Pedersen, Brent S., and Aaron R. Quinlan.
  ["Who’s who? Detecting and resolving sample anomalies in human DNA sequencing studies with peddy."](https://www.sciencedirect.com/science/article/pii/S0002929717300174)
  The American Journal of Human Genetics 100.3 (2017): 406-413.
- Lee, Sejoon, et al.
  ["NGSCheckMate: software for validating sample identity in next-generation sequencing studies within and across data types."](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5499645/)
  Nucleic acids research 45.11 (2017): e103-e103.
- Schröder, Jan, Vincent Corbin, and Anthony T. Papenfuss.
  ["HYSYS: have you swapped your samples?."](https://academic.oup.com/bioinformatics/article/33/4/596/2624551)
  Bioinformatics 33.4 (2017): 596-598.
- Messerschmidt, Clemens, Manuel Holtgrewe, and Dieter Beule.
  ["HLA-MA: simple yet powerful matching of samples using HLA typing results."](https://academic.oup.com/bioinformatics/article/33/14/2241/3064507)
  Bioinformatics 33.14 (2017): 2241-2242.
