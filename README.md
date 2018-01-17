This repository contains details about the data availability and the computational methods developed for the *submitted* manuscript (add link).


## Data

All the sequencing data have been deposited in the [ArrayExpress database](https://www.ebi.ac.uk/arrayexpress/) at EMBL-EBI under accession number E-MTAB-6367 (currently reviewer access only):

 https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6367
 
 
## Code

- **Pilot screen Pool8 A375**:
  - [Preliminary exploration](scripts/20160503_exploration.md)
  - [Data processing 1](scripts/20160616_processing.md)
  - [Analysis 1](scripts/20160623_analysis.md)
  - [Data processing 2](scripts/20170301_processing.md)
  - [Analysis 2](scripts/20170301_analysis.md) 

- **Genome-wide screen A375**:
  - *1st sequencing round*:
    - Pools 1, 2, 5 and 10: [data processing](scripts/20160627_processing.md) and [analysis](scripts/20160627_analysis.md)
    - Pools 3, 4, 6 and 7: [data processing](scripts/20160902_processing.md) and [analysis](scripts/20160902_analysis.md)
    - Pools 8, 9, 11 and 12: [data processing](scripts/20161011_processing.md) and [analysis](scripts/20161011_analysis.md)
  - *2nd sequencing round*:
    - Pools 1-6 and 10-12: [data processing](scripts/20161108_processing.md) and [analysis](scripts/20161108_analysis.md)
  - *3rd sequencing round*:
    - Pools 3, 11 and 12: [data processing](scripts/20170105_processing.md) and [analysis](scripts/20170105_analysis.md)
  - *Integration* of all sequencing runs: creating the [integrated dataset](scripts/20170120_integration_finaldataset.md) and additional [plots](scripts/20170328_integration_finaldataset_additional_plots.md)

- **Focussed screens A375 and HT1080**:
  - [Data processing](scripts/20170404_processing.md)
  - [Analysis](scripts/20170404_analysis.md)

- **Functional analyses**:
  - [GO and KEGG](scripts/20170306_functional_analysis.md)
  - COSMIC: [preliminary](scripts/20170314_cosmic.md) and [advanced](scripts/20170821_cosmic_tcga_icgc.md)
 
- **Literature and database searches**:
  - [UniprotKB and GO](scripts/20170727_literaturesearch.md#go-and-uniprot-2)
  - [MeSH and PolySearch2](scripts/20170727_literaturesearch.md#polysearch2-2)

- **Submission to ArrayExpress**:
  - [Script](scripts/20171014_data.md)


## Sofware requirements

- [FastQC v0.11.3](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [FASTX-Toolkit v0.0.14](http://hannonlab.cshl.edu/fastx_toolkit/)
- [Bowtie 2 v2.2.6](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- Unix [tools](http://opengroup.org/unix) (cat, cut, awk, sort and uniq)
- [Python programming language v2.7.10](https://www.python.org/). Libraries used:
  - os
  - csv
  - Biopython
- [R programming language v3.2.1](https://cran.r-project.org/). Libraries used:
  - [edgeR](http://bioconductor.org/packages/release/bioc/html/edgeR.html)
  - data.table
  - ggplot2
  - VennDiagram
  - reshape
  - Biostrings
- [PolySearch2](http://polysearch.cs.ualberta.ca/)
