# RNAseq-fusion-nf

## Fusion-genes discovery from RNAseq data using STAR-Fusion

<!---[![CircleCI](https://circleci.com/gh/IARCbioinfo/template-nf.svg?style=svg)](https://circleci.com/gh/IARCbioinfo/template-nf)--->
<!---[![Docker Hub](https://img.shields.io/badge/docker-ready-blue.svg)](https://hub.docker.com/r/iarcbioinfo/template-nf/)--->
<!---[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/1404)--->
<!---[![DOI](https://zenodo.org/badge/94193130.svg)](https://zenodo.org/badge/latestdoi/94193130)--->

![Workflow representation](RNAseq-fusion-nf.png)

## Description
Performs Fusion-gene discovery from junction reads identified by STAR during alignment. See citation of STAR-Fusion below.

## Dependencies

1. This pipeline is based on [nextflow](https://www.nextflow.io). As we have several nextflow pipelines, we have centralized the common information in the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository. Please read it carefully as it contains essential information for the installation, basic usage and configuration of nextflow and our pipelines.
2. External software:
- [STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion/wiki); note that you can use the bioconda receipe 

You can avoid installing all the external software by only installing Docker. See the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository for more information.

In addition, STAR-Fusion requires a [CTAT bundle](https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/) with reference genome and annotations.

## Input
  | Type      | Description     |
  |-----------|---------------|
  | input_folder    | Folder containing fastq files and STAR junction files |

  Specify the test files location

## Parameters

  * #### Mandatory
| Name      | Example value | Description     |
|-----------|---------------|-----------------|
| --CTAT_folder   | CTAT | Folder with STAR-Fusion bundle (CTAT) |

  * #### Optional
| Name      | Default value | Description     |
|-----------|---------------|-----------------|
| --output_folder   |         . | Output folder |
| --fastq_ext       |   fq.gz    |            Extension of fastq files |
| --suffix1         |  _1 |   Suffix of 1st element of fastq files pair |
| --suffix2         |  _2 |     Suffix of 2nd element of fastq files pair |
| --junction_suffix | Chimeric.SJ.out.junction |       Suffix of STAR chimeric junction files |
| --cpu             |  2 |         Number of cpu used by bwa mem and sambamba |
| --mem             |  2 |   Size of memory used for mapping (in GB)|


  * #### Flags

Flags are special parameters without value.

| Name      | Description     |
|-----------|-----------------|
| --help    | Display help |


## Usage
  ```
  nextflow run iarcbioinfo/RNAseq-fusion-nf --input_folder input --CTAT_folder CTAT --output_folder output
  ```

## Output
  | Type      | Description     |
  |-----------|---------------|
  | output   | Folder with fusion genes file |


<!--- ## Detailed description (optional section) --->

## Directed Acyclic Graph
[![DAG](dag.png)](http://htmlpreview.github.io/?https://github.com/IARCbioinfo/RNAseq-fusion-nf/blob/master/dag.html)

## Contributions

  | Name      | Email | Description     |
  |-----------|---------------|-----------------|
  | Nicolas Alcala    |  alcalan@fellows.iarc.fr | Developer to contact for support |

## References
Haas, B., Dobin, A., Stransky, N., Li, B., Yang, X., Tickle, T., ... & Sun, J. (2017). STAR-Fusion: fast and accurate fusion transcript detection from RNA-Seq. BioRxiv, 120295.



<!--- ## FAQ (optional)--->
