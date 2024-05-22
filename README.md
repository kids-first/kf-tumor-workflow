# Kids First DRC Tumor Only Pipeline

This repository contains tools and workflows for processing of tumor-only samples.
It is currently in beta phase.
Much of the components have been borrowed from the Kids First Somatic Workflow.
It can also be used to process PDX data by first pre-processing reads using the Xenome tool, explained more here in documentation.

<p align="center">
  <img src="docs/kids_first_logo.svg" alt="Kids First repository logo" width="660px" />
</p>
<p align="center">
  <a href="https://github.com/kids-first/kf-tumor-workflow/blob/main/LICENSE"><img src="https://img.shields.io/github/license/kids-first/kf-tumor-workflow.svg?style=for-the-badge"></a>
</p>

## Import info on cloning the git repo
This repository takes advantage of the git submodule feature.
The Single Nucleotide Variant annotation workflow is maintained in our [Annotation Tools Repository](https://github.com/kids-first/kf-annotation-tools).
Therefore, in order to get the code for a submodule you can either:
1. Clone the repository recursively with `git clone --recursive`
2. After cloning, run: `git submodule init && git submodule update`
More info on how this worked [here](https://git-scm.com/book/en/v2/Git-Tools-Submodules)

## Main workflow
The wrapper workflow that runs most of the tools is [found here](./workflows/kfdrc_tumor_only_dna_wf.cwl).

## Tools run
### Single Nucleotide Variant (SNV)
 - Mutect2 from GATK 4.2.2.0
 - Annotation using the [Kids First DRC Somatic SNV Annotation Workflow](https://github.com/kids-first/kf-annotation-tools/blob/master/workflows/kfdrc-somatic-snv-annot-workflow.cwl)
### Copy Number Variant (CNV)
 - ControlFREEC v11.6
### Structural Variant (SV)
 - Manta v1.4.0

## Inputs
Most inputs have recommended values that should auto import both files and parameters
### Recommended file/param defaults:
 - `reference_fasta`: Homo_sapiens_assembly38.fasta
 - `reference_fai`: Homo_sapiens_assembly38.fasta.fai
 - `reference_dict`: Homo_sapiens_assembly38.dict
 - `mutect2_af_only_gnomad_vcf`: af-only-gnomad.hg38.vcf.gz
 - `mutect2_exac_common_vcf`: small_exac_common_3.hg38.vcf.gz
 - `gem_mappability_file`: hg38_canonical_150.mappability # will need note on file generation
 - `cfree_chr_len`: hs38_chr.len
 - `b_allele`: dbSNP_v153_ucsc-compatible.converted.vt.decomp.norm.common_snps.vcf.gz # will need note on file generation
 - `hg38_strelka_bed`: hg38_strelka.bed.gz
 - `vep_cache`: homo_sapiens_merged_vep_105_indexed_GRCh38.tar.gz
 - `genomic_hotspots`: tert.bed # bed file with TERT gene promoter region
 - `protein_snv_hotspots`: protein_snv_cancer_hotspots_v2.ENS105_liftover.tsv
 - `protein_indel_hotspots`: protein_indel_cancer_hotspots_v2.ENS105_liftover.tsv
 - `echtvar_anno_zips`: gnom ad.v3.1.1.custom.echtvar.zip
### Necessary for user to define:
 - `input_tumor_aligned`: <input BAM or CRAM file, indexed>
 - `input_tumor_name`: sample name, should match what is in BAM/CRAM
 - `panel_of_normals`: Mutect2 Panel of Normals
 - `output_basename`: A file name prefix for all output files
 - `wgs_or_wxs`: Choose whether input is Whole Genome Sequencing (WGS) or Whole Exome Sequencing or Panel (WXS)
 - `wgs_calling_interval_list`: if WGS, recommend wgs_canonical_calling_regions.hg38.bed
# TODO: ADD BLACKLIST REGIONS?
 - `padded_capture_regions`: if WXS, recommend 100bp padded intervals of capture kit used
 - `i_flag`: for CNV calling, whether to intersect b allele file. Set to `N` for WGS or to skip.
 - `cfree_sex`: for CNV calling, set to XX for female, XY for male
 - `cfree_ploidy`: Array of ploidy possibilities for ControlFREEC to try. Recommend 2-4
 - `unpadded_capture_regions`: if WXS, for CNV, capture regions with NO padding
 - `gatk_filter_name`: `["GNOMAD_AF_HIGH"]` Highly recommended for SNV annotations
# TODO: ADD AF, AD, lower gnomadAF - `gatk_filter_expression`: `["gnomad_3_1_1_AF != '.' && gnomad_3_1_1_AF > 0.001 && gnomad_3_1_1_FILTER == 'PASS'"]` Highly recommended for SNV annotations
 - `output_basename`: String value to use as basename for outputs

## Output Files
### Mutect2
 - `mutect2_protected_outputs`: VCF with SNV, MNV, and INDEL variant calls and of pipeline soft FILTER-added values in MAF and  VCF format with annotation, VCF index, and MAF format output
 - `mutect2_public_outputs`: Protected outputs, except MAF and VCF have had entries with soft FILTER values removed
 - `mutect2_bam`: BAM generated will be written as BAM. Really for debugging purposes only
### ControlFREEC CNV
 - `ctrlfreec_pval`: Copy number call with GT (if BAF provided) and p values. Most people want this
 - `ctrlfreec_config`: Config file used to run
 - `ctrlfreec_pngs`: Visualization of CN and BAF
 - `ctrlfreec_bam_ratio`: Calls as log2 ratio
 - `ctrlfreec_bam_seg`: Custom made microarray-style SEG file
 - `ctrlfreec_baf`: b allele frequency file
 - `ctrlfreec_info`: Calculated information, like ploidy, if a range was given 
### Manta SV
 - `manta_pass_vcf`: VCF file with SV calls that PASS
 - `manta_prepass_vcf`: VCF file with all SV calls
 - `annotsv_annotated_calls`: Manta calls annotated with AnnotSV
 - `annotsv_unannotated_calls`: Manta calls not annotated with AnnotSV
