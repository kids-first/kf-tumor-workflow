# Kids First Tumor Only Pipeline

This repo contains tools and workflows for processing of tumor-only samples.
It is currently in beta phase.
Much of the components have been borrowed from the Kids First Somatic Workflow.
It can also be used to process PDX data by first pre-processing reads using the Xenome tool, explained more here in documentation.

<p align="center">
  <img src="docs/kids_first_logo.svg" alt="Kids First repository logo" width="660px" />
</p>
<p align="center">
  <a href="https://github.com/kids-first/kf-tumor-workflow/blob/main/LICENSE"><img src="https://img.shields.io/github/license/kids-first/kf-template-repo.svg?style=for-the-badge"></a>
</p>

## Import info on cloning the git repo
This repo takes advantage of the git submodule feature.
The SNV annotation workflow is maintained in a different repo.
Therefore, in order to get the rest of the code after cloning, you need to run: `git submodule init` and `git submodule update`.
Currently this workflow uses tools from `v4.3.6` of the somatic workflow.
If that is updated, submodule should be as well.
More info on how this worked [here](https://git-scm.com/book/en/v2/Git-Tools-Submodules)

## Main workflow
The wrapper workflow that runs most of the tools is `workflows/kfdrc_tumor_only_dna_wf.cwl`.

## Tools run
SNV
 - Mutect2 from GATK 4.2.2.0
 - Annotation run using same as [Kids First DRC Somatic Variant Annotation Workflow](https://github.com/kids-first/kf-somatic-workflow/blob/v4.3.0/docs/kfdrc_annotation_wf.md)
CNV
 - ControlFreeC v11.6
SV
 - Manta v1.4.0

## Inputs
Most inputs have recommended values that should auto import for both files and parameters
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
 - `bcftools_annot_vcf`: gnomad_3.1.1.vwb_subset.vcf.gz # An export from the Variant Workbench of WGS gnomAD v3.1.1

### Necessary for user to define:
 - `input_tumor_aligned`: <input bam or cram file, indexed>
 - `input_tumor_name`: sample name, should match what is in bam/cram
 - `output_basename`: A file name prefix for all output files
 - `wgs_or_wxs`: Choose whether input is `WGS` or `WXS`. `WXS` works for both whole exome and panel
 - `wgs_calling_interval_list`: if WGS, recommend wgs_canonical_calling_regions.hg38.bed
 - `padded_capture_regions`: if WXS, recommend 100bp padded intervals of capture kit used
 - `i_flag`: for CNV calling, whether to intersect b allele file. Set to `N` for WGS or to skip.
 - `cfree_sex`: for CNV calling, set to XX for female, XY for male
 - `cfree_ploidy`: Array of ploidy possibilities for ControlFreeC to try. Recommend 2-4
 - `unpadded_capture_regions`: if WXS, for CNV, capture regions with NO padding
 - `gatk_filter_name`: ["GNOMAD_AF_HIGH"] # if annotating SNVs, highly recommended
 - `gatk_filter_expression`: ["gnomad_3_1_1_AF > 0.001"] # if annotating SNVs, highly recommended
 - `output_basename`: String value to use as basename for outputs

## Output Files

### Mutect2
 - `mutect2_protected_outputs`: VCF with SNV, MNV, and INDEL variant calls and of pipeline soft FILTER-added values in MAF and  VCF format with annotation, VCF index, and MAF format output
 - `mutect2_public_outputs`: Protected outputs, except MAF and VCF have had entries with soft FILTER values removed
 - `mutect2_bam`: BAM generated will be written as BAM. Really for debugging purposes only
### ControlFreeC CNV
 - `ctrlfreec_pval`: Copy number call with GT (if BAF provided) and p values. Most people want this
 - `ctrlfreec_config`: Config file used to run
 - `ctrlfreec_pngs`: Visualization of CN and BAF
 - `ctrlfreec_bam_ratio`: Calls as log2 ratio
 - `ctrlfreec_bam_seg`: Custom made microarray-style seg file
 - `ctrlfreec_baf`: b allele frequency file
 - `ctrlfreec_info`: Calculated inforamtion, like ploidy, if a range was given 
### Manta SV
 - `manta_pass_vcf`: VCF file with SV calls that PASS
 - `manta_prepass_vcf`: VCF file with all SV calls
 - `annotsv_annotated_calls`: Manta calls annotated with AnnotSV
 - `annotsv_unannotated_calls`: Manta calls not annotated with AnnotSV