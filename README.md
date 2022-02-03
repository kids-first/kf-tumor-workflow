<p align="center">
  <img src="docs/kids_first_logo.svg" alt="Kids First repository logo" width="660px" />
</p>
<p align="center">
  <a href="https://github.com/kids-first/kf-tumor-workflow/blob/main/LICENSE"><img src="https://img.shields.io/github/license/kids-first/kf-template-repo.svg?style=for-the-badge"></a>
</p>

# Kids First Tumor Only Pipeline

This repo contains tools and workflows for processing of tumor-only samples.
It is currently in beta phase.
Much of the components have been borrowed from the Kids First Somatic Workflow.
It can also be used to process PDX data by first pre-processing reads using the Xenome tool, explained more here in documentation.

## Main workflow
The wrapper workflow that runs most of the tools is `workflows/kfdrc_tumor_only_dna_wf.cwl`.

### Tools run
SNV
 - Mutect2 from GATK 4.2.2.0
CNV
 - ControlFreeC v11.6
SV
 - Manta v1.4.0

### Inputs
Most inputs have recommended values that should auto import for both files and parameters
#### Recommended file/param defaults:
```
reference_fasta: Homo_sapiens_assembly38.fasta
reference_fai: Homo_sapiens_assembly38.fasta.fai
reference_dict: Homo_sapiens_assembly38.dict
mutect2_af_only_gnomad_vcf: af-only-gnomad.hg38.vcf.gz
mutect2_exac_common_vcf: small_exac_common_3.hg38.vcf.gz
gem_mappability_file: hg38_canonical_150.mappability # will need note on file generation
cfree_chr_len: hs38_chr.len
b_allele: dbSNP_v153_ucsc-compatible.converted.vt.decomp.norm.common_snps.vcf.gz # will need note on file generation
hg38_strelka_bed: hg38_strelka.bed.gz
# If annotating SNVs, recommended
vep_cache: homo_sapiens_vep_93_GRCh38.tar.gz
genomic_hotspots: tert.bed # bed file with TERT gene promoter region
protein_snv_hotspots: protein_snv_cancer_hotspots_v2.tsv
protein_indel_hotspots: protein_indel_cancer_hotspots_v2.tsv
bcftools_annot_vcf: af-only-gnomad.hg38.vcf.gz # yes, same as mutect2 input
```
#### Necessary for user to define:
```
input_tumor_aligned: <input bam or cram file, indexed>
input_tumor_name: sample name, should match what is in bam/cram
output_basename: A file name prefix for all output files
wgs_or_wxs: Choose whether input is `WGS` or `WXS`. `WXS` works for both whole exome and panel
wgs_calling_interval_list: if WGS, recommend wgs_canonical_calling_regions.hg38.bed
padded_capture_regions: if WXS, recommend 100bp padded intervals of capture kit used
i_flag: for CNV calling, whether to intersect b allele file. Set to `N` for WGS or to skip.
cfree_sex: for CNV calling, set to XX for female, XY for male
unpadded_capture_regions: if WXS, for CNV, capture regions with NO padding
gatk_filter_name: ["GNOMAD_AF_HIGH"] # if annotating SNVs, highly recommended
gatk_filter_expression: ["AF > 0.001"] # if annotating SNVs, highly recommended
```

#### Comprehensive:
```yaml
  # Shared
  reference_fasta: {type: 'File', doc: "Reference fasta", "sbg:suggestedValue": {
      class: File, path: 60639014357c3a53540ca7a3, name: Homo_sapiens_assembly38.fasta},
    "sbg:fileTypes": "FASTA, FA"}
  reference_fai: {type: 'File?', doc: "FAI index for reference_fasta", "sbg:suggestedValue": {
      class: File, path: 60639016357c3a53540ca7af, name: Homo_sapiens_assembly38.fasta.fai},
    "sbg:fileTypes": "FAI"}
  reference_dict: {type: 'File?', doc: "DICT index for reference_fasta", "sbg:suggestedValue": {
      class: File, path: 60639019357c3a53540ca7e7, name: Homo_sapiens_assembly38.dict},
    "sbg:fileTypes": "DICT"}
  input_tumor_aligned: {type: 'File', secondaryFiles: [{pattern: ".bai", required: false},
      {pattern: "^.bai", required: false}, {pattern: ".crai", required: false}, {
        pattern: "^.crai", required: false}], doc: "BAM/SAM/CRAM file containing tumor\
      \ reads", "sbg:fileTypes": "BAM, CRAM, SAM"}
  input_tumor_name: {type: 'string', doc: "BAM sample name of tumor. May be URL-encoded\
      \ as output by GetSampleName with -encode argument."}
  output_basename: {type: 'string', doc: "String value to use as basename for outputs"}

  # GATK Files
  panel_of_normals: {type: 'File?', secondaryFiles: ['.tbi'], doc: "VCF file (and\
      \ index) of sites observed in normal. A panel of normals can be a useful (optional)\
      \ input to help filter out commonly seen sequencing noise that may appear as\
      \ low allele-fraction somatic variants.", "sbg:fileTypes": "VCF, VCF.GZ"}
  mutect2_af_only_gnomad_vcf: {type: 'File', secondaryFiles: ['.tbi'], doc: "Population\
      \ vcf (and index) of germline sequencing containing allele fractions in VCF\
      \ format.", "sbg:suggestedValue": {class: File, path: 5f50018fe4b054958bc8d2e3,
      name: af-only-gnomad.hg38.vcf.gz, secondaryFiles: [{class: File, path: 5f50018fe4b054958bc8d2e5,
          name: af-only-gnomad.hg38.vcf.gz.tbi}]}, "sbg:fileTypes": "VCF, VCF.GZ"}
  mutect2_alleles_vcf: {type: 'File?', secondaryFiles: ['.tbi'], doc: "VCF file (and\
      \ index) containing set of alleles for which to force genotyping regardless\
      \ of evidence", "sbg:fileTypes": "VCF, VCF.GZ"}
  mutect2_exac_common_vcf: {type: 'File?', secondaryFiles: ['.tbi'], doc: "Exac Common\
      \ VCF (and index) used for calculating contamination values used in filtering\
      \ the Mutect VCF. If do not wish to perfom this filtering, remove this input.",
    "sbg:suggestedValue": {class: File, path: 5f500135e4b0370371c051ad, name: small_exac_common_3.hg38.vcf.gz,
      secondaryFiles: [{class: File, path: 5f500135e4b0370371c051af, name: small_exac_common_3.hg38.vcf.gz.tbi}]},
    "sbg:fileTypes": "VCF, VCF.GZ"}
  bwa_mem_index_image: {type: 'File?', doc: "BWA-mem index image used for FilterAlignmentArtifacts.\
      \ WARNING!!! Experimental Software!!! Do not provide unless you are absolutely\
      \ certain you want to run FilterAlignmentArtifacts.", "sbg:fileTypes": "IMG"}

  # GATK Arguments
  tool_name: {type: 'string?', doc: "Tool name for this workflow", default: "mutect2"}
  make_bamout: {type: 'boolean?', default: false, doc: "Should Mutect2 create a BAM\
      \ output? Turning on this mode may result in serious performance cost for the\
      \ caller."}
  run_orientation_bias_mixture_model_filter: {type: 'boolean?', default: true, doc: "Should\
      \ Orientation Bias Mixture Model Filter be applied to Mutect2 outputs?"}
  mutect2_extra_args: {type: 'string?', doc: "Any additional arguments for Mutect2.\
      \ See GATK documentation for all available options."}
  filtermutectcalls_extra_args: {type: 'string?', doc: "Any additional arguments for\
      \ FilterMutectCalls. See GATK documentation for all available options."}
  select_vars_mode: {type: ['null', {type: enum, name: select_vars_mode, symbols: [
          "gatk", "grep"]}], default: "gatk", doc: "Choose 'gatk' for SelectVariants\
      \ tool, or 'grep' for grep expression"}
  scatter_count: {type: 'int?', default: 50, doc: "The number of files into which\
      \ to scatter the resulting list by locus; in some situations, fewer intervals\
      \ may be emitted"}

  # WGS or WXS
  wgs_or_wxs: {type: {type: enum, name: wgs_or_wxs, symbols: ["WGS", "WXS"]}, doc: "Select\
      \ input type WGS or WXS"}
  wgs_calling_interval_list: {type: 'File?', doc: "GATK intervals list-style, or bed\
      \ file.  Recommend canocical chromosomes with N regions removed", "sbg:suggestedValue": {
      class: File, path: 5f500135e4b0370371c051b6, name: wgs_canonical_calling_regions.hg38.bed},
    "sbg:fileTypes": "BED, INTERVALS, INTERVAL_LIST, LIST"}
  padded_capture_regions: {type: 'File?', doc: "Recommend 100bp pad, for somatic variant.\
      \ Used as calling intervals for exomes.", "sbg:fileTypes": "BED, INTERVALS,\
      \ INTERVAL_LIST, LIST"}

  # ControlFreeC CNV
  mate_copynumber_file_sample: {type: 'File?', doc: "Tumor cpn file from previous run. If used, will override bam use"}
  gem_mappability_file: {type: 'File?', doc: "GEM mappability file to make read count adjustments with"}
  min_subclone_presence: {type: 'float?', doc: "Use if you want to detect subclones. Recommend 0.2 for WGS, 0.3 for WXS"}
  cfree_chr_len: { type: File, doc: "file with chromosome lengths" }
  cfree_ploidy: { type: 'int[]', doc: "Array of ploidy possibilities for ControlFreeC to try" }
  cfree_mate_orientation_sample: { type: ['null', { type: enum, name: mate_orientation_sample, symbols: ["0", "FR", "RF", "FF"] }], default: "FR", doc: "0 (for single ends), RF (Illumina mate-pairs), FR (Illumina paired-ends), FF (SOLiD mate-pairs)" }

  # Optional with Multiple Defaults (handled in choose_defaults)
  i_flag: { type: 'string?', doc: "Flag to intersect germline calls on padded regions. Use N if you want to skip this or have a WGS run" }

  # Optional
  b_allele: { type: 'File?', doc: "germline calls, needed for BAF.  GATK HC VQSR input recommended. Tool will prefilter for germline and pass if expression given" }
  b_allele_index: { type: 'File?', doc: "Tabix index for b_allele" }
  cfree_coeff_var: { type: 'float?', default: 0.05, doc: "Coefficient of variation to set window size. Default 0.05 recommended" }
  cfree_sex: { type: ['null', { type: enum, name: sex, symbols: ["XX", "XY"] }], doc: "If known, XX for female, XY for male" }

  # Manta SV
  hg38_strelka_bed: { type: File, doc: "Bgzipped interval bed file. Recommned padding 100bp for WXS; Recommend canonical chromosomes for WGS" }
  hg38_strelka_tbi: { type: 'File?', doc: "Tabix index for hg38_strelka_bed" }

  # WXS only Fields
  unpadded_capture_regions: { type: 'File?', doc: "Capture regions with NO padding for cnv calling" }


  # Annotation Arguments
  run_annotation: {type: 'boolean?', default: true, doc: "Should annotation be run\
      \ on the final VCF?"}
  vep_cache: {type: 'File?', doc: "tar gzipped cache from ensembl/local converted\
      \ cache", "sbg:suggestedValue": {class: File, path: 607713829360f10e3982a425,
      name: homo_sapiens_vep_93_GRCh38.tar.gz}, "sbg:fileTypes": "TAR,\
      \ TAR.GZ, TGZ"}
  vep_ref_build: {type: 'string?', default: "GRCh38", doc: "Genome ref build used,\
      \ should line up with cache"}
  genomic_hotspots: {type: 'File[]?', doc: "Tab-delimited BED formatted file(s) containing\
      \ hg38 genomic positions corresponding to hotspots", "sbg:suggestedValue": [{
        class: File, path: 607713829360f10e3982a423, name: tert.bed}], "sbg:fileTypes": "BED,\
      \ TSV"}
  protein_snv_hotspots: {type: 'File[]?', doc: "Column-name-containing, tab-delimited\
      \ file(s) containing protein names and amino acid positions corresponding to\
      \ hotspots", "sbg:suggestedValue": [{class: File, path: 607713829360f10e3982a426,
        name: protein_snv_cancer_hotspots_v2.tsv}], "sbg:fileTypes": "BED,\
      \ TSV"}
  protein_indel_hotspots: {type: 'File[]?', doc: "Column-name-containing, tab-delimited\
      \ file(s) containing protein names and amino acid position ranges corresponding\
      \ to hotspots", "sbg:suggestedValue": [{class: File, path: 607713829360f10e3982a424,
        name: protein_indel_cancer_hotspots_v2.tsv}], "sbg:fileTypes": "BED,\
      \ TSV"}
  retain_info: {type: 'string?', doc: "csv string with INFO fields that you want to\
      \ keep", default: "MBQ,TLOD,HotSpotAllele"}
  retain_fmt: {type: 'string?', doc: "csv string with FORMAT fields that you want\
      \ to keep"}
  add_common_fields: {type: 'boolean?', doc: "Set to true if input is a strelka2 vcf\
      \ that hasn't had common fields added", default: false}
  bcftools_annot_columns: {type: 'string?', doc: "csv string of columns from annotation\
      \ to port into the input vcf, i.e INFO/AF", default: "INFO/AF"}
  bcftools_annot_vcf: {type: 'File?', secondaryFiles: [.tbi], doc: "bgzipped annotation\
      \ vcf file", "sbg:suggestedValue": {class: File, path: 5f50018fe4b054958bc8d2e3,
      name: af-only-gnomad.hg38.vcf.gz, secondaryFiles: [{class: File, path: 5f50018fe4b054958bc8d2e5,
          name: af-only-gnomad.hg38.vcf.gz.tbi}]}, "sbg:fileTypes": "VCF, VCF.GZ"}
  bcftools_public_filter: {type: 'string?', doc: "Will hard filter final result to\
      \ create a public version", default: FILTER="PASS"|INFO/HotSpotAllele=1}
  gatk_filter_name: {type: 'string[]?', doc: "Array of names for each filter tag to\
      \ add, recommend: [\"NORM_DP_LOW\", \"GNOMAD_AF_HIGH\"]"}
  gatk_filter_expression: {type: 'string[]?', doc: "Array of filter expressions to\
      \ establish criteria to tag variants with. See https://gatk.broadinstitute.org/hc/en-us/articles/360036730071-VariantFiltration,\
      \ recommend: \"vc.getGenotype('\" + inputs.input_normal_name + \"').getDP()\
      \ <= 7\"), \"AF > 0.001\"]"}
  disable_hotspot_annotation: {type: 'boolean?', doc: "Disable Hotspot Annotation\
      \ and skip this task.", default: false}
  maf_center: {type: 'string?', doc: "Sequencing center of variant called", default: "."}

  # Resource Control
  mutect_cores: {type: 'int?', doc: "CPUs to allocate to GATK Mutect2"}
  mutect_memory: {type: 'int?', doc: "GB of memory to allocate to GATK Mutect2 (hard-capped)"}
  getpileup_memory: {type: 'int?', doc: "GB of memory to allocate to GATK GetPileupSummaries\
      \ (hard-capped)"}
  learnorientation_memory: {type: 'int?', doc: "GB of memory to allocate to GATK LearnReadOrientationModel\
      \ (hard-capped)"}
  filtermutectcalls_memory: {type: 'int?', doc: "GB of memory to allocate to GATK\
      \ FilterMutectCalls (hard-capped)"}
  filteralignmentartifacts_cores: {type: 'int?', doc: "CPUs to allocate to GATK FilterAlignmentArtifacts"}
  filteralignmentartifacts_memory: {type: 'int?', doc: "GB of memory to allocate to\
      \ GATK FilterAlignmentArtifacts (hard-capped)"}
  cfree_threads: { type: 'int?', default: 16, doc: "For ControlFreeC. Recommend 16 max, as I/O gets saturated after that losing any advantage" }
  manta_memory: {type: 'int?', doc: "GB of memory to allocate to Manta; defaults to 10 (soft-capped)"}
  manta_cores: {type: 'int?', doc: "Number of cores to allocate to Manta; defaults to 18"}
  ```

### Output Files

```yaml
  # Mutect2
  mutect2_prepass_vcf: {type: 'File', outputSource: run_mutect2/mutect2_filtered_vcf,
    doc: "VCF with SNV, MNV, and INDEL variant calls."}
  mutect2_protected_outputs: {type: 'File[]?', outputSource: run_mutect2/mutect2_protected_outputs,
    doc: "Array of files containing MAF format of PASS hits, PASS VCF with annotation\
      \ pipeline soft FILTER-added values, and VCF index"}
  mutect2_public_outputs: {type: 'File[]?', outputSource: run_mutect2/mutect2_public_outputs,
    doc: "Protected outputs, except MAF and VCF have had entries with soft FILTER\
      \ values removed"}
  mutect2_bam: {type: 'File?', outputSource: run_mutect2/mutect2_bam, doc: "BAM generated\
      \ by Mutect2, if desired. The assembled haplotypes and locally realigned reads\
      \ will be written as BAM. Really for debugging purposes only."}
  # ControlFreeC CNV
  ctrlfreec_pval: { type: File, outputSource: run_controlfreec/ctrlfreec_pval, doc: 'Copy number call with GT (if BAF provided) and p values. Most people want this' }
  ctrlfreec_config: { type: File, outputSource: run_controlfreec/ctrlfreec_config, doc: 'Config file used to run' }
  ctrlfreec_pngs: { type: 'File[]', outputSource: run_controlfreec/ctrlfreec_pngs, doc: 'Visualization of CN and BAF' }
  ctrlfreec_bam_ratio: { type: File, outputSource: run_controlfreec/ctrlfreec_bam_ratio, doc: 'Calls as log2 ratio' }
  ctrlfreec_bam_seg: { type: File, outputSource: run_controlfreec/ctrlfreec_bam_seg, doc: 'Custom made microarray-style seg file' }
  ctrlfreec_baf: { type: File, outputSource: run_controlfreec/ctrlfreec_baf, doc: 'b allele frequency file' }
  ctrlfreec_info: { type: File, outputSource: run_controlfreec/ctrlfreec_info, doc: 'Calculated inforamtion, like ploidy, if a range was given' }
  # Manta SV
  manta_pass_vcf: { type: File, outputSource: run_manta/manta_pass_vcf, doc: 'VCF file with SV calls that PASS' }
  manta_prepass_vcf: { type: File, outputSource: run_manta/manta_prepass_vcf, doc: 'VCF file with all SV calls' }
```