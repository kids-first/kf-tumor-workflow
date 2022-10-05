cwlVersion: v1.2
class: Workflow
id: kfdrc_mutect2_production_wf
label: KFDRC Mutect2 Production Workflow
doc: |
  # Mutect2 Subworkflow

  ## Minimal Run
  - Mutect
  - MergeMutectStats
  - MergeVcfs
  - FilterMutectCalls

  ## Mutect Filter Support (creates contamination and segmentation tables for FilterMutectCalls)
  __Requires `exac_common_vcf` input__

  - GetPileupSummaries
  - GatherPileupSummaries
  - CalculateContamination

  ## Orientation Bias Mixture Model Filter (creates artifact prior tables for FilterMutectCalls)
  __Requires `run_orientation_bias_mixture_model_filter` set to `true`__

  - Additional f1r2 outputs from Mutect
  - LearnReadOrientationModel

  ## Mutect BAM Creation
  __Requires `make_bamout` set to `true`__

  - Additional BAM outputs from Mutect
  - GatherBamFiles
  - SortSam
  - BuildBamIndex

  The assembled haplotypes and locally realigned reads will be written as BAM to
  this file if requested. Really for debugging purposes only. Note that the output
  here does not include uninformative reads so that not every input read is emitted
  to the bam. Turning on this mode may result in serious performance cost for the
  caller. It's really only appropriate to use in specific areas where you want to
  better understand why the caller is making specific calls.

  ## Filter Alignment Artifacts
  __Requires `bwa_mem_index_image` input__

  - FilterAlignmentArtifacts

  ## Annotation
  __Requires `run_annotation` set to `true`__

  - SelectVariants (PASS)
  - vt normalize VCF
  - bcftools strip annotations
  - add strelka2 standard fields (Requires add_common_fields set to true; Not used for Mutect VCFs)
  - VEP annotate
  - bcftools gnomad annotate
  - VariantFiltration
  - hotspot annotation
  - vcf2maf public
  - bcftools hardfilter
  - vcf2maf protected
  - rename protected files
  - rename public files


  # Mutect2 Production Workflow

  ## Default Functionality
  - Minimal Run
  - Mutect Filter Support
  - Orientation Bias Mixture Model Filter
  - Annotation

  Any segments not used by default can still be used activated bt the user.
  We do not create BAMs for storage cost reasons.
  FilterAlignmentArtifacts is still an experimental tool and not recommended for production pipelines! The current tool still has several issues as of GATK 4.2.0.0.

  ## Somatic vs Tumor-Only
  The production workflow is capable of running in both somatic and tumor-only modes.

  To run the workflow in somatic mode, the following inputs are necessary:

  - `input_tumor_aligned`
  - `input_normal_aligned`

  To run the workflow in tumor-only mode, the following inputs are necessary:

  - `input_tumor_aligned`
  - `panel_of_normals`

  ## WGS vs WXS
  The production workflow is capable of running both WGS and WXS inputs.

  ### Running the tool in WGS mode:

  __Requires `WGS` for `wgs_or_wxs`__

  - Sets `--BREAK_BANDS_AT_MULTIPLES_OF` to `80000000` in IntervalListTools

  ### Running the tool in WXS mode:

  __Requires `WXS` for `wgs_or_wxs`__

  - Sets `--BREAK_BANDS_AT_MULTIPLES_OF` to `0` in IntervalListTools
  - Additional `--disable-adaptive-pruning` flag set in Mutect
requirements:
- class: ScatterFeatureRequirement
- class: MultipleInputFeatureRequirement
- class: SubworkflowFeatureRequirement
- class: InlineJavascriptRequirement

inputs:
  # GATK Files
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
  input_normal_aligned: {type: 'File?', secondaryFiles: [{pattern: ".bai", required: false},
      {pattern: "^.bai", required: false}, {pattern: ".crai", required: false}, {
        pattern: "^.crai", required: false}], doc: "BAM/SAM/CRAM file containing normal\
      \ reads", "sbg:fileTypes": "BAM, CRAM, SAM"}
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
  input_tumor_name: {type: 'string', doc: "BAM sample name of tumor. May be URL-encoded\
      \ as output by GetSampleName with -encode argument."}
  input_normal_name: {type: 'string?', doc: "BAM sample name of normal(s), if any.\
      \ May be URL-encoded as output by GetSampleName with -encode argument."}
  tool_name: {type: 'string?', doc: "Tool name for this workflow", default: "mutect2"}
  output_basename: {type: 'string', doc: "String value to use as basename for outputs"}
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

  # VEP params
  vep_cache: {type: 'File', doc: "tar gzipped cache from ensembl/local converted cache",  "sbg:suggestedValue": {class: File, path: 6332f8e47535110eb79c794f,
      name: homo_sapiens_merged_vep_105_indexed_GRCh38.tar.gz}}
  vep_ram: {type: 'int?', default: 32, doc: "In GB, may need to increase this value depending on the size/complexity of input"}
  vep_cores: {type: 'int?', default: 16, doc: "Number of cores to use. May need to increase for really large inputs"}
  vep_buffer_size: {type: 'int?', default: 1000, doc: "Increase or decrease to balance speed and memory usage"}
  dbnsfp: { type: 'File?', secondaryFiles: [.tbi,^.readme.txt], doc: "VEP-formatted plugin file, index, and readme file containing dbNSFP annotations" }
  dbnsfp_fields: { type: 'string?', doc: "csv string with desired fields to annotate. Use ALL to grab all"}
  merged: { type: 'boolean?', doc: "Set to true if merged cache used", default: true }
  cadd_indels: { type: 'File?', secondaryFiles: [.tbi], doc: "VEP-formatted plugin file and index containing CADD indel annotations" }
  cadd_snvs: { type: 'File?', secondaryFiles: [.tbi], doc: "VEP-formatted plugin file and index containing CADD SNV annotations" }
  run_cache_existing: { type: 'boolean?', doc: "Run the check_existing flag for cache" }
  run_cache_af: { type: 'boolean?', doc: "Run the allele frequency flags for cache" }

  # annotation vars
  genomic_hotspots: { type: 'File[]?', doc: "Tab-delimited BED formatted file(s) containing hg38 genomic positions corresponding to hotspots", "sbg:suggestedValue": [{class: File, path: 607713829360f10e3982a423, name: tert.bed}] }
  protein_snv_hotspots: { type: 'File[]?', doc: "Column-name-containing, tab-delimited file(s) containing protein names and amino acid positions corresponding to hotspots", "sbg:suggestedValue": [{class: File, path: 607713829360f10e3982a426, name: protein_snv_cancer_hotspots_v2.tsv}] }
  protein_indel_hotspots: { type: 'File[]?', doc: "Column-name-containing, tab-delimited file(s) containing protein names and amino acid position ranges corresponding to hotspots", "sbg:suggestedValue": [{class: File, path: 607713829360f10e3982a424, name: protein_indel_cancer_hotspots_v2.tsv}] }
  retain_info: {type: 'string?', doc: "csv string with INFO fields that you want to keep", default: "gnomad_3_1_1_AC,gnomad_3_1_1_AN,gnomad_3_1_1_AF,gnomad_3_1_1_nhomalt,gnomad_3_1_1_AC_popmax,gnomad_3_1_1_AN_popmax,gnomad_3_1_1_AF_popmax,gnomad_3_1_1_nhomalt_popmax,gnomad_3_1_1_AC_controls_and_biobanks,gnomad_3_1_1_AN_controls_and_biobanks,gnomad_3_1_1_AF_controls_and_biobanks,gnomad_3_1_1_AF_non_cancer,gnomad_3_1_1_primate_ai_score,gnomad_3_1_1_splice_ai_consequence,gnomad_3_1_1_AC,gnomad_3_1_1_AN,gnomad_3_1_1_AF,gnomad_3_1_1_nhomalt,gnomad_3_1_1_AC_popmax,gnomad_3_1_1_AN_popmax,gnomad_3_1_1_AF_popmax,gnomad_3_1_1_nhomalt_popmax,gnomad_3_1_1_AC_controls_and_biobanks,gnomad_3_1_1_AN_controls_and_biobanks,gnomad_3_1_1_AF_controls_and_biobanks,gnomad_3_1_1_AF_non_cancer,gnomad_3_1_1_primate_ai_score,gnomad_3_1_1_splice_ai_consequence,MBQ,TLOD,HotSpotAllele"}
  retain_fmt: {type: 'string?', doc: "csv string with FORMAT fields that you want to keep"}
  retain_ann: { type: 'string?', doc: "csv string of annotations (within the VEP CSQ/ANN) to retain as extra columns in MAF", default: "HGVSg" }
  bcftools_annot_columns: {type: 'string?', doc: "csv string of columns from annotation to port into the input vcf, i.e INFO/AF", default: "INFO/gnomad_3_1_1_AC:=INFO/AC,INFO/gnomad_3_1_1_AN:=INFO/AN,INFO/gnomad_3_1_1_AF:=INFO/AF,INFO/gnomad_3_1_1_nhomalt:=INFO/nhomalt,INFO/gnomad_3_1_1_AC_popmax:=INFO/AC_popmax,INFO/gnomad_3_1_1_AN_popmax:=INFO/AN_popmax,INFO/gnomad_3_1_1_AF_popmax:=INFO/AF_popmax,INFO/gnomad_3_1_1_nhomalt_popmax:=INFO/nhomalt_popmax,INFO/gnomad_3_1_1_AC_controls_and_biobanks:=INFO/AC_controls_and_biobanks,INFO/gnomad_3_1_1_AN_controls_and_biobanks:=INFO/AN_controls_and_biobanks,INFO/gnomad_3_1_1_AF_controls_and_biobanks:=INFO/AF_controls_and_biobanks,INFO/gnomad_3_1_1_AF_non_cancer:=INFO/AF_non_cancer,INFO/gnomad_3_1_1_primate_ai_score:=INFO/primate_ai_score,INFO/gnomad_3_1_1_splice_ai_consequence:=INFO/splice_ai_consequence"}
  bcftools_strip_columns: {type: 'string?', doc: "csv string of columns to strip if needed to avoid conflict, i.e INFO/AF"}
  bcftools_annot_vcf: {type: 'File', doc: "bgzipped annotation vcf file", "sbg:suggestedValue": {
      class: File, path: 6324ef5ad01163633daa00d8, name: gnomad_3.1.1.vwb_subset.vcf.gz}}
  bcftools_annot_vcf_index: {type: 'File?', doc: "index of bcftools_annot_vcf", "sbg:suggestedValue": {
      class: File, path: 6324ef5ad01163633daa00d7, name: gnomad_3.1.1.vwb_subset.vcf.gz.tbi}}
  bcftools_public_filter: {type: 'string?', doc: "Will hard filter final result to create a public version", default: FILTER="PASS"|INFO/HotSpotAllele=1}
  gatk_filter_name: {type: 'string[]', doc: "Array of names for each filter tag to add, recommend: [\"NORM_DP_LOW\", \"GNOMAD_AF_HIGH\"]"}
  gatk_filter_expression: {type: 'string[]', doc: "Array of filter expressions to establish criteria to tag variants with. See https://gatk.broadinstitute.org/hc/en-us/articles/360036730071-VariantFiltration, recommend: \"vc.getGenotype('\" + inputs.input_normal_name + \"').getDP() <= 7\"), \"gnomad_3_1_1_AF > 0.001\"]"}
  disable_hotspot_annotation: { type: 'boolean?', doc: "Disable Hotspot Annotation and skip this task.", default: false }
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

outputs:
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

steps:
  choose_defaults:
    run: ../tools/mode_defaults.cwl
    in:
      input_mode: wgs_or_wxs
    out: [out_exome_flag, out_cnvkit_wgs_mode, out_i_flag, out_lancet_padding, out_lancet_window,
      out_vardict_padding]

  prepare_reference:
    run: ../subworkflows/prepare_reference.cwl
    in:
      input_fasta: reference_fasta
      input_fai: reference_fai
      input_dict: reference_dict
    out: [indexed_fasta, reference_dict]
  
  index_bcftools_annot_vcf:
    run: ../tools/tabix_index.cwl
    in:
      input_file: bcftools_annot_vcf
      input_index: bcftools_annot_vcf_index
    out: [output]

  select_interval_list:
    run: ../tools/mode_selector.cwl
    in:
      input_mode: wgs_or_wxs
      wgs_input: wgs_calling_interval_list
      wxs_input: padded_capture_regions
    out: [output]

  gatk_intervallisttools:
    run: ../tools/gatk_intervallisttool.cwl
    in:
      interval_list: select_interval_list/output
      reference_dict: prepare_reference/reference_dict
      exome_flag: choose_defaults/out_exome_flag
      scatter_ct: scatter_count
      bands:
        valueFrom: ${return 80000000}
    out: [output]

  run_mutect2:
    hints:
    - class: 'sbg:AWSInstanceType'
      value: c5.9xlarge
    run: ../subworkflows/kfdrc_mutect2_sub_wf.cwl
    in:
      indexed_reference_fasta: prepare_reference/indexed_fasta
      reference_dict: prepare_reference/reference_dict
      bed_invtl_split: gatk_intervallisttools/output
      af_only_gnomad_vcf: mutect2_af_only_gnomad_vcf
      alleles: mutect2_alleles_vcf
      exac_common_vcf: mutect2_exac_common_vcf
      input_tumor_aligned: input_tumor_aligned
      input_tumor_name: input_tumor_name
      input_normal_aligned: input_normal_aligned
      input_normal_name: input_normal_name
      panel_of_normals: panel_of_normals
      disable_adaptive_pruning:
        source: wgs_or_wxs
        valueFrom: $(self == 'WXS')
      mutect2_extra_args: mutect2_extra_args
      filtermutectcalls_extra_args: filtermutectcalls_extra_args
      select_vars_mode: select_vars_mode
      retain_info: retain_info
      retain_fmt: retain_fmt
      retain_ann: retain_ann
      bcftools_annot_columns: bcftools_annot_columns
      bcftools_strip_columns: bcftools_strip_columns
      bcftools_annot_vcf: index_bcftools_annot_vcf/output
      bcftools_public_filter: bcftools_public_filter
      dbnsfp: dbnsfp
      dbnsfp_fields: dbnsfp_fields
      merged: merged
      cadd_indels: cadd_indels
      cadd_snvs: cadd_snvs
      run_cache_af: run_cache_af
      run_cache_existing: run_cache_existing
      gatk_filter_name: gatk_filter_name
      gatk_filter_expression: gatk_filter_expression
      vep_cache: vep_cache
      vep_ram: vep_ram
      vep_cores: vep_cores
      vep_buffer_size: vep_buffer_size
      disable_hotspot_annotation: disable_hotspot_annotation
      genomic_hotspots: genomic_hotspots
      protein_snv_hotspots: protein_snv_hotspots
      protein_indel_hotspots: protein_indel_hotspots
      maf_center: maf_center
      bwa_mem_index_image: bwa_mem_index_image
      make_bamout: make_bamout
      run_orientation_bias_mixture_model_filter: run_orientation_bias_mixture_model_filter
      mutect_cores: mutect_cores
      mutect_memory: mutect_memory
      getpileup_memory: getpileup_memory
      learnorientation_memory: learnorientation_memory
      filtermutectcalls_memory: filtermutectcalls_memory
      filteralignmentartifacts_cores: filteralignmentartifacts_cores
      filteralignmentartifacts_memory: filteralignmentartifacts_memory
      tool_name: tool_name
      output_basename: output_basename
    out: [mutect2_filtered_stats, mutect2_filtered_vcf, mutect2_protected_outputs,
      mutect2_public_outputs, mutect2_bam]

$namespaces:
  sbg: https://sevenbridges.com
hints:
- class: 'sbg:maxNumberOfParallelInstances'
  value: 4
"sbg:license": Apache License 2.0
"sbg:publisher": KFDRC
"sbg:categories":
- BAM
- CRAM
- GATK
- GERMLINE
- MUTECT2
- SOMATIC
- TUMORONLY
- VCF
