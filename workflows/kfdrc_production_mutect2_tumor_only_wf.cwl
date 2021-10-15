cwlVersion: v1.2
class: Workflow
id: kfdrc_production_mutect2_tumor_only_wf
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement
  - class: InlineJavascriptRequirement

inputs:
  # GATK Files
  reference_fasta: { type: 'File', doc: "Reference fasta" }
  reference_fai: { type: 'File?', doc: "FAI index for reference_fasta" }
  reference_dict: { type: 'File?', doc: "DICT index for reference_fasta" }
  input_tumor_aligned: { type: 'File', secondaryFiles: [{ pattern: ".bai", required: false },{ pattern: "^.bai", required: false },{ pattern: ".crai", required: false },{ pattern: "^.crai", required: false }], doc: "BAM/SAM/CRAM file containing tumor reads" }
  input_normal_aligned: { type: 'File?', secondaryFiles: [{ pattern: ".bai", required: false },{ pattern: "^.bai", required: false },{ pattern: ".crai", required: false },{ pattern: "^.crai", required: false }], doc: "BAM/SAM/CRAM file containing normal reads"}
  panel_of_normals: { type: 'File?', secondaryFiles: ['.tbi'], doc: "VCF file (and index) of sites observed in normal. A panel of normals can be a useful (optional) input to help filter out commonly seen sequencing noise that may appear as low allele-fraction somatic variants." }
  mutect2_af_only_gnomad_vcf: { type: 'File', secondaryFiles: ['.tbi'], doc: "Population vcf (and index) of germline sequencing containing allele fractions in VCF format." }
  mutect2_alleles_vcf: { type: 'File?', secondaryFiles: ['.tbi'], doc: "VCF file (and index) containing set of alleles for which to force genotyping regardless of evidence" }
  mutect2_exac_common_vcf: { type: 'File?', secondaryFiles: ['.tbi'], doc: "Exac Common VCF (and index) used for calculating contamination values used in filtering the Mutect VCF. If do not wish to perfom this filtering, remove this input." }
  bwa_mem_index_image: { type: 'File?', doc: "BWA-mem index image used for FilterAlignmentArtifacts. WARNING!!! Experimental Software!!! Do not provide unless you are absolutely certain you want to run FilterAlignmentArtifacts." }

  # GATK Arguments
  input_tumor_name: { type: 'string', doc: "BAM sample name of tumor. May be URL-encoded as output by GetSampleName with -encode argument." }
  input_normal_name: { type: 'string?', doc: "BAM sample name of normal(s), if any. May be URL-encoded as output by GetSampleName with -encode argument." }
  tool_name: { type: 'string?', doc: "Tool name for this workflow", default: "mutect2" }
  output_basename: { type: 'string', doc: "String value to use as basename for outputs" }
  make_bamout: { type: 'boolean?', default: false, doc: "Should Mutect2 create a BAM output?" }
  run_orientation_bias_mixture_model_filter: { type: 'boolean?', default: true, doc: "Should Orientation Bias Mixture Model Filter be applied to Mutect2 outputs?" }
  mutect2_extra_args: { type: 'string?', doc: "Any additional arguments for Mutect2. See GATK documentation for all available options." }
  filtermutectcalls_extra_args: { type: 'string?', doc: "Any additional arguments for FilterMutectCalls. See GATK documentation for all available options." }
  select_vars_mode: { type: ['null', { type: enum, name: select_vars_mode, symbols: ["gatk", "grep"] }], default: "gatk", doc: "Choose 'gatk' for SelectVariants tool, or 'grep' for grep expression" }
  scatter_count: {type: 'int?', default: 50, doc: "The number of files into which to scatter the resulting list by locus; in some situations, fewer intervals may be emitted"}

  # WGS or WXS
  wgs_or_wxs: { type: { type: enum, name: wgs_or_wxs, symbols: ["WGS", "WXS"] }, doc: "Select input type WGS or WXS" }
  wgs_calling_interval_list: { type: 'File?', doc: "GATK intervals list-style, or bed file.  Recommend canocical chromosomes with N regions removed" }
  padded_capture_regions: { type: 'File?', doc: "Recommend 100bp pad, for somatic variant. Used as calling intervals for exomes." }

  # Annotation Arguments
  run_annotation: {type: 'boolean?', default: true, doc: "Should annotation be run on the final VCF?"}
  vep_cache: { type: 'File?', doc: "tar gzipped cache from ensembl/local converted cache" }
  vep_ref_build: { type: 'string?', default: "GRCh38", doc: "Genome ref build used, should line up with cache" }
  genomic_hotspots: { type: 'File[]?', doc: "Tab-delimited BED formatted file(s) containing hg38 genomic positions corresponding to hotspots", "sbg:suggestedValue": [{class: File, path: 607713829360f10e3982a423, name: tert.bed}] }
  protein_snv_hotspots: { type: 'File[]?', doc: "Column-name-containing, tab-delimited file(s) containing protein names and amino acid positions corresponding to hotspots", "sbg:suggestedValue": [{class: File, path: 607713829360f10e3982a426, name: protein_snv_cancer_hotspots_v2.tsv}] }
  protein_indel_hotspots: { type: 'File[]?', doc: "Column-name-containing, tab-delimited file(s) containing protein names and amino acid position ranges corresponding to hotspots", "sbg:suggestedValue": [{class: File, path: 607713829360f10e3982a424, name: protein_indel_cancer_hotspots_v2.tsv}] }
  retain_info: {type: 'string?', doc: "csv string with INFO fields that you want to keep", default: "MBQ,TLOD,HotSpotAllele"}
  retain_fmt: {type: 'string?', doc: "csv string with FORMAT fields that you want to keep"}
  add_common_fields: {type: 'boolean?', doc: "Set to true if input is a strelka2 vcf that hasn't had common fields added", default: false}
  bcftools_annot_columns: {type: 'string?', doc: "csv string of columns from annotation to port into the input vcf, i.e INFO/AF", default: "INFO/AF"}
  bcftools_annot_vcf: {type: 'File?', secondaryFiles: [.tbi], doc: "bgzipped annotation vcf file", "sbg:suggestedValue": {class: File, path: 5f50018fe4b054958bc8d2e3, name: af-only-gnomad.hg38.vcf.gz} }
  bcftools_public_filter: {type: 'string?', doc: "Will hard filter final result to create a public version", default: FILTER="PASS"|INFO/HotSpotAllele=1}
  gatk_filter_name: {type: 'string[]?', doc: "Array of names for each filter tag to add, recommend: [\"NORM_DP_LOW\", \"GNOMAD_AF_HIGH\"]"}
  gatk_filter_expression: {type: 'string[]?', doc: "Array of filter expressions to establish criteria to tag variants with. See https://gatk.broadinstitute.org/hc/en-us/articles/360036730071-VariantFiltration, recommend: \"vc.getGenotype('\" + inputs.input_normal_name + \"').getDP() <= 7\"), \"AF > 0.001\"]"}
  disable_hotspot_annotation: { type: 'boolean?', doc: "Disable Hotspot Annotation and skip this task.", default: false }
  maf_center: {type: 'string?', doc: "Sequencing center of variant called", default: "."}

  # Resource Control
  mutect_cores: { type: 'int?', doc: "CPUs to allocate to GATK Mutect2" }
  mutect_memory: { type: 'int?', doc: "GB of memory to allocate to GATK Mutect2 (hard-capped)" }
  getpileup_memory: { type: 'int?', doc: "GB of memory to allocate to GATK GetPileupSummaries (hard-capped)" }
  learnorientation_memory: { type: 'int?', doc: "GB of memory to allocate to GATK LearnReadOrientationModel (hard-capped)" }
  filtermutectcalls_memory: { type: 'int?', doc: "GB of memory to allocate to GATK FilterMutectCalls (hard-capped)" }
  filteralignmentartifacts_cores: { type: 'int?', doc: "CPUs to allocate to GATK FilterAlignmentArtifacts" }
  filteralignmentartifacts_memory: { type: 'int?', doc: "GB of memory to allocate to GATK FilterAlignmentArtifacts (hard-capped)" }

outputs:
  mutect2_prepass_vcf: { type: 'File', outputSource: run_mutect2/mutect2_filtered_vcf }
  mutect2_protected_outputs: { type: 'File[]?', outputSource: run_mutect2/mutect2_protected_outputs }
  mutect2_public_outputs: { type: 'File[]?', outputSource: run_mutect2/mutect2_public_outputs }
  mutect2_bam: { type: 'File?', outputSource: run_mutect2/mutect2_bam }

steps:
  choose_defaults:
    run: ../tools/mode_defaults.cwl
    in:
      input_mode: wgs_or_wxs
    out: [out_exome_flag,out_cnvkit_wgs_mode,out_i_flag,out_lancet_padding,out_lancet_window,out_vardict_padding]

  prepare_reference:
    run: ../subworkflows/prepare_reference.cwl
    in:
      input_fasta: reference_fasta
      input_fai: reference_fai
      input_dict: reference_dict
    out: [indexed_fasta,reference_dict]

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
      vep_cache: vep_cache
      vep_ref_build: vep_ref_build
      select_vars_mode: select_vars_mode
      genomic_hotspots: genomic_hotspots
      protein_snv_hotspots: protein_snv_hotspots
      protein_indel_hotspots: protein_indel_hotspots
      retain_info: retain_info
      retain_fmt: retain_fmt
      add_common_fields: add_common_fields
      bcftools_annot_columns: bcftools_annot_columns
      bcftools_annot_vcf: bcftools_annot_vcf
      bcftools_public_filter: bcftools_public_filter
      gatk_filter_name: gatk_filter_name
      gatk_filter_expression: gatk_filter_expression
      maf_center: maf_center
      bwa_mem_index_image: bwa_mem_index_image
      make_bamout: make_bamout
      run_orientation_bias_mixture_model_filter: run_orientation_bias_mixture_model_filter
      run_annotation: run_annotation
      disable_hotspot_annotation: disable_hotspot_annotation
      mutect_cores: mutect_cores
      mutect_memory: mutect_memory
      getpileup_memory: getpileup_memory
      learnorientation_memory: learnorientation_memory
      filtermutectcalls_memory: filtermutectcalls_memory
      filteralignmentartifacts_cores: filteralignmentartifacts_cores
      filteralignmentartifacts_memory: filteralignmentartifacts_memory
      tool_name: tool_name
      output_basename: output_basename
    out:
      [mutect2_filtered_stats, mutect2_filtered_vcf, mutect2_protected_outputs, mutect2_public_outputs, mutect2_bam]

$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 4
