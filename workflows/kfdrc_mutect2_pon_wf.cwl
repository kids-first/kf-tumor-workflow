cwlVersion: v1.2
class: Workflow
id: kfdrc_mutect2_create_pon_wf
label: KFDRC Mutect2 Create Panel Workflow
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


  # Create Somatic Panel of Normals Workflow

  ## Workflow Outline
  For each `input_normal_aligned, input_normal_name` pair:
  - Mutect2 Minimal Run

  Using the collected outputs from above:
  - SplitIntevals
  - GenomicsDBImport
  - CreateSomaticPanelOfNormals
  - MergeVcfs

  ## Workflow Description
  The Create Somatic Panel of Normals workflow is a CWL port of the [GATK Best Pratice WDL](https://github.com/broadinstitute/gatk/blob/master/scripts/mutect2_wdl/mutect2_pon.wdl).
requirements:
- class: ScatterFeatureRequirement
- class: MultipleInputFeatureRequirement
- class: SubworkflowFeatureRequirement
- class: InlineJavascriptRequirement

inputs:
  # GATK Files
  reference_fasta: {type: 'File', doc: "Reference fasta", "sbg:suggestedValue": {
      class: File, path: 60639014357c3a53540ca7a3, name: Homo_sapiens_assembly38.fasta},
    "sbg:fileTypes": "FASTA,FA"}
  reference_fai: {type: 'File?', doc: "FAI index for reference_fasta", "sbg:suggestedValue": {
      class: File, path: 60639016357c3a53540ca7af, name: Homo_sapiens_assembly38.fasta.fai},
    "sbg:fileTypes": "FAI"}
  reference_dict: {type: 'File?', doc: "DICT index for reference_fasta", "sbg:suggestedValue": {
      class: File, path: 60639019357c3a53540ca7e7, name: Homo_sapiens_assembly38.dict},
    "sbg:fileTypes": "DICT"}
  input_normal_aligned: {type: 'File[]', secondaryFiles: [{pattern: ".bai", required: false},
      {pattern: "^.bai", required: false}, {pattern: ".crai", required: false}, {
        pattern: "^.crai", required: false}], doc: "BAM/SAM/CRAM files containing\
      \ normal reads", "sbg:fileTypes": "BAM,CRAM,SAM"}
  mutect2_af_only_gnomad_vcf: {type: 'File', secondaryFiles: ['.tbi'], doc: "Population\
      \ vcf (and index) of germline sequencing containing allele fractions in VCF\
      \ format. Used in the panel of normal creation to subtract germline variants\
      \ from the panel of normals, keeping only technical artifacts.", "sbg:suggestedValue": {
      class: File, path: 5f50018fe4b054958bc8d2e3, name: af-only-gnomad.hg38.vcf.gz,
      secondaryFiles: [{class: File, path: 5f50018fe4b054958bc8d2e5, name: af-only-gnomad.hg38.vcf.gz.tbi}]},
    "sbg:fileTypes": "VCF,VCF.GZ"}

  # GATK Arguments
  input_normal_name: {type: 'string[]', doc: "BAM sample name of normal(s), if any.\
      \ May be URL-encoded as output by GetSampleName with -encode argument."}
  tool_name: {type: 'string?', doc: "Tool name for this workflow", default: "mutect2"}
  output_basename: {type: 'string', doc: "String value to use as basename for outputs"}
  mutect2_extra_args: {type: 'string?', doc: "Any additional arguments for Mutect2.\
      \ See GATK documentation for all available options."}
  mutect_scatter_count: {type: 'int?', default: 50, doc: "The number of files into\
      \ which to scatter the resulting list by locus; in some situations, fewer intervals\
      \ may be emitted"}
  filtermutectcalls_extra_args: {type: 'string?', doc: "Any additional arguments for\
      \ FilterMutectCalls. See GATK documentation for all available options."}
  pon_scatter_count: {type: 'int?', default: 24, doc: "The number of files into which\
      \ the intervals will be scattered for Panel of Normal creation"}
  min_contig_size: {type: 'int?', default: 1000000, doc: "Minimum contig size to keep\
      \ if getting intervals from the reference"}
  genomicsdbimport_extra_args: {type: 'string?', doc: "Any additional arguments for\
      \ GenomicsDBImport. See GATK documentation for all available options."}
  createpanel_extra_args: {type: 'string?', doc: "Any additional arguments for CreateSomaticPanelOfNormals.\
      \ See GATK documentation for all available options."}

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

  # Resource Control
  mutect_cores: {type: 'int?', doc: "CPUs to allocate to GATK Mutect2"}
  mutect_memory: {type: 'int?', doc: "GB of memory to allocate to GATK Mutect2 (hard-capped)"}
  filtermutectcalls_memory: {type: 'int?', doc: "GB of memory to allocate to GATK\
      \ FilterMutectCalls (hard-capped)"}
  genomicsdbimport_cores: {type: 'int?', doc: "CPUs to allocate to GATK GenomicsDBImport"}
  genomicsdbimport_memory: {type: 'int?', doc: "GB of memory to allocate to GATK GenomicsDBImport"}
  createpanel_cores: {type: 'int?', doc: "CPUs to allocate to GATK CreateSomaticPanelOfNormals"}
  createpanel_memory: {type: 'int?', doc: "GB of memory to allocate to GATK CreateSomaticPanelOfNormals"}

outputs:
  mutect2_panel_of_normals: {type: 'File', outputSource: gatk_mergevcfs/merged_vcf,
    doc: "The panel of normals (PoN) containing germline and artifactual sites for\
      \ use with Mutect2."}
  mutect2_prepass_vcfs: {type: 'File[]', outputSource: run_mutect2/mutect2_filtered_vcf,
    doc: "Collect of VCFs (one for each input aligned file with SNV, MNV, and INDEL\
      \ variant calls. Contains all soft FILTER values generated by variant caller."}

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
      scatter_ct: mutect_scatter_count
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
      input_tumor_aligned: input_normal_aligned
      input_tumor_name: input_normal_name
      disable_adaptive_pruning:
        source: wgs_or_wxs
        valueFrom: $(self == 'WXS')
      mutect2_extra_args:
        source: mutect2_extra_args
        valueFrom: ${ return [self,"--max-mnp-distance 0"].join(" ")}
      filtermutectcalls_extra_args: filtermutectcalls_extra_args
      mutect_cores: mutect_cores
      mutect_memory: mutect_memory
      filtermutectcalls_memory: filtermutectcalls_memory
      tool_name: tool_name
      output_basename:
        valueFrom: $(inputs.input_tumor_aligned.nameroot)
    scatter: [input_tumor_aligned, input_tumor_name]
    scatterMethod: dotproduct
    out: [mutect2_filtered_stats, mutect2_filtered_vcf]

  gatk_splitintervals_pon:
    run: ../tools/gatk_splitintervals.cwl
    in:
      interval_list: select_interval_list/output
      scatter_ct: pon_scatter_count
      reference_fasta: prepare_reference/indexed_fasta
      extra_args:
        source: min_contig_size
        valueFrom: ${ return "--dont-mix-contigs --min-contig-size " + self }
    out: [output]

  gatk_genomicsdbimport:
    run: ../tools/gatk_genomicsdbimport.cwl
    in:
      intervals: gatk_splitintervals_pon/output
      input_vcfs: run_mutect2/mutect2_filtered_vcf
      genomicsdb_path:
        valueFrom: ${ return "pon_db_" + inputs.intervals.nameroot  }
      reference_fasta: prepare_reference/indexed_fasta
      extra_args: genomicsdbimport_extra_args
      cores: genomicsdbimport_cores
      max_memory: genomicsdbimport_memory
    scatter: [intervals]
    out: [output]

  gatk_createpanelofnormals:
    hints:
    - class: 'sbg:AWSInstanceType'
      value: c5.9xlarge
    run: ../tools/gatk_createsomaticpanelofnormals.cwl
    in:
      input_genomicsdb: gatk_genomicsdbimport/output
      output_filename:
        valueFrom: $(inputs.input_genomicsdb.basename).vcf.gz
      reference_fasta: prepare_reference/indexed_fasta
      germline_resource: mutect2_af_only_gnomad_vcf
      extra_args: createpanel_extra_args
      cores: createpanel_cores
      max_memory: createpanel_memory
    scatter: [input_genomicsdb]
    out: [output]

  gatk_mergevcfs:
    run: ../tools/gatk_mergevcfs.cwl
    in:
      input_vcfs: gatk_createpanelofnormals/output
      reference_dict: prepare_reference/reference_dict
      tool_name: tool_name
      output_basename: output_basename
    out: [merged_vcf]

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
- MUTECT2
- TUMORONLY
- VCF
