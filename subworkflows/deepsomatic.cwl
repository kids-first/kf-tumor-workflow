cwlVersion: v1.2
class: Workflow
id: deepsomatic 
doc: "DeepSomatic Workflow"

requirements:
- class: ScatterFeatureRequirement
- class: StepInputExpressionRequirement
- class: MultipleInputFeatureRequirement
- class: InlineJavascriptRequirement
- class: SubworkflowFeatureRequirement

inputs:
  call_variants_extra_args: { type: 'string?', doc: "Any advanced or extra arguments to hand to call_variants" }
  customized_model: { type: 'File?', doc: "Optional. A path to a model checkpoint to load for the `call_variants` step. If not set, the default for each --model_type will be used" }
  make_examples_extra_args: { type: 'string?', doc: "Any advanced or extra arguments to hand to make_examples" }
  model_type:
    type:
      - type: enum
        name: model_type
        symbols: ["WGS", "WES", "PACBIO", "ONT", "FFPE_WGS", "FFPE_WES", "WGS_TUMOR_ONLY", "PACBIO_TUMOR_ONLY", "ONT_TUMOR_ONLY"]
    doc: |
      Type of model to use for variant calling. Set this flag to use the default model
      associated with each type, and it will set necessary flags corresponding to each model.
      If you want to use a customized model, add --customized_model flag in addition to this flag.
  num_shards: { type: 'int?', default: 36, doc: "Number of shards to create." }
  output_basename: { type: 'string', doc: "String to use as basename for outputs." }
  pon_filtering: { type: 'string?', doc: "Optional. Path to VCF file (in the docker) with Panel of Normals (PON) data.If set, the output VCF will be filtered: any variants that appear in PON will be marked with a PON filter, and PASS filter value will be removed." }
  pon_filtering_custom: { type: 'File?', doc: "Optional. Custom VCF with Panel of Normals (PON) data.If set, the output VCF will be filtered: any variants that appear in PON will be marked with a PON filter, and PASS filter value will be removed." }
  process_somatic: { type: 'boolean?', doc: "Optional. Is specified the input is treated as somatic." }
  reads_normal: { type: 'File[]?', secondaryFiles: [{ pattern: '.bai', required: false }, { pattern: '^.bai', required: false }, { pattern: '.crai', required: false }, { pattern: '^.crai', required: false }], doc: "Reads from the normal matched sample. Aligned, sorted, indexed BAM/CRAM file. Should be aligned to a reference genome compatible with --ref. Can provide multiple BAMs/CRAMs (comma-separated)." }
  reads_tumor: { type: 'File[]', secondaryFiles: [{ pattern: '.bai', required: false }, { pattern: '^.bai', required: false }, { pattern: '.crai', required: false }, { pattern: '^.crai', required: false }], doc: "Reads from the tumor sample. Aligned, sorted, indexed BAM/CRAM file. Should be aligned to a reference genome compatible with --ref. Can provide multiple BAMs/CRAMs (comma-separated)." }
  indexed_reference_fasta: { type: 'File', secondaryFiles: [{ pattern: ".fai", required: true }], doc: "Genome reference to use. Must have an associated FAI index as well. Supports text or gzipped references. Should match the reference used to align the provided reads." }
  regions_file: { type: 'File[]?', doc: "Optional. Space-separated list of regions we want to process." }
  report_title: { type: 'string?', doc: "Title for the VCF stats report (HTML).If not provided, the title will be the sample name." }
  sample_name_normal: { type: 'string?', doc: "Sample name to use instead of the sample name from the input reads_normal BAM (SM tag in the header). This flag is used for both make_examples_somatic and postprocess_variants." }
  sample_name_tumor: { type: 'string?', doc: "Sample name to use instead of the sample name from the input reads_tumor BAM (SM tag in the header). This flag is used for both make_examples_somatic and postprocess_variants." }
  use_candidate_partition: { type: 'boolean?', doc: "If set, make_examples is run over partitions that contain an equal number of candidates. Default value is False.NOTE: This is experimental. So far, we have seen runtime improvements on PACBIO data." }
  use_default_pon_filtering: { type: 'boolean?', doc: "If true then default PON filtering will be used in tumor-only models." }
  vcf_stats_report: { type: 'boolean?', doc: "Output a visual report (HTML) of statistics about the output VCF." }

outputs:
  deepvariant_all_vcf: { type: 'File', outputSource: bcftools_index/output }
  deepvariant_pass_vcf: { type: 'File', outputSource: bcftools_filter_index/output }
  deepvariant_all_vcf_stats: { type: 'File', outputSource: bcftools_stats_all/stats }
  deepvariant_pass_vcf_stats: { type: 'File', outputSource: bcftools_stats_pass/stats }
  bcftools_vcf: { type: 'File?', outputSource: canine_annotation_module/bcftools_vcf }
  tumor_only_vcf: { type: 'File?', outputSource: canine_annotation_module/tumor_only_vcf }
  snpeff_all_vcf: { type: 'File?', outputSource: canine_annotation_module/snpeff_all_vcf }
  snpeff_canon_vcf: { type: 'File?', outputSource: canine_annotation_module/snpeff_canon_vcf }
  vep_all_vcf: { type: 'File?', outputSource: canine_annotation_module/vep_all_vcf }
  vep_con_vcf: { type: 'File?', outputSource: canine_annotation_module/vep_con_vcf }

steps:
  make_task_array:
    run:
      cwlVersion: v1.2
      class: CommandLineTool
      doc: "Make an array of input length with values beginning at input index"
      requirements: [class: InlineJavascriptRequirement]
      baseCommand: [echo, done]
      inputs:
        index: { type: 'int?', default: 1 }
        length: { type: 'int' }
      outputs:
        array: { type: 'int[]', outputBinding: { outputEval: $(Array(inputs.length).fill(inputs.index).map(function(x, y) { return  x + y })) } }
    in:
      length: num_shards
    out: [array]
  deepsomatic_make_examples_candidate_sweep:
    run: ../tools/deepsomatic_make_examples.cwl
    scatter: [task]
    when: $(inputs.use_candidate_partition)
    in:
      use_candidate_partition: use_candidate_partition
    out: []
  deepsomatic_make_examples_calling:
    run: ../tools/deepsomatic_make_examples.cwl
    scatter: [task]
    in:
    out: []
  deepsomatic_call_variants:
    hints:
      - class: 'sbg:AWSInstanceType'
        value: p3.2xlarge
    run: ../tools/deepsomatic_call_variants.cwl
    in:
    out: [] 
  deepsomatic_postprocess_variants:
    run: ../tools/deepsomatic_postprocess_variants.cwl
    in:
    out: []
  deepsomatic_vcf_stats_report:
    run: ../tools/deepsomatic_vcf_stats_report.cwl
    when: $(inputs.vcf_stats_report)
    in:
      vcf_stats_report: vcf_stats_report


  expr_make_int_array:
    run: ../tools/expr_make_int_array.cwl
    in:
      length: num_shards
    out: [output]

  deepvariant_make_examples:
    run: ../tools/deepvariant_make_examples.cwl
    scatter: [task]
    in:
      mode:
        valueFrom: "calling"
      task: expr_make_int_array/output
      task_total: num_shards
      ref: indexed_reference_fasta
      reads:
        source: input_reads
        valueFrom: $([self])
      examples_outname:
        source: output_basename
        valueFrom: $(self).ex.tfrecord@$(inputs.task_total).gz
      gvcf_outname: 
        source: output_basename
        valueFrom: $(self).gvcf.tfrecord@$(inputs.task_total).gz
    out: [examples, gvcf]

  deepvariant_call_variants:
    hints:
      - class: 'sbg:AWSInstanceType'
        value: p3.2xlarge
    run: ../tools/deepvariant_call_variants.cwl
    in:
      sample_type: sample_type 
      examples: deepvariant_make_examples/examples
      examples_name:
        source: num_shards
        valueFrom: |
          $(inputs.examples[0].basename.split('tfrecord')[0])tfrecord@$(self).gz
      outfile:
        source: output_basename
        valueFrom: $(self).cvo.tfrecord.gz
    out: [output]

  deepvariant_postprocess_variants:
    hints:
      - class: 'sbg:AWSInstanceType'
        value: p3.2xlarge
    run: ../tools/deepvariant_postprocess_variants.cwl
    in:
      ref: indexed_reference_fasta
      infile:
        source: deepvariant_call_variants/output
        valueFrom: $([self])
      outfile:
        source: output_basename
        valueFrom: $(self).deepvariant.all.vcf.gz
      nonvariant_site_tfrecord_path:
        source: num_shards
        valueFrom: |
          $(inputs.nonvariant_site_tfrecords[0].basename.split('tfrecord')[0])tfrecord@$(self).gz
      nonvariant_site_tfrecords: deepvariant_make_examples/gvcf
      gvcf_outfile:
        source: output_basename
        valueFrom: $(self).deepvariant.all.g.vcf.gz
    out: [output, gvcf]


$namespaces:
  sbg: https://sevenbridges.com

hints:
- class: "sbg:maxNumberOfParallelInstances"
  value: 2
