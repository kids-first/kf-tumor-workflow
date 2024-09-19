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
  num_shards: { type: 'int?', default: 32, doc: "Number of shards to create." }
  output_gvcf: { type: 'boolean?', doc: "Set to true to return a gVCF output." }
  output_basename: { type: 'string', doc: "String to use as basename for outputs." }
  pon_filtering_custom: { type: 'File?', doc: "Optional. Custom VCF with Panel of Normals (PON) data.If set, the output VCF will be filtered: any variants that appear in PON will be marked with a PON filter, and PASS filter value will be removed." }
  process_somatic: { type: 'boolean?', default: true, doc: "Optional. Is specified the input is treated as somatic." }
  reads_normal: { type: 'File?', secondaryFiles: [{ pattern: '.bai', required: false }, { pattern: '^.bai', required: false }, { pattern: '.crai', required: false }, { pattern: '^.crai', required: false }], doc: "Reads from the normal matched sample. Aligned, sorted, indexed BAM/CRAM file. Should be aligned to a reference genome compatible with --ref. Can provide multiple BAMs/CRAMs (comma-separated)." }
  reads_tumor: { type: 'File', secondaryFiles: [{ pattern: '.bai', required: false }, { pattern: '^.bai', required: false }, { pattern: '.crai', required: false }, { pattern: '^.crai', required: false }], doc: "Reads from the tumor sample. Aligned, sorted, indexed BAM/CRAM file. Should be aligned to a reference genome compatible with --ref. Can provide multiple BAMs/CRAMs (comma-separated)." }
  indexed_reference_fasta: { type: 'File', secondaryFiles: [{ pattern: ".fai", required: true }], doc: "Genome reference to use. Must have an associated FAI index as well. Supports text or gzipped references. Should match the reference used to align the provided reads." }
  regions_file: { type: 'File[]?', doc: "Optional. Space-separated list of regions we want to process." }
  report_title: { type: 'string?', doc: "Title for the VCF stats report (HTML).If not provided, the title will be the sample name." }
  sample_name_normal: { type: 'string?', doc: "Sample name to use instead of the sample name from the input reads_normal BAM (SM tag in the header). This flag is used for both make_examples_somatic and postprocess_variants." }
  sample_name_tumor: { type: 'string?', doc: "Sample name to use instead of the sample name from the input reads_tumor BAM (SM tag in the header). This flag is used for both make_examples_somatic and postprocess_variants." }
  use_candidate_partition: { type: 'boolean?', doc: "If set, make_examples is run over partitions that contain an equal number of candidates. Default value is False.NOTE: This is experimental. So far, we have seen runtime improvements on PACBIO data." }
  use_default_pon_filtering: { type: 'boolean?', doc: "If true then default PON filtering will be used in tumor-only models." }
  vcf_stats_report: { type: 'boolean?', doc: "Output a visual report (HTML) of statistics about the output VCF." }

outputs:
  deepsomatic_all_vcf: { type: 'File', outputSource: deepsomatic_postprocess_variants/output }
  deepsomatic_all_vcf_stats: { type: 'File?', outputSource: deepsomatic_vcf_stats_report/stats_html }
  deepsomatic_gvcf: { type: 'File?', outputSource: deepsomatic_postprocess_variants/gvcf }

steps:
  make_task_array:
    run:
      cwlVersion: v1.2
      class: CommandLineTool
      doc: "Make an array of input length with values beginning at input index"
      requirements: [class: InlineJavascriptRequirement]
      baseCommand: [echo, done]
      inputs:
        index: { type: 'int?', default: 0 }
        length: { type: 'int' }
      outputs:
        array:
          type: 'int[]'
          outputBinding:
            outputEval: |
              $(Array(inputs.length).fill(inputs.index).map(function(x, y) { return  x + y }))
    in:
      length: num_shards
    out: [array]
  deepsomatic_make_examples_candidate_sweep:
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c6i.8xlarge
    run: ../tools/deepsomatic_make_examples.cwl
    when: $(inputs.use_candidate_partition == true)
    scatter: [task]
    in:
      use_candidate_partition: use_candidate_partition
      task: make_task_array/array
      task_total: num_shards
      mode:
        valueFrom: "candidate_sweep"
      output_gvcf: output_gvcf
      output_candidate_positions:
        valueFrom: $(1 == 1)
      ref: indexed_reference_fasta
      reads_tumor:
        source: reads_tumor
        valueFrom: $([self])
      reads_normal:
        source: reads_normal
        valueFrom: |
          $(self ? [self] : null)
      checkpoint:
        source: model_type
        valueFrom: |
          ${
            if (inputs.checkpoint_custom != null) { return null };
            var MODEL_TYPE_MAP = {
              'WGS': '/opt/models/deepsomatic/wgs',
              'WES': '/opt/models/deepsomatic/wes',
              'PACBIO': '/opt/models/deepsomatic/pacbio',
              'ONT': '/opt/models/deepsomatic/ont',
              'FFPE_WGS': '/opt/models/deepsomatic/ffpe_wgs',
              'FFPE_WES': '/opt/models/deepsomatic/ffpe_wes',
              'WGS_TUMOR_ONLY': '/opt/models/deepsomatic/wgs_tumor_only',
              'PACBIO_TUMOR_ONLY': '/opt/models/deepsomatic/pacbio_tumor_only',
              'ONT_TUMOR_ONLY': '/opt/models/deepsomatic/ont_tumor_only',
            };
            return MODEL_TYPE_MAP[self];
          }
      checkpoint_custom: customized_model
      regions_file: regions_file
      sample_name_tumor: sample_name_tumor
      sample_name_normal: sample_name_normal
      population_vcf_files:
        source: [use_default_pon_filtering, model_type]
        valueFrom: |
          ${
            if (self[0] != true) { return null };
            var PON_MAP = {
              'WGS_TUMOR_ONLY': '/opt/models/deepsomatic/pons/AF_ilmn_PON_DeepVariant.GRCh38.AF0.05.vcf.gz',
              'PACBIO_TUMOR_ONLY': '/opt/models/deepsomatic/pons/AF_pacbio_PON_CoLoRSdb.GRCh38.AF0.05.vcf.gz',
              'ONT_TUMOR_ONLY': '/opt/models/deepsomatic/pons/AF_pacbio_PON_CoLoRSdb.GRCh38.AF0.05.vcf.gz',
            };
            return PON_MAP[self[1]];
          }
      population_vcf_files_custom:
        source: [use_default_pon_filtering, pon_filtering_custom]
        valueFrom: |
          $(self[0] != true ? self[1] : null)
      extra_args:
        source: [make_examples_extra_args, model_type]
        valueFrom: |
          ${
            var MAKE_EXAMPLE_EA_MAP = {
              'WGS': '--vsc_min_fraction_indels 0.05 --vsc_min_fraction_snps 0.029 --vsc_max_fraction_indels_for_non_target_sample 0.5 --vsc_max_fraction_snps_for_non_target_sample 0.5',
              'WES': '--vsc_min_fraction_indels 0.05 --vsc_min_fraction_snps 0.029 --vsc_max_fraction_indels_for_non_target_sample 0.5 --vsc_max_fraction_snps_for_non_target_sample 0.5',
              'FFPE_WGS': '--vsc_min_fraction_indels 0.05 --vsc_min_fraction_snps 0.029 --vsc_max_fraction_indels_for_non_target_sample 0.5 --vsc_max_fraction_snps_for_non_target_sample 0.5',
              'FFPE_WES': '--vsc_min_fraction_indels 0.05 --vsc_min_fraction_snps 0.029 --vsc_max_fraction_indels_for_non_target_sample 0.5 --vsc_max_fraction_snps_for_non_target_sample 0.5',
              'WGS_TUMOR_ONLY': '--vsc_min_fraction_indels 0.07 --vsc_min_fraction_snps 0.05 --vsc_max_fraction_indels_for_non_target_sample 0.5 --vsc_max_fraction_snps_for_non_target_sample 0.5',
              'PACBIO': '--alt_aligned_pileup diff_channels --max_reads_for_dynamic_bases_per_region 1500 --max_reads_per_partition 0 --min_mapping_quality 5 --parse_sam_aux_fields --partition_size 25000 --phase_reads --pileup_image_width 99 --norealign_reads --sort_by_haplotypes --track_ref_reads --trim_reads_for_pileup --vsc_max_fraction_indels_for_non_target_sample 0.5 --vsc_max_fraction_snps_for_non_target_sample 0.5 --vsc_min_count_snps 1 --vsc_min_fraction_indels 0.1 --vsc_min_fraction_snps 0.02',
              'PACBIO_TUMOR_ONLY': '--alt_aligned_pileup diff_channels --max_reads_for_dynamic_bases_per_region 1500 --max_reads_per_partition 0 --min_mapping_quality 5 --parse_sam_aux_fields --partition_size 25000 --phase_reads --pileup_image_width 99 --norealign_reads --sort_by_haplotypes --track_ref_reads --trim_reads_for_pileup --vsc_max_fraction_indels_for_non_target_sample 0.5 --vsc_max_fraction_snps_for_non_target_sample 0.5 --vsc_min_count_snps 1 --vsc_min_fraction_indels 0.1 --vsc_min_fraction_snps 0.02',
              'ONT': '--alt_aligned_pileup diff_channels --max_reads_for_dynamic_bases_per_region 1500 --max_reads_per_partition 0 --min_mapping_quality 5 --parse_sam_aux_fields --partition_size 25000 --phase_reads --pileup_image_width 99 --norealign_reads --sort_by_haplotypes --track_ref_reads --trim_reads_for_pileup --vsc_max_fraction_indels_for_non_target_sample 0.6 --vsc_max_fraction_snps_for_non_target_sample 0.6 --vsc_min_fraction_indels 0.1 --vsc_min_fraction_snps 0.05',
              'ONT_TUMOR_ONLY': '--alt_aligned_pileup diff_channels --max_reads_for_dynamic_bases_per_region 1500 --max_reads_per_partition 0 --min_mapping_quality 5 --parse_sam_aux_fields --partition_size 25000 --phase_reads --pileup_image_width 99 --norealign_reads --sort_by_haplotypes --track_ref_reads --trim_reads_for_pileup --vsc_max_fraction_indels_for_non_target_sample 0.6 --vsc_max_fraction_snps_for_non_target_sample 0.6 --vsc_min_fraction_indels 0.1 --vsc_min_fraction_snps 0.05',
            };
            if (self[0] != null) {
              return self[0] + " " + MAKE_EXAMPLE_EA_MAP[self[1]];
            } else { 
              return MAKE_EXAMPLE_EA_MAP[self[1]];
            }
          }
    out: [examples, candidates, gvcf]
  deepsomatic_make_examples_calling:
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c6i.8xlarge
    run: ../tools/deepsomatic_make_examples.cwl
    scatter: [task]
    in:
      task: make_task_array/array
      task_total: num_shards
      mode:
        valueFrom: "calling"
      output_gvcf: output_gvcf
      candidate_positions:
        source: [use_candidate_partition, deepsomatic_make_examples_candidate_sweep/candidates]
        valueFrom: |
          $(self[0] != true ? null : self[1])
      ref: indexed_reference_fasta
      reads_tumor:
        source: reads_tumor
        valueFrom: $([self])
      reads_normal:
        source: reads_normal
        valueFrom: |
          $(self ? [self] : null)
      checkpoint:
        source: model_type
        valueFrom: |
          ${
            if (inputs.checkpoint_custom != null) { return null };
            var MODEL_TYPE_MAP = {
              'WGS': '/opt/models/deepsomatic/wgs',
              'WES': '/opt/models/deepsomatic/wes',
              'PACBIO': '/opt/models/deepsomatic/pacbio',
              'ONT': '/opt/models/deepsomatic/ont',
              'FFPE_WGS': '/opt/models/deepsomatic/ffpe_wgs',
              'FFPE_WES': '/opt/models/deepsomatic/ffpe_wes',
              'WGS_TUMOR_ONLY': '/opt/models/deepsomatic/wgs_tumor_only',
              'PACBIO_TUMOR_ONLY': '/opt/models/deepsomatic/pacbio_tumor_only',
              'ONT_TUMOR_ONLY': '/opt/models/deepsomatic/ont_tumor_only',
            };
            return MODEL_TYPE_MAP[self];
          }
      checkpoint_custom: customized_model
      regions_file: regions_file
      sample_name_tumor: sample_name_tumor
      sample_name_normal: sample_name_normal
      population_vcf_files:
        source: [use_default_pon_filtering, model_type]
        valueFrom: |
          ${
            if (self[0] != true) { return null };
            var PON_MAP = {
              'WGS_TUMOR_ONLY': '/opt/models/deepsomatic/pons/AF_ilmn_PON_DeepVariant.GRCh38.AF0.05.vcf.gz',
              'PACBIO_TUMOR_ONLY': '/opt/models/deepsomatic/pons/AF_pacbio_PON_CoLoRSdb.GRCh38.AF0.05.vcf.gz',
              'ONT_TUMOR_ONLY': '/opt/models/deepsomatic/pons/AF_pacbio_PON_CoLoRSdb.GRCh38.AF0.05.vcf.gz',
            };
            return PON_MAP[self[1]];
          }
      population_vcf_files_custom:
        source: [use_default_pon_filtering, pon_filtering_custom]
        valueFrom: |
          $(self[0] != true ? self[1] : null)
      extra_args:
        source: [make_examples_extra_args, model_type]
        valueFrom: |
          ${
            var MAKE_EXAMPLE_EA_MAP = {
              'WGS': '--vsc_min_fraction_indels 0.05 --vsc_min_fraction_snps 0.029 --vsc_max_fraction_indels_for_non_target_sample 0.5 --vsc_max_fraction_snps_for_non_target_sample 0.5',
              'WES': '--vsc_min_fraction_indels 0.05 --vsc_min_fraction_snps 0.029 --vsc_max_fraction_indels_for_non_target_sample 0.5 --vsc_max_fraction_snps_for_non_target_sample 0.5',
              'FFPE_WGS': '--vsc_min_fraction_indels 0.05 --vsc_min_fraction_snps 0.029 --vsc_max_fraction_indels_for_non_target_sample 0.5 --vsc_max_fraction_snps_for_non_target_sample 0.5',
              'FFPE_WES': '--vsc_min_fraction_indels 0.05 --vsc_min_fraction_snps 0.029 --vsc_max_fraction_indels_for_non_target_sample 0.5 --vsc_max_fraction_snps_for_non_target_sample 0.5',
              'WGS_TUMOR_ONLY': '--vsc_min_fraction_indels 0.07 --vsc_min_fraction_snps 0.05 --vsc_max_fraction_indels_for_non_target_sample 0.5 --vsc_max_fraction_snps_for_non_target_sample 0.5',
              'PACBIO': '--alt_aligned_pileup diff_channels --max_reads_for_dynamic_bases_per_region 1500 --max_reads_per_partition 0 --min_mapping_quality 5 --parse_sam_aux_fields --partition_size 25000 --phase_reads --pileup_image_width 99 --norealign_reads --sort_by_haplotypes --track_ref_reads --trim_reads_for_pileup --vsc_max_fraction_indels_for_non_target_sample 0.5 --vsc_max_fraction_snps_for_non_target_sample 0.5 --vsc_min_count_snps 1 --vsc_min_fraction_indels 0.1 --vsc_min_fraction_snps 0.02',
              'PACBIO_TUMOR_ONLY': '--alt_aligned_pileup diff_channels --max_reads_for_dynamic_bases_per_region 1500 --max_reads_per_partition 0 --min_mapping_quality 5 --parse_sam_aux_fields --partition_size 25000 --phase_reads --pileup_image_width 99 --norealign_reads --sort_by_haplotypes --track_ref_reads --trim_reads_for_pileup --vsc_max_fraction_indels_for_non_target_sample 0.5 --vsc_max_fraction_snps_for_non_target_sample 0.5 --vsc_min_count_snps 1 --vsc_min_fraction_indels 0.1 --vsc_min_fraction_snps 0.02',
              'ONT': '--alt_aligned_pileup diff_channels --max_reads_for_dynamic_bases_per_region 1500 --max_reads_per_partition 0 --min_mapping_quality 5 --parse_sam_aux_fields --partition_size 25000 --phase_reads --pileup_image_width 99 --norealign_reads --sort_by_haplotypes --track_ref_reads --trim_reads_for_pileup --vsc_max_fraction_indels_for_non_target_sample 0.6 --vsc_max_fraction_snps_for_non_target_sample 0.6 --vsc_min_fraction_indels 0.1 --vsc_min_fraction_snps 0.05',
              'ONT_TUMOR_ONLY': '--alt_aligned_pileup diff_channels --max_reads_for_dynamic_bases_per_region 1500 --max_reads_per_partition 0 --min_mapping_quality 5 --parse_sam_aux_fields --partition_size 25000 --phase_reads --pileup_image_width 99 --norealign_reads --sort_by_haplotypes --track_ref_reads --trim_reads_for_pileup --vsc_max_fraction_indels_for_non_target_sample 0.6 --vsc_max_fraction_snps_for_non_target_sample 0.6 --vsc_min_fraction_indels 0.1 --vsc_min_fraction_snps 0.05',
            };
            if (self[0] != null) {
              return self[0] + " " + MAKE_EXAMPLE_EA_MAP[self[1]];
            } else { 
              return MAKE_EXAMPLE_EA_MAP[self[1]];
            }
          }
      ram:
        source: model_type
        valueFrom: |
          $(self.search("PACBIO") != -1 ? 3 : null)
    out: [examples, candidates, gvcf]
  deepsomatic_call_variants:
    hints:
      - class: 'sbg:AWSInstanceType'
        value: g5.2xlarge
    run: ../tools/deepsomatic_call_variants.cwl
    in:
      checkpoint:
        source: model_type
        valueFrom: |
          ${
            if (inputs.checkpoint_custom != null) { return null };
            var MODEL_TYPE_MAP = {
              'WGS': '/opt/models/deepsomatic/wgs',
              'WES': '/opt/models/deepsomatic/wes',
              'PACBIO': '/opt/models/deepsomatic/pacbio',
              'ONT': '/opt/models/deepsomatic/ont',
              'FFPE_WGS': '/opt/models/deepsomatic/ffpe_wgs',
              'FFPE_WES': '/opt/models/deepsomatic/ffpe_wes',
              'WGS_TUMOR_ONLY': '/opt/models/deepsomatic/wgs_tumor_only',
              'PACBIO_TUMOR_ONLY': '/opt/models/deepsomatic/pacbio_tumor_only',
              'ONT_TUMOR_ONLY': '/opt/models/deepsomatic/ont_tumor_only',
            };
            return MODEL_TYPE_MAP[self];
          }
      checkpoint_custom: customized_model
      examples: deepsomatic_make_examples_calling/examples
      examples_name:
        source: num_shards
        valueFrom: |
          $(inputs.examples[0].basename.split('tfrecord')[0])tfrecord@$(self).gz
      outfile:
        source: output_basename
        valueFrom: $(self).cvo.tfrecord.gz
      extra_args: call_variants_extra_args
    out: [output] 
  deepsomatic_postprocess_variants:
    hints:
      - class: 'sbg:AWSInstanceType'
        value: g5.2xlarge
    run: ../tools/deepsomatic_postprocess_variants.cwl
    in:
      gvcf_outfile:
        source: [output_gvcf, output_basename]
        valueFrom: |
          $(self[0] == true ? self[1] + ".g.vcf.gz" : null)
      infile: deepsomatic_call_variants/output
      infile_path:
        source: output_basename
        valueFrom: $(self).cvo.tfrecord.gz
      nonvariant_site_tfrecords: deepsomatic_make_examples_calling/gvcf
      nonvariant_site_tfrecord_path:
        source: num_shards
        valueFrom: |
          $(inputs.nonvariant_site_tfrecords[0].basename.split('tfrecord')[0])tfrecord@$(self).gz
      outfile:
        source: output_basename
        valueFrom: $(self).deepsomatic.all.vcf.gz
      pon_filtering:
        source: [use_default_pon_filtering, model_type]
        valueFrom: |
          ${
            if (self[0] != true) { return null };
            var PON_MAP = {
              'WGS_TUMOR_ONLY': '/opt/models/deepsomatic/pons/PON_dbsnp138_gnomad_ILMN1000g_pon.vcf.gz',
              'PACBIO_TUMOR_ONLY': '/opt/models/deepsomatic/pons/PON_dbsnp138_gnomad_PB1000g_pon.vcf.gz',
              'ONT_TUMOR_ONLY': '/opt/models/deepsomatic/pons/PON_dbsnp138_gnomad_PB1000g_pon.vcf.gz',
            };
            return PON_MAP[self[1]];
          }
      pon_filtering_custom:
        source: [use_default_pon_filtering, pon_filtering_custom]
        valueFrom: |
          $(self[0] != true ? self[1] : null)
      ref: indexed_reference_fasta
      extra_args:
        source: process_somatic
        valueFrom: |
          $(self ? "--process_somatic" : null)
    out: [output, gvcf]
  deepsomatic_vcf_stats_report:
    run: ../tools/deepsomatic_vcf_stats_report.cwl
    when: $(inputs.vcf_stats_report)
    in:
      vcf_stats_report: vcf_stats_report
      input_vcf: deepsomatic_postprocess_variants/output
      outfile_base: output_basename
      title: report_title
    out: [stats_html]

$namespaces:
  sbg: https://sevenbridges.com

hints:
- class: "sbg:maxNumberOfParallelInstances"
  value: 2
