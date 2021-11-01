cwlVersion: v1.2
class: Workflow
id: mutect2_filter_support
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: InlineJavascriptRequirement

inputs:
  indexed_reference_fasta: {type: 'File', secondaryFiles: [.fai, ^.dict]}
  reference_dict: File
  wgs_calling_interval_list: "File[]"
  input_tumor_aligned:
    type: File
    secondaryFiles: [{ pattern: ".bai", required: false },{ pattern: "^.crai", required: false }]
    doc: "tumor BAM or CRAM"
  input_normal_aligned:
    type: 'File?'
    secondaryFiles: [{ pattern: ".bai", required: false },{ pattern: "^.crai", required: false }]
    doc: "normal BAM or CRAM"

  exac_common_vcf: {type: 'File', secondaryFiles: [.tbi]}
  getpileup_memory: {type: 'int?'}
  learnorientation_memory: {type: 'int?'}
  output_basename: string
  tool_name: string

outputs:
  contamination_table: {type: 'File', outputSource: gatk_calculate_contamination/contamination_table}
  segmentation_table: {type: 'File', outputSource: gatk_calculate_contamination/segmentation_table}
  
steps:
  gatk_get_tumor_pileup_summaries:
    run: ../tools/gatk_getpileupsummaries.cwl
    in:
      aligned_reads: input_tumor_aligned
      reference: indexed_reference_fasta
      interval_list: wgs_calling_interval_list
      exac_common_vcf: exac_common_vcf
      max_memory: getpileup_memory
    scatter: [interval_list]
    out: [pileup_table]

  gatk_get_normal_pileup_summaries:
    run: ../tools/gatk_getpileupsummaries.cwl
    when: $(inputs.aligned_reads != null)
    in:
      aligned_reads: input_normal_aligned
      reference: indexed_reference_fasta
      interval_list: wgs_calling_interval_list
      exac_common_vcf: exac_common_vcf
      max_memory: getpileup_memory
    scatter: [interval_list]
    out: [pileup_table]

  gatk_gather_tumor_pileup_summaries:
    run: ../tools/gatk_gatherpileupsummaries.cwl
    in:
      input_tables: gatk_get_tumor_pileup_summaries/pileup_table
      output_basename: output_basename
      reference_dict: reference_dict
      tool_name: tool_name
    out: [merged_table]

  gatk_gather_normal_pileup_summaries:
    run: ../tools/gatk_gatherpileupsummaries.cwl
    when: $(inputs.input_tables.some(function(e) { return e != null; }))
    in:
      input_tables: gatk_get_normal_pileup_summaries/pileup_table
      output_basename: output_basename
      reference_dict: reference_dict
      tool_name: tool_name
    out: [merged_table]
  
  gatk_calculate_contamination:
    run: ../tools/gatk_calculatecontamination.cwl
    in:
      tumor_pileup: gatk_gather_tumor_pileup_summaries/merged_table
      normal_pileup: gatk_gather_normal_pileup_summaries/merged_table
      output_contam_table:
        source: output_basename
        valueFrom: $(self).contamination.table
      output_seg_table:
        source: output_basename
        valueFrom: $(self).segmentation.table
    out: [contamination_table, segmentation_table]

$namespaces:
  sbg: https://sevenbridges.com
