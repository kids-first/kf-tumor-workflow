cwlVersion: v1.0
class: CommandLineTool
id: gatk_calulcate_contamination
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/gatk:4.1.1.0'
  - class: ResourceRequirement
    ramMin: $(inputs.max_memory * 1000)
    coresMin: $(inputs.cores)
baseCommand: [/gatk, CalculateContamination]
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      --java-options "-Xmx${return Math.floor(inputs.max_memory*1000/1.074-1)}m"

inputs:
  tumor_pileup: { type: 'File', inputBinding: { prefix: "--input", position: 2 }, doc: "The input table from the tumor pileup" }
  normal_pileup: { type: 'File?', inputBinding: { prefix: "--matched-normal", position: 2 }, doc: "The matched normal input table from the normal pileup" }
  output_contam_table: { type: 'string', inputBinding: { prefix: "--output", position: 2 }, doc: "The output table" }
  output_seg_table: { type: 'string?', inputBinding: { prefix: "--tumor-segmentation", position: 2 }, doc: "The output table containing segmentation of the tumor by minor allele fraction" }

  extra_args: { type: 'string?', inputBinding: { position: 3 }, doc: "Any additional arguments for this tool. See GATK Documentation for complete list of options. Example input: --interval-merging-rule OVERLAPPING_ONLY" }

  cores: { type: 'int?', default: 2, doc: "Minimum reserved number of CPU cores for the task." }
  max_memory: { type: 'int?', default: 4, doc: "GB of RAM to allocate to the task." }
outputs:
  contamination_table:
    type: File
    outputBinding:
      glob: $(inputs.output_contam_table)
  segmentation_table:
    type: File
    outputBinding:
      glob: $(inputs.output_seg_table)

