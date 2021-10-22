cwlVersion: v1.0
class: CommandLineTool
id: gatk_gatherpileupsummaries
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: 'broadinstitute/gatk:4.2.2.0'
  - class: ResourceRequirement
    ramMin: $(inputs.max_memory * 1000)
    coresMin: $(inputs.cores)
baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      gatk --java-options "-Xmx${return Math.floor(inputs.max_memory*1000/1.074-1)}m" GatherPileupSummaries
      --sequence-dictionary $(inputs.reference_dict.path)
      -O $(inputs.output_basename).$(inputs.tool_name).merged.pileup.table

inputs:
  input_tables:
    type:
      type: array
      items: File
      inputBinding:
        prefix: -I
    inputBinding:
      position: 1
  reference_dict: File
  tool_name: string
  output_basename: string
  max_memory: { type: 'int?', default: 4, doc: "GB of memory to allocate ot this task" }
  cores: { type: 'int?', default: 2, doc: "CPUs to allocate ot this task" }
outputs:
  merged_table:
    type: File
    outputBinding:
      glob: '*.merged.pileup.table'
