cwlVersion: v1.0
class: CommandLineTool
id: gatk_mergemutectstats
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
      gatk --java-options "-Xmx${return Math.floor(inputs.max_memory*1000/1.074-1)}m" MergeMutectStats
      -O $(inputs.output_basename).Mutect2.merged.stats

inputs:
  input_stats:
    type:
      type: array
      items: File
      inputBinding:
        prefix: --stats
    inputBinding:
      position: 1
  output_basename: string
  max_memory: { type: 'int?', default: 4, doc: "Maximum memory to allocate to this task" }
  cores: { type: 'int?', default: 2, doc: "CPUs to allocate to this task" }
outputs:
  merged_stats:
    type: File
    outputBinding:
      glob: '*.merged.stats'
