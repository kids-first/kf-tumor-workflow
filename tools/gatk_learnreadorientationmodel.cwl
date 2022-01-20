cwlVersion: v1.0
class: CommandLineTool
id: gatk_learnreadorientationmodel
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
      gatk --java-options "-Xmx${return Math.floor(inputs.max_memory*1000/1.074-1)}m" LearnReadOrientationModel
      -O $(inputs.output_basename).$(inputs.tool_name).f1r2_bias.tar.gz

inputs:
  input_tgz:
    type:
      type: array
      items: File
      inputBinding:
        prefix: -I
    inputBinding:
      position: 1
  tool_name: string
  enable_tool: { type: 'boolean?', doc: "Should this tool be run? This option may only be used in a workflow." }
  output_basename: string
  max_memory: {type: 'int?', default: 4, doc: "Maximum memory in GB for GATK LearnReadOrientationModel"}
  cores: { type: 'int?', default: 2, doc: "CPUs to allocate to this task" }
outputs:
  f1r2_bias:
    type: File
    outputBinding:
      glob: '*.f1r2_bias.tar.gz'
