cwlVersion: v1.0
class: CommandLineTool
id: gatk_mergevcfs
doc: "Merge input vcfs"
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/gatk:4.1.1.0'
  - class: ResourceRequirement
    ramMin: $(inputs.max_memory * 1000)
    coresMin: $(inputs.cores)
baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      /gatk --java-options "-Xmx${return Math.floor(inputs.max_memory*1000/1.074-1)}m" MergeVcfs
      --TMP_DIR=./TMP
      --CREATE_INDEX=true
      --SEQUENCE_DICTIONARY=$(inputs.reference_dict.path)
      ${
        var cmd = "--OUTPUT=" + inputs.output_basename + "." + inputs.tool_name + ".merged.vcf.gz "
        if (typeof inputs.silent_flag !== 'undefined' && inputs.silent_flag == 1){
          cmd += "--VALIDATION_STRINGENCY SILENT"
        }
        return cmd
      }

inputs:
  input_vcfs:
    type:
      type: array
      items: File
      inputBinding:
        prefix: -I
    secondaryFiles: [.tbi]
    inputBinding:
      position: 1
  reference_dict: File
  tool_name: string
  output_basename: string
  silent_flag: { type: 'int?' }
  max_memory: { type: 'int?', default: 4, doc: "Maximum memory in GB to allocate to this task" }
  cores: { type: 'int?', default: 2, doc: "CPUs to allocate to this task" }
outputs:
  merged_vcf:
    type: File
    outputBinding:
      glob: '*.merged.vcf.gz'
    secondaryFiles: [.tbi]
