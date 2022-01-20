cwlVersion: v1.0
class: CommandLineTool
id: gatk_variantfiltration
doc: "Simple tool add a FILTER tag based on criteria"
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: $(inputs.max_memory * 1000)
    coresMin: $(inputs.cores)
  - class: DockerRequirement
    dockerPull: 'broadinstitute/gatk:4.2.2.0'

baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      gatk --java-options "-Xmx${return Math.floor(inputs.max_memory*1000/1.074-1)}m" VariantFiltration
      -R $(inputs.reference.path)
      -V $(inputs.input_vcf.path)
      -O $(inputs.output_basename).$(inputs.tool_name).gatk.soft_filtered.vcf.gz
      ${
        var args = "";
        for (var i = 0; i < inputs.filter_name.length; i++){
          args += "--filter-name \"" + inputs.filter_name[i] + "\" --filter-expression \"" + inputs.filter_expression[i] + "\" ";
        }
        return args
      }

inputs:
  input_vcf: {type: 'File', secondaryFiles: ['.tbi']}
  reference: {type: 'File', secondaryFiles: [^.dict, .fai]}
  filter_name: {type: 'string[]', doc: "Array of names for each filter tag to add"}
  filter_expression: {type: 'string[]', doc: "Array of filter expressions to establish criteria to tag variants with. See https://gatk.broadinstitute.org/hc/en-us/articles/360036730071-VariantFiltration for clues"}
  output_basename: string
  tool_name: string
  max_memory: { type: 'int?', default: 8, doc: "Maximum memory in GB to allocate to this task" }
  cores: { type: 'int?', default: 4, doc: "CPUs to allocate to this task" }

outputs:
  gatk_soft_filtered_vcf:
    type: File
    outputBinding:
      glob: '*.vcf.gz'
    secondaryFiles: ['.tbi']
