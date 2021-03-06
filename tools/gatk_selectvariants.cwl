cwlVersion: v1.0
class: CommandLineTool
id: gatk_selectvariants
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: 'broadinstitute/gatk:4.2.2.0'
  - class: ResourceRequirement
    ramMin: $(inputs.max_memory * 1000)
    coresMin: $(inputs.cores)
baseCommand: ["/bin/bash", "-c"]
arguments:
  - position: 0
    shellQuote: true
    valueFrom: >-
      set -eo pipefail

      ${
        var run_mode = inputs.mode;
        if (run_mode == 'grep' || run_mode == 'gatk'){
          var in_vcf = inputs.input_vcf.path;
          var out_vcf = inputs.output_basename + '.' + inputs.tool_name + '.PASS.vcf.gz';
          var cmd = 'gatk SelectVariants --java-options "-Xmx' + Math.floor(inputs.max_memory*1000/1.074-1) + 'm" -V ' + in_vcf +  ' -O ' + out_vcf + ' --exclude-filtered TRUE';
          if (run_mode == 'grep'){
            cmd = 'zcat ' + in_vcf + ' | grep -E "^#|PASS" | bgzip > ' + out_vcf + '; tabix ' + out_vcf;
          }
          return cmd;
        }
        else{
          throw new Error(run_mode + ' is a not a valid mode.  Choices are gatk or grep.');
        }
      }

inputs:
  input_vcf: {type: File, secondaryFiles: [.tbi]}
  output_basename: string
  tool_name: string
  mode: {type: ['null', {type: enum, name: select_vars_mode, symbols: ["gatk", "grep"]}], doc: "Choose 'gatk' for SelectVariants tool, or 'grep' for grep expression", default: "gatk"}
  enable_tool: {type: 'boolean?'}
  max_memory: { type: 'int?', default: 8, doc: "Maximum memory in GB to allocate to this task" }
  cores: { type: 'int?', default: 4, doc: "CPUs to allocate to this task" }

outputs:
  pass_vcf:
    type: File
    outputBinding:
      glob: '*.PASS.vcf.gz'
    secondaryFiles: ['.tbi']
