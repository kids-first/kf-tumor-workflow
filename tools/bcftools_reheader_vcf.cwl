cwlVersion: v1.0
class: CommandLineTool
id: bcftools_reheader_vcf
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/bvcftools:latest'
  - class: ResourceRequirement
    ramMin: 1000
    coresMin: 1
  - class: InlineJavascriptRequirement
baseCommand: []
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      ${
        if (!inputs.input_normal_name && !inputs.input_tumor_name){
          var err = '>&2 echo "Need to give at least one of normal name or tumor name"; exit 1;'
          return err;
        }
        var cmd = "";
        if (inputs.input_normal_name){
          cmd += 'echo ' + inputs.input_normal_name + ' >> sample_list.txt;'
        }
        if (inputs.input_tumor_name){
          cmd += 'echo ' + inputs.input_tumor_name + ' >> sample_list.txt;'
        }
        return cmd;
      }
      bcftools reheader -s sample_list.txt $(inputs.input_vcf.path) > $(inputs.input_vcf.nameroot.replace(".vcf", ".reheadered.vcf.gz")) &&
      tabix $(inputs.input_vcf.nameroot.replace(".vcf", ".reheadered.vcf.gz"))

inputs:
  input_vcf: File
  input_normal_name: 'string?'
  input_tumor_name: 'string?'
outputs:
  reheadered_vcf:
    type: File
    outputBinding:
      glob: '*.vcf.gz'
    secondaryFiles: [.tbi]
  sample_list:
    type: File
    outputBinding:
      glob: 'sample_list.txt'
