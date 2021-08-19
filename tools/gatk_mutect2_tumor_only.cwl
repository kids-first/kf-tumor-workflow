cwlVersion: v1.2
class: CommandLineTool
id: gatk_mutect2_tumor_only
label: GATK Mutect2 Tumor Only Mode

requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/gatk:4.1.1.0'
  - class: ResourceRequirement
    ramMin: ${return inputs.ram * 1000}
    coresMin: ${return inputs.cpu}

baseCommand: [/gatk, Mutect2]

arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      --java-options ${var str = "-Xmx"; return str.concat(inputs.ram, "g")}
      -O ${return (inputs.input_bcram).basename.concat(".vcf.gz")}

inputs:
  reference: {type: "File", secondaryFiles: [^.dict, .fai], inputBinding: {prefix: -R}}
  input_bcram: {type: "File", secondaryFiles: ["^.bai?", ".bai?", "^.crai?", ".crai?"], inputBinding: {prefix: -I}, doc: "Input bam / cram file"}
  mnp_distance: {type: "int?", inputBinding: {prefix: --max-mnp-distance}, doc: "Max distance to merge substitutions into one MNP"}
  panel_of_normals: {type: "File?", secondaryFiles: [".tbi"], inputBinding: {prefix: --panel-of-normals}, doc: "Panel of normal vcf file"}
  ram: {type: "int?", default: 2, doc: "In GB"}
  cpu: {type: "int?", default: 4, doc: "Number of CPUs to request"}

outputs:
  mutect2_output_vcf:
    type: File
    outputBinding:
      glob: '*.vcf.gz'
