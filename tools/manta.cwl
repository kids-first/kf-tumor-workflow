cwlVersion: v1.0
class: CommandLineTool
id: kfdrc-manta-sv
label: Manta sv caller
doc: 'Calls structural variants. Tool designed to pick correct run mode based on if tumor, normal, or both crams are given'
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: $(inputs.ram * 1000)
    coresMin: $(inputs.cores)
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/manta:1.4.0'

baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      /manta-1.4.0.centos6_x86_64/bin/configManta.py --runDir=./
  - position: 9
    shellQuote: false
    valueFrom: >-
      $(inputs.input_normal_aligned != null ? inputs.input_tumor_aligned != null ? "--normalBam " + inputs.input_normal_aligned.path : "--bam " + inputs.input_normal_aligned.path : "")
  - position: 10
    shellQuote: false
    prefix: "&&"
    valueFrom: >-
      ./runWorkflow.py -m local
  - position: 20
    shellQuote: false
    valueFrom: >-
      && mv results/variants/diploidSV.vcf.gz $([inputs.output_basename, inputs.tool_name, "diploidSV.vcf.gz"].join(".")) || :
      && mv results/variants/diploidSV.vcf.gz.tbi $([inputs.output_basename, inputs.tool_name, "diploidSV.vcf.gz.tbi"].join(".")) || :
      && mv results/variants/tumorSV.vcf.gz $([inputs.output_basename, inputs.tool_name, "tumorSV.vcf.gz"].join(".")) || :
      && mv results/variants/tumorSV.vcf.gz.tbi $([inputs.output_basename, inputs.tool_name, "tumorSV.vcf.gz.tbi"].join(".")) || :
      && mv results/variants/somaticSV.vcf.gz $([inputs.output_basename, inputs.tool_name, "somaticSV.vcf.gz"].join(".")) || :
      && mv results/variants/somaticSV.vcf.gz.tbi $([inputs.output_basename, inputs.tool_name, "somaticSV.vcf.gz.tbi"].join(".")) || :
      && mv results/variants/candidateSmallIndels.vcf.gz $([inputs.output_basename, inputs.tool_name, "candidateSmallIndels.vcf.gz"].join("."))
      && mv results/variants/candidateSmallIndels.vcf.gz.tbi $([inputs.output_basename, inputs.tool_name, "candidateSmallIndels.vcf.gz.tbi"].join("."))

inputs:
    reference: {type: File, secondaryFiles: [^.dict, .fai], inputBinding: {position: 2, prefix: "--ref"}}
    hg38_strelka_bed: {type: File, secondaryFiles: [.tbi], inputBinding: {position: 2, prefix: "--callRegions"}}
    input_tumor_aligned: { type: 'File?', secondaryFiles: ["^.bai?", ".bai?", "^.crai?", ".crai?"], inputBinding: {position: 9, prefix: "--tumorBam"},  doc: "tumor BAM or CRAM" }
    input_normal_aligned: { type: 'File?', secondaryFiles: ["^.bai?", ".bai?", "^.crai?", ".crai?"], doc: "normal BAM or CRAM" }
    cores: {type: ['null', int], default: 16, inputBinding: {position: 12, prefix: "-j"}}
    ram: {type: "int?", default: 10, doc: "GB of RAM an instance must have to run the task"}
    output_basename: string
    tool_name: { type: 'string?', default: "manta", doc: "Tool name to use in outputs." }
outputs:
  output_sv:
    type: File
    outputBinding:
      glob: '*SV.vcf.gz'
    secondaryFiles: [.tbi]
  small_indels:
    type: File
    outputBinding:
      glob: '*SmallIndels.vcf.gz'
    secondaryFiles: [.tbi]
