cwlVersion: v1.2
class: Workflow
id: gatk_mutect2_generate_pon
label: GATK Mutect2 Generate Panel of Normals
doc: |-
  # GATK Mutect2 Generate Panel of Normals Workflow
  Workflow for generating a Mutect2 Panel of Normals for use in the Tumor Only Workflow

  ![data service logo](https://github.com/d3b-center/d3b-research-workflows/raw/master/doc/kfdrc-logo-sm.png)

requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement

inputs:
  reference: {type: "File", secondaryFiles: [^.dict, .fai]}
  input_bcram: {type: "File[]", secondaryFiles: ["^.bai?", ".bai?", "^.crai?", ".crai?"], inputBinding: {prefix: -I}, doc: "Input bam / cram file"}
  interval_list: {type: "File?", doc: "GATK intervals list-style, or bed file.  Recommend canocical chromosomes with N regions removed"}
  output_basename: {type: string, doc: "Output panel of normal vcf basename"}
  cpus: { type: 'int?', default: 4, doc: "CPUs to allocate to call task"}
  ram: { type: 'int?', default: 8, doc: "RAM to allocate to call task in gb"}

outputs:
  panel_of_normal_vcf: {type: 'File', outputSource: create_pon/panel_of_normal_vcf}

steps:

  mutect2_tumor_only:
    run: ../tools/gatk_mutect2_tumor_only.cwl
    scatter: input_bcram
    in:
      reference: reference
      input_bcram: input_bcram
      mnp_distance: 0
      ram: ram
      cpus: cpus
    out: [mutect2_output_vcf]

  create_pon:
    run: ../tools/rnaseqc.cwl
    in:
      reference: reference
      input_vcfs: mutect2_tumor_only/mutect2_output_vcf
      interval_list: interval_list
      output_basename: output_basename
      ram: ram
      cpus: cpus
    out: [panel_of_normal_vcf]

$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 2
