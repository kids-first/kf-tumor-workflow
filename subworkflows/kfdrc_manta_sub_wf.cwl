cwlVersion: v1.2
class: Workflow
id: kfdrc_manta_tumor_only_sub_wf
label: KFDRC Manta Tumor Only
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement

inputs:
  indexed_reference_fasta: {type: File, secondaryFiles: [.fai, ^.dict]}
  hg38_strelka_bed: {type: File, secondaryFiles: ['.tbi']}
  input_tumor_aligned: { type: File, secondaryFiles: [{ pattern: ".bai", required: false },{ pattern: "^.bai", required: false },{ pattern: ".crai", required: false },{ pattern: "^.crai", required: false }], doc: "tumor BAM or CRAM" }
  input_tumor_name: string
  output_basename: string
  manta_memory: {type: "int?"}
  manta_cores: {type: "int?"}
  select_vars_mode: {type: ['null', {type: enum, name: select_vars_mode, symbols: ["gatk", "grep"]}], doc: "Choose 'gatk' for SelectVariants tool, or 'grep' for grep expression", default: "gatk"}
  annotsv_annotations_dir_tgz: {type: 'File', doc: "TAR.GZ'd Directory containing annotations for AnnotSV"}

outputs:
  manta_prepass_vcf: {type: File, outputSource: rename_manta_samples/reheadered_vcf}
  manta_pass_vcf: {type: File, outputSource: gatk_selectvariants_manta/pass_vcf}
  manta_small_indels: {type: File, outputSource: manta/small_indels}
  annotsv_annotated_calls: {type: 'File?', outputSource: annotsv/annotated_calls}
  annotsv_unannotated_calls: {type: 'File?', outputSource: annotsv/unannotated_calls}

steps:
  manta:
    run: ../tools/manta.cwl
    in:
      input_tumor_aligned: input_tumor_aligned
      output_basename: output_basename
      ram: manta_memory
      cores: manta_cores
      reference: indexed_reference_fasta
      hg38_strelka_bed: hg38_strelka_bed
    out: [output_sv, small_indels]

  rename_manta_samples:
    run: ../tools/bcftools_reheader_vcf.cwl
    in:
      input_vcf: manta/output_sv
      input_tumor_name: input_tumor_name
    out: [reheadered_vcf]

  gatk_selectvariants_manta:
    run: ../tools/gatk_selectvariants.cwl
    label: GATK Select Manta PASS
    in:
      input_vcf: rename_manta_samples/reheadered_vcf
      output_basename: output_basename
      tool_name:
        valueFrom: $("manta")
      mode: select_vars_mode
    out: [pass_vcf]

  annotsv:
    run: ../tools/annotsv.cwl
    in:
      annotations_dir_tgz: annotsv_annotations_dir_tgz
      sv_input_file: gatk_selectvariants_manta/pass_vcf
    out: [annotated_calls, unannotated_calls]
