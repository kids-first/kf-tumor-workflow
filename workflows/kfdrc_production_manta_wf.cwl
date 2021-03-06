cwlVersion: v1.2
class: Workflow
id: kfdrc_production_manta_wf
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement
inputs:
  # Required
  reference_fasta: { type: File }
  reference_fai: { type: 'File?' }
  reference_dict: { type: 'File?' }
  input_tumor_aligned: { type: File, secondaryFiles: ["^.bai?", ".bai?", "^.crai?", ".crai?"], doc: "tumor BAM or CRAM" }
  input_tumor_name: string
  hg38_strelka_bed: { type: File, doc: "Bgzipped interval bed file. Recommned padding 100bp for WXS; Recommend canonical chromosomes for WGS" }
  hg38_strelka_tbi: { type: 'File?', doc: "Tabix index for hg38_strelka_bed" }
  output_basename: { type: string, doc: "String value to use as basename for outputs" }

  # Optional with One Default
  select_vars_mode: { type: ['null', { type: enum, name: select_vars_mode, symbols: ["gatk", "grep"] }], default: "gatk", doc: "Choose 'gatk' for SelectVariants tool, or 'grep' for grep expression" }
  manta_memory: {type: 'int?', doc: "GB of memory to allocate to Manta; defaults to 10 (soft-capped)"}
  manta_cores: {type: 'int?', doc: "Number of cores to allocate to Manta; defaults to 18"}

outputs:
  manta_pass_vcf: { type: File, outputSource: run_manta/manta_pass_vcf }
  manta_prepass_vcf: { type: File, outputSource: run_manta/manta_prepass_vcf }
  manta_small_indels: { type: File, outputSource: run_manta/manta_small_indels }

steps:
  prepare_reference:
    run: ../subworkflows/prepare_reference.cwl
    in:
      input_fasta: reference_fasta
      input_fai: reference_fai
      input_dict: reference_dict
    out: [indexed_fasta,reference_dict]

  index_strelka_bed:
    run: ../tools/tabix_index.cwl
    in:
      input_file: hg38_strelka_bed
      input_index: hg38_strelka_tbi
    out: [output]

  run_manta:
    run: ../subworkflows/kfdrc_manta_sub_wf.cwl
    in:
      indexed_reference_fasta: prepare_reference/indexed_fasta
      hg38_strelka_bed: index_strelka_bed/output
      input_tumor_aligned: input_tumor_aligned
      input_tumor_name: input_tumor_name
      output_basename: output_basename
      manta_memory: manta_memory
      manta_cores: manta_cores
    out:
      [manta_prepass_vcf, manta_pass_vcf, manta_small_indels]

$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 2
