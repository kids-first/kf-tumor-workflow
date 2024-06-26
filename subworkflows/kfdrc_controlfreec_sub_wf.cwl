cwlVersion: v1.2
class: Workflow
id: kfdrc_controlfreec_tumor_only_sub_wf
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement

inputs:
  wgs_or_wxs: {type: {type: enum, name: wgs_or_wxs, symbols: ["WGS", "WXS"]}, doc: "Select if this run is WGS or WXS"}
  input_tumor_aligned: {type: File, secondaryFiles: ['^.bai']}
  input_tumor_name: {type: string, doc: "Sample name to put into the converted seg file"}
  threads: {type: int, doc: "Number of threads to run controlfreec.  Going above 16 is not recommended, there is no apparent added value"}
  output_basename: {type: string, doc: "basename to use for naming output files."}
  ploidy: {type: 'int[]', doc: "Array of ploidy possibilities for ControlFreeC to try"}
  mate_copynumber_file_sample: {type: 'File?', doc: "Tumor cpn file from previous run. If used, will override bam use"}
  gem_mappability_file: {type: 'File?', doc: "GEM mappability file to make read count adjustments with"}
  min_subclone_presence: {type: 'int?', doc: "Tool default 100 (meaning \"do not look for subclones\"). Suggested: 20 (or 0.2) for WGS and 30 (or 0.3) for WES."}
  mate_orientation_sample: {type: ['null', {type: enum, name: mate_orientation_sample, symbols: ["0", "FR", "RF", "FF"]}], default: "FR", doc: "0 (for single ends), RF (Illumina mate-pairs), FR (Illumina paired-ends), FF (SOLiD mate-pairs)"}
  calling_regions: { type: File, doc: "BED containing the regions over which Control-FREEC should be run. For WXS, these are the bait intervals." }
  indexed_reference_fasta: {type: File, secondaryFiles: [.fai]}
  b_allele: {type: ['null', File], doc: "germline calls, needed for BAF.  VarDict input recommended.  Tool will prefilter for germline and pass if expression given"}
  coeff_var: {type: float, default: 0.05, doc: "Coefficient of variantion to set window size.  Default 0.05 recommended"}
  cfree_sex: {type: ['null', {type: enum, name: sex, symbols: ["XX", "XY"] }], doc: "If known, XX for female, XY for male"}
  tool_name: { type: 'string?', doc: "Tool name to use in outputs." }

outputs:
  ctrlfreec_cnvs: {type: File, outputSource: rename_outputs/ctrlfreec_cnvs}
  ctrlfreec_pval: {type: File, outputSource: rename_outputs/ctrlfreec_pval}
  ctrlfreec_config: {type: File, outputSource: rename_outputs/ctrlfreec_config}
  ctrlfreec_pngs: {type: 'File[]', outputSource: rename_outputs/ctrlfreec_pngs}
  ctrlfreec_bam_ratio: {type: File, outputSource: rename_outputs/ctrlfreec_bam_ratio}
  ctrlfreec_bam_seg: {type: File, outputSource: convert_ratio_to_seg/ctrlfreec_ratio2seg}
  ctrlfreec_baf: {type: 'File?', outputSource: rename_outputs/ctrlfreec_baf}
  ctrlfreec_info: {type: File, outputSource: rename_outputs/ctrlfreec_info}

steps:
  awk_chrlen_builder:
    run: ../tools/awk_chrlen_builder.cwl
    in:
      input_intervals: calling_regions
      reference_fai:
        source: indexed_reference_fasta
        valueFrom: |
          $(self.secondaryFiles.filter(function(e) { return e.basename.search(/.fai$/) != -1 })[0])
    out: [chrlen]

  samtools_mpileup:
    run: ../tools/samtools_mpileup.cwl
    in:
      input_reads: input_tumor_aligned
      reference: indexed_reference_fasta
      snp_vcf: b_allele
      calling_regions: awk_chrlen_builder/chrlen
    out: [pileup]

  control_free_c:
    run: ../tools/control-freec-11-6-sbg.cwl
    in:
      mate_copynumber_file_sample: mate_copynumber_file_sample
      gem_mappability_file: gem_mappability_file
      min_subclone_presence: min_subclone_presence
      mate_file_sample: input_tumor_aligned
      mate_orientation_sample: mate_orientation_sample
      mini_pileup_sample: samtools_mpileup/pileup
      chr_len: awk_chrlen_builder/chrlen
      ploidy: ploidy
      capture_regions:
        source: [wgs_or_wxs, calling_regions]
        valueFrom: |
          $(self[0] == 'WXS' ? self[1] : null)
      max_threads: threads
      reference: indexed_reference_fasta
      snp_file: b_allele
      coeff_var: coeff_var
      sex: cfree_sex
    out: [cnvs, cnvs_pvalue, config_script, pngs, ratio, sample_BAF, info_txt]

  rename_outputs:
    run: ../tools/ubuntu_rename_outputs.cwl
    in:
      input_files: [control_free_c/cnvs, control_free_c/cnvs_pvalue, control_free_c/config_script, control_free_c/ratio, control_free_c/sample_BAF, control_free_c/info_txt]
      input_pngs: control_free_c/pngs
      tool_name: tool_name
      output_basename: output_basename
    out: [ctrlfreec_cnvs, ctrlfreec_pval, ctrlfreec_config, ctrlfreec_pngs, ctrlfreec_bam_ratio, ctrlfreec_baf, ctrlfreec_info]

  convert_ratio_to_seg:
    run: ../tools/ubuntu_ratio2seg.cwl
    in:
      reference_fai:
        source: indexed_reference_fasta
        valueFrom: |
          $(self.secondaryFiles.filter(function(e) { return e.basename.search(/.fai$/) != -1 })[0])
      ctrlfreec_ratio: control_free_c/ratio
      sample_name: input_tumor_name
      output_basename: output_basename
      tool_name: tool_name
    out: [ctrlfreec_ratio2seg]
