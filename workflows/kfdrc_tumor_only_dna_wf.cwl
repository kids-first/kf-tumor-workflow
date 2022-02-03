cwlVersion: v1.2
class: Workflow
id: kfdrc_tomor_only_dna_wf
label: KFDRC DNA Tumor Only Beta Production Workflow
doc: "Beta workflow for tumor-only samples"
requirements:
- class: ScatterFeatureRequirement
- class: MultipleInputFeatureRequirement
- class: SubworkflowFeatureRequirement
- class: InlineJavascriptRequirement

inputs:
  # Shared
  reference_fasta: {type: 'File', doc: "Reference fasta", "sbg:suggestedValue": {
      class: File, path: 60639014357c3a53540ca7a3, name: Homo_sapiens_assembly38.fasta},
    "sbg:fileTypes": "FASTA, FA"}
  reference_fai: {type: 'File?', doc: "FAI index for reference_fasta", "sbg:suggestedValue": {
      class: File, path: 60639016357c3a53540ca7af, name: Homo_sapiens_assembly38.fasta.fai},
    "sbg:fileTypes": "FAI"}
  reference_dict: {type: 'File?', doc: "DICT index for reference_fasta", "sbg:suggestedValue": {
      class: File, path: 60639019357c3a53540ca7e7, name: Homo_sapiens_assembly38.dict},
    "sbg:fileTypes": "DICT"}
  input_tumor_aligned: {type: 'File', secondaryFiles: [{pattern: ".bai", required: false},
      {pattern: "^.bai", required: false}, {pattern: ".crai", required: false}, {
        pattern: "^.crai", required: false}], doc: "BAM/SAM/CRAM file containing tumor\
      \ reads", "sbg:fileTypes": "BAM, CRAM, SAM"}
  input_tumor_name: {type: 'string', doc: "BAM sample name of tumor. May be URL-encoded\
      \ as output by GetSampleName with -encode argument."}
  output_basename: {type: 'string', doc: "String value to use as basename for outputs"}

  # GATK Files
  panel_of_normals: {type: 'File?', secondaryFiles: ['.tbi'], doc: "VCF file (and\
      \ index) of sites observed in normal. A panel of normals can be a useful (optional)\
      \ input to help filter out commonly seen sequencing noise that may appear as\
      \ low allele-fraction somatic variants.", "sbg:fileTypes": "VCF, VCF.GZ"}
  mutect2_af_only_gnomad_vcf: {type: 'File', secondaryFiles: ['.tbi'], doc: "Population\
      \ vcf (and index) of germline sequencing containing allele fractions in VCF\
      \ format.", "sbg:suggestedValue": {class: File, path: 5f50018fe4b054958bc8d2e3,
      name: af-only-gnomad.hg38.vcf.gz, secondaryFiles: [{class: File, path: 5f50018fe4b054958bc8d2e5,
          name: af-only-gnomad.hg38.vcf.gz.tbi}]}, "sbg:fileTypes": "VCF, VCF.GZ"}
  mutect2_alleles_vcf: {type: 'File?', secondaryFiles: ['.tbi'], doc: "VCF file (and\
      \ index) containing set of alleles for which to force genotyping regardless\
      \ of evidence", "sbg:fileTypes": "VCF, VCF.GZ"}
  mutect2_exac_common_vcf: {type: 'File?', secondaryFiles: ['.tbi'], doc: "Exac Common\
      \ VCF (and index) used for calculating contamination values used in filtering\
      \ the Mutect VCF. If do not wish to perfom this filtering, remove this input.",
    "sbg:suggestedValue": {class: File, path: 5f500135e4b0370371c051ad, name: small_exac_common_3.hg38.vcf.gz,
      secondaryFiles: [{class: File, path: 5f500135e4b0370371c051af, name: small_exac_common_3.hg38.vcf.gz.tbi}]},
    "sbg:fileTypes": "VCF, VCF.GZ"}
  bwa_mem_index_image: {type: 'File?', doc: "BWA-mem index image used for FilterAlignmentArtifacts.\
      \ WARNING!!! Experimental Software!!! Do not provide unless you are absolutely\
      \ certain you want to run FilterAlignmentArtifacts.", "sbg:fileTypes": "IMG"}

  # GATK Arguments
  tool_name: {type: 'string?', doc: "Tool name for this workflow", default: "mutect2"}
  make_bamout: {type: 'boolean?', default: false, doc: "Should Mutect2 create a BAM\
      \ output? Turning on this mode may result in serious performance cost for the\
      \ caller."}
  run_orientation_bias_mixture_model_filter: {type: 'boolean?', default: true, doc: "Should\
      \ Orientation Bias Mixture Model Filter be applied to Mutect2 outputs?"}
  mutect2_extra_args: {type: 'string?', doc: "Any additional arguments for Mutect2.\
      \ See GATK documentation for all available options."}
  filtermutectcalls_extra_args: {type: 'string?', doc: "Any additional arguments for\
      \ FilterMutectCalls. See GATK documentation for all available options."}
  select_vars_mode: {type: ['null', {type: enum, name: select_vars_mode, symbols: [
          "gatk", "grep"]}], default: "gatk", doc: "Choose 'gatk' for SelectVariants\
      \ tool, or 'grep' for grep expression"}
  scatter_count: {type: 'int?', default: 50, doc: "The number of files into which\
      \ to scatter the resulting list by locus; in some situations, fewer intervals\
      \ may be emitted"}

  # WGS or WXS
  wgs_or_wxs: {type: {type: enum, name: wgs_or_wxs, symbols: ["WGS", "WXS"]}, doc: "Select\
      \ input type WGS or WXS"}
  wgs_calling_interval_list: {type: 'File?', doc: "GATK intervals list-style, or bed\
      \ file.  Recommend canocical chromosomes with N regions removed", "sbg:suggestedValue": {
      class: File, path: 5f500135e4b0370371c051b6, name: wgs_canonical_calling_regions.hg38.bed},
    "sbg:fileTypes": "BED, INTERVALS, INTERVAL_LIST, LIST"}
  padded_capture_regions: {type: 'File?', doc: "Recommend 100bp pad, for somatic variant.\
      \ Used as calling intervals for exomes.", "sbg:fileTypes": "BED, INTERVALS,\
      \ INTERVAL_LIST, LIST"}

  # ControlFreeC CNV
  mate_copynumber_file_sample: {type: 'File?', doc: "Tumor cpn file from previous run. If used, will override bam use"}
  gem_mappability_file: {type: 'File?', doc: "GEM mappability file to make read count adjustments with"}
  min_subclone_presence: {type: 'float?', doc: "Use if you want to detect subclones. Recommend 0.2 for WGS, 0.3 for WXS"}
  cfree_chr_len: { type: File, doc: "file with chromosome lengths" }
  cfree_ploidy: { type: 'int[]', doc: "Array of ploidy possibilities for ControlFreeC to try" }
  cfree_mate_orientation_sample: { type: ['null', { type: enum, name: mate_orientation_sample, symbols: ["0", "FR", "RF", "FF"] }], default: "FR", doc: "0 (for single ends), RF (Illumina mate-pairs), FR (Illumina paired-ends), FF (SOLiD mate-pairs)" }

  # Optional with Multiple Defaults (handled in choose_defaults)
  i_flag: { type: 'string?', doc: "Flag to intersect germline calls on padded regions. Use N if you want to skip this or have a WGS run" }

  # Optional
  b_allele: { type: 'File?', doc: "germline calls, needed for BAF.  GATK HC VQSR input recommended. Tool will prefilter for germline and pass if expression given" }
  b_allele_index: { type: 'File?', doc: "Tabix index for b_allele" }
  cfree_coeff_var: { type: 'float?', default: 0.05, doc: "Coefficient of variation to set window size. Default 0.05 recommended" }
  cfree_sex: { type: ['null', { type: enum, name: sex, symbols: ["XX", "XY"] }], doc: "If known, XX for female, XY for male" }

  # Manta SV
  hg38_strelka_bed: { type: File, doc: "Bgzipped interval bed file. Recommned padding 100bp for WXS; Recommend canonical chromosomes for WGS" }
  hg38_strelka_tbi: { type: 'File?', doc: "Tabix index for hg38_strelka_bed" }

  # WXS only Fields
  unpadded_capture_regions: { type: 'File?', doc: "Capture regions with NO padding for cnv calling" }


  # Annotation Arguments
  run_annotation: {type: 'boolean?', default: true, doc: "Should annotation be run\
      \ on the final VCF?"}
  vep_cache: {type: 'File?', doc: "tar gzipped cache from ensembl/local converted\
      \ cache", "sbg:suggestedValue": {class: File, path: 607713829360f10e3982a425,
      name: homo_sapiens_vep_93_GRCh38.tar.gz}, "sbg:fileTypes": "TAR,\
      \ TAR.GZ, TGZ"}
  vep_ref_build: {type: 'string?', default: "GRCh38", doc: "Genome ref build used,\
      \ should line up with cache"}
  genomic_hotspots: {type: 'File[]?', doc: "Tab-delimited BED formatted file(s) containing\
      \ hg38 genomic positions corresponding to hotspots", "sbg:suggestedValue": [{
        class: File, path: 607713829360f10e3982a423, name: tert.bed}], "sbg:fileTypes": "BED,\
      \ TSV"}
  protein_snv_hotspots: {type: 'File[]?', doc: "Column-name-containing, tab-delimited\
      \ file(s) containing protein names and amino acid positions corresponding to\
      \ hotspots", "sbg:suggestedValue": [{class: File, path: 607713829360f10e3982a426,
        name: protein_snv_cancer_hotspots_v2.tsv}], "sbg:fileTypes": "BED,\
      \ TSV"}
  protein_indel_hotspots: {type: 'File[]?', doc: "Column-name-containing, tab-delimited\
      \ file(s) containing protein names and amino acid position ranges corresponding\
      \ to hotspots", "sbg:suggestedValue": [{class: File, path: 607713829360f10e3982a424,
        name: protein_indel_cancer_hotspots_v2.tsv}], "sbg:fileTypes": "BED,\
      \ TSV"}
  retain_info: {type: 'string?', doc: "csv string with INFO fields that you want to\
      \ keep", default: "MBQ,TLOD,HotSpotAllele"}
  retain_fmt: {type: 'string?', doc: "csv string with FORMAT fields that you want\
      \ to keep"}
  add_common_fields: {type: 'boolean?', doc: "Set to true if input is a strelka2 vcf\
      \ that hasn't had common fields added", default: false}
  bcftools_annot_columns: {type: 'string?', doc: "csv string of columns from annotation\
      \ to port into the input vcf, i.e INFO/AF", default: "INFO/AF"}
  bcftools_annot_vcf: {type: 'File?', secondaryFiles: [.tbi], doc: "bgzipped annotation\
      \ vcf file", "sbg:suggestedValue": {class: File, path: 5f50018fe4b054958bc8d2e3,
      name: af-only-gnomad.hg38.vcf.gz, secondaryFiles: [{class: File, path: 5f50018fe4b054958bc8d2e5,
          name: af-only-gnomad.hg38.vcf.gz.tbi}]}, "sbg:fileTypes": "VCF, VCF.GZ"}
  bcftools_public_filter: {type: 'string?', doc: "Will hard filter final result to\
      \ create a public version", default: FILTER="PASS"|INFO/HotSpotAllele=1}
  gatk_filter_name: {type: 'string[]?', doc: "Array of names for each filter tag to\
      \ add, recommend: [\"NORM_DP_LOW\", \"GNOMAD_AF_HIGH\"]"}
  gatk_filter_expression: {type: 'string[]?', doc: "Array of filter expressions to\
      \ establish criteria to tag variants with. See https://gatk.broadinstitute.org/hc/en-us/articles/360036730071-VariantFiltration,\
      \ recommend: \"vc.getGenotype('\" + inputs.input_normal_name + \"').getDP()\
      \ <= 7\"), \"AF > 0.001\"]"}
  disable_hotspot_annotation: {type: 'boolean?', doc: "Disable Hotspot Annotation\
      \ and skip this task.", default: false}
  maf_center: {type: 'string?', doc: "Sequencing center of variant called", default: "."}

  # Resource Control
  mutect_cores: {type: 'int?', doc: "CPUs to allocate to GATK Mutect2"}
  mutect_memory: {type: 'int?', doc: "GB of memory to allocate to GATK Mutect2 (hard-capped)"}
  getpileup_memory: {type: 'int?', doc: "GB of memory to allocate to GATK GetPileupSummaries\
      \ (hard-capped)"}
  learnorientation_memory: {type: 'int?', doc: "GB of memory to allocate to GATK LearnReadOrientationModel\
      \ (hard-capped)"}
  filtermutectcalls_memory: {type: 'int?', doc: "GB of memory to allocate to GATK\
      \ FilterMutectCalls (hard-capped)"}
  filteralignmentartifacts_cores: {type: 'int?', doc: "CPUs to allocate to GATK FilterAlignmentArtifacts"}
  filteralignmentartifacts_memory: {type: 'int?', doc: "GB of memory to allocate to\
      \ GATK FilterAlignmentArtifacts (hard-capped)"}
  cfree_threads: { type: 'int?', default: 16, doc: "For ControlFreeC. Recommend 16 max, as I/O gets saturated after that losing any advantage" }
  manta_memory: {type: 'int?', doc: "GB of memory to allocate to Manta; defaults to 10 (soft-capped)"}
  manta_cores: {type: 'int?', doc: "Number of cores to allocate to Manta; defaults to 18"}


outputs:
  # Mutect2
  mutect2_prepass_vcf: {type: 'File', outputSource: run_mutect2/mutect2_filtered_vcf,
    doc: "VCF with SNV, MNV, and INDEL variant calls."}
  mutect2_protected_outputs: {type: 'File[]?', outputSource: run_mutect2/mutect2_protected_outputs,
    doc: "Array of files containing MAF format of PASS hits, PASS VCF with annotation\
      \ pipeline soft FILTER-added values, and VCF index"}
  mutect2_public_outputs: {type: 'File[]?', outputSource: run_mutect2/mutect2_public_outputs,
    doc: "Protected outputs, except MAF and VCF have had entries with soft FILTER\
      \ values removed"}
  mutect2_bam: {type: 'File?', outputSource: run_mutect2/mutect2_bam, doc: "BAM generated\
      \ by Mutect2, if desired. The assembled haplotypes and locally realigned reads\
      \ will be written as BAM. Really for debugging purposes only."}
  # ControlFreeC CNV
  ctrlfreec_pval: { type: File, outputSource: run_controlfreec/ctrlfreec_pval, doc: 'Copy number call with GT (if BAF provided) and p values. Most people want this' }
  ctrlfreec_config: { type: File, outputSource: run_controlfreec/ctrlfreec_config, doc: 'Config file used to run' }
  ctrlfreec_pngs: { type: 'File[]', outputSource: run_controlfreec/ctrlfreec_pngs, doc: 'Visualization of CN and BAF' }
  ctrlfreec_bam_ratio: { type: File, outputSource: run_controlfreec/ctrlfreec_bam_ratio, doc: 'Calls as log2 ratio' }
  ctrlfreec_bam_seg: { type: File, outputSource: run_controlfreec/ctrlfreec_bam_seg, doc: 'Custom made microarray-style seg file' }
  ctrlfreec_baf: { type: File, outputSource: run_controlfreec/ctrlfreec_baf, doc: 'b allele frequency file' }
  ctrlfreec_info: { type: File, outputSource: run_controlfreec/ctrlfreec_info, doc: 'Calculated inforamtion, like ploidy, if a range was given' }
  # Manta SV
  manta_pass_vcf: { type: File, outputSource: run_manta/manta_pass_vcf, doc: 'VCF file with SV calls that PASS' }
  manta_prepass_vcf: { type: File, outputSource: run_manta/manta_prepass_vcf, 'VCF file with all SV calls' }


steps:
  choose_defaults:
    run: ../tools/mode_defaults.cwl
    in:
      input_mode: wgs_or_wxs
      i_flag: i_flag
    out: [out_exome_flag, out_cnvkit_wgs_mode, out_i_flag, out_lancet_padding, out_lancet_window,
      out_vardict_padding]

  prepare_reference:
    run: ../subworkflows/prepare_reference.cwl
    in:
      input_fasta: reference_fasta
      input_fai: reference_fai
      input_dict: reference_dict
    out: [indexed_fasta, reference_dict]

  select_interval_list:
    run: ../tools/mode_selector.cwl
    in:
      input_mode: wgs_or_wxs
      wgs_input: wgs_calling_interval_list
      wxs_input: padded_capture_regions
    out: [output]

  gatk_intervallisttools:
    run: ../tools/gatk_intervallisttool.cwl
    in:
      interval_list: select_interval_list/output
      reference_dict: prepare_reference/reference_dict
      exome_flag: choose_defaults/out_exome_flag
      scatter_ct: scatter_count
      bands:
        valueFrom: ${return 80000000}
    out: [output]

  run_mutect2:
    hints:
    - class: "sbg:AWSInstanceType"
      value: c5.9xlarge
    run: ../subworkflows/kfdrc_mutect2_sub_wf.cwl
    in:
      indexed_reference_fasta: prepare_reference/indexed_fasta
      reference_dict: prepare_reference/reference_dict
      bed_invtl_split: gatk_intervallisttools/output
      af_only_gnomad_vcf: mutect2_af_only_gnomad_vcf
      alleles: mutect2_alleles_vcf
      exac_common_vcf: mutect2_exac_common_vcf
      input_tumor_aligned: input_tumor_aligned
      input_tumor_name: input_tumor_name
      input_normal_aligned: input_normal_aligned
      input_normal_name: input_normal_name
      panel_of_normals: panel_of_normals
      disable_adaptive_pruning:
        source: wgs_or_wxs
        valueFrom: $(self == 'WXS')
      mutect2_extra_args: mutect2_extra_args
      filtermutectcalls_extra_args: filtermutectcalls_extra_args
      vep_cache: vep_cache
      vep_ref_build: vep_ref_build
      select_vars_mode: select_vars_mode
      genomic_hotspots: genomic_hotspots
      protein_snv_hotspots: protein_snv_hotspots
      protein_indel_hotspots: protein_indel_hotspots
      retain_info: retain_info
      retain_fmt: retain_fmt
      add_common_fields: add_common_fields
      bcftools_annot_columns: bcftools_annot_columns
      bcftools_annot_vcf: bcftools_annot_vcf
      bcftools_public_filter: bcftools_public_filter
      gatk_filter_name: gatk_filter_name
      gatk_filter_expression: gatk_filter_expression
      maf_center: maf_center
      bwa_mem_index_image: bwa_mem_index_image
      make_bamout: make_bamout
      run_orientation_bias_mixture_model_filter: run_orientation_bias_mixture_model_filter
      run_annotation: run_annotation
      disable_hotspot_annotation: disable_hotspot_annotation
      mutect_cores: mutect_cores
      mutect_memory: mutect_memory
      getpileup_memory: getpileup_memory
      learnorientation_memory: learnorientation_memory
      filtermutectcalls_memory: filtermutectcalls_memory
      filteralignmentartifacts_cores: filteralignmentartifacts_cores
      filteralignmentartifacts_memory: filteralignmentartifacts_memory
      tool_name: tool_name
      output_basename: output_basename
    out: [mutect2_filtered_stats, mutect2_filtered_vcf, mutect2_protected_outputs,
      mutect2_public_outputs, mutect2_bam]

  index_b_allele:
    run: ../tools/tabix_index.cwl
    in:
      input_file: b_allele
      input_index: b_allele_index
    out: [output]

  bedtools_intersect_germline:
    run: ../tools/bedtools_intersect.cwl
    in:
      input_vcf: index_b_allele/output
      output_basename: output_basename
      input_bed_file: unpadded_capture_regions
      flag: choose_defaults/out_i_flag
    out:
      [intersected_vcf]

  gatk_filter_germline:
    run: ../tools/gatk_filter_germline_variant.cwl
    in:
      input_vcf: bedtools_intersect_germline/intersected_vcf
      reference_fasta: prepare_reference/indexed_fasta
      output_basename: output_basename
    out:
      [filtered_vcf, filtered_pass_vcf]

  samtools_cram2bam_plus_calmd_tumor:
    run: ../tools/samtools_calmd.cwl
    in:
      input_reads: input_tumor_aligned
      threads:
        valueFrom: ${return 16;}
      reference: prepare_reference/indexed_fasta
    out: [bam_file]

  run_controlfreec:
    run: ../subworkflows/kfdrc_controlfreec_sub_wf.cwl
    in:
      mate_copynumber_file_sample: mate_copynumber_file_sample
      gem_mappability_file: gem_mappability_file
      min_subclone_presence: min_subclone_presence
      input_tumor_aligned: samtools_cram2bam_plus_calmd_tumor/bam_file
      input_tumor_name: input_tumor_name
      threads: cfree_threads
      output_basename: output_basename
      ploidy: cfree_ploidy
      mate_orientation_sample: cfree_mate_orientation_sample
      capture_regions: unpadded_capture_regions
      indexed_reference_fasta: prepare_reference/indexed_fasta
      reference_fai: reference_fai
      b_allele: gatk_filter_germline/filtered_pass_vcf
      chr_len: cfree_chr_len
      coeff_var: cfree_coeff_var
      cfree_sex: cfree_sex
    out:
      [ctrlfreec_cnvs, ctrlfreec_pval, ctrlfreec_config, ctrlfreec_pngs, ctrlfreec_bam_ratio, ctrlfreec_bam_seg, ctrlfreec_baf, ctrlfreec_info]

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
- class: "sbg:maxNumberOfParallelInstances"
  value: 4
"sbg:license": Apache License 2.0
"sbg:publisher": KFDRC
"sbg:categories":
- BAM
- CRAM
- GATK
- MUTECT2
- CONTROLFREEC
- MANTA
- SOMATIC
- TUMORONLY
- VCF
