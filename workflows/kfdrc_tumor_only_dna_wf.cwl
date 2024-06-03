cwlVersion: v1.2
class: Workflow
id: kfdrc_tomor_only_dna_wf
label: KFDRC DNA Tumor Only Beta Production Workflow
doc: |
  # Kids First DRC Tumor Only Pipeline

  This repository contains tools and workflows for processing of tumor-only samples.
  It is currently in beta phase.
  Much of the components have been borrowed from the Kids First Somatic Workflow.
  It can also be used to process PDX data by first pre-processing reads using the Xenome tool, explained more here in documentation.

  <p align="center">
    <img src="docs/kids_first_logo.svg" alt="Kids First repository logo" width="660px" />
  </p>
  <p align="center">
    <a href="https://github.com/kids-first/kf-tumor-workflow/blob/main/LICENSE"><img src="https://img.shields.io/github/license/kids-first/kf-tumor-workflow.svg?style=for-the-badge"></a>
  </p>

  ## Import info on cloning the git repo
  This repository takes advantage of the git submodule feature.
  The Single Nucleotide Variant annotation workflow is maintained in our [Annotation Tools Repository](https://github.com/kids-first/kf-annotation-tools).
  Therefore, in order to get the code for a submodule you can either:
  - Clone the repository recursively with `git clone --recursive`
  - After cloning, run: `git submodule init && git submodule update`
  More info on how this worked [here](https://git-scm.com/book/en/v2/Git-Tools-Submodules)

  ## Main workflow
  The wrapper workflow that runs most of the tools is [found here](./workflows/kfdrc_tumor_only_dna_wf.cwl).

  ## Tools run
  ### Single Nucleotide Variant (SNV)
   - Mutect2 from GATK 4.2.2.0
   - Annotation using the [Kids First DRC Somatic SNV Annotation Workflow](https://github.com/kids-first/kf-annotation-tools/blob/master/workflows/kfdrc-somatic-snv-annot-workflow.cwl)
  ### Copy Number Variant (CNV)
   - ControlFREEC v11.6
  ### Structural Variant (SV)
   - Manta v1.4.0

  ## Inputs
  Most inputs have recommended values that should auto import both files and parameters
  ### Recommended file/param defaults:
   - `indexed_reference_fasta`: FAI and DICT indexed Homo_sapiens_assembly38.fasta
   - `mutect2_af_only_gnomad_vcf`: af-only-gnomad.hg38.vcf.gz
   - `mutect2_exac_common_vcf`: small_exac_common_3.hg38.vcf.gz
   - `gem_mappability_file`: hg38_canonical_150.mappability # will need note on file generation
   - `b_allele`: dbSNP_v153_ucsc-compatible.converted.vt.decomp.norm.common_snps.vcf.gz # will need note on file generation
   - `vep_cache`: homo_sapiens_merged_vep_105_indexed_GRCh38.tar.gz
   - `genomic_hotspots`: tert.bed # bed file with TERT gene promoter region
   - `protein_snv_hotspots`: protein_snv_cancer_hotspots_v2.ENS105_liftover.tsv
   - `protein_indel_hotspots`: protein_indel_cancer_hotspots_v2.ENS105_liftover.tsv
   - `echtvar_anno_zips`: gnomad.v3.1.1.custom.echtvar.zip
  ### Necessary for user to define:
   - `input_tumor_aligned`: Indexed BAM/CRAM/SAM file
   - `input_tumor_name`: sample name, should match read group sample name in `input_tumor_aligned`
   - `panel_of_normals`: Mutect2 Panel of Normals
   - `wgs_or_wxs`: Choose whether input is Whole Genome Sequencing (WGS) or Whole Exome Sequencing or Panel (WXS)
   - `calling_regions`:
      - For WGS: wgs_canonical_calling_regions.hg38.bed
      - For WXS: Unpadded experimental bait capture regions
   - `blacklist_regions`:
      - For WGS: hg38-blacklist.v2.bed.gz
      - For WXS: none
   - `cnv_blacklist_regions`:
      - For WGS: somatic-hg38_CNV_and_centromere_blacklist.hg38liftover.bed
      - For WXS: none
   - `i_flag`: for CNV calling, whether to intersect b allele file. Set to `N` skip
   - `cfree_sex`: for CNV calling, set to XX for female, XY for male
   - `cfree_ploidy`: Array of ploidy possibilities for ControlFREEC to try. Recommend [2,3,4]
   - `filtermutectcalls_extra_args`: "--min-allele-fraction 0.01"
   - `gatk_filter_name`: `["GNOMAD_AF_HIGH", "ALT_DEPTH_LOW"]`
   - `gatk_filter_expression`: `["gnomad_3_1_1_AF != '.' && gnomad_3_1_1_AF > 0.001 && gnomad_3_1_1_FILTER == 'PASS'", "vc.getGenotype('<input_tumor_name>').getAD().1 < 1"]`
   - `output_basename`: String value to use as basename for outputs

  ## Output Files
  ### Mutect2
   - `mutect2_protected_outputs`: VCF with SNV, MNV, and INDEL variant calls and of pipeline soft FILTER-added values in MAF and  VCF format with annotation, VCF index, and MAF format output
   - `mutect2_public_outputs`: Protected outputs, except MAF and VCF have had entries with soft FILTER values removed
   - `mutect2_bam`: BAM generated will be written as BAM. Useful for debugging
  ### ControlFREEC CNV
   - `ctrlfreec_pval`: Copy number call with GT (if BAF provided) and p values. Most people want this
   - `ctrlfreec_config`: Config file used to run
   - `ctrlfreec_pngs`: Visualization of CN and BAF
   - `ctrlfreec_bam_ratio`: Calls as log2 ratio
   - `ctrlfreec_bam_seg`: Custom made microarray-style SEG file
   - `ctrlfreec_baf`: b allele frequency file
   - `ctrlfreec_info`: Calculated information, like ploidy, if a range was given
  ### Manta SV
   - `manta_pass_vcf`: VCF file with SV calls that PASS
   - `manta_prepass_vcf`: VCF file with all SV calls
   - `annotsv_annotated_calls`: Manta calls annotated with AnnotSV
   - `annotsv_unannotated_calls`: Manta calls not annotated with AnnotSV
requirements:
- class: ScatterFeatureRequirement
- class: MultipleInputFeatureRequirement
- class: SubworkflowFeatureRequirement
- class: InlineJavascriptRequirement
- class: StepInputExpressionRequirement
inputs:
  # Runtime Switches
  wgs_or_wxs: {type: {type: enum, name: wgs_or_wxs, symbols: ["WGS", "WXS"]}, doc: "Select input type WGS or WXS"}
  run_mutect2: {type: 'boolean?', default: true, doc: "Set to false to disable Mutect2"}
  run_manta: {type: 'boolean?', default: true, doc: "Set to false to disable Manta"}
  run_controlfreec: {type: 'boolean?', default: true, doc: "Set to false to disable Control-FREEC"}
  run_calmd_bam: {type: 'boolean?', default: false, doc: "Set to true to enable samtools calmd running for BAM inputs"}
  make_bamout: {type: 'boolean?', default: false, doc: "Should Mutect2 create a BAM output? Turning on this mode may result in serious
      performance cost for the caller."}
  run_orientation_bias_mixture_model_filter: {type: 'boolean?', default: true, doc: "Should Orientation Bias Mixture Model Filter
      be applied to Mutect2 outputs?"}
  i_flag: {type: ['null', {type: enum, name: i_flag, symbols: ["Y", "N"]}], doc: "Flag to intersect germline calls on padded regions.
      Use N if you want to skip this or have a WGS run"}
  disable_hotspot_annotation: {type: 'boolean?', default: false, doc: "Set to true to disable Hotspot Annotation"}

  # Shared
  indexed_reference_fasta: {type: 'File', secondaryFiles: [{pattern: ".fai", required: true}, {pattern: "^.dict", required: true}],
    doc: "Reference fasta to which the input reads were aligned. Must include FAI and DICT indicies", "sbg:suggestedValue": {class: File,
      path: 60639014357c3a53540ca7a3, name: Homo_sapiens_assembly38.fasta, secondaryFiles: [{class: File, path: 60639016357c3a53540ca7af,
          name: Homo_sapiens_assembly38.fasta.fai}, {class: File, path: 60639019357c3a53540ca7e7, name: Homo_sapiens_assembly38.dict}]}, "sbg:fileTypes": "FASTA, FA"}
  calling_regions: {type: 'File', doc: "BED or INTERVALLIST file containing a set of genomic regions over which the callers will be
      run. For WGS, this should be the wgs_calling_regions.interval_list. For WXS, the user must provide the appropriate regions for
      their analysis.", "sbg:fileTypes": "BED, INTERVALLIST, INTERVAL_LIST"}
  blacklist_regions: {type: 'File?', doc: "BED or INTERVALLIST file containing a set of genomic regions to remove from the calling
      regions for SNV and SV calling.", "sbg:suggestedValue": {class: File, path: 665df995a193b420129c7830, name: hg38-blacklist.v2.bed.gz}, "sbg:fileTypes": "BED, BED.GZ, INTERVALLIST, INTERVAL_LIST"}
  cnv_blacklist_regions: {type: 'File?', doc: "BED or INTERVALLIST file containing a set of genomic regions to remove from the calling
      regions for CNV calling only!", "sbg:suggestedValue": {class: File, path: 665df995a193b420129c782f, name: somatic-hg38_CNV_and_centromere_blacklist.hg38liftover.bed},
      "sbg:fileTypes": "BED, INTERVALLIST, INTERVAL_LIST"}
  input_tumor_aligned: {type: 'File', secondaryFiles: [{pattern: ".bai", required: false}, {pattern: "^.bai", required: false}, {
        pattern: ".crai", required: false}, {pattern: "^.crai", required: false}], doc: "BAM/SAM/CRAM file containing tumor reads",
    "sbg:fileTypes": "BAM, CRAM, SAM"}
  input_tumor_name: {type: 'string', doc: "BAM sample name of tumor. May be URL-encoded as output by GetSampleName with -encode argument."}
  old_tumor_name: {type: 'string?', doc: "If `SM:` sample name in the align file is different than `input_tumor_name`, you **must**
      provide it here"}
  output_basename: {type: 'string', doc: "String value to use as basename for outputs"}
  panel_of_normals: {type: 'File?', secondaryFiles: ['.tbi'], doc: "VCF file (and index) of sites observed in normal. A panel of normals
      can be a useful (optional) input to help filter out commonly seen sequencing noise that may appear as low allele-fraction somatic
      variants.", "sbg:suggestedValue": {class: File, path: 65fb2dfef70c1f39338b2aa8, name: 1000g_pon.hg38.vcf.gz, secondaryFiles: [{class: File,
         path: 665df995a193b420129c782e,  name: 1000g_pon.hg38.vcf.gz.tbi}]}, "sbg:fileTypes": "VCF, VCF.GZ"}
  mutect2_af_only_gnomad_vcf: {type: 'File', secondaryFiles: ['.tbi'], doc: "Population vcf (and index) of germline sequencing containing
      allele fractions in VCF format.", "sbg:suggestedValue": {class: File, path: 5f50018fe4b054958bc8d2e3, name: af-only-gnomad.hg38.vcf.gz,
      secondaryFiles: [{class: File, path: 5f50018fe4b054958bc8d2e5, name: af-only-gnomad.hg38.vcf.gz.tbi}]}, "sbg:fileTypes": "VCF,
      VCF.GZ"}
  mutect2_alleles_vcf: {type: 'File?', secondaryFiles: ['.tbi'], doc: "VCF file (and index) containing set of alleles for which to
      force genotyping regardless of evidence", "sbg:fileTypes": "VCF, VCF.GZ"}
  mutect2_exac_common_vcf: {type: 'File?', secondaryFiles: ['.tbi'], doc: "Exac Common VCF (and index) used for calculating contamination
      values used in filtering the Mutect VCF. If do not wish to perfom this filtering, remove this input.", "sbg:suggestedValue": {
      class: File, path: 5f500135e4b0370371c051ad, name: small_exac_common_3.hg38.vcf.gz, secondaryFiles: [{class: File, path: 5f500135e4b0370371c051af,
          name: small_exac_common_3.hg38.vcf.gz.tbi}]}, "sbg:fileTypes": "VCF, VCF.GZ"}
  bwa_mem_index_image: {type: 'File?', doc: "BWA-mem index image used for FilterAlignmentArtifacts. WARNING!!! Experimental Software!!!
      Do not provide unless you are absolutely certain you want to run FilterAlignmentArtifacts.", "sbg:fileTypes": "IMG"}
  mutect2_extra_args: {type: 'string?', doc: "Any additional arguments for Mutect2. See GATK documentation for all available options."}
  filtermutectcalls_extra_args: {type: 'string?', default: "--min-allele-fraction 0.01", doc: "Any additional arguments for FilterMutectCalls.
      See GATK documentation for all available options. Recommend: --min-allele-fraction 0.01"}
  select_vars_mode: {type: ['null', {type: enum, name: select_vars_mode, symbols: ["gatk", "grep"]}], default: "gatk", doc: "Choose
      'gatk' for SelectVariants tool, or 'grep' for grep expression"}
  mate_copynumber_file_sample: {type: 'File?', doc: "Tumor cpn file from previous run. If used, will override bam use"}
  gem_mappability_file: {type: 'File?', doc: "GEM mappability file to make read count adjustments with", "sbg:suggestedValue": {class: File,
      path: 663d2bcc27374715fccd8c68, name: hg38_canonical_150.mappability}}
  min_subclone_presence: {type: 'int?', doc: "Tool default 100 (meaning \"do not look for subclones\"). Suggested: 20 (or 0.2) for
      WGS and 30 (or 0.3) for WES."}
  cfree_ploidy: {type: 'int[]?', default: [2, 3, 4], doc: "Array of ploidy possibilities for ControlFreeC to try"}
  cfree_mate_orientation_sample: {type: ['null', {type: enum, name: mate_orientation_sample, symbols: ["0", "FR", "RF", "FF"]}], default: "FR",
    doc: "0 (for single ends), RF (Illumina mate-pairs), FR (Illumina paired-ends), FF (SOLiD mate-pairs)"}
  b_allele: {type: 'File?', secondaryFiles: [{pattern: ".tbi", required: true}], doc: "germline calls, needed for BAF.  GATK HC VQSR
      input recommended. Tool will prefilter for germline and pass if expression given", "sbg:suggestedValue": {class: File, path: 663e7c69f180a9177dd4212b,
      name: dbSNP_v153_ucsc-compatible.converted.vt.decomp.norm.common_snps.vcf.gz, secondaryFiles: [{class: File, path: 6509b6a37f6417197fc158fc,
          name: dbSNP_v153_ucsc-compatible.converted.vt.decomp.norm.common_snps.vcf.gz.tbi}]}}
  cfree_coeff_var: {type: 'float?', default: 0.05, doc: "Coefficient of variation to set window size. Default 0.05 recommended"}
  cfree_sex: {type: ['null', {type: enum, name: sex, symbols: ["XX", "XY"]}], doc: "If known, XX for female, XY for male"}

  # Annotation
  annotsv_annotations_dir_tgz: {type: 'File?', doc: "TAR.GZ'd Directory containing annotations for AnnotSV", "sbg:fileTypes": "TAR,
      TAR.GZ, TGZ", "sbg:suggestedValue": {class: File, path: 6328ab26d01163633dabcc2e, name: annotsv_311_plus_ens105_annotations_dir.tgz}}
  vep_cache: {type: 'File', doc: "tar gzipped cache from ensembl/local converted cache", "sbg:suggestedValue": {class: File, path: 6332f8e47535110eb79c794f,
      name: homo_sapiens_merged_vep_105_indexed_GRCh38.tar.gz}}
  vep_buffer_size: {type: 'int?', default: 5000, doc: "Increase or decrease to balance speed and memory usage"}
  dbnsfp: {type: 'File?', secondaryFiles: [.tbi, ^.readme.txt], doc: "VEP-formatted plugin file, index, and readme file containing
      dbNSFP annotations"}
  dbnsfp_fields: {type: 'string?', doc: "csv string with desired fields to annotate. Use ALL to grab all"}
  merged: {type: 'boolean?', default: true, doc: "Set to true if merged cache used"}
  cadd_indels: {type: 'File?', secondaryFiles: [.tbi], doc: "VEP-formatted plugin file and index containing CADD indel annotations"}
  cadd_snvs: {type: 'File?', secondaryFiles: [.tbi], doc: "VEP-formatted plugin file and index containing CADD SNV annotations"}
  run_cache_existing: {type: 'boolean?', doc: "Run the check_existing flag for cache"}
  run_cache_af: {type: 'boolean?', doc: "Run the allele frequency flags for cache"}
  genomic_hotspots: {type: 'File[]?', doc: "Tab-delimited BED formatted file(s) containing hg38 genomic positions corresponding to
      hotspots", "sbg:suggestedValue": [{class: File, path: 607713829360f10e3982a423, name: tert.bed}]}
  protein_snv_hotspots: {type: 'File[]?', doc: "Column-name-containing, tab-delimited file(s) containing protein names and amino acid
      positions corresponding to hotspots", "sbg:suggestedValue": [{class: File, path: 663d2bcc27374715fccd8c6a, name: protein_snv_cancer_hotspots_v2.ENS105_liftover.tsv}]}
  protein_indel_hotspots: {type: 'File[]?', doc: "Column-name-containing, tab-delimited file(s) containing protein names and amino
      acid position ranges corresponding to hotspots", "sbg:suggestedValue": [{class: File, path: 663d2bcc27374715fccd8c6f, name: protein_indel_cancer_hotspots_v2.ENS105_liftover.tsv}]}
  mutect2_retain_info: {type: 'string?', doc: "csv string with INFO fields that you want to keep", default: "gnomad_3_1_1_AC,gnomad_3_1_1_AN,gnomad_3_1_1_AF,gnomad_3_1_1_nhomalt,gnomad_3_1_1_AC_popmax,gnomad_3_1_1_AN_popmax,gnomad_3_1_1_AF_popmax,gnomad_3_1_1_nhomalt_popmax,gnomad_3_1_1_AC_controls_and_biobanks,gnomad_3_1_1_AN_controls_and_biobanks,gnomad_3_1_1_AF_controls_and_biobanks,gnomad_3_1_1_AF_non_cancer,gnomad_3_1_1_primate_ai_score,gnomad_3_1_1_splice_ai_consequence,gnomad_3_1_1_AF_non_cancer_afr,gnomad_3_1_1_AF_non_cancer_ami,gnomad_3_1_1_AF_non_cancer_asj,gnomad_3_1_1_AF_non_cancer_eas,gnomad_3_1_1_AF_non_cancer_fin,gnomad_3_1_1_AF_non_cancer_mid,gnomad_3_1_1_AF_non_cancer_nfe,gnomad_3_1_1_AF_non_cancer_oth,gnomad_3_1_1_AF_non_cancer_raw,gnomad_3_1_1_AF_non_cancer_sas,gnomad_3_1_1_AF_non_cancer_amr,gnomad_3_1_1_AF_non_cancer_popmax,gnomad_3_1_1_AF_non_cancer_all_popmax,gnomad_3_1_1_FILTER,MQ,MQ0,QSI,HotSpotAllele"}
  mutect2_retain_fmt: {type: 'string?', doc: "csv string with FORMAT fields that you want to keep"}
  mutect2_retain_ann: {type: 'string?', doc: "csv string of annotations (within the VEP CSQ/ANN) to retain as extra columns in MAF",
    default: "HGVSg"}
  bcftools_strip_columns: {type: 'string?', doc: "csv string of columns to strip if needed to avoid conflict, i.e INFO/AF"}
  bcftools_public_filter: {type: 'string?', doc: "Will hard filter final result to create a public version", default: FILTER="PASS"|INFO/HotSpotAllele=1}
  echtvar_anno_zips: {type: 'File[]?', doc: "Annotation ZIP files for echtvar anno", "sbg:suggestedValue": [{class: File, path: 65c64d847dab7758206248c6,
        name: gnomad.v3.1.1.custom.echtvar.zip}]}
  gatk_filter_name: {type: 'string[]', doc: "Array of names for each filter tag to add, recommend: [\"GNOMAD_AF_HIGH\", \"ALT_DEPTH_LOW\"\
      ]"}
  gatk_filter_expression: {type: 'string[]', doc: "Array of filter expressions to establish criteria to tag variants with. See https://gatk.broadinstitute.org/hc/en-us/articles/360036730071-VariantFiltration,
      recommend: [\"gnomad_3_1_1_AF != '.' && gnomad_3_1_1_AF > 0.001 && gnomad_3_1_1_FILTER == 'PASS'\", \"vc.getGenotype('<input_tumor_name>').getAD().1
      < 1\"]"}
  maf_center: {type: 'string?', doc: "Sequencing center of variant called", default: "."}
  custom_enst: {type: 'File?', doc: "Use a file with ens tx IDs for each gene to override VEP PICK", "sbg:suggestedValue": {class: File,
      path: 6480c8a61dfc710d24a3a368, name: kf_isoform_override.tsv}}

  # Resource Control
  vep_ram: {type: 'int?', default: 32, doc: "In GB, may need to increase this value depending on the size/complexity of input"}
  vep_cores: {type: 'int?', default: 16, doc: "Number of cores to use. May need to increase for really large inputs"}
  mutect_cores: {type: 'int?', doc: "CPUs to allocate to GATK Mutect2"}
  mutect_memory: {type: 'int?', doc: "GB of memory to allocate to GATK Mutect2 (hard-capped)"}
  getpileup_memory: {type: 'int?', doc: "GB of memory to allocate to GATK GetPileupSummaries (hard-capped)"}
  learnorientation_memory: {type: 'int?', doc: "GB of memory to allocate to GATK LearnReadOrientationModel (hard-capped)"}
  filtermutectcalls_memory: {type: 'int?', doc: "GB of memory to allocate to GATK FilterMutectCalls (hard-capped)"}
  filteralignmentartifacts_cores: {type: 'int?', doc: "CPUs to allocate to GATK FilterAlignmentArtifacts"}
  filteralignmentartifacts_memory: {type: 'int?', doc: "GB of memory to allocate to GATK FilterAlignmentArtifacts (hard-capped)"}
  cfree_threads: {type: 'int?', default: 16, doc: "For ControlFreeC. Recommend 16 max, as I/O gets saturated after that losing any
      advantage"}
  manta_memory: {type: 'int?', doc: "GB of memory to allocate to Manta; defaults to 10 (soft-capped)"}
  manta_cores: {type: 'int?', doc: "Number of cores to allocate to Manta; defaults to 18"}
outputs:
  # Mutect2
  mutect2_protected_outputs: {type: 'File[]?', outputSource: mutect2/mutect2_protected_outputs, doc: "VCF with SNV, MNV, and INDEL
      variant calls and of pipeline soft FILTER-added values in MAF and  VCF format with annotation, VCF index, and MAF format output"}
  mutect2_public_outputs: {type: 'File[]?', outputSource: mutect2/mutect2_public_outputs, doc: "Protected outputs, except MAF and
      VCF have had entries with soft FILTER values removed"}
  mutect2_bam: {type: 'File?', outputSource: mutect2/mutect2_bam, doc: "BAM generated by Mutect2, if desired. The assembled haplotypes
      and locally realigned reads will be written as BAM. Really for debugging purposes only."}
  ctrlfreec_pval: {type: 'File?', outputSource: controlfreec/ctrlfreec_pval, doc: 'Copy number call with GT (if BAF provided) and
      p values. Most people want this'}
  ctrlfreec_config: {type: 'File?', outputSource: controlfreec/ctrlfreec_config, doc: 'Config file used to run'}
  ctrlfreec_pngs: {type: 'File[]?', outputSource: controlfreec/ctrlfreec_pngs, doc: 'Visualization of CN and BAF'}
  ctrlfreec_bam_ratio: {type: 'File?', outputSource: controlfreec/ctrlfreec_bam_ratio, doc: 'Calls as log2 ratio'}
  ctrlfreec_bam_seg: {type: 'File?', outputSource: controlfreec/ctrlfreec_bam_seg, doc: 'Custom made microarray-style seg file'}
  ctrlfreec_baf: {type: 'File?', outputSource: controlfreec/ctrlfreec_baf, doc: 'b allele frequency file'}
  ctrlfreec_info: {type: 'File?', outputSource: controlfreec/ctrlfreec_info, doc: 'Calculated inforamtion, like ploidy, if a range
      was given'}
  manta_pass_vcf: {type: 'File?', outputSource: manta/manta_pass_vcf, doc: 'VCF file with SV calls that PASS'}
  manta_prepass_vcf: {type: 'File?', outputSource: manta/manta_prepass_vcf, doc: 'VCF file with all SV calls'}
  annotsv_annotated_calls: {type: 'File?', outputSource: manta/annotsv_annotated_calls, doc: 'Manta calls annotated with AnnotSV'}
  annotsv_unannotated_calls: {type: 'File?', outputSource: manta/annotsv_unannotated_calls, doc: 'Manta calls not annotated with AnnotSV'}
steps:
  prepare_regions_padded:
    run: ../subworkflows/prepare_regions.cwl
    when: $(inputs.wgs_or_wxs == 'WXS' && inputs.run_tool.some(function(e) { return e }))
    in:
      wgs_or_wxs: wgs_or_wxs
      run_tool:
        source: [run_mutect2, run_manta]
      reference_dict:
        source: indexed_reference_fasta
        valueFrom: |
          $(self.secondaryFiles.filter(function(e) { return e.basename.search(/.dict$/) != -1 })[0])
      calling_regions: calling_regions
      calling_padding:
        valueFrom: $(100)
      blacklist_regions: blacklist_regions
      break_bands_at_multiples_of:
        valueFrom: $(0)
      scatter_count:
        valueFrom: $(50)
    out: [prescatter_intervallist, prescatter_bed, prescatter_bedgz, scattered_intervallists, scattered_beds]
  prepare_regions_unpadded:
    run: ../subworkflows/prepare_regions.cwl
    when: $(inputs.wgs_or_wxs == 'WGS' && inputs.run_tool.some(function(e) { return e }))
    in:
      wgs_or_wxs: wgs_or_wxs
      run_tool:
        source: [run_mutect2, run_manta]
      reference_dict:
        source: indexed_reference_fasta
        valueFrom: |
          $(self.secondaryFiles.filter(function(e) { return e.basename.search(/.dict$/) != -1 })[0])
      calling_regions: calling_regions
      blacklist_regions: blacklist_regions
      break_bands_at_multiples_of:
        valueFrom: $(80000000)
      scatter_count:
        valueFrom: $(50)
    out: [prescatter_intervallist, prescatter_bed, prescatter_bedgz, scattered_intervallists, scattered_beds]
  prepare_regions_unpadded_cnv:
    run: ../subworkflows/prepare_regions.cwl
    when: $(inputs.run_tool)
    in:
      run_tool: run_controlfreec
      reference_dict:
        source: indexed_reference_fasta
        valueFrom: |
          $(self.secondaryFiles.filter(function(e) { return e.basename.search(/.dict$/) != -1 })[0])
      calling_regions: calling_regions
      blacklist_regions: cnv_blacklist_regions
      scatter_count:
        valueFrom: $(0)
    out: [prescatter_intervallist, prescatter_bed, prescatter_bedgz, scattered_intervallists, scattered_beds]
  mutect2:
    run: ../subworkflows/kfdrc_mutect2_sub_wf.cwl
    when: $(inputs.run_mutect2)
    in:
      run_mutect2: run_mutect2
      indexed_reference_fasta: indexed_reference_fasta
      reference_dict:
        source: indexed_reference_fasta
        valueFrom: |
          $(self.secondaryFiles.filter(function(e) { return e.basename.search(/.dict$/) != -1 })[0])
      bed_invtl_split:
        source: [prepare_regions_padded/scattered_beds, prepare_regions_unpadded/scattered_beds]
        valueFrom: |
          $(self.some(function(e) { return e != null }) ? self.filter(function(e) { return e != null })[0] : null)
      af_only_gnomad_vcf: mutect2_af_only_gnomad_vcf
      alleles: mutect2_alleles_vcf
      exac_common_vcf: mutect2_exac_common_vcf
      input_tumor_aligned: input_tumor_aligned
      input_tumor_name: input_tumor_name
      old_tumor_name: old_tumor_name
      panel_of_normals: panel_of_normals
      disable_adaptive_pruning:
        source: wgs_or_wxs
        valueFrom: $(self == 'WXS')
      mutect2_extra_args: mutect2_extra_args
      filtermutectcalls_extra_args: filtermutectcalls_extra_args
      select_vars_mode: select_vars_mode
      retain_info: mutect2_retain_info
      retain_fmt: mutect2_retain_fmt
      retain_ann: mutect2_retain_ann
      bcftools_strip_columns: bcftools_strip_columns
      bcftools_public_filter: bcftools_public_filter
      echtvar_anno_zips: echtvar_anno_zips
      dbnsfp: dbnsfp
      dbnsfp_fields: dbnsfp_fields
      merged: merged
      cadd_indels: cadd_indels
      cadd_snvs: cadd_snvs
      run_cache_af: run_cache_af
      run_cache_existing: run_cache_existing
      gatk_filter_name: gatk_filter_name
      gatk_filter_expression: gatk_filter_expression
      vep_cache: vep_cache
      vep_ram: vep_ram
      vep_cores: vep_cores
      vep_buffer_size: vep_buffer_size
      disable_hotspot_annotation: disable_hotspot_annotation
      genomic_hotspots: genomic_hotspots
      protein_snv_hotspots: protein_snv_hotspots
      protein_indel_hotspots: protein_indel_hotspots
      maf_center: maf_center
      custom_enst: custom_enst
      bwa_mem_index_image: bwa_mem_index_image
      make_bamout: make_bamout
      run_orientation_bias_mixture_model_filter: run_orientation_bias_mixture_model_filter
      mutect_cores: mutect_cores
      mutect_memory: mutect_memory
      getpileup_memory: getpileup_memory
      learnorientation_memory: learnorientation_memory
      filtermutectcalls_memory: filtermutectcalls_memory
      filteralignmentartifacts_cores: filteralignmentartifacts_cores
      filteralignmentartifacts_memory: filteralignmentartifacts_memory
      tool_name:
        valueFrom: "mutect2.tumor_only"
      output_basename: output_basename
    out: [mutect2_filtered_stats, mutect2_protected_outputs, mutect2_public_outputs, mutect2_bam]
  bedtools_intersect_germline:
    run: ../tools/bedtools_intersect.cwl
    when: $(inputs.run_tool)
    in:
      run_tool: run_controlfreec
      input_vcf: b_allele
      output_basename: output_basename
      input_bed_file: prepare_regions_unpadded_cnv/prescatter_bed
      flag:
        source: [i_flag, wgs_or_wxs]
        valueFrom: |
          $(self[0] != null ? self[0] : self[1] == "WXS" ? "Y" : "N")
    out: [intersected_vcf]
  gatk_filter_germline:
    run: ../tools/gatk_filter_germline_variant.cwl
    when: $(inputs.run_tool)
    in:
      run_tool: run_controlfreec
      input_vcf: bedtools_intersect_germline/intersected_vcf
      reference_fasta: indexed_reference_fasta
      output_basename: output_basename
    out: [filtered_vcf, filtered_pass_vcf]
  samtools_cram2bam_plus_calmd_tumor:
    run: ../tools/samtools_calmd.cwl
    when: $(inputs.run_tool == true && (inputs.input_reads.basename.search(/.cram$/) != -1 || inputs.run_anyway))
    in:
      run_tool: run_controlfreec
      run_anyway: run_calmd_bam
      input_reads: input_tumor_aligned
      threads:
        valueFrom: $(16)
      reference: indexed_reference_fasta
    out: [bam_file]
  controlfreec:
    run: ../subworkflows/kfdrc_controlfreec_sub_wf.cwl
    when: $(inputs.run_controlfreec)
    in:
      run_controlfreec: run_controlfreec
      wgs_or_wxs: wgs_or_wxs
      mate_copynumber_file_sample: mate_copynumber_file_sample
      gem_mappability_file: gem_mappability_file
      min_subclone_presence: min_subclone_presence
      input_tumor_aligned:
        source: [samtools_cram2bam_plus_calmd_tumor/bam_file, input_tumor_aligned]
        pickValue: first_non_null
      input_tumor_name: input_tumor_name
      threads: cfree_threads
      output_basename: output_basename
      ploidy: cfree_ploidy
      mate_orientation_sample: cfree_mate_orientation_sample
      calling_regions: prepare_regions_unpadded_cnv/prescatter_bed
      indexed_reference_fasta: indexed_reference_fasta
      b_allele: gatk_filter_germline/filtered_pass_vcf
      coeff_var: cfree_coeff_var
      cfree_sex: cfree_sex
    out: [ctrlfreec_cnvs, ctrlfreec_pval, ctrlfreec_config, ctrlfreec_pngs, ctrlfreec_bam_ratio, ctrlfreec_bam_seg, ctrlfreec_baf,
      ctrlfreec_info]
  manta:
    run: ../subworkflows/kfdrc_manta_sub_wf.cwl
    when: $(inputs.run_manta)
    in:
      run_manta: run_manta
      indexed_reference_fasta: indexed_reference_fasta
      hg38_strelka_bed:
        source: [prepare_regions_padded/prescatter_bedgz, prepare_regions_unpadded/prescatter_bedgz]
        valueFrom: |
          $(self.some(function(e) { return e != null }) ? self.filter(function(e) { return e != null })[0] : null)
      input_tumor_aligned: input_tumor_aligned
      input_tumor_name: input_tumor_name
      annotsv_annotations_dir_tgz: annotsv_annotations_dir_tgz
      output_basename: output_basename
      manta_memory: manta_memory
      manta_cores: manta_cores
    out: [manta_prepass_vcf, manta_pass_vcf, manta_small_indels, annotsv_annotated_calls, annotsv_unannotated_calls]
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

"sbg:links":
- id: 'https://github.com/kids-first/kf-tumor-workflow/tree/v0.4.1-beta'
  label: github-release
