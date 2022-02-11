cwlVersion: v1.2
class: Workflow
id: kfdrc_mutect2_sub_wf
label: KFDRC Mutect2 Subworkflow
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement
  - class: InlineJavascriptRequirement

inputs:
  # MultiStep
  indexed_reference_fasta: { type: 'File', secondaryFiles: [.fai, ^.dict] }
  reference_dict: { type: 'File' }
  bed_invtl_split: { type: 'File[]', doc: "Bed file intervals passed on from and outside pre-processing step" }
  input_tumor_aligned: { type: 'File', secondaryFiles: [{ pattern: ".bai", required: false },{ pattern: "^.bai", required: false },{ pattern: ".crai", required: false },{ pattern: "^.crai", required: false }] }
  input_tumor_name: { type: 'string' }
  input_normal_aligned: { type: 'File?', secondaryFiles: [{ pattern: ".bai", required: false },{ pattern: "^.bai", required: false },{ pattern: ".crai", required: false },{ pattern: "^.crai", required: false }] }
  input_normal_name: { type: 'string?' }
  tool_name: { type: 'string?', doc: "String to describe what tool was run as part of file name", default: "mutect2" }
  output_basename: { type: 'string' }

  # Mutect2 Files
  af_only_gnomad_vcf: { type: 'File', secondaryFiles: ['.tbi'] }
  alleles: { type: 'File?', secondaryFiles: ['.tbi'] }
  panel_of_normals: { type: 'File?', secondaryFiles: ['.tbi'] }

  # Mutect2 Arguments (optional)
  disable_adaptive_pruning: { type: 'boolean?' }
  make_bamout: { type: 'boolean?' }
  run_orientation_bias_mixture_model_filter: { type: 'boolean?' }
  mutect2_extra_args: { type: 'string?' }

  # Filtration Support (optional)
  exac_common_vcf: { type: 'File?', secondaryFiles: ['.tbi'] }

  # FilterMutectCalls Arguments (optional)
  filtermutectcalls_extra_args: { type: 'string?' }

  # FilterAlignmentArtifacts Files (optional)
  bwa_mem_index_image: { type: 'File?' }

  # SelectVariants
  select_vars_mode: { type: ['null', { type: enum, name: select_vars_mode, symbols: ["gatk", "grep"] }], doc: "Choose 'gatk' for SelectVariants tool, or 'grep' for grep expression", default: "gatk" }

  # Annotation Files (optional)
  vep_cache: { type: 'File?', doc: "tar gzipped cache from ensembl/local converted cache" }
  genomic_hotspots: { type: 'File[]?', doc: "Tab-delimited BED formatted file(s) containing hg38 genomic positions corresponding to hotspots" }
  protein_snv_hotspots: { type: 'File[]?', doc: "Column-name-containing, tab-delimited file(s) containing protein names and amino acid positions corresponding to hotspots" }
  protein_indel_hotspots: { type: 'File[]?', doc: "Column-name-containing, tab-delimited file(s) containing protein names and amino acid position ranges corresponding to hotspots" }
  bcftools_annot_vcf: { type: 'File?', secondaryFiles: ['.tbi'], doc: "bgzipped annotation vcf file" }

  # Annotation Arguments (optional)
  run_annotation: { type: 'boolean?' }
  vep_ref_build: { type: 'string?', doc: "Genome ref build used, should line up with cache.", default: "GRCh38" }
  retain_info: { type: 'string?', doc: "csv string with INFO fields that you want to keep", default: "MBQ,TLOD,HotSpotAllele" }
  retain_fmt: { type: 'string?', doc: "csv string with FORMAT fields that you want to keep" }
  add_common_fields: { type: 'boolean?', doc: "Set to true if input is a strelka2 vcf that hasn't had common fields added", default: false }
  bcftools_annot_columns: { type: 'string?', doc: "csv string of columns from annotation to port into the input vcf, i.e INFO/AF", default: "INFO/AF" }
  bcftools_public_filter: { type: 'string?', doc: "Will hard filter final result to create a public version", default: FILTER="PASS"|INFO/HotSpotAllele=1 }
  gatk_filter_name: { type: 'string[]?', doc: "Array of names for each filter tag to add, recommend: [\"NORM_DP_LOW\", \"GNOMAD_AF_HIGH\"]" }
  gatk_filter_expression: { type: 'string[]?', doc: "Array of filter expressions to establish criteria to tag variants with. See https://gatk.broadinstitute.org/hc/en-us/articles/360036730071-VariantFiltration, recommend: \"vc.getGenotype('\" + inputs.input_normal_name + \"').getDP() <= 7\"), \"AF > 0.001\"]" }
  disable_hotspot_annotation: { type: 'boolean?', doc: "Disable Hotspot Annotation and skip this task.", default: false }
  maf_center: { type: 'string?', doc: "Sequencing center of variant called", default: "." }

  # Resource Control
  mutect_cores: { type: 'int?' }
  mutect_memory: { type: 'int?' }
  getpileup_memory: { type: 'int?' }
  learnorientation_memory: { type: 'int?' }
  filtermutectcalls_memory: { type: 'int?' }
  filteralignmentartifacts_cores: { type: 'int?' }
  filteralignmentartifacts_memory: { type: 'int?' }

outputs:
  mutect2_filtered_stats: { type: 'File', outputSource: filter_mutect2_vcf/stats_table }
  mutect2_filtered_vcf: { type: 'File', outputSource: pickvalue_workaround/output }
  mutect2_protected_outputs: { type: 'File[]?', outputSource: rename_protected/renamed_files }
  mutect2_public_outputs: { type: 'File[]?', outputSource: rename_public/renamed_files }
  mutect2_bam: { type: 'File?', outputSource: gatk_gathersortindexbams/output }

steps:
  mutect2:
    run: ../tools/gatk_mutect2.cwl
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c5.9xlarge
    in:
      reference: indexed_reference_fasta
      input_tumor_aligned: input_tumor_aligned
      input_normal_aligned: input_normal_aligned
      input_tumor_name: input_tumor_name
      input_normal_name: input_normal_name
      interval_list: bed_invtl_split
      germline_resource_vcf: af_only_gnomad_vcf
      panel_of_normals: panel_of_normals
      alleles: alleles
      output_vcf_name: { valueFrom: $(inputs.input_tumor_aligned.nameroot).$(inputs.interval_list.nameroot).Mutect2.vcf.gz }
      output_f1r2_name:
        source: run_orientation_bias_mixture_model_filter
        valueFrom: '$(self ? inputs.input_tumor_aligned.nameroot+"."+inputs.interval_list.nameroot+".f1r2_counts.tar.gz" : null)'
      output_bam_name:
        source: make_bamout
        valueFrom: '$(self ? inputs.input_tumor_aligned.nameroot + "." + inputs.interval_list.nameroot + ".bam" : null)'
      disable_read_filter: { valueFrom: '${return ["MateOnSameContigOrNoMappedMateReadFilter"]}' }
      disable_adaptive_pruning: disable_adaptive_pruning
      extra_args: mutect2_extra_args
      cores: mutect_cores
      max_memory: mutect_memory
    scatter: [interval_list]
    out: [mutect2_vcf, f1r2_counts, mutect_stats, mutect2_bam]

  mutect2_filter_support:
    run: ../subworkflows/kfdrc_mutect2_filter_support_subwf.cwl
    when: $(inputs.exac_common_vcf != null)
    in:
      indexed_reference_fasta: indexed_reference_fasta
      reference_dict: reference_dict
      wgs_calling_interval_list: bed_invtl_split
      input_tumor_aligned: input_tumor_aligned
      input_normal_aligned: input_normal_aligned
      exac_common_vcf: exac_common_vcf
      tool_name: tool_name
      output_basename: output_basename
      getpileup_memory: getpileup_memory
    out: [contamination_table, segmentation_table]

  gatk_gathersortindexbams:
    run: ../tools/gatk_gathersortindexbams.cwl
    when: '$(inputs.enable_tool ? true : false)'
    in:
      reference: indexed_reference_fasta
      input_bams: mutect2/mutect2_bam
      enable_tool: make_bamout
      output_basename: output_basename
    out: [output]

  gatk_learn_orientation_bias:
    run: ../tools/gatk_learnreadorientationmodel.cwl
    when: '$(inputs.enable_tool ? true : false)'
    in:
      input_tgz: mutect2/f1r2_counts
      tool_name: tool_name
      enable_tool: run_orientation_bias_mixture_model_filter
      output_basename: output_basename
      max_memory: learnorientation_memory
    out: [f1r2_bias]

  merge_mutect2_vcf:
    run: ../tools/gatk_mergevcfs.cwl
    in:
      input_vcfs: mutect2/mutect2_vcf
      output_basename: output_basename
      reference_dict: reference_dict
      tool_name: tool_name
    out: [merged_vcf]

  merge_mutect2_stats:
    run: ../tools/gatk_mergemutectstats.cwl
    in:
      input_stats: mutect2/mutect_stats
      output_basename: output_basename
    out: [merged_stats]

  filter_mutect2_vcf:
    run: ../tools/gatk_filtermutectcalls.cwl
    in:
      output_vcf_name:
        source: output_basename
        valueFrom: $(self).mutect2_filtered.vcf.gz
      output_filtering_stats:
        source: output_basename
        valueFrom: $(self).mutect2_filtered.txt
      reference: indexed_reference_fasta
      mutect_vcf: merge_mutect2_vcf/merged_vcf
      mutect_stats: merge_mutect2_stats/merged_stats
      ob_priors: gatk_learn_orientation_bias/f1r2_bias
      contamination_table: mutect2_filter_support/contamination_table
      segmentation_table: mutect2_filter_support/segmentation_table
      extra_args: filtermutectcalls_extra_args
      max_memory: filtermutectcalls_memory
    out: [stats_table, filtered_vcf]

  gatk_filteralignmentartifacts:
    run: ../tools/gatk_filteralignmentartifacts.cwl
    when: $(inputs.bwa_mem_index_image != null)
    in:
      input_vcf: filter_mutect2_vcf/filtered_vcf
      reference: indexed_reference_fasta
      input_reads: input_tumor_aligned
      bwa_mem_index_image: bwa_mem_index_image
      output_vcf_name:
        source: output_basename
        valueFrom: $(self).mutect2_filtered.artifact_filtered.vcf.gz
      cores: filteralignmentartifacts_cores
      max_memory: filteralignmentartifacts_memory
    out: [output]

  pickvalue_workaround:
    run: ../tools/expression_pickvalue_workaround.cwl
    in:
      input_file:
        source: [gatk_filteralignmentartifacts/output, filter_mutect2_vcf/filtered_vcf]
        pickValue: first_non_null
    out: [output]

  gatk_selectvariants_mutect2:
    run: ../tools/gatk_selectvariants.cwl
    when: '$(inputs.enable_tool ? true : false)'
    in:
      input_vcf:
        source: [gatk_filteralignmentartifacts/output, filter_mutect2_vcf/filtered_vcf]
        pickValue: first_non_null
      output_basename: output_basename
      tool_name: tool_name
      mode: select_vars_mode
      enable_tool: run_annotation
    out: [pass_vcf]

  annotate:
    run: ../subworkflows/kfdrc_annot_vcf_sub_wf.cwl
    when: '$(inputs.enable_workflow ? true : false)'
    in:
      indexed_reference_fasta: indexed_reference_fasta
      input_vcf: gatk_selectvariants_mutect2/pass_vcf
      input_tumor_name: input_tumor_name
      input_normal_name: input_normal_name
      add_common_fields: add_common_fields
      retain_info: retain_info
      retain_fmt: retain_fmt
      bcftools_annot_columns: bcftools_annot_columns
      bcftools_annot_vcf: bcftools_annot_vcf
      bcftools_public_filter: bcftools_public_filter
      gatk_filter_name: gatk_filter_name
      gatk_filter_expression: gatk_filter_expression
      vep_cache: vep_cache
      vep_ref_build: vep_ref_build
      disable_hotspot_annotation: disable_hotspot_annotation
      genomic_hotspots: genomic_hotspots
      protein_snv_hotspots: protein_snv_hotspots
      protein_indel_hotspots: protein_indel_hotspots
      maf_center: maf_center
      output_basename: output_basename
      tool_name: tool_name
      enable_workflow: run_annotation
    out: [annotated_protected_vcf, annotated_protected_maf, annotated_public_vcf, annotated_public_maf]

  rename_protected:
    run: ../tools/generic_rename_outputs.cwl
    when: '$(inputs.enable_tool ? true : false)'
    in:
      input_files:
        source: [annotate/annotated_protected_vcf, annotate/annotated_protected_maf]
        valueFrom: "${return [self[0],self[0].secondaryFiles[0],self[1]]}"
      rename_to:
        source: [output_basename, tool_name]
        valueFrom: "${var pro_vcf=self[0] + '.' + self[1] + '.norm.annot.protected.vcf.gz'; \
        var pro_tbi=self[0] + '.' + self[1] + '.norm.annot.protected.vcf.gz.tbi'; \
        var pro_maf=self[0] + '.' + self[1] + '.norm.annot.protected.maf'; \
        return [pro_vcf, pro_tbi, pro_maf];}"
      enable_tool: run_annotation
    out: [renamed_files]

  rename_public:
    run: ../tools/generic_rename_outputs.cwl
    when: '$(inputs.enable_tool ? true : false)'
    in:
      input_files:
        source: [annotate/annotated_public_vcf, annotate/annotated_public_maf]
        valueFrom: "${return [self[0],self[0].secondaryFiles[0],self[1]]}"
      rename_to:
        source: [output_basename, tool_name]
        valueFrom: "${var pub_vcf=self[0] + '.' + self[1] + '.norm.annot.public.vcf.gz'; \
        var pub_tbi=self[0] + '.' + self[1] + '.norm.annot.public.vcf.gz.tbi'; \
        var pub_maf=self[0] + '.' + self[1] + '.norm.annot.public.maf'; \
        return [pub_vcf, pub_tbi, pub_maf];}"
      enable_tool: run_annotation
    out: [renamed_files]

$namespaces:
  sbg: https://sevenbridges.com
