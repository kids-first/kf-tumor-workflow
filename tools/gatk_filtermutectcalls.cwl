cwlVersion: v1.0
class: CommandLineTool
id: gatk_filtermutectcalls
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: 'broadinstitute/gatk:4.2.0.0'
  - class: ResourceRequirement
    ramMin: ${ return inputs.max_memory * 1000 }
    coresMin: $(inputs.cores)
baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      gatk --java-options "-Xmx${return Math.floor(inputs.max_memory*1000/1.074-1)}m" FilterMutectCalls

inputs:
  output_vcf_name: { type: 'string', inputBinding: { prefix: "--output", position: 2 }, doc: "Name of the output filtered VCF file" }
  output_filtering_stats: { type: 'string', inputBinding: { prefix: "--filtering-stats", position: 2 }, doc: "Name of the output filtering stats file" }
  reference: { type: 'File', secondaryFiles: ['.fai','^.dict'], inputBinding: { prefix: "--reference", position: 2 }, doc: "Reference sequence file with fai and dict indices" }
  mutect_vcf: { type: 'File', secondaryFiles: ['.tbi'], inputBinding: { prefix: "--variant", position: 2 }, doc: "A VCF file from Mutect2 containing variants" }
  mutect_stats: { type: 'File', inputBinding: { prefix: "--stats", position: 2 }, doc: "The Mutect stats file output by Mutect2" }

  ob_priors: { type: 'File?', inputBinding: { prefix: "--orientation-bias-artifact-priors", position: 2 }, doc: "One or more .tar.gz files containing tables of prior artifact probabilities for the read orientation filter model, one table per tumor sample" }
  contamination_table: { type: 'File?', inputBinding: { prefix: "--contamination-table", position: 2 }, doc: "Tables containing contamination information." }
  segmentation_table: { type: 'File?', inputBinding: { prefix: "--tumor-segmentation", position: 2 }, doc: "Tables containing tumor segments' minor allele fractions for germline hets emitted by CalculateContamination" }

  extra_args: { type: 'string?', inputBinding: { position: 3, shellQuote: false }, doc: "Any additional arguments for this tool. See GATK Documentation for complete list of options. Example input: --interval-merging-rule OVERLAPPING_ONLY" }

  cores: {type: 'int?', default: 2}
  max_memory: {type: 'int?', default: 4, doc: "GB of memory to allocate to the task"}

outputs:
  stats_table:
    type: File
    outputBinding:
      glob: $(inputs.output_filtering_stats)

  filtered_vcf:
    type: File
    outputBinding:
      glob: $(inputs.output_vcf_name)
    secondaryFiles: ['.tbi']
