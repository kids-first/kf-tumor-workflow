cwlVersion: v1.1
class: CommandLineTool
id: gatk_createsomaticpanelofnormals
doc: |
  Make a panel of normals for use with Mutect2
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/gatk:4.1.1.0'
  - class: ResourceRequirement
    ramMin: $(inputs.max_memory * 1000)
    coresMin: $(inputs.cores)

baseCommand: []
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      /gatk --java-options "-Xmx${return Math.floor(inputs.max_memory*1000/1.074-1)}m" CreateSomaticPanelOfNormals

inputs:
  input_vcf: { type: 'File?', inputBinding: { prefix: "--variant", position: 2 }, secondaryFiles: [{ pattern: ".tbi", required: false}], doc: "A VCF file containing variants." }
  input_genomicsdb: { type: 'Directory?', inputBinding: { prefix: "--variant gendb://", separate: false, shellQuote: false, position: 2 }, doc: "GenomicsDB containing variants." }
  output_filename: { type: 'string', inputBinding: { prefix: "--output", position: 2 }, doc: "Output file with only the sample name in it." }
  reference_fasta: { type: 'File', secondaryFiles: [{ pattern: ".fai", required: true },{ pattern: "^.dict", required: true }], inputBinding: { prefix: "--reference", position: 2 }, doc: "Reference sequence fasta and index" }
  germline_resource: { type: 'File', secondaryFiles: ['.tbi'], inputBinding: {prefix: "--germline-resource", position: 2 }, doc: "Population vcf of germline sequencing containing allele fractions." }
  extra_args: { type: 'string?', inputBinding: { position: 3, shellQuote: false }, doc: "Any additional arguments for this tool. See GATK Documentation for complete list of options. Example input: --interval-merging-rule OVERLAPPING_ONLY" }
  cores: { type: 'int?', default: 2, doc: "Minimum reserved number of CPU cores for the task." }
  max_memory: { type: 'int?', default: 8, doc: "GB of RAM to allocate to the task." }
outputs:
  output:
    type: File
    secondaryFiles: [{ pattern: ".tbi", required: false },{ pattern: ".idx", required: false }]
    outputBinding:
      glob: $(inputs.output_filename)
