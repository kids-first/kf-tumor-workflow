cwlVersion: v1.1
class: CommandLineTool
id: gatk_bwamemindeximagecreator
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: 'broadinstitute/gatk:4.2.2.0'
  - class: ResourceRequirement
    ramMin: $(inputs.max_memory * 1000)
    coresMin: $(inputs.cores)
baseCommand: []
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      gatk --java-options "-Xmx${return Math.floor(inputs.max_memory*1000/1.074-1)}m" BwaMemIndexImageCreator

inputs:
  reference: { type: 'File', secondaryFiles: [^.dict, .fai], inputBinding: { prefix: "--input", position: 2 }, doc: "Refernce fasta, with fai and dict indices" }

  output_name: { type: 'string', inputBinding: { prefix: "--output", position: 2 }, doc: "Name of vcf.gz file into which variants should be written" }

  extra_args: { type: 'string?', inputBinding: { position: 3, shellQuote: false }, doc: "Any additional arguments for this tool. See GATK Documentation for complete list of options. Example input: --interval-merging-rule OVERLAPPING_ONLY" }

  cores: { type: 'int?', default: 1, doc: "Minimum reserved number of CPU cores for the task." }
  max_memory: { type: 'int?', default: 2, doc: "GB of RAM to allocate to the task." }

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.output_name)
