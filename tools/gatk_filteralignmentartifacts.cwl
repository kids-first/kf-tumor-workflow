cwlVersion: v1.1
class: CommandLineTool
id: gatk_filteralignmentartifacts
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/gatk:4.2.0.0R'
  - class: ResourceRequirement
    ramMin: $(inputs.max_memory * 1000)
    coresMin: $(inputs.cores)
baseCommand: []
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      /gatk --java-options "-Xmx${return Math.floor(inputs.max_memory*1000/1.074-1)}m" FilterAlignmentArtifacts

inputs:
  reference: { type: 'File', secondaryFiles: [^.dict, .fai], inputBinding: { prefix: "--reference", position: 2 }, doc: "Refernce fasta, with fai and dict indices" }
  input_reads: { type: 'File', secondaryFiles: [{ pattern: ".bai", required: false },{ pattern: "^.bai", required: false },{ pattern: ".crai", required: false },{ pattern: "^.crai", required: false }], inputBinding: { prefix: "--input", position: 2 }, doc: "BAM/SAM/CRAM file containing reads" }
  bwa_mem_index_image: { type: 'File', inputBinding: { prefix: "--bwa-mem-index-image", position: 2 }, doc: "BWA-mem index image" }
  input_vcf: { type: 'File', secondaryFiles: [.tbi], inputBinding: { prefix: "--variant", position: 2 }, doc: "VCF files containing variants" }

  output_vcf_name: { type: 'string', inputBinding: { prefix: "--output", position: 2 }, doc: "Name of vcf.gz file into which variants should be written" }

  extra_args: { type: 'string?', inputBinding: { position: 3, shellQuote: false }, doc: "Any additional arguments for this tool. See GATK Documentation for complete list of options. Example input: --interval-merging-rule OVERLAPPING_ONLY" }

  cores: { type: 'int?', default: 2, doc: "Minimum reserved number of CPU cores for the task." }
  max_memory: { type: 'int?', default: 10, doc: "GB of RAM to allocate to the task." }

outputs:
  output:
    type: File
    outputBinding:
      glob: '*.vcf.gz'
    secondaryFiles: [.tbi]
