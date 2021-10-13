cwlVersion: v1.1
class: CommandLineTool
id: gatk_gathersortindexbams
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/gatk:4.1.1.0'
  - class: ResourceRequirement
    ramMin: $(inputs.max_memory * 1000)
    coresMin: $(inputs.cores)
baseCommand: ["/bin/bash","-c"]
arguments:
  - position: 1
    shellQuote: true
    valueFrom: >-
      set -e

      /gatk --java-options "-Xmx${return Math.floor(inputs.max_memory*1000/1.074-1)}m" GatherBamFiles
      -I $(inputs.input_bams.map(function(e){return e.path;}).join(" -I "))
      -O unsorted.out.bam
      -R $(inputs.reference.path) 

      /gatk --java-options "-Xmx${return Math.floor(inputs.max_memory*1000/1.074-1)}m" SortSam
      -I unsorted.out.bam
      -O $(inputs.output_basename).sorted.bam
      --SORT_ORDER coordinate
      -VALIDATION_STRINGENCY LENIENT

      /gatk --java-options "-Xmx${return Math.floor(inputs.max_memory*1000/1.074-1)}m" BuildBamIndex
      -I $(inputs.output_basename).sorted.bam
      -VALIDATION_STRINGENCY LENIENT

inputs:
  reference: { type: 'File', secondaryFiles: [^.dict, .fai], doc: "Refernce fasta, with fai and dict indices" }
  input_bams: { type: 'File[]', secondaryFiles: [{ pattern: ".bai", required: false },{ pattern: "^.bai", required: false }], doc: "List of BAM files to merge then sort and index." }
  output_basename: { type: 'string?', doc: "Basename for output BAM" }
  enable_tool: { type: 'boolean?', doc: "Should this tool be run? This option may only be used in a workflow." }

  cores: { type: 'int?', default: 1, doc: "Minimum reserved number of CPU cores for the task." }
  max_memory: { type: 'int?', default: 2, doc: "GB of RAM to allocate to the task." }

outputs:
  output:
    type: File
    outputBinding:
      glob: '*.sorted.bam'
    secondaryFiles: [^.bai]
