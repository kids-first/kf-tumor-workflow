cwlVersion: v1.2
class: CommandLineTool
id: xenome_classify
doc: "Tool to classify reads and output only host using a prebuilt index"
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: ${ return inputs.ram * 1000 }
    coresMin: $(inputs.cores)
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/gossamer:1.0.0'
baseCommand: [tar, -xzvf]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      $(inputs.xenome_index.path) --strip-components 1
  - position: 3
    shellQuote: false
    valueFrom: >-
      && xenome classify
      --pairs
      -M ${ return Math.floor(inputs.ram * (4/5)) }
      --output-filename-prefix $(inputs.output_basename)
  - position: 5
    shellQuote: false
    valueFrom: >-
      > $(inputs.output_basename).Xenome.Classification.Stats.txt

inputs:
  xenome_index: {type: File, doc: "Xenome index made form host and graft fasta"}
  host_name: {type: "string?", inputBinding: {prefix: --host-name, position: 4}, doc: 'name to use describing model organism receiving graft', default: "mouse"}
  graft_name: {type: "string?", inputBinding: {prefix: --graft-name, position: 4}, doc: 'name to use describing organism that grafted tissue came from', default: "human"}
  cores: {type: "int?", inputBinding: {prefix: -T, position: 4}, doc: "Num cores to use", default: 8}
  ram: {type: "int?", doc: "Mem to use in GB", default: 8}
  idx_prefix: {type: string, inputBinding: {prefix: -P, position: 4}, doc: "String prefix of index files when decompressed"}
  fastq_reads: 
    type:
      type: array
      items: File
      inputBinding:
        prefix: -i
    inputBinding:
      position: 4
  output_basename: string

outputs:
  graft_fastqs:
    type: 'File[]'
    outputBinding:
      glob: '*$(inputs.graft_name)*.fastq'
  host_fastqs:
    type: 'File[]'
    outputBinding:
      glob: '*$(inputs.host_name)*.fastq'
  output_stats:
    type: File
    outputBinding:
      glob: '*.Xenome.Classification.Stats.txt'
