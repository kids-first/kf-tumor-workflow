cwlVersion: v1.2
class: CommandLineTool
id: xenome_index
doc: "Tool to create an index for classifying host and graft reads"
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: ${ return inputs.ram * 1000 }
    coresMin: $(inputs.cores)
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/gossamer:1.0.0'
baseCommand: [mkdir]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      $(inputs.output_basename)
      && xenome index
      -T $(inputs.cores)
      -M ${ return Math.floor(inputs.ram * (4/5)) }
      -H $(inputs.host_fasta.path)
      -G $(inputs.graft_fasta.path)
      -v
      --prefix $(inputs.output_basename)/$(inputs.output_basename)
      && tar -I pigz -cf $(inputs.output_basename).tgz $(inputs.output_basename)

inputs:
  host_fasta: {type: File, doc: "Fasta file of organism hosting foreign tissue"}
  graft_fasta: {type: File, doc: "Fasta file of organism that tissue was transferred from"}
  output_basename: string
  cores: {type: "int?", doc: "Num cores to use", default: 16}
  ram: {type: "int?", doc: "Mem to use in GB", default: 64}


outputs:
  output_xenome:
    type: File
    outputBinding:
      glob: $(inputs.output_basename).tgz
