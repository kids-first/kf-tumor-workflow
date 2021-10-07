cwlVersion: v1.0
class: CommandLineTool
id: gatk_splitintervals
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/gatk:4.1.1.0'
  - class: ResourceRequirement
    ramMin: $(inputs.ram * 1000)
    coresMin: $(inputs.cpu)

baseCommand: []
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      /gatk --java-options "-Xmx${return Math.floor(inputs.max_memory*1000/1.074-1)}m" SplitIntervals

inputs:
  interval_list: {type: 'File?', inputBinding: {prefix: "--intervals", position: 2}, doc: "One or more genomic intervals as a list in a file."}
  scatter_ct: {type: 'int?', default: 50, inputBinding: {prefix: "--scatter-count", position: 2}, doc: "number of output interval files to split into"}
  reference_fasta: {type: 'File', inputBinding: {prefix: "--reference", position: 2}, doc: "Reference sequence"}
  output_directory: {type: 'string?', default: '.', inputBinding: {prefix: "--output", position: 2}, doc: "The directory into which to write the scattered interval sub-directories."}
  extra_args: {type: 'string?', inputBinding: {position: 3}, doc: "Any additional arguments for SplitIntervals. See GATK Documentation for complete list of options."}
outputs:
  output:
    type: 'File[]'
    outputBinding:
      glob: '$(inputs.output_directory)/*.interval_list'

$namespaces:
  sbg: https://sevenbridges.com
