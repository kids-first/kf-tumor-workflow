cwlVersion: v1.1
class: CommandLineTool
id: gatk_splitintervals
doc: |
  This tool takes in intervals via the standard arguments of IntervalArgumentCollection
  and splits them into interval files for scattering. The resulting files contain
  equal number of bases.

  Standard GATK engine arguments include -L and -XL, interval padding, and interval
  set rule etc. For example, for the -L argument, the tool accepts GATK-style intervals
  (.list or .intervals), BED files and VCF files. See --subdivision-mode parameter
  for more options.
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/gatk:4.1.1.0'
  - class: ResourceRequirement
    ramMin: $(inputs.ram * 1000)
    coresMin: $(inputs.cores)

baseCommand: []
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      /gatk --java-options "-Xmx${return Math.floor(inputs.max_memory*1000/1.074-1)}m" SplitIntervals

inputs:
  interval_list: {type: 'File?', inputBinding: {prefix: "--intervals", position: 2}, doc: "One or more genomic intervals as a list in a file."}
  scatter_ct: {type: 'int?', default: 50, inputBinding: {prefix: "--scatter-count", position: 2}, doc: "number of output interval files to split into"}
  reference_fasta: {type: 'File', secondaryFiles: [{pattern: ".fai", required: true},{pattern: "^.dict", required: true}], inputBinding: {prefix: "--reference", position: 2}, doc: "Reference sequence fasta and index"}
  output_directory: {type: 'string?', default: 'intervals', inputBinding: {prefix: "--output", position: 2}, doc: "The directory into which to write the scattered interval sub-directories."}
  extra_args: {type: 'string?', inputBinding: {position: 3}, doc: "Any additional arguments for this tool. See GATK Documentation for complete list of options. Example input: --interval-merging-rule OVERLAPPING_ONLY"}
  cores: { type: 'int?', default: 2, doc: "Minimum reserved number of CPU cores for the task." }
  max_memory: { type: 'int?', default: 4, doc: "GB of RAM to allocate to the task." }
outputs:
  split_intervals:
    type: 'File[]'
    outputBinding:
      glob: $(inputs.output_directory)/*.interval_list
