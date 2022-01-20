cwlVersion: v1.1
class: CommandLineTool
id: gatk_genomicsdbimport
doc: |
  Import VCFs to GenomicsDB
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
      gatk --java-options "-Xmx${return Math.floor(inputs.max_memory*1000/1.074-1)}m" GenomicsDBImport

inputs:
  reference_fasta: { type: 'File', secondaryFiles: [{ pattern: ".fai", required: true },{ pattern: "^.dict", required: true }], inputBinding: { prefix: "--reference", position: 2 }, doc: "Reference sequence fasta and index" }
  input_vcfs: { type: { type: array, items: File, inputBinding: { prefix: "--variant" }}, inputBinding: { position: 2 }, doc: "VCF files to be imported to GenomicsDB. Each file must containdata for only a single sample." }
  intervals: { type: 'File?', inputBinding: { prefix: "--intervals", position: 2 }, doc: "One or more genomic intervals over which to operate" }
  genomicsdb_path: { type: 'string?', inputBinding: { prefix: "--genomicsdb-workspace-path", position: 2 }, doc: "Workspace for GenomicsDB. Can be a POSIX file system absolute or relative path" }
  extra_args: { type: 'string?', inputBinding: { position: 3, shellQuote: false }, doc: "Any additional arguments for this tool. See GATK Documentation for complete list of options. Example input: --interval-merging-rule OVERLAPPING_ONLY" }
  cores: { type: 'int?', default: 1, doc: "Minimum reserved number of CPU cores for the task." }
  max_memory: { type: 'int?', default: 2, doc: "GB of RAM to allocate to the task." }
outputs:
  output:
    type: 'Directory'
    outputBinding:
      glob: $(inputs.genomicsdb_path)
