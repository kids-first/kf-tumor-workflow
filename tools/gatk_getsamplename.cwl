cwlVersion: v1.1
class: CommandLineTool
id: gatk_getsamplename
doc: |
  Emit a single sample name from the bam header into an output file. The sample name
  is that in the read group (RG) sample (SM) field
  
  Note: If the bam has zero or more than one sample names in the header, this tool
  will error, by design. This tool has not been tested extensively. Most options 
  supported by the GATK are irrelevant for this tool.
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
      /gatk --java-options "-Xmx${return Math.floor(inputs.max_memory*1000/1.074-1)}m" GetSampleName

inputs:
  override_samplename: {type: 'string?', doc: "If the samplename is known ahead of time or different than the one in the header, provide it here. When used in a workflow, this value can be used to skip this task."}
  input_reads: {type: 'File', inputBinding: {prefix: "--input", position: 2}, doc: "BAM/SAM/CRAM file containing reads."}
  output_filename: {type: 'string?', default: "sample_name.txt", inputBinding: {prefix: "--output", position: 2}, doc: "Output file with only the sample name in it."}
  reference_fasta: {type: 'File', secondaryFiles: [{pattern: ".fai", required: true},{pattern: "^.dict", required: true}], inputBinding: {prefix: "--reference", position: 2}, doc: "Reference sequence fasta and index"}
  extra_args: {type: 'string?', inputBinding: {position: 3}, doc: "Any additional arguments for this tool. See GATK Documentation for complete list of options. Example input: --interval-merging-rule OVERLAPPING_ONLY"}
  cores: { type: 'int?', default: 1, doc: "Minimum reserved number of CPU cores for the task." }
  max_memory: { type: 'int?', default: 2, doc: "GB of RAM to allocate to the task." }
outputs:
  output:
    type: 'string'
    outputBinding:
      glob: $(inputs.output_filename)
      loadContents: true
      outputEval: $(self[0].contents)
