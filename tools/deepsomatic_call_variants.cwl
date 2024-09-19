cwlVersion: v1.2
class: CommandLineTool
id: deepsomatic_call_variants
doc: "Deepsomatic Call Variants"
requirements:
  - class: InlineJavascriptRequirement
  - class: ShellCommandRequirement
  - class: ResourceRequirement
    ramMin: $(inputs.ram*1000)
    coresMin: $(inputs.cpu)
  - class: DockerRequirement
    dockerPull: 'google/deepsomatic:1.7.0-gpu'
  - class: InitialWorkDirRequirement
    listing: $(inputs.examples)
baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      call_variants
inputs:
  checkpoint: { type: 'string?', inputBinding: { position: 2, prefix: "--checkpoint"}, doc: "Path to the TensorFlow model checkpoint." }
  checkpoint_custom: { type: 'File?', inputBinding: { position: 2, prefix: "--checkpoint"}, doc: "Custom TensorFlow model checkpoint." }
  examples: { type: 'File[]', secondaryFiles: [{pattern: ".example_info.json", required: true}], doc: "Required. tf.Example protos containing DeepVariant candidate variants in TFRecord format, as emitted by make_examples. Can be a comma-separated list of files, and the file names can contain wildcard characters." }
  examples_name: { type: 'string', inputBinding: { position: 2, prefix: "--examples" }, doc: "Short name convention for all example files (e.g. filename@10.gz)" }
  outfile: { type: 'string', inputBinding: { position: 2, prefix: "--outfile"}, doc: "Required. Destination path where we will write output candidate variants with additional likelihood information in TFRecord format of CallVariantsOutput protos." }

  # Catchall for Remaining Options
  extra_args: { type: 'string?', inputBinding: { position: 2, shellQuote: false }, doc: "Any additional args for deepsomatic make_call_variants" }

  cpu:
    type: 'int?'
    default: 8
    doc: "Number of CPUs to allocate to this task."
  ram:
    type: 'int?'
    default: 32
    doc: "Maximum GB of RAM to allocate for this tool."
outputs:
  output: 
    type: File[]
    outputBinding:
      glob: "*cvo-*.tfrecord.gz"
