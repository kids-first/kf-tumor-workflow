cwlVersion: v1.2
class: CommandLineTool
id: deepsomatic_vcf_stats_report
doc: "Deepsomatic VCF Stats Report"
requirements:
  - class: InlineJavascriptRequirement
  - class: ShellCommandRequirement
  - class: ResourceRequirement
    ramMin: $(inputs.max_memory*1000)
    coresMin: $(inputs.cpu)
  - class: DockerRequirement
    dockerPull: 'google/deepsomatic:1.7.0'
baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      vcf_stats_report
inputs:
  input_vcf: { type: 'File', inputBinding: { position: 2, prefix: "--input_vcf"}, doc: "Path to the input VCF." }
  num_records: { type: 'int?', inputBinding: { position: 2, prefix: "--num_records"}, doc: "Maximum number of VCF lines to read. If negative, read the whole VCF." }
  outfile_base: { type: 'string', inputBinding: { position: 2, prefix: "--outfile_base"}, doc: "Base path for output. The HTML report will be created at <outfile_base>.visual_report.html." }
  title: { type: 'string?', inputBinding: { position: 2, prefix: "--title"}, doc: "Title to show at the top of the HTML report. Default: Use the sample name from the VCF." }

  cpu:
    type: 'int?'
    default: 1
    doc: "Number of CPUs to allocate to this task."
  ram:
    type: 'int?'
    default: 2
    doc: "Maximum GB of RAM to allocate for this tool."
outputs:
  stats_html: 
    type: File
    outputBinding:
      glob: "*.html" 
