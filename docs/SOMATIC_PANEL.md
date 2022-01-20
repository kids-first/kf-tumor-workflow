# Mutect2 Subworkflow

## Minimal Run
- Mutect
- MergeMutectStats
- MergeVcfs
- FilterMutectCalls

## Mutect Filter Support (creates contamination and segmentation tables for FilterMutectCalls)
__Requires `exac_common_vcf` input__

- GetPileupSummaries
- GatherPileupSummaries
- CalculateContamination

## Orientation Bias Mixture Model Filter (creates artifact prior tables for FilterMutectCalls)
__Requires `run_orientation_bias_mixture_model_filter` set to `true`__

- Additional f1r2 outputs from Mutect
- LearnReadOrientationModel

## Mutect BAM Creation
__Requires `make_bamout` set to `true`__

- Additional BAM outputs from Mutect
- GatherBamFiles
- SortSam
- BuildBamIndex

## Filter Alignment Artifacts
__Requires `bwa_mem_index_image` input__

- FilterAlignmentArtifacts

## Annotation
__Requires `run_annotation` set to `true`__

- SelectVariants (PASS)
- vt normalize VCF
- bcftools strip annotations
- add strelka2 standard fields (Requires add_common_fields set to true; Not used for Mutect VCFs)
- VEP annotate
- bcftools gnomad annotate
- VariantFiltration
- hotspot annotation
- vcf2maf public
- bcftools hardfilter
- vcf2maf protected
- rename protected files
- rename public files


# Create Somatic Panel of Normals Workflow

## Workflow Outline
For each `input_normal_aligned, input_normal_name` pair:
- Mutect2 Minimal Run

Using the collected outputs from above:
- SplitIntevals
- GenomicsDBImport
- CreateSomaticPanelOfNormals
- MergeVcfs

## Workflow Description
The Create Somatic Panel of Normals workflow is a CWL port of the [GATK Best Pratice WDL](https://github.com/broadinstitute/gatk/blob/master/scripts/mutect2_wdl/mutect2_pon.wdl).
