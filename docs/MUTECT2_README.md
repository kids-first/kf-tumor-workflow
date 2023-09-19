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

The assembled haplotypes and locally realigned reads will be written as BAM to
this file if requested. Really for debugging purposes only. Note that the output
here does not include uninformative reads so that not every input read is emitted
to the bam. Turning on this mode may result in serious performance cost for the
caller. It's really only appropriate to use in specific areas where you want to
better understand why the caller is making specific calls.

## Filter Alignment Artifacts
__Requires `bwa_mem_index_image` input__

- FilterAlignmentArtifacts

## Annotation
__Requires `run_annotation` set to `true`__

- SelectVariants (PASS|HotSpotAllele=1)
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


# Mutect2 Production Workflow

## Default Functionality
- Minimal Run
- Mutect Filter Support
- Orientation Bias Mixture Model Filter
- Annotation

Any segments not used by default can still be used activated bt the user.
We do not create BAMs for storage cost reasons.
FilterAlignmentArtifacts is still an experimental tool and not recommended for production pipelines! The current tool still has several issues as of GATK 4.2.0.0.

## Somatic vs Tumor-Only
The production workflow is capable of running in both somatic and tumor-only modes.

To run the workflow in somatic mode, the following inputs are necessary:

- `input_tumor_aligned`
- `input_normal_aligned`

To run the workflow in tumor-only mode, the following inputs are necessary:

- `input_tumor_aligned`
- `panel_of_normals`

## WGS vs WXS
The production workflow is capable of running both WGS and WXS inputs.

### Running the tool in WGS mode:

__Requires `WGS` for `wgs_or_wxs`__

- Sets `--BREAK_BANDS_AT_MULTIPLES_OF` to `80000000` in IntervalListTools

### Running the tool in WXS mode:

__Requires `WXS` for `wgs_or_wxs`__

- Sets `--BREAK_BANDS_AT_MULTIPLES_OF` to `0` in IntervalListTools
- Additional `--disable-adaptive-pruning` flag set in Mutect
