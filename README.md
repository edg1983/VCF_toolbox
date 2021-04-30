# VCF toolbox

A set of small utilities in nim for fast operations on VCF files.
At the moment this includes

- update_info  :   Simple manipulation of INFO fields (select, drop, fill)
- annotate_multiallele:   Add an INFO flag to annotate multi-allelic variants
- deepvar_fix  :   Various fixes and for multisample deepvar decomposed VCF created by GLnexus

## udpate_info
- Allows to select or drop specific INFO fields.
- Given an INFO tag already presetn in the VCF, the fill option allows to fill it with a default value when absent.

## annotate_multiallelic
Process a **NOT** decomposed VCF and annotate multiallelic variant with a flag
- MULTIALLELIC_INDEL: multiallelic var involving at least one indel
- MULTIALLELIC_SNV: multiallelic var involving only SNVs

After decomposing, this can be useful to keep track of SNVs that were actually part of multiallelic vars involving and INDEL (likely artifacts)

## deepvar_fix 
Process standard VCF to:
- fix missing GQ values (.) setting them to 0
- set GT to ./. when DP is zero
- annotate allele balance (AB) from allele depth (AD) for the ALT allele
- convert half genotypes (1/., ./1, ./0, 0/.) to standard genotypes



**NOTES**

- Multiallelic annotation is meant to be applied to VCF before decomposition
- All other tools expect the input VCF to be decomposed so that only 1 ALT is present per variant
