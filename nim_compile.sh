#Compile the single executable
#nim compile -d:nimEmulateOverflowChecks AB_annotator.nim

#Compile a static build with no dependencies
singularity exec --bind /gpfs3/well/gel/HICF2/software/BRC_tools/VCF_toolbox:/gpfs3/well/gel/HICF2/software/BRC_tools/VCF_toolbox --bind /gpfs3/well/gel/HICF2/software/BRC_tools/VCF_toolbox:/load/ /well/gel/HICF2/software/singularity/musl-hts-nim.sif /usr/local/bin/nsb  -s /gpfs3/well/gel/HICF2/software/BRC_tools/VCF_toolbox/src/VCF_toolbox.nim --nimble-file /gpfs3/well/gel/HICF2/software/BRC_tools/VCF_toolbox/VCF_toolbox.nimble --  -d:danger -d:release
