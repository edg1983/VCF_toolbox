# Annotate multiallele
# Author: Edoardo Giacopuzzi
# Process standard VCF before decompose and annotate multi-allelic vars
# Two tags are added to info
# 1. MULTIALLELIC_INDEL, when at least 1 allele is an indel this is likely an artifact
# 2. MULTIALLELIC_SNV, only SNVs alleles, depending on VCF these call may be good

import argparse
import hts
import times
import system
from ./utils import progress_counter
# import sequtils

proc main* (dropfirst:bool=false) =
    # Parse arguments
    var p = newParser("VCF_fixes"):
        option("-v", "--vcf", help="path to VCF/BCF")
        option("-o", "--out_vcf", help="output VCF")

    var argv = commandLineParams()
    if len(argv) > 0 and argv[0] == "annotate_multiallele":
        argv = argv[1..argv.high]
    if len(argv) == 0: argv = @["--help"]
    var opts = p.parse(argv)

    # Quit if help is used
    if opts.help:
        quit "END of help message", QuitSuccess

    echo "Input VCF: ", opts.vcf
    echo "Output: ", opts.out_vcf
    echo "############################################"
    echo "Start VCF processing"

    # Open VCF for reading
    var v:VCF
    if not open(v, opts.vcf):
        quit "Couldn't open input VCF", QuitFailure

    # Open new VCF for writing
    var wtr:VCF
    doAssert(open(wtr, opts.out_vcf, mode="w", threads=5))
    wtr.header = v.header
    doAssert wtr.header.add_info("MULTIALLELIC_INDEL", "0", "Flag", "Variant is part of a multi-allelic variant including at least one indel") == Status.OK
    doAssert wtr.header.add_info("MULTIALLELIC_SNV", "0", "Flag", "Variant is part of a multi-allelic variant including only SNVs") == Status.OK
    doAssert(wtr.write_header())

    var
        multiallelic_vars: int32
        start_time = cpuTime()
        t0 = cpuTime()
        n = 0
        interval = 5000
        multiallele_tag: string

    for rec in v:
        # Variants counter
        n = n + 1
        progress_counter(n, interval, t0)
        
        multiallele_tag = "MULTIALLELIC_SNV"
        
        # Get REF and ALT alleles
        if len(rec.ALT) > 1:
            multiallelic_vars = multiallelic_vars + 1
            for a in rec.ALT:
                if len(a) != len(rec.REF):
                    multiallele_tag = "MULTIALLELIC_INDEL"     
            
            doAssert rec.info.set(multiallele_tag,true) == Status.OK
        
        doAssert wtr.write_variant(rec)

    close(v)
    close(wtr)
    var elapsed_time = cpuTime() - start_time

    echo "Finished ", n, " variants in ", elapsed_time
    echo "Annotated multiallelic records: ", multiallelic_vars
