# Deepvariant fixes
# Author: Edoardo Giacopuzzi
# Several small-fixes to cohort VCF
# - fix missing GQ values (.) setting them to 0
# - annotate allele balance (AB) from allele depth (AD) for the ALT allele
# - convert half genotypes 1/. or ./1 to standard het genotypes
# - set GT to ./. when DP is zero
# Expects input VCF to be decomposed so that only 1 ALT is present

import argparse
import hts
import times
import math
import system
from ./utils import progress_counter
# import sequtils

proc getADvalues(ad_values: seq[int32], geno: seq[int32], DP: int32): seq[int32] =
    # var ad_values = ad_string.split(",")
    var allele_counts = ad_values
    if allele_counts[0] < 0: allele_counts[0] = 0
    
    if (geno[0] == 4 and geno[1] == 0) or (geno[0] == 0 and geno[1] == 4):
        allele_counts[0] = DP - allele_counts[1]
    else: discard

    return allele_counts

proc calculateAB(ad_values: seq[int32]): float32 =
    var ab_value: float32
    if sum(ad_values) > 0:
        ab_value = ad_values[1] / sum(ad_values)
    else:
        ab_value = 0
    return round(ab_value,2)

proc main* (dropfirst:bool=false) =

    # Parse arguments
    var p = newParser("VCF_fixes"):
        option("-v", "--vcf", help="path to VCF/BCF")
        option("-o", "--out_vcf", help="output VCF")
        flag("-a", "--computeAB", help="fill in AB tag for ALT allele")
        flag("-g", "--fixGQ", help="Fix missing GQ setting it to zero")
        flag("-m", "--setMissing", help="Set GT to missing when DP equals zero")
        flag("-t", "--fixhalfGT", help="Fix half GT setting 1/. and ./1 to het calls")
        option("-d", "--DPfield", help="Name of the format field containing read depth", default="DP")
        option("-f", "--ADfield", help="Name of the format field containing allele detpth", default="AD")

    var argv = commandLineParams()
    if len(argv) > 0 and argv[0] == "deepvar_fix":
        argv = argv[1..argv.high]
    if len(argv) == 0: argv = @["--help"]
    var opts = p.parse(argv)

    # Quit if help is used
    if opts.help:
        quit "END of help message", QuitSuccess

    # Quit if no operation is set
    if not (opts.computeAB or opts.fixGQ or opts.setMissing or opts.fixhalfGT):
        quit "No operation required. Please set at least one. See help for options", QuitFailure    

    echo "Input VCF: ", opts.vcf
    echo "Output: ", opts.out_vcf
    echo "AD field: ", opts.ADfield, "; DP field: ", opts.DPfield
    echo "Operations: "
    echo "\tFix GQ: ", opts.fixGQ
    echo "\tSet GT to missing on DP zero: ", opts.setMissing
    echo "\tUpdate half GT: ", opts.fixhalfGT 
    echo "\tFill AB value:", opts.computeAB
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
    if opts.computeAB == true:
        doAssert wtr.header.add_format("AB", "A", "Float", "Allele balance for the ALT allele") == Status.OK
    doAssert(wtr.write_header())

    var
        changed_GQ: int64
        changed_GT: int64
        changed_halfGT: int64
        start_time = cpuTime()
        t0 = cpuTime()
        n = 0
        interval = 5000
        dps = new_seq[int32](len(v.samples))
        gts = new_seq[int32](len(v.samples))
        gqs = new_seq[int32](len(v.samples))
        ads = new_seq[int32](len(v.samples))
        abs = new_seq[float32](len(v.samples))

    for rec in v:
        # Variants counter
        n = n + 1
        progress_counter(n, interval, t0)
    
        # Get relevant values from VCF record
        doAssert rec.format.get(opts.ADfield, ads) == Status.OK
        doAssert rec.format.get(opts.DPfield, dps) == Status.OK
        doAssert rec.format.get("GT", gts) == Status.OK
        doAssert rec.format.get("GQ", gqs) == Status.OK
        # var genos = rec.format.genotypes(gts)

        # Apply requested operation to sample data
        for i in 0..high(dps):

            # Fix GQ values
            if opts.fixGQ:
                if gqs[i] < 0: 
                    gqs[i] = 0
                    changed_GQ = changed_GQ + 1

            # Set GT to missing for samples with DP == 0
            if opts.setMissing:
                if dps[i] == 0: 
                    gts[i*2] = 0
                    gts[i*2+1] = 0
                    changed_GT = changed_GT + 1
            
            # Compute AB
            if opts.computeAB:
                var ad_counts = getADvalues(@[ads[i*2],ads[i*2+1]], @[gts[i*2],gts[i*2+1]], dps[i]) 
                abs[i] = calculateAB(ad_counts)
                #echo "DP: ", dps[i], "\tAD: ", ads[i], "\tGT: ", genos[i], "\tAB: ", abs[i]
            
            # Fix half GT to het genotypes
            if opts.fixhalfGT:
                # gts are values for the single alleles in each GT 
                # 0 = miss, 2 = 0, 4 = 1, 6 = 2
                if gts[i*2] == 4 and gts[i*2+1] == 0:
                    gts[i*2] = 4
                    gts[i*2+1] = 2
                elif gts[i*2] == 0 and gts[i*2+1] == 4:
                    gts[i*2] = 2
                    gts[i*2+1] = 4

        # Update variant record and write to output file
        if opts.computeAB: 
            doAssert rec.format.set(opts.ADfield, ads) == Status.OK
            doAssert rec.format.set("AB", abs) == Status.OK
        if opts.fixhalfGT or opts.setMissing:
            doAssert rec.format.set("GT", gts) == Status.OK
        if opts.fixGQ:
            doAssert rec.format.set("GQ", gqs) == Status.OK
        
        doAssert wtr.write_variant(rec)

    close(v)
    close(wtr)
    var elapsed_time = cpuTime() - start_time

    echo "Finished ", n, " variants in ", elapsed_time, " secs"
    if opts.computeAB: echo "AB annotated for all variants"
    if opts.fixGQ: echo "Invalid GQ values corrected: ", changed_GQ
    if opts.setMissing: echo "GT with DP0 set to miss: ", changed_GT
    if opts.fixhalfGT: echo "Updated half GT: ", changed_halfGT 
