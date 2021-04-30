# Extract variants base on vids
# Author: Edoardo Giacopuzzi
# Process a decomposed VCF and lookup for variants in input vid file
# vids are string that can be parsed using the expression provided

import argparse
import hts
import times
import system
import re
from ./utils import progress_counter_small
# import sequtils

type
    VAR = object
        chrom: string
        pos: string
        ref_allele: string
        alt_allele: string

proc main* (dropfirst:bool=false) =
    # Parse arguments
    var p = newParser("Subset VCF by vids"):
        option("-i", "--vcf", help="path to VCF/BCF")
        option("-o", "--out_vcf", help="output VCF")
        option("-v", "--ids", help="file with list of vids one per line")
        option("-f", "--format", help="regex format to parse vids", default=r"(chr[0-9XYM]+):([0-9]+)(\w+)>(\w+)")

    var argv = commandLineParams()
    if len(argv) > 0 and argv[0] == "extract_vars":
        argv = argv[1..argv.high]
    if len(argv) == 0: argv = @["--help"]
    var opts = p.parse(argv)

    # Quit if help is used
    if opts.help:
        quit "END of help message", QuitSuccess

    echo "Input VCF: ", opts.vcf
    echo "Output: ", opts.out_vcf
    echo "Variant ids: ", opts.ids
    echo "Variant ID format: ", opts.format
    echo "############################################"
    
    # Compile regexp from opts
    let vid_regexp = re(opts.format)

    # Load vids and make them in a list of var objects
    echo "Reading VIDs from input file"
    var
        f: File
        line : string
        matches: array[4, string]
        n = 0
        skipped = 0
        vars_of_interest: seq[VAR]
    
    try:
        f = open(opts.ids)
    except IOError:
        quit "Couldn't open vids file", QuitFailure

    while f.read_line(line):
        if match(line, vid_regexp, matches): 
            n = n + 1
            vars_of_interest.add(VAR(
                chrom: matches[0], 
                pos: matches[1],
                ref_allele: matches[2],
                alt_allele: matches[3]))
        else:
            skipped = skipped + 1
    f.close()

    echo n, " vids correcly loaded from ", opts.ids
    if skipped > 0: echo "WARN - ", skipped, " can not be parsed based on ", opts.format
    if n == 0: quit "No variants can be loaded from your list", QuitFailure

    # Open VCF for reading
    echo "Start VCF processing"
    var v:VCF
    if not open(v, opts.vcf):
        quit "Couldn't open input VCF", QuitFailure

    # Open new VCF for writing
    var wtr:VCF
    doAssert(open(wtr, opts.out_vcf, mode="w", threads=5))
    wtr.header = v.header
    doAssert(wtr.write_header())

    var
        start_time = cpuTime()
        t0 = cpuTime()
        written_vars = 0
        interval = 50
    
    n = 0
    for vid in vars_of_interest:
        n = n + 1
        progress_counter_small(n, interval, t0)
        for rec in v.query(vid.chrom & ":" & vid.pos & "-" & vid.pos):
            if rec.REF == vid.ref_allele and rec.ALT[0] == vid.alt_allele:
                doAssert wtr.write_variant(rec)
                written_vars = written_vars + 1
            
    close(v)
    close(wtr)
    var elapsed_time = cpuTime() - start_time

    echo "Finished search for ", n, " variants in ", elapsed_time
    echo "Written variants: ", written_vars
