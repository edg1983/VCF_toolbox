# Manipulate INFO tags
# Author: Edoardo Giacopuzzi
# Subsetting INFO fields in VCF or filling missing INFO tags

import argparse
import hts
import times
import strutils
import ./utils
import system
import sequtils
# import sequtils

type
    Fill_item = tuple
        tag: string
        value: string

proc any_empty_element(s: seq[string]): bool {.inline.} =
    result = any(s, proc (x: string): bool = return x == "")

proc is_empty[T](s: seq[T]): bool {.inline.} =
    if s.len == 1 and s[0] == "":
        return true
    else:
        return false

proc main* (dropfirst:bool=false) =
    var p = newParser("update_info"):
        option("-v", "--vcf", help="path to VCF/BCF")
        option("-o", "--out_vcf", help="output VCF")
        option("-d", "--drop", help="comma separated list of INFO fields to drop")
        option("-s", "--select", help="comma separated list of INFO fields to select")
        option("-f", "--fill", multiple=true, help="(tag,value). INFO tag to fill with value when missing. Can be specified multiple times")

    var argv = commandLineParams()
    if len(argv) > 0 and argv[0] == "update_info":
        argv = argv[1..argv.high]
    if len(argv) == 0: argv = @["--help"]
    var opts = p.parse(argv)

    # Quit if help is used
    if opts.help:
        quit "END of help message", QuitSuccess

    var
        drop_fields: seq[string]
        select_fields: seq[string]
        fill_items: seq[Fill_item]        

    drop_fields = split(opts.drop, ",")
    select_fields = split(opts.select, ",")
    
    # Check at least one operation is required and drop / select are not set together
    if is_empty(drop_fields) and is_empty(select_fields) and is_empty(opts.fill):
        quit "At least one of --drop, --select or --fill must be set", QuitFailure
    
    if not is_empty(drop_fields) and not is_empty(select_fields):
        quit "Only one of --drop or --select can be set", QuitFailure
    
    if not is_empty(drop_fields) and any_empty_element(drop_fields):
        quit "Parsing tags from --drop resulted in at least one empty element", QuitFailure
    if not is_empty(select_fields) and any_empty_element(select_fields):
        quit "Parsing tags from --select resulted in at least one empty element", QuitFailure

    if not is_empty(opts.fill):
        for x in opts.fill:
            let x_values = split(x, ",")
            if any_empty_element(x_values):
                quit "Parsing tags from --fill resulted in at least one empty element", QuitFailure
            fill_items.add((x_values[0], x_values[1])) 
        
    for x in fill_items:
        if x.tag in drop_fields:
            quit x.tag & "field specified by --fill is to be removed by --drop", QuitFailure

    echo "Input VCF: ", opts.vcf
    echo "Output: ", opts.out_vcf
    if not is_empty(drop_fields): echo "Drop fields: ", opts.drop
    if not is_empty(select_fields): echo "Select fields: ", opts.select
    if fill_items.len > 0: echo "Fill fields: ", join(opts.fill, "; ")
    echo "############################################"
    echo "Start VCF processing"
    
    # Open VCF for reading
    var v:VCF
    if not open(v, opts.vcf):
        quit "couldn't open input VCF", QuitFailure
    
    for x in fill_items:
        var itype: string
        try:
            itype = v.header.get(x.tag, BCF_HEADER_TYPE.BCF_HL_INFO)["Type"]
        except KeyError:
            quit x.tag & "field specified by --fill but not found in VCF header", QuitFailure
        
        if check_info_datatype(x.value, itype):
            echo x.tag," will be filled with value ", x.value," of type ",itype
        else:
            quit "Can not convert " & x.value & " to type " & itype, QuitFailure

    # Open new VCF for writing
    var wtr:VCF
    doAssert(open(wtr, opts.out_vcf, mode="w", threads=8))
    wtr.header = v.header

    if len(drop_fields) > 0: 
        for d in drop_fields:
            doAssert wtr.header.remove_info(d) == Status.OK

    doAssert(wtr.write_header())

    var
        drop_var = 0
        select_var = 0
        fill_var = 0
        start_time = cpuTime()
        t0 = cpuTime()
        n = 0
        interval = 5000

    for rec in v:
        # Variants counter
        n = n + 1
        progress_counter(n, interval, t0)
        var info_fields: seq[string]

        for f in rec.info.fields(): info_fields.add(f.name)

        # Get relevant values from VCF record
        if not is_empty(drop_fields):
            for d in drop_fields:
                if d in info_fields:
                    doAssert rec.info.delete(d) == Status.OK
                    drop_var = drop_var + 1   
        
        if not is_empty(select_fields):
            for f in info_fields:
                if f notin select_fields:
                    doAssert rec.info.delete(f) == Status.OK
                    select_var = select_var + 1
        
        if fill_items.len > 0:
            for x in fill_items:
                if x.tag notin info_fields:
                    var fill_value = x.value
                    doAssert rec.info.set(x.tag,fill_value) == Status.OK
                    fill_var = fill_var + 1

        doAssert wtr.write_variant(rec)

    close(v)
    close(wtr)
    var elapsed_time = cpuTime() - start_time

    echo "Finished ", n, " variants in ", elapsed_time
    if not is_empty(drop_fields): echo "Vars with field(s) removed due to drop: ", drop_var
    if not is_empty(select_fields): echo "Vars with field(s) removed due to select: ", select_var
    if len(fill_items) > 0: echo "Vars with field(s) annotated due to fill: ", fill_var
