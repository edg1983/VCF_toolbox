import argparse
import os
import tables
import strformat
import ./modules/manipulate_info
import ./modules/multiallele_anno
import ./modules/deepvar_fix
import ./modules/extract_vars

proc main*() =
  type pair = object
    f: proc(dropfirst:bool)
    description: string

  var dispatcher = {
    "update_info": pair(f:manipulate_info.main, description:"Simple manipulation of INFO fields (select, drop, fill)"),
    "annotate_multiallele": pair(f:multiallele_anno.main, description:"Add an INFO flag to annotate multi-allelic variants"),
    "deepvar_fix": pair(f:deepvar_fix.main, description:"Various fixes and for multisample deepvar decomposed VCF created by GLnexus"),
    "extract_vars": pair(f:extract_vars.main, description:"Extract variants based on a list of VIDs")
    }.toOrderedTable

  var args = commandLineParams()
  
  if len(args) == 0 or not (args[0] in dispatcher):
    stderr.write_line "\nCommands: "
    for k, v in dispatcher:
      echo &"  {k:<13}:   {v.description}"
    if len(args) > 0 and (args[0] notin dispatcher) and args[0] notin @["-h", "-help"]:
      echo &"unknown program '{args[0]}'"
    quit ""

  dispatcher[args[0]].f(false)

when isMainModule:
  main()