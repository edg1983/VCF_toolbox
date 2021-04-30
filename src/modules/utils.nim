import times
import strformat
import strutils
import parseutils
import system
import math

proc check_info_datatype* (value:string, new_type:string): bool =
    ## Check if a string can be converted to proper VCF INFO datatype
    ## Return true/false if conversion pass / fail
    ## Raise ValueError if new_type is not among allowed types for INFO data
    case new_type
    of "String":
        return true
    of "Char":
        return true
    of "Integer":
        var x: int
        if parseInt(value,x,0) == value.len:
            return true
        else:
            return false
    of "Float":
        var x: float
        if parseFloat(value,x,0) == value.len:
            return true
        else:
            return false
    of "Flag":
        try:
            var x: bool
            x = parseBool(value)
            return true
        except ValueError:
            return false
    else:
        raise newException(ValueError, "Unknown datatype " & new_type & " defined in VCF header")

proc progress_counter* (n:int, interval: var int, t0: var float) {.discardable.} =
    case n
        of 20000: interval = 10000
        of 50000: interval = 25000
        of 150000: interval = 50000
        of 500000: interval = 100000
        of 1000000: interval = 500000
        else: discard
    
    if floorMod(n, interval) == 0:
            var elapsed_time = cpuTime() - t0
            t0 = cpuTime()
            echo n, " variants processed. Last batch: ", interval, " in ", fmt"{elapsed_time:1.2f}", " seconds"

proc progress_counter_small* (n:int, interval: var int, t0: var float) {.discardable.} =
    case n
        of 500: interval = 250
        of 1000: interval = 1000
        of 10000: interval = 5000
        else: discard
    
    if floorMod(n, interval) == 0:
            var elapsed_time = cpuTime() - t0
            t0 = cpuTime()
            echo n, " variants processed. Last batch: ", interval, " in ", fmt"{elapsed_time:1.2f}", " seconds"