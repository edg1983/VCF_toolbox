# Package
version     = "0.1.0"
author      = "Edoardo Giacopuzzi"
description = "Simple operations on VCF files"
license     = "MIT"

# Deps
requires "nim >= 0.10.0", "hts >= 0.3.4"
requires "argparse 0.10.1"

srcDir = "src"
bin = @["VCF_toolbox"]

skipDirs = @["test"]
