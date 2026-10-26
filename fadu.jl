#!/usr/bin/env julia

# Launcher for FADU. Loads the FADU package from this directory, so the
# command `julia fadu.jl [options]` keeps working from a clone of the repository.
pushfirst!(LOAD_PATH, @__DIR__)
using FADU
exit(FADU.main(ARGS))
