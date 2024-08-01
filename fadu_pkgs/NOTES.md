# Notes for package author

I am sure there is a better way to do this, but this is how I have done things for the time being.

## Updating the Manifest.toml file

1. `docker run -it --exec /bin/bash FADU`
2. `cd /tmp/.julia/environments/v1.7/`
3. `rm Manifest.toml` - This needs to be removed so a future step will instantiate a new one
4. `mkdir src; cd src; touch FADU.jl` - Since the Project.toml lists FADU as a package, this needs to be done or else none of the package manager steps will work
5. `cd ..`
6. Run `julia` to bring up the REPL
7. Hit `]` to bring up the Pkg manager.  Run `instantiate` which will generate an updated Manifest.toml file that can be copied out of the container.