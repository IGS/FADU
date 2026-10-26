# Notes for package author

## Updating the Manifest.toml file

Project.toml is a plain Julia environment (no name/uuid/version), so the Manifest can be regenerated directly from the repository root with the Julia version you want to target (currently v1.12):

1. `rm fadu_pkgs/Manifest.toml` - This needs to be removed so the next step will resolve a new one
2. `julia --project=fadu_pkgs -e 'using Pkg; Pkg.resolve(); Pkg.instantiate(); Pkg.precompile()'`
3. Run `bash test_data/tiny/run_test.sh` to make sure the output has not changed
