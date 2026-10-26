# Notes for package author

## Updating the Manifest.toml file

Manifest.toml pins the dependency versions used for releases, Docker and bioconda. Regenerate it from the repository root with the Julia version you want to target (currently v1.12):

1. `rm Manifest.toml` - This needs to be removed so the next step will resolve a new one
2. `julia --project=. -e 'using Pkg; Pkg.resolve(); Pkg.instantiate(); Pkg.precompile()'`
3. Run `julia --project=. -e 'using Pkg; Pkg.test()'` to make sure the output has not changed

## Making a release

1. Bump `version` in Project.toml and `VERSION_NUMBER` in src/FADU.jl (a test fails if they differ)
2. Add a CHANGELOG.md entry and update "Current FADU release" in README.md and docs/index.md
3. Commit, then tag with `git tag -a vX.Y.Z -m "vX.Y.Z"` and push the tag
