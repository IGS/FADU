using FADU
using Logging
using Test

const DATA_DIR = joinpath(@__DIR__, "data")
const TINY_ARGS = ["-g", joinpath(DATA_DIR, "tiny.gff3"), "-b", joinpath(DATA_DIR, "tiny.bam"),
                   "-s", "reverse", "-f", "CDS", "-a", "ID"]

# Run FADU quietly and return the exit code and the lines of the counts file
function run_fadu(extra_args::Vector{String}=String[])
    mktempdir() do out_dir
        exit_code = with_logger(NullLogger()) do
            FADU.main([TINY_ARGS; "-o"; out_dir; extra_args])
        end
        return exit_code, readlines(joinpath(out_dir, "tiny.counts.txt"))
    end
end

@testset "FADU" begin
    @testset "version matches Project.toml" begin
        @test FADU.VERSION_NUMBER == string(pkgversion(FADU))
    end

    expected = readlines(joinpath(DATA_DIR, "expected.counts.txt"))

    @testset "tiny dataset matches expected output" begin
        exit_code, lines = run_fadu()
        @test exit_code == 0
        @test lines == expected
    end

    @testset "--no_output_header omits header line" begin
        exit_code, lines = run_fadu(["--no_output_header"])
        @test exit_code == 0
        @test lines == expected[2:end]
    end

    @testset "--remove_multimapped runs" begin
        exit_code, lines = run_fadu(["--remove_multimapped"])
        @test exit_code == 0
        @test length(lines) == length(expected)
    end

    @testset "missing BAM file returns non-zero exit code" begin
        mktempdir() do out_dir
            args = ["-g", joinpath(DATA_DIR, "tiny.gff3"), "-b", joinpath(DATA_DIR, "missing.bam"), "-o", out_dir]
            @test with_logger(NullLogger()) do
                FADU.main(args)
            end == 1
        end
    end
end
