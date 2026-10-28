using FADU
using Logging
using Test

const DATA_DIR = joinpath(@__DIR__, "data")
const TINY_ARGS = ["-g", joinpath(DATA_DIR, "tiny.gff3"), "-b", joinpath(DATA_DIR, "tiny.bam"),
                   "-f", "CDS", "-a", "ID"]
const EXCLUDE_BED = joinpath(DATA_DIR, "tiny.exclude.bed")
const EXCLUDED_FEATURE = "cds-WP_001339389.1"    # Entirely covered by tiny.exclude.bed

# Run FADU quietly and return the exit code and the lines of the counts file
function run_fadu(extra_args::Vector{String}=String[]; stranded::String="reverse")
    mktempdir() do out_dir
        exit_code = with_logger(NullLogger()) do
            FADU.main([TINY_ARGS; "-s"; stranded; "-o"; out_dir; extra_args])
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

    @testset "--exclude_regions_file (stranded=$stranded)" for stranded in ["no", "reverse"]
        _, baseline = run_fadu(; stranded)
        exit_code, lines = run_fadu(["--exclude_regions_file", EXCLUDE_BED]; stranded)
        @test exit_code == 0
        # The excluded feature loses all of its unique coordinates and counts...
        excluded = split(only(filter(l -> startswith(l, EXCLUDED_FEATURE * "\t"), lines)), '\t')
        @test excluded[2] == "0"
        @test excluded[4] == "0.00"
        # ...while every other feature keeps its unique coordinates
        uniq_lens(ls) = Dict(f[1] => f[2] for f in split.(ls[2:end], '\t') if f[1] != EXCLUDED_FEATURE)
        @test uniq_lens(lines) == uniq_lens(baseline)
    end

    @testset "BAM index lookup" begin
        mktempdir() do dir
            # Extensionless BAM must not be mistaken for its own index
            bam = joinpath(dir, "sample")
            cp(joinpath(DATA_DIR, "tiny.bam"), bam)
            @test FADU.find_bam_index(bam) === nothing
            cp(joinpath(DATA_DIR, "tiny.bam.bai"), bam * ".bai")
            @test FADU.find_bam_index(bam) == bam * ".bai"

            # Only the final ".bam" extension is replaced
            sub = mkpath(joinpath(dir, "run.bam_files"))
            bam = joinpath(sub, "x.bam")
            cp(joinpath(DATA_DIR, "tiny.bam"), bam)
            cp(joinpath(DATA_DIR, "tiny.bam.bai"), joinpath(sub, "x.bai"))
            @test FADU.find_bam_index(bam) == joinpath(sub, "x.bai")
        end
    end

    @testset "extensionless BAM runs" begin
        mktempdir() do dir
            bam = joinpath(dir, "sample")
            cp(joinpath(DATA_DIR, "tiny.bam"), bam)
            cp(joinpath(DATA_DIR, "tiny.bam.bai"), bam * ".bai")
            args = ["-g", joinpath(DATA_DIR, "tiny.gff3"), "-b", bam, "-o", dir, "-s", "reverse", "-f", "CDS", "-a", "ID"]
            @test with_logger(NullLogger()) do
                FADU.main(args)
            end == 0
            @test readlines(joinpath(dir, "sample.counts.txt")) == expected
        end
    end

    @testset "missing BAM index returns non-zero exit code" begin
        mktempdir() do dir
            bam = joinpath(dir, "tiny.bam")
            cp(joinpath(DATA_DIR, "tiny.bam"), bam)
            args = ["-g", joinpath(DATA_DIR, "tiny.gff3"), "-b", bam, "-o", dir]
            @test with_logger(NullLogger()) do
                FADU.main(args)
            end == 1
        end
    end

    @testset "--log_level is case-insensitive" begin
        exit_code, lines = run_fadu(["--log_level", "warning"])
        @test exit_code == 0
        @test lines == expected
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
