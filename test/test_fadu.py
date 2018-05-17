import os
import pytest
import sys
sys.path.append("..")

TEST_BAM = os.path.abspath("test.bam")
TEST_GFF = os.path.abspath("test.gff")

class BamTest():
    pass

class GffTest():
    pass

# Tests to add
# 1) KeyError where properly paired reads have "N" CIGAR string areas where the depth at the position ends up being 0, but one read fully overlaps the other read
# 2) OSError on Mac when reading Samtools depth
