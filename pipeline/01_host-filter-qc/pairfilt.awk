# Drops an interleaved read pair (8 lines) when either mate has an empty sequence
# line. Removes the zero-length reads that samtools bam2fq sometimes emits, which crash the
# next minimap2 step, without desynchronizing pairs. Pass-through on clean input.
# Usage: ... | samtools fastq ... - | awk -f pairfilt.awk > out.fq

{ b[(NR - 1) % 8 + 1] = $0 }

NR % 8 == 0 {
    if (length(b[2]) > 0 && length(b[6]) > 0)
        for (i = 1; i <= 8; i++) print b[i]
    else
        dropped++
}

END {
    if (NR % 8 != 0)
        print "pairfilt.awk: WARNING — input not a multiple of 8 lines (NR=" NR "); pairing may be broken upstream" > "/dev/stderr"
    print "pairfilt.awk: dropped " dropped + 0 " pair(s) with an empty read" > "/dev/stderr"
}
