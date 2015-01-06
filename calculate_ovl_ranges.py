#!/usr/bin/env python
import math
import sys

def get_chunks_per_file(num_total_chunks):
    """This returns the largest n such that the nth triangular number
    is <= num_total_chunks.

    Probably!
    """

    return int(math.floor(.5*(math.sqrt(8*num_total_chunks + 1) - 1)))

def format_ovl_range(ovl_range):
    """Print the ovl_range as parameters for celera assembler's
    overlapInCore.
    """

    return ("-h {hs}-{he} -r {rs}-{re}"
            .format(hs=ovl_range[0][0]+1, he=ovl_range[0][1],
                    rs=ovl_range[1][0]+1, re=ovl_range[1][1]))

def get_ranges(num_reads, num_total_chunks):
    """Get read ranges to be passed to overlapInCore."""

    chunks_per_file = get_chunks_per_file(num_total_chunks)

    bounds = [x*(num_reads/chunks_per_file) for x in xrange(chunks_per_file)]
    bounds += [num_reads]
    starts = bounds[:-1]
    ends = bounds[1:]

    intervals = zip(starts, ends)
    
    ovl_ranges = []

    for hash_int in intervals:
        for ref_int in intervals:
            if ref_int[0] >= hash_int[1]:
                continue

            ovl_ranges.append((hash_int, ref_int))

    return ovl_ranges

def main(num_reads, max_num_chunks):
    ovl_ranges = get_ranges(num_reads, max_num_chunks)
    for ovl_range in ovl_ranges:
        print format_ovl_range(ovl_range)

if __name__ == '__main__':
    main(int(sys.argv[1]), int(sys.argv[2]))
