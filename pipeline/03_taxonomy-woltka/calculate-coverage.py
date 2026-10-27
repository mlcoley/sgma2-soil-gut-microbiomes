#!/usr/bin/env python
# copied from https://github.com/qiita-spots/qp-woltka/blob/222d68d84e0000fce63c5b4ead567e04ae4dbaf2/qp_woltka/util.py
# BSD 3-clause; license text in THIRD_PARTY_LICENSES at the repository root.
from collections import defaultdict
import sys

def _merge_ranges(files):
    # the lines below are borrowed from zebra filter but they are sligthly
    # modified; mainly the autocompress parameter was deleted so it always
    # autocompresses
    class SortedRangeList:
        def __init__(self):
            self.ranges = []

        def add_range(self, start, end):
            self.ranges.append((start, end))
            self.compress()

        def compress(self):
            # Sort ranges by start index
            self.ranges.sort(key=lambda r: r[0])

            new_ranges, start_val, end_val = [], None, None
            for r in self.ranges:
                if end_val is None:
                    # case 1: no active range, start active range.
                    start_val = r[0]
                    end_val = r[1]
                elif end_val >= r[0] - 1:
                    # case 2: active range continues through this range
                    # extend active range
                    end_val = max(end_val, r[1])
                else:  # if end_val < r[0] - 1:
                    # case 3: active range ends before this range begins
                    # write new range out, then start new active range
                    new_range = (start_val, end_val)
                    new_ranges.append(new_range)
                    start_val = r[0]
                    end_val = r[1]

            if end_val is not None:
                new_range = (start_val, end_val)
                new_ranges.append(new_range)

            self.ranges = new_ranges

        def compute_length(self):
            total = 0
            for r in self.ranges:
                total += r[1] - r[0] + 1
            return total

    merger = defaultdict(SortedRangeList)
    for file in files:
        with open(file, 'r') as f:
            for range_line in f:
                mems = range_line.split()
                gotu = mems.pop(0)
                for srange, erange in zip(*[iter(mems)]*2):
                    merger[gotu].add_range(int(srange), int(erange))

    return merger

def coverage_percentage(files, length_map):
    with open(length_map, 'r') as f:
        length_map = dict()
        for line in f:
            line = line.split()
            length_map[line[0]] = int(line[1])

    merger = _merge_ranges(files)

    lines = []
    for gotu in merger:
        length = merger[gotu].compute_length()
        coverage = float(length)/length_map[gotu]*100
        if coverage > 100:
            raise ValueError(f'{gotu} yielded {coverage}; please contact '
                             'qiita.help@gmail.com')
        lines.append('%s\t%.6f\n' % (gotu, coverage))

    return lines

if __name__ == '__main__':
    with open(sys.argv[3], 'w') as f:
        for line in coverage_percentage([sys.argv[2]], sys.argv[1]):
            f.write(line)
