from itertools import groupby
import sys
import math

def fastaread(fasta_name):
    f = open(fasta_name)
    faiter = (x[1] for x in groupby(f, lambda line: line.startswith(">")))
    for header in faiter:
        header = next(header)[1:].strip()
        seq = "".join(s.strip() for s in next(faiter))
        yield header, seq


def center_header(header, start_offset, length):
    chrom, location = header.split(':')
    added_info = ''
    if '|' in location:
        location, added_info = location.split('|', 1)
        added_info = '|' + added_info
    start, _ = location.split('-')
    start = int(start) + start_offset
    end = int(start) + length
    location = str(start) + '-' + str(end)
    return chrom + ':' + location + added_info


def middle_subseqs(path_in, target_length):
    """
    :param path_in:
    :return: the middle target_length letters from the sequences in the input file
    """
    faiter = fastaread(path_in)
    for header, seq in faiter:
        l = len(seq)
        seq = seq.upper()
        if l < target_length:
            sys.stderr.write('target sequence length is longer than a sequence in the file.\n')
            # exit(1)
        else:
            start_idx = math.floor((l - target_length) // 2)
            header = center_header(header, start_idx, target_length)
            yield header, seq[start_idx:start_idx + target_length]