import numpy as np
import sys

bases = ['A', 'C', 'G', 'T']


def get_motif_pwm(homer_motif_path):
    with open(homer_motif_path, 'r') as f:
        motif_lines = []
        for line in f.readlines():
            if line.strip() != '':
                motif_lines.append(line.strip())
        header = motif_lines[0]
        pwm = np.zeros((4, len(motif_lines) - 1))
        for i in range(1, len(motif_lines)):
            probs = list(map(lambda str_p: float(str_p), motif_lines[i].split('\t')))
            pwm[:, i-1] = probs
    return header, pwm


def sample_motif(pwm):
    motif = []
    cols = np.shape(pwm)[1]
    for c in range(cols):
        bases_idx = np.random.multinomial(1, pwm[:, c])
        base = bases[np.argmax(bases_idx)]
        motif.append(base)
    return motif


def replace_at(seq, motif, start_idx):
    """
    Replaces the letters in 'seq', starting at 'start_idx', with 'motif'
    :param seq:
    :param motif:
    :param start_idx:
    :return:
    """
    motif_length = len(motif)
    if start_idx + motif_length > len(seq):
        sys.stderr.write('Starting the motif at the given index surpasses the sequence bounds.\n')
        return None
    else:
        seq_arr = list(seq)
        seq_arr[start_idx: start_idx + motif_length] = motif
        return ''.join(seq_arr)