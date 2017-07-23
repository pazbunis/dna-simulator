from os import listdir
from os.path import join
import motifs as mf
import numpy as np
import utils
from dnautils import idxs_to_seq
import random


# read motif_pwms
motif_dir = 'motif_pwms/'
motif_files = sorted(listdir(motif_dir))
motifs_and_headers = []
for motif_name in motif_files:
    motif_path = join(motif_dir, motif_name)
    h, pwm = mf.get_motif_pwm(motif_path)
    motifs_and_headers.append((h, pwm))

# read background sequences
seq_path = 'NEnhancers_more_than_15Kb_from_Enhancer_peaks.tfa'
seq_length = 500
seqs_and_headers = []
c = 0
for header, seq in utils.middle_subseqs(seq_path, seq_length):
    seqs_and_headers.append((header, seq))
    c+=1
    if c > 100000:
        break
# uniform negatives:
#for i in range(100000):
#    header = 'random:30000-30500'
#    seq = idxs_to_seq(np.random.randint(0, 4, 500))
#    seqs_and_headers.append((header, seq))

random.shuffle(seqs_and_headers)
pos_samples = [seqs_and_headers[i] for i in range(len(seqs_and_headers)) if i%2 == 0]
neg_samples = [seqs_and_headers[i] for i in range(len(seqs_and_headers)) if i%2 != 0]

# embed motif in center
h1, pwm1 = motifs_and_headers[1]
#h2, pwm2 = motifs_and_headers[1]
for i in range(len(pos_samples)):
    offset_idx = random.randint(-50, 50)
    header, seq = pos_samples[i]
    replaced_seq = mf.replace_at(seq, mf.sample_motif(pwm1), 250 + offset_idx)
#    replaced_seq = mf.replace_at(replaced_seq, mf.sample_motif(pwm2), 150 + offset_idx)
    pos_samples[i] = (header, replaced_seq)


# write files
pos_out_path = 'single.motif.50offset.50K.pos.fasta'
neg_out_path = 'single.motif.50offset.50K.neg.fasta'
with open(pos_out_path, 'w') as f:
    for header, seq in pos_samples:
        f.write('>' + header + '\n')
        f.write(seq + '\n')

with open(neg_out_path, 'w') as f:
    for header, seq in neg_samples:
        f.write('>' + header + '\n')
        f.write(seq + '\n')
