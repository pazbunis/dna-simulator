from os import listdir
from os.path import join
import motifs as mf
import utils
import random


# read motifs
motif_dir = '/cs/grad/pazbu/paz/dev/projects/dna-simulator/motifs/'
motif_files = sorted(listdir(motif_dir))
motifs_and_headers = []
for motif_name in motif_files:
    motif_path = join(motif_dir, motif_name)
    h, pwm = mf.get_motif_pwm(motif_path)
    motifs_and_headers.append((h, pwm))

# read background sequences
seq_path = '/cs/grad/pazbu/paz/dev/projects/data/enhancers/NEnhancers_more_than_15Kb_from_Enhancer_peaks.tfa'
seq_length = 500
seqs_and_headers = []
for header, seq in utils.middle_subseqs(seq_path, seq_length):
    seqs_and_headers.append((header, seq))

idxs = list(range(len(seqs_and_headers)))
random.shuffle(idxs)
pos_samples = [seqs_and_headers[i] for i in range(len(seqs_and_headers)) if i%2 == 0]
neg_samples = [seqs_and_headers[i] for i in range(len(seqs_and_headers)) if i%2 != 0]

# embed motif in center
h, pwm = motifs_and_headers[0]
for i in range(len(pos_samples)):
    header, seq = pos_samples[i]
    replaced_seq = mf.replace_at(seq, mf.sample_motif(pwm), 250)
    pos_samples[i] = (header, replaced_seq)


# write files
pos_out_path = 'pos.fasta'
neg_out_path = 'neg.fasta'
with open(pos_out_path, 'w') as f:
    for header, seq in pos_samples:
        f.write('>' + header + '\n')
        f.write(seq + '\n')

with open(neg_out_path, 'w') as f:
    for header, seq in neg_samples:
        f.write('>' + header + '\n')
        f.write(seq + '\n')
