from os import listdir
from os.path import join
import motifs as mf
import numpy as np
from dnautils import idxs_to_seq


class Simulator:
    def __init__(self, motifs_dir, bg_path, seq_length):
        self.motifs_dir = motifs_dir
        self.bg_path = bg_path
        self.seq_length = seq_length

    def generate_sample(self, n_samples):
        headers, pwms = self.load_motifs()
        bg_seqs = self.generate_background_seqs(n_samples)
        embedded_seqs = self.embed_motifs(bg_seqs, pwms)
        return embedded_seqs

    def load_motifs(self):
        motif_paths = listdir(self.motifs_dir)
        headers = []
        pwms = []
        for motif_path in motif_paths:
            header, pwm = mf.get_motif_pwm(join(self.motifs_dir, motif_path))
            headers.append(header)
            pwms.append(pwm)
        return headers, pwms

    def generate_background_seqs(self, n_samples):
        # uniform negatives:
        print("generating bg sequences...")
        seqs = []
        for i in range(n_samples):
            seq = idxs_to_seq(np.random.randint(0, 4, self.seq_length))
            seqs.append(seq)
        return seqs

    def embed_motifs(self, bg_seqs, pwms):
        pass

    def save(self, seqs, labels, path):
        np.save(join(path, 'seqs.npy'), seqs)
        np.save(join(path, 'labels.npy'), labels)

# # read motif_pwms
# motif_dir = '/cs/grad/pazbu/paz/dev/projects/dna-simulator/JASPAR_motifs/npy'
# motif_files = sorted(listdir(motif_dir))
# motifs_and_headers = []
# for motif_name in motif_files:
#     motif_path = join(motif_dir, motif_name)
#     pwm = np.load(motif_path)
#     motifs_and_headers.append((motif_name, pwm))
#
# # read background sequences
# # seq_path = '/cs/grad/pazbu/paz/dev/projects/data/enhancers/NEnhancers_more_than_15Kb_from_Enhancer_peaks.tfa'
# seq_length = 1000
# seqs_and_headers = []
# # for header, seq in utils.middle_subseqs(seq_path, seq_length):
# #     seqs_and_headers.append((header, seq))
#
# # uniform negatives:
# for i in range(50000):
#     header = 'random_motif'
#     seq = idxs_to_seq(np.random.randint(0, 4, seq_length))
#     seqs_and_headers.append((header, seq))
#
#
# # seqs_and_headers = seqs_and_headers[:20000]
#
# idxs = list(range(len(seqs_and_headers)))
# random.shuffle(idxs)
# pos_samples = [seqs_and_headers[i] for i in range(len(seqs_and_headers)) if i%2 == 0]
# neg_samples = [seqs_and_headers[i] for i in range(len(seqs_and_headers)) if i%2 != 0]
#
# # embed motif in center
# h1, pwm1 = motifs_and_headers[1]
# #h2, pwm2 = motifs_and_headers[1]
# for i in range(len(pos_samples)):
#     idx = np.random.randint(0,28)
#     offset_idx = random.randint(-50, 50)
#     header, seq = pos_samples[i]
#     replaced_seq = mf.replace_at(seq, mf.sample_motif(pwm1), 250 + offset_idx)
#     #replaced_seq = mf.replace_at(replaced_seq, mf.sample_motif(pwm2), 150 + offset_idx)
#     pos_samples[i] = (header, replaced_seq)
#
#
# # write files
# pos_out_path = 'single_motif_50offset_revcomp.10KP.pos.fasta'
# neg_out_path = 'single_motif_50offset_revcomp.10KP.neg.fasta'
# with open(pos_out_path, 'w') as f:
#     for header, seq in pos_samples:
#         f.write('>' + header + '\n')
#         f.write(seq + '\n')
#
# with open(neg_out_path, 'w') as f:
#     for header, seq in neg_samples:
#         f.write('>' + header + '\n')
#         f.write(seq + '\n')
#
#
