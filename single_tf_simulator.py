import numpy as np
from simulator import Simulator
import motifs as mf
import utils

class SingleTFSimulator(Simulator):
    """Creates a binary data set where positive sequences are defined by a single motif"""
    def embed_motifs(self, bg_seqs, pwms):
        print("embedding motifs...")
        embedded_seqs = []
        pwm = pwms[0]
        locs = []
        for seq in bg_seqs:
            mf_instance = mf.sample_motif(pwm, random_rev_comp=False)
            loc = np.random.randint(0, 500)
            locs.append(loc)
            embedded_seqs.append(mf.replace_at(seq, mf_instance, loc))
        utils.plot_hist(locs)
        return embedded_seqs

    def generate_sample(self, n_samples):
        headers, pwms = self.load_motifs()
        bg_seqs = self.generate_background_seqs(n_samples)
        pos_seqs = bg_seqs[:n_samples//2]
        neg_seqs = bg_seqs[n_samples//2:]
        embedded_seqs = self.embed_motifs(pos_seqs, pwms)
        labels = np.zeros((n_samples,1))
        labels[:n_samples//2] = 1
        return embedded_seqs+neg_seqs, labels

if __name__ == "__main__":
    motifs_path = '/cs/grad/pazbu/paz/dev/projects/dna-simulator/homer_motifs/single_motif'
    output_path = '/cs/grad/pazbu/paz/dev/projects/dna-simulator/datasets/single_motif'
    seq_length = 1000
    bs = SingleTFSimulator(motifs_path, None, seq_length)
    seqs, labels = bs.generate_sample(20000)
    bs.save(seqs, labels, output_path)
