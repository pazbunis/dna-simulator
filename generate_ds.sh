#!/bin/csh
source /cs/grad/pazbu/paz/dev/python3.5+tf1+gpu/bin/activate.csh
python /cs/grad/pazbu/paz/dev/projects/dna-simulator/single_tf_simulator.py
python /cs/grad/pazbu/proj/seq2dataset/seq2simple_simulations.py

