#!/usr/bin/env python3

#from __future__ import print_function

import argparse
import os
import os.path
import numpy as np
import random
import subprocess


def gen_samples(args):
    
    profile = np.genfromtxt(args.profile, delimiter=',').astype(int)
            
    try:
        os.makedirs(args.out_dir)
    except:
        print("Output dir already exists!")
        exit()

    G = profile.shape[0] 
    S = profile.shape[1] - 1

    for cur_s in range(1, S+1):

        print("Making {} from {} samples".format(cur_s, S))

        sample_dir  = os.path.join(args.out_dir, "sample" + str(cur_s))

        if not os.path.exists(sample_dir):
            os.makedirs(sample_dir)

        for cur_g in range(G):
            isolate_name = profile[cur_g, 0]
            seed = random.randint(1, 100)

            for order in [1, 2]:

                input_file = os.path.join(args.input_dir, 'SRR182451{}_{}.fastq.gz'.format(isolate_name, order))
                
                os.system("seqtk sample -s{seed} {input_file} {num_reads} >> {output}".format(
                    seed=seed, 
                    input_file=input_file, 
                    num_reads=profile[cur_g, cur_s], 
                    output=os.path.join(sample_dir, 'R{}.fastq'.format(order))))
                
                #TODO: why is it not working?
                #reads_file = open(os.path.join(curdir, 'R{}.fastq'.format(order)), 'a')
                #subprocess.check_call(["~/tools/seqtk/seqtk sample", "-s", str(seed), input_file, str(num_reads)], shell=True, stdout=reads_file)
    


parser = argparse.ArgumentParser(description="Metagenomic sampler from isolates experiments")

parser.add_argument("profile", type=str, help="File with coverage profiles (first column should be the num of an experiment)")
parser.add_argument("--input-dir", "-i", type=str, required=True, help="Directory with isolates reads")
parser.add_argument("--out-dir", "-o", type=str, required=True, help="Output directory to create")

args = parser.parse_args()
gen_samples(args)
