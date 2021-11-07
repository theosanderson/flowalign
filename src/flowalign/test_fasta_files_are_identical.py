
import argparse
import os
import sys
from Bio import SeqIO
#import md5
import hashlib
import gzip
import tqdm
import re
def strip_ending_ns(seq):
    """
    Strip the ending Ns from a sequence.
    """

    r_right = re.compile('[^N]N+$')
    m_right = re.search(r_right, seq)
    if m_right:
        seq = seq[0:m_right.start()+1]
    return seq
    



# define argparser
parser = argparse.ArgumentParser(description='Test if two fasta files are identical.')
parser.add_argument('fasta_file_1', help='First fasta file.')
parser.add_argument('fasta_file_2', help='Second fasta file.')
args = parser.parse_args()

hashed_seqs = {}

handle1 = gzip.open(args.fasta_file_1, 'rt') if args.fasta_file_1.endswith('.gz') else open(args.fasta_file_1, 'r')
# read fasta 1 and hash sequences:
for seq_record in tqdm.tqdm(SeqIO.parse(handle1, "fasta"),desc= "Reading first file"):
    seq_id = seq_record.id
    seq = str(seq_record.seq)
    seq = strip_ending_ns(seq)
    hashed_seqs[seq_id] = hashlib.md5(seq.encode('utf-8')).hexdigest()

handle2 = gzip.open(args.fasta_file_2, 'rt') if args.fasta_file_2.endswith('.gz') else open(args.fasta_file_2, 'r')
# check fasta 2 matches

not_found = 0
matches = 0
mismatches =0
for seq_record in tqdm.tqdm(SeqIO.parse(handle2, "fasta"),desc= "Reading second file"):
    seq_id = seq_record.id
    
    if seq_id not in hashed_seqs:
        #print("Sequence id {} not found in fasta 1.".format(seq_id))
        not_found += 1
        
    else:
        seq = str(seq_record.seq)
        seq = strip_ending_ns(seq)
 
        if hashed_seqs[seq_id] != hashlib.md5(seq.encode('utf-8')).hexdigest():
            print("Sequence id {} does not match.".format(seq_id))
            mismatches += 1
        else:
            matches += 1

print(f"Results:")
print(f"{matches} matches.")
print(f"{mismatches} mismatches.")
print(f"{not_found} sequences not found in fasta 1.")
        
