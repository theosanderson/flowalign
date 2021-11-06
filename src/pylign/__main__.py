import multiprocessing
import mappy as mp
import tqdm
from . import functions_based_on_sam_2_fasta
from . import helpers
from . import pyligner
import sys

#Create argparser
import argparse


yield_aligned = pyligner.yield_aligned
    
    
def main():
    parser = argparse.ArgumentParser(description='Alignment')
    parser.add_argument('input', help='Input fasta to align')
    parser.add_argument('--reference',
                        help='Input reference sequence',
                        required=True)
    parser.add_argument('--output', help='Output aligned fasta')
    parser.add_argument('--threads',
                        help='Number of threads',
                        required=False,
                        default=multiprocessing.cpu_count(), type=int)

    args = parser.parse_args()

    aligned = yield_aligned(args.input, args.reference, args.threads)

    if args.output:
        output_file = open(args.output, "wt")
    else:
        output_file = sys.stdout
    for name, result in tqdm.tqdm(aligned):
        output_file.write(f">{name}\n{result}\n")

if __name__ == "__main__":
    main()
    