import multiprocessing
import mappy as mp
import tqdm
from . import functions_based_on_sam_2_fasta
from . import helpers
import sys

#Create argparser
import argparse

class Aligner(object):
    def __init__(self, reference_filename):
        self.mappy_aligner = helpers.get_aligner(reference_filename)
        if len(self.mappy_aligner) > 1:
            raise ValueError("More than one sequence found for reference file: {}".format(reference_filename))
        # TODO finish this

    
    

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
                        default=multiprocessing.cpu_count())

    args = parser.parse_args()

    aligner = helpers.get_aligner(args.reference)

    if len(aligner.seq_names) > 1:
        print(
            "Multiple reference sequences found. Please provide a single reference sequence.",
            file=sys.stderr)

        sys.exit(1)

    rlen = len(aligner.seq(aligner.seq_names[0], start=0, end=0x7fffffff))


    def do_alignment(the_tuple):
        name, seq, qual = the_tuple  # This is the tuple returned by fast_x_read
        seq = seq.replace("-", "")
        hits = aligner.map(seq)
        result = functions_based_on_sam_2_fasta.get_seq_from_query(
            hits, seq, rlen, True, name)
        return name, result


    pool = multiprocessing.Pool(args.threads)
    reader = mp.fastx_read(args.input)
    aligned = pool.imap(do_alignment, reader, chunksize=10)

    if args.output:
        output_file = open(args.output, "wt")
    else:
        output_file = sys.stdout
    for name, result in tqdm.tqdm(aligned):
        output_file.write(f">{name}\n{result}\n")

if __name__ == "__main__":
    main()