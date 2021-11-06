from . import helpers
from . import functions_based_on_sam_2_fasta
import multiprocessing
import mappy as mp
import sys

def yield_aligned(input_filename,reference_filename, threads):
    aligner = helpers.get_aligner(reference_filename)

    if len(aligner.seq_names) > 1:
        raise ValueError("Only one reference sequence allowed")

    rlen = len(aligner.seq(aligner.seq_names[0], start=0, end=0x7fffffff))

    global do_alignment # THIS IS OBVIOUSLY A HORRIBLE WAY TO DO THINGS, OH WELL..
                        # (something is necessary because of how multiprocessing works)
    def do_alignment(the_tuple):
        name, seq, qual = the_tuple  # This is the tuple returned by fast_x_read
        seq = seq.replace("-", "")
        #strip the sequence of Ns at start and end
        seq = helpers.strip_starting_and_ending_ns(seq)
        hits = aligner.map(seq)
        
        result = functions_based_on_sam_2_fasta.get_seq_from_query(
            hits, seq, rlen, True, name)
        return name, result
    reader = mp.fastx_read(input_filename)
   
    if threads>1:
        pool = multiprocessing.Pool(threads)
        
        aligned = pool.imap(do_alignment, reader, chunksize=20)
    else:
        aligned = map(do_alignment, reader)
    return aligned
