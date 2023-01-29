from . import helpers
from . import functions_based_on_sam_2_fasta
import multiprocessing
import mappy as mp
import sys

def fasta_parser(stream):
    """
    Parses a fasta file and yields tuples of (name, sequence)
    """
    seq=""
    name = ""
    for line in stream:
        if line.startswith(">"):
            if name and seq:
              yield (name, seq,None)
            name = line.strip()[1:]
            seq = ""
        else:
            seq += line.strip()
            
    yield (name, seq,None)

def yield_aligned(input, reference, threads = multiprocessing.cpu_count() ):
    aligner = helpers.get_aligner(reference)

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
        if result is None:
            result = "N" * rlen
        return name, result

    # if input type is string:
    if isinstance(input, str):
        reader = mp.fastx_read(input)
    else:
        reader = fasta_parser(input)
   
    if threads>1:
        pool = multiprocessing.Pool(threads)
        
        aligned = pool.imap(do_alignment, reader, chunksize=20)
    else:
        aligned = map(do_alignment, reader)
    return aligned
