import mappy as mp
import re

def get_aligner(reference_filename):

    #asm5:
    # mo->a = 1, mo->b = 19, mo->q = 39, mo->q2 = 81, mo->e = 3, mo->e2 = 1, mo->zdrop = mo->zdrop_inv = 200;
    # #asm20:
    # mo->a = 1, mo->b = 4, mo->q = 6, mo->q2 = 26, mo->e = 2, mo->e2 = 1, mo->zdrop = mo->zdrop_inv = 200;
	# 		io->w = 10;

 
    opt_a = 1
    opt_b = 19
    opt_q = 39
    opt_q2 =81
    opt_e = 3
    opt_e2 = 1
    opt_sc_ambi = 0


    MM_F_SAM_HIT_ONLY = 0x40000000
    MM_F_NO_PRINT_2ND = 0x4000

    scoring = (opt_a, opt_b, opt_q, opt_e, opt_q2, opt_e2, opt_sc_ambi)
    aligner = mp.Aligner(reference_filename,
                         preset="asm5",
                         extra_flags=MM_F_SAM_HIT_ONLY | MM_F_NO_PRINT_2ND,
                         scoring=scoring)  # load or build index
    if not aligner:
        raise Exception("ERROR: failed to load/build index")
    return aligner

def strip_starting_and_ending_ns(seq):
    """
    Strip the starting and ending Ns from a sequence.
    """
    r_left = re.compile('^N+[^N]')
    m_left = re.search(r_left, seq)
    if m_left:
        seq = seq[m_left.end()-1:]
    r_right = re.compile('[^N]N+$')
    m_right = re.search(r_right, seq)
    if m_right:
        seq = seq[0:m_right.start()+1]
    return seq
    
