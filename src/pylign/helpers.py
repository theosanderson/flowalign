import mappy as mp


def get_aligner(reference_filename):
    opt_a = 1
    opt_b = 19
    opt_q = 30
    opt_q2 = 81
    opt_e = 3
    opt_e2 = 1
    opt_sc_ambi = 0

    MM_F_SAM_HIT_ONLY = 0x40000000
    MM_F_NO_PRINT_2ND = 0x4000

    scoring = (opt_a, opt_b, opt_q, opt_e, opt_q2, opt_e2, opt_sc_ambi)
    aligner = mp.Aligner(reference_filename,
                         preset="asm5",
                         extra_flags=MM_F_SAM_HIT_ONLY + MM_F_NO_PRINT_2ND,
                         scoring=scoring)  # load or build index
    if not aligner:
        raise Exception("ERROR: failed to load/build index")
    return aligner
