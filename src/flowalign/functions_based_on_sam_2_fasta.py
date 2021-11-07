"""This code is entirely (with the exception of any errors in adaptation)
based on datafunk's sam_2_fasta created by Ben Jackson.

https://github.com/cov-ert/datafunk/blob/master/datafunk/sam_2_fasta.py

"""

import re, sys

lambda_dict = {
    'M': (lambda query_start, ref_start, length, seq:
          (query_start + length, ref_start + length, seq[
              query_start:query_start + length])),
    'I': (lambda query_start, ref_start, length, seq:
          (query_start + length, ref_start, '')),
    'D': (lambda query_start, ref_start, length, seq:
          (query_start, ref_start + length, '-' * length)),
    'N': (lambda query_start, ref_start, length, seq:
          (query_start, ref_start + length, '-' * length)),
    'S': (lambda query_start, ref_start, length, seq:
          (query_start + length, ref_start, '')),
    'H': (lambda query_start, ref_start, length, seq:
          (query_start, ref_start, '')),
    'P': (lambda query_start, ref_start, length, seq:
          (query_start, ref_start, '')),
    '=': (lambda query_start, ref_start, length, seq:
          (query_start + length, ref_start + length, seq[
              query_start:query_start + length])),
    'X': (lambda query_start, ref_start, length, seq:
          (query_start + length, ref_start + length, seq[
              query_start:query_start + length]))
}


def split_sam_cigar_operation(one_operation):
    """
    o is a tuple with the format(operation, size)
    e.g. ('M', 2377) or ('D', 1)
    """
    type = one_operation[-1:]
    size = int(one_operation[:-1])
    o = (type, size)
    return (o)


def split_sam_cigar(cigar):
    """
    l is a list of strings (e.g. ['1000M', '4I'])
    """
    r = re.compile('\d{1,}[A-Z]{1}')
    l = re.findall(r, cigar)
    return (l)


def get_sam_cigar_operations(cigar):
    """
    operations is a list of tuples that correspond to
    operations to apply in order:
    [('M', 10000),('I', 3),('M', 19763)]
    """
    operations_raw = split_sam_cigar(cigar)
    operations = [split_sam_cigar_operation(x) for x in operations_raw]
    return (operations)


def get_one_string(hit, qseq, rlen):
    """
    Transform one line of the SAM alignment into sample sequence in unpadded
    reference coordinates (insertions relative to the reference are omitted,
    but are logged to a global dict if log_inserts = True).
    """

    # CIGAR STRING
    CIGAR = hit.cigar_str

    POS = hit.r_st 

    # Query seq:
    SEQ = qseq[hit.q_st:hit.q_en] #this is where mappy works differently to
                                  #minimap2 which would do softclipping with "S"


    if POS < 0:
        return (None)

    # parse the CIGAR string to get the operations:
    operations = get_sam_cigar_operations(CIGAR)

    # left-pad the new sequence with gaps if required
    new_seq = '*' * POS

    # then build the sequence:
    qstart = 0
    rstart = POS
    for o in operations:
        operation = o[0]
        size = o[1]

        # based on this CIGAR operation, call the relavent lambda function
        # from the dict of lambda functions, returns sequence to be appended
        # and the next set of coordinates
        new_qstart, new_rstart, extension = lambda_dict[operation](qstart,
                                                                   rstart,
                                                                   size, SEQ)

        new_seq = new_seq + extension

        qstart = new_qstart
        rstart = new_rstart

    rightpad = '*' * (rlen - len(new_seq))

    new_seq = new_seq + rightpad

    return (new_seq)


def check_and_get_flattened_site(site, QNAME):
    """
    A per-site check that there isn't any ambiguity between
    alignments within a single sequence
    """
    site = set(site)

    check = sum([x.isalpha() for x in site])
    if check > 1:
        sys.stderr.write('ambiguous overlapping alignment: ' + QNAME + '\n')
        return ('N')

    # because {A, C, G, T} > {-} > {*}, we can use max()
    base = max(site)
    return (base)


def swap_in_gaps_Ns(seq, pad):
    """
    replace internal runs of '*'s with 'N's
    and external runs of '*'s with '-'s
    """
    r_internal = re.compile('[A-Z]\*+[A-Z]')
    for x in re.findall(r_internal, seq):
        seq = seq.replace(x, x[0] + x[1:-1].replace('*', 'N') + x[-1])

    r_left = re.compile('^\*+[A-Z]')
    m_left = re.search(r_left, seq)
    if m_left:
        g_left = m_left.group()
        if pad:
            seq = seq.replace(g_left,
                              g_left[:-1].replace('*', 'N') + g_left[-1])
        else:
            seq = seq.replace(g_left,
                              g_left[:-1].replace('*', '-') + g_left[-1])

    r_right = re.compile('[A-Z]\*+$')
    m_right = re.search(r_right, seq)
    if m_right:
        g_right = m_right.group()
        if pad:
            seq = seq.replace(g_right,
                              g_right[0] + g_right[1:].replace('*', 'N'))
        else:
            seq = seq.replace(g_right,
                              g_right[0] + g_right[1:].replace('*', '-'))

    return (seq)


def get_seq_from_query(hits,
                       query_seq,
                       rlen,
                       pad,
                       query_name,
                       check_carefully=True):

    block_lines_sites_list = [
        get_one_string(hit, query_seq, rlen) for hit in hits
    ]
    block_lines_sites_list = [x for x in block_lines_sites_list if x]

    if len(block_lines_sites_list) == 1:
        seq_flat_no_internal_gaps = swap_in_gaps_Ns(block_lines_sites_list[0],
                                                    pad=pad)
        return (seq_flat_no_internal_gaps)

    elif len(block_lines_sites_list) > 1:
        # # as an alternative to check_and_get_flattened_site() we can flatten
        # # the site with no checks (about three times as fast):
        # flattened_site_list = [max(x) for x in zip(*[list(x) for x in block_lines_sites_list])]

        if check_carefully:
            flattened_site_list = [
                check_and_get_flattened_site(x, query_name)
                for x in zip(*[list(x) for x in block_lines_sites_list])
            ]
        else:
            flattened_site_list = [
                max(x) for x in zip(*[list(x) for x in block_lines_sites_list])
            ]
        seq_flat = ''.join(flattened_site_list)

        # replace central '*'s with 'N's, and external '*'s with '-'s
        seq_flat_no_internal_gaps = swap_in_gaps_Ns(seq_flat, pad=pad)
        return (seq_flat_no_internal_gaps)

    else:
        return (None)
