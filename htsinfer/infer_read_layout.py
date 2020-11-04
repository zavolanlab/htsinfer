"""Infer read layout from sample data."""


def infer():
    """Main function coordinating the execution of all other functions.
    Should be imported/called from main app and return results to it.
    """
    pass
    # implement me


def find_overlaps(motif, read, min_overlap, full_contain=False):
    '''
    Note:
        motif < read
    Args:
        motif = motif sequence. must be string and big character
        read= = read sequence. must be string and big character
        min_overlap = min threshold for match. must be integer
        full_contain=boolean parameter. If True, only report matches that are
                    fully contained in the read
    returns:
        returns a list of tuple(a,b) where:
        a is overlap start position of read. (integer)
        b is fraction of motif that overlaps (real number)
    '''
    # check the type of arguments
    if not isinstance(motif, str):
        raise TypeError('Incorrect argument type: motif')
    if not isinstance(read, str):
        raise TypeError('Incorrect argument type: read')
    if not isinstance(min_overlap, int):
        raise TypeError('Incorrect argument type: min_overlap')
    if not isinstance(full_contain, bool):
        raise TypeError('Incorrect argument type: full_contain')
    # check value of arguments
    if len(motif) == 0 or len(read) == 0:
        raise ValueError('length of motif or read must be bigger than 0')
    if min_overlap < 1:
        raise ValueError('min_overlap must be longer than 1')
    if len(read) < len(motif):
        raise ValueError('read should be >= than motif')
    if not motif.isupper() or not read.isupper():
        raise ValueError('motif or read should be all big characters')
    # compute partial overlaps of the motif at the start of the read
    # 1st case
    partial_overlaps_start = []
    if len(read) >= min_overlap:
        # iterate from min_overlap to make it faster
        partial_overlaps_start = [(0, ov / len(motif))
                                  for ov in range(min_overlap, len(motif))
                                  if read[0:ov] == motif[len(motif)-ov:]]
    # compute matches of the motif inside the read
    # 2nd case
    # you compare the full motif length here
    full_overlaps = []
    full_overlaps = [(i, 1) for i in range(0, len(read) - len(motif)+1)
                     if read[i:i+len(motif)] == motif]
    # compute partial overlaps of the motif at the end of the read
    # 3rd case
    partial_overlaps_end = []
    if len(read) >= min_overlap:
        partial_overlaps_end = [(len(read)-ov, ov / len(motif))
                                for ov in range(min_overlap, len(motif))
                                if read[len(read)-ov:] == motif[0:ov]]
    # print('success')
    # return the list of overlaps
    if full_contain:
        return full_overlaps
    return partial_overlaps_start + full_overlaps + partial_overlaps_end
