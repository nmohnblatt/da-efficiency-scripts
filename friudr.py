#!/usr/bin/env python
import math
from schemes import *

GRINDING = 20
RO_QUERIES = 60
STATISTICAL_SECURITY = 40
FRI_SOUNDNESS = STATISTICAL_SECURITY + RO_QUERIES - GRINDING

REVIEWER = False
REVIEWERPROXDELTA = 0.25


def sizeMerkleOpening(numleafs, tuplesize, fsize):
    '''
    assume a Merkle tree that represents numleafs tuples of elements
    where each element has size fsize and one tuple contains tuplesize
    many elements. Each leaf of the tree contains one such tuple.
    size of one opening (i.e., Merkle path + element)
    of a Merkle tree that represents n elements
    and each element has size fsize.
    '''
    tupleItself = tuplesize * fsize
    sibling = tuplesize * fsize
    treedepth = math.ceil(math.log2(numleafs))
    copath = (treedepth - 1) * HASH_SIZE
    return tupleItself + sibling + copath


def friAuthSize(domainsize, rate, fsize, batchsize, fanin, basedimension):
    '''
    size of the information needed to open one position
    in the FRI base layer, including all Merkle paths
    domainsize is the number of elements in the base layer
    fsize is the field size, i.e., size of one element
    '''

    size = 0
    # batching phase if FRI_BATCH_SIZE > 1:
    # we need to open one symbol ( = FRI_BATCH_SIZE field elements)
    # of the interleaved code, so one leaf of the Merkle tree
    # representing the batch
    if batchsize > 1:
        size += sizeMerkleOpening(domainsize, batchsize, fsize)

    # now assume that we have already opened the batching
    # i.e., it remains to open everything from oracle G_0
    # to oracle G_r. In every oracle, the FRI verifier
    # queries on a set of fanin positions.
    # We put this entire set into the same Merkle leaf.
    # We do not put the final oracle in a Merkle tree.
    # Instead, we put it in plain into the commitment.
    # This makes sense as we have to open it entirely.
    ncurr = domainsize
    while ncurr * rate > basedimension:
        numleafs = ncurr // fanin
        size += sizeMerkleOpening(numleafs, fanin, fsize)
        ncurr = numleafs
    return size


def friNumRounds(mink, fanin, basedimension):
    '''
    Function to compute the number of rounds, given the number of field
    elements needed to represent the data, the fanin of the folding, 
    and the dimension at which we stop folding
    '''

    # if we do no round, we can represent basedimension many elements
    # if we do one round, we can represent basedimension * fanin many elements
    # if we do r rounds, we can represent basedimension * (fanin**r) many elements
    dimension = basedimension
    rnd = 0
    while dimension < mink:
        dimension *= fanin
        rnd += 1
    return rnd


def friNumRepetitions(rate, domainsize, fsize, batchsize, fanin):
    '''
    Function to compute the number of repetitions of the query phase
    '''

    # first make sure that the soundness error
    # induced by LuckySet (e.g., distortion) is small
    maxbf = max(fanin, batchsize)
    logeps1 = 1 + math.ceil(math.log2(domainsize * (maxbf - 1))) - fsize
    assert (logeps1 <= - FRI_SOUNDNESS)

    # now determine number of repetitions such that the
    # soundness error related to the query phase is small
    # recall: the soundness error for L repetitions is
    # (1-delta^*/F)^L, and we need to get it below 2^{-FRI_SOUNDNESS}
    deltastar = 0.5 * (1.0 - rate)

    if REVIEWER:
        deltastar = REVIEWERPROXDELTA

    base = 1.0 - deltastar
    logbase = math.log2(base)
    assert (logbase < 0)
    L = - FRI_SOUNDNESS / logbase
    return math.ceil(L)


def makeFRIUDRScheme(datasize, invrate=4, fsize=128, verbose=False):
    # determine k. Should be "compatible" with the fan-in
    # we need k to be at least ceil(datasize / fsize)
    minfe = math.ceil(datasize / fsize)
    if verbose:
        print("Need at least dimension minfe = " + str(minfe) +
              " field elements to represent the data.")

    # call algorithm to find good batchsize, fanin, and base dimension
    (batchsize, fanin, basedimension) = friGoodParameters(minfe, fsize, invrate)

    mink = math.ceil(minfe / batchsize)
    if verbose:
        print("With batch size B = " + str(batchsize) +
              ", we need at least dimension mink = " + str(mink))
        print("Use fanin F = " + str(fanin) +
              " and base dimension = " + str(basedimension))
    # now determine the number of rounds to get at least
    # dimension mink in the base layer
    r = friNumRounds(mink, fanin, basedimension)
    if verbose:
        print("Need " + str(r) + " rounds.")
    # with that, we get the actual k and n
    k = basedimension * (fanin ** r)
    n = invrate * k
    rate = 1.0 / invrate
    if verbose:
        print("Need dimension k = " + str(k) +
              " and evaluation domain size n = " + str(n) + ".")
    # determine the number of repetitions we need
    # to get good soundness guarantees
    L = friNumRepetitions(rate, n, fsize, batchsize, fanin)
    if verbose:
        print("Need " + str(L) + " repetitions of the query phase.")

    # determine the size of one opening
    authsize = friAuthSize(n, rate, fsize, batchsize, fanin, basedimension)

    # now compile the scheme
    rs = makeRSCode(fsize, k, n)
    if REVIEWER:
        deltastar = 0.5 * (1.0 - rate)
        rs.reception = n-math.ceil((deltastar - REVIEWERPROXDELTA)*n)
        rs.samples = samples_from_reception(40, rs.reception, rs.codeword_len)

    # we include all openings for the final layer in the commitment, and no Merkle root for it
    # if we do batching, we need one root more
    final = basedimension * fsize
    openings = L * authsize
    roots = r * HASH_SIZE + (batchsize > 1) * HASH_SIZE

    opening_overhead = authsize-batchsize*fsize

    if REVIEWER and batchsize > 1:
        opening_overhead = sizeMerkleOpening(n, batchsize, fsize)

    if REVIEWER and batchsize == 1:
        ncurr = n
        numleafs = ncurr // fanin
        opening_overhead = sizeMerkleOpening(numleafs, fanin, fsize)

    return Scheme(
        com_size=roots + final + openings,
        code=rs.interleave(batchsize),
        opening_overhead=opening_overhead
    )

# --------------------------------------------------------------------------#
#                           OPTIMIZATION SECTION                           #
# --------------------------------------------------------------------------#


def friGoodBatchsize(minfe, fsize, invrate, basedimension, fanin):
    '''
    given the minimum number of field elements we need to represent (minfe),
    the field size (fsize), the inverse rate (invrate), the basedimension, and
    the fanin, this function computes a good batchsize. Good means that the
    batchsize minimizes (in a certain range) the size of a single opening
    '''
    batchsizerange = range(1, 257)
    batchsize = 1
    mink = math.ceil(minfe / batchsize)
    r = friNumRounds(mink, fanin, basedimension)
    minauthsize = friAuthSize(basedimension * (fanin**r) * invrate,
                              1.0 / invrate, fsize, batchsize, fanin, basedimension)
    for b in batchsizerange:
        mink = math.ceil(minfe / b)
        r = friNumRounds(mink, fanin, basedimension)
        currauthsize = friAuthSize(
            basedimension * (fanin**r) * invrate, 1.0 / invrate, fsize, b, fanin, basedimension)
        if currauthsize <= minauthsize:
            batchsize = b
            minauthsize = currauthsize
    return batchsize


def friGoodParameters(minfe, fsize, invrate):
    '''
    given the minimum number of field elements we need to represent (minfe),
    the field size (fsize) and the inverse rate (invrate), this function
    computes (batchsize, fanin, basedimension) for FRI that works reasonably
    well. This is for sure not always the optimal setting, especially if a
    specific metric should be optimized, e.g., communication per query
    '''

    # overall idea is to minimize the gap between the dimension on the largest layer
    # and the dimension we would actually need to represent minfe elements. That is,
    # we minimimize gap = basedimension * fanin^rounds - minfe / batchsize, ensuring
    # that gap >= 0. To do so, we try a few reasonable fanins and base dimensions
    faninrange = [4, 8, 16]
    basedimensionrange = [2, 4, 6, 8, 16, 32, 64, 128]

    # start minimazation loop. Iterate over all combinations (fanin, basedimension)
    optfanin = 0
    optbasedimension = 0
    optbatchsize = 0
    mingap = -1
    for fanin in faninrange:
        for basedimension in basedimensionrange:
            # if we want to compute the gap for the pair (fanin, basedimension),
            # we need to know a suitable batchsize first. To find it, we want
            # to minimize the size of an opening, i.e., minimize friAuthSize
            batchsize = friGoodBatchsize(
                minfe, fsize, invrate, basedimension, fanin)
            mink = math.ceil(minfe / batchsize)

            # determine the number of rounds that we need now
            r = friNumRounds(mink, fanin, basedimension)
            # compute gap for this fanin, basedimension, and batchsize
            gap = basedimension * (fanin**r) - mink
            # update if it is better
            if mingap == -1 or (gap >= 0 and gap <= mingap):
                mingap = gap
                optfanin = fanin
                optbasedimension = basedimension
                optbatchsize = batchsize

    return (optbatchsize, optfanin, optbasedimension)
