#!/usr/bin/env python

import math
import sys
import csv
import os

from schemes import *
from friudr import *

# the graphs will be for data sizes 
# i*DATASIZEUNIT for every i in DATASIZERANGE
DATASIZEUNIT = 8000*1000  # Megabytes
DATASIZERANGE = range(1, 156, 15)


def writeCSV(path, d):
    with open(path, mode="w") as outfile:
        writer = csv.writer(outfile, delimiter=',')
        for x in d:
            writer.writerow([x, d[x]])


def writeScheme(name, makeScheme):
    '''
    Writes the graphs for a given scheme into a csv file
    The scheme should be specified by a function makeScheme
    that takes as input the datasize.
    '''

    commitment = {}
    commpq = {}
    commtotal = {}
    encoding = {}

    for s in DATASIZERANGE:
        datasize = s*DATASIZEUNIT
        scheme = makeScheme(datasize)
        commitment[s] = scheme.com_size / 8000000  # MB
        commpq[s] = scheme.comm_per_query() / 8000  # KB
        commtotal[s] = scheme.total_comm() / 8000000000  # GB
        encoding[s] = scheme.encoding_size() / 8000000000  # GB

    if not os.path.exists("./csvdata/"):
        os.makedirs("./csvdata")

    writeCSV("./csvdata/"+name+"_com.csv", commitment)
    writeCSV("./csvdata/"+name+"_comm_pq.csv", commpq)
    writeCSV("./csvdata/"+name+"_comm_total.csv", commtotal)
    writeCSV("./csvdata/"+name+"_encoding.csv", encoding)


############################################
writeScheme("rs", makeKZGScheme)
writeScheme("tensor", makeTensorScheme)
writeScheme("hash", makeHashBasedScheme)
writeScheme("homhash", makeHomHashBasedScheme)
writeScheme("fri", makeFRIUDRScheme)
