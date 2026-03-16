#!/usr/bin/env python

import math
import sys
from tabulate import tabulate

from schemes import *
from friudr import *


def makeRow(name, scheme, tex):
    comsize = '{:.2f}'.format(round(scheme.com_size/8000.0, 2))
    encodingsize = '{:.2f}'.format(
        round(scheme.encoding_size() / 8000000.0, 2))
    commpqsize = '{:.2f}'.format(round(scheme.comm_per_query() / 8000.0, 2))
    reception = scheme.reception()
    encodinglength = scheme.encoding_length()
    samples = scheme.samples()
    commsize = '{:.2f}'.format(round(scheme.total_comm() / 8000000.0, 2))
    if tex:
        row = ["\Inst"+name, comsize, encodingsize, commpqsize, commsize]
    else:
        row = [name, comsize, encodingsize, commpqsize,
               (reception, encodinglength), samples, commsize]
    return row

#####################################################################


opts = [opt for opt in sys.argv[1:] if opt.startswith("-")]
args = [arg for arg in sys.argv[1:] if not arg.startswith("-")]

if len(args) == 0:
    print("Missing Argument: Datasize in Megabytes.")
    print("Hint: To print the table in LaTeX code, add the option -l.")
    sys.exit(-1)

datasize = int(args[0])*8000000


# Print to LaTeX
tex = "-l" in opts

if tex:
    table = [["Name", "|com|", "|Encoding|", "Comm. p. Q.", "Comm Total"]]
else:
    table = [["Name", "|com| [KB]", "|Encoding| [MB]",
              "Comm. p. Q. [KB]", "Reception", "Samples", "Comm Total [MB]"]]


scheme = makeNaiveScheme(datasize)
table.append(makeRow("Naive", scheme, tex))

scheme = makeMerkleScheme(datasize)
table.append(makeRow("Merkle", scheme, tex))

scheme = makeKZGScheme(datasize)
table.append(makeRow("RS", scheme, tex))

scheme = makeTensorScheme(datasize)
table.append(makeRow("Tensor", scheme, tex))

scheme = makeHashBasedScheme(datasize)
table.append(makeRow("Hash", scheme, tex))

scheme = makeHomHashBasedScheme(datasize)
table.append(makeRow("HomHash", scheme, tex))

scheme = makeFRIUDRScheme(datasize)
table.append(makeRow("FRI (UDR)", scheme, tex))


if tex:
    print(tabulate(table, headers='firstrow',
          tablefmt='latex_raw', disable_numparse=True))
else:
    print(tabulate(table, headers='firstrow', tablefmt='fancy_grid'))
