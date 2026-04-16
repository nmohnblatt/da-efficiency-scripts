#!/usr/bin/env python

import math
import sys
import matplotlib.pyplot as plt

from schemes import *
from friudr import makeFRIUDRScheme
from frijohnson import makeFRIJohnsonScheme


def plot_commitment_vs_invrate(datasize_mb=32, invrates=None):
    """
    Plot the commitment size of FRI UDR and FRI Johnson as a function of
    the inverse rate.

    Parameters
    ----------
    datasize_mb : int
        Data size in megabytes.
    invrates : list[int] or None
        List of inverse rates to evaluate. Defaults to [2, 4, 8, 16, 32, 64].
    """
    if invrates is None:
        invrates = [2, 4, 8, 16, 32, 64]

    datasize = datasize_mb * 8_000_000  # convert MB to bits

    com_udr = []
    com_johnson = []

    for invrate in invrates:
        scheme_udr = makeFRIUDRScheme(datasize, invrate=invrate)
        com_udr.append(scheme_udr.com_size / 8_000)  # bits to KB

        scheme_johnson = makeFRIJohnsonScheme(datasize, invrate=invrate)
        com_johnson.append(scheme_johnson.com_size / 8_000)  # bits to KB

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(invrates, com_udr, marker='o', label='FRI (Unique Decoding)')
    ax.plot(invrates, com_johnson, marker='s', label='FRI (Johnson)')

    ax.set_xlabel('Inverse Rate')
    ax.set_ylabel('Commitment Size [KB]')
    ax.set_title(f'Commitment Size vs Inverse Rate (data = {datasize_mb} MB)')
    ax.legend()
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.set_xticks(range(0, max(invrates) + 1, 2))

    plt.tight_layout()
    plt.savefig('commitment_vs_invrate.png', dpi=150)
    plt.show()


if __name__ == '__main__':
    datasize_mb = int(sys.argv[1]) if len(sys.argv) > 1 else 32
    plot_commitment_vs_invrate(datasize_mb)
