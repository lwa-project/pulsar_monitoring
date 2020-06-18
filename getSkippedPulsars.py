#!/usr/bin/env python

"""
Script to help manage the LK009 pulsar database.
"""

from __future__ import print_function

import os
import sys
import time
import ephem
import numpy
import shutil
import argparse
from datetime import datetime

from lsl.common.mcs import datetime_to_mjdmpm as dt2mjd

from runLK009 import _CATALOG_FILENAME, Pulsar


def main(args):
    # Get the current MJD
    mjd, _ = dt2mjd(datetime.utcnow())
    
    # Load in the target list
    ## Read
    fh = open(_CATALOG_FILENAME, 'r')
    lines = fh.readlines()
    fh.close()
    ## Parse
    bdys = []
    for line in lines:
        line = line.strip().rstrip()
        if line[0] == '#':
            continue
        elif len(line) < 3:
            continue
        bdy = Pulsar.from_line(line)
        bdys.append(bdy)
        
    # Search for missed pulsars
    missed = []
    for bdy in bdys:
        if bdy.cadence < 0:
            continue 
        if mjd - bdy.last_mjd > 1.9*bdy.cadence:
            missed.append(bdy)
            
    # Sort by "level of egregiousness" and report
    missed.sort(key=lambda x: (mjd-x.last_mjd)/x.cadence)
    missed.reverse()
    for bdy in missed:
        print("%s was last observed %i days (%i cycles) ago on %i" % (bdy.name,
                                                                      mjd-bdy.last_mjd,
                                                                      (mjd-bdy.last_mjd)/bdy.cadence,
                                                                      bdy.last_mjd))
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='List pulsars that have not been observed for more than one cadence period')
    args = parser.parse_args()
    main(args)
