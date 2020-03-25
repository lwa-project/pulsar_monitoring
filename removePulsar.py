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

from runLK009 import _CATALOG_FILENAME, Pulsar


def main(args):
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
        
    # Find the entry to remove
    idx = None
    for i,bdy in enumerate(bdys):
        if bdy.name == args.name:
            idx = i
            break
    if idx is None:
        raise RuntimeError("Cannot find pulsar '%s' in the database" % args.name)
        
    # Prompt
    if args.force:
        yn = 'y'
    else:
        yn = raw_input('remove %s (entry #%i)? ' % (args.name, idx))
    if yn.lower() in ('y', 'yea', 'yes'):
        ## Remove
        del bdys[idx]
       
        ## Backup the catalog and write out the new version
        shutil.copy(_CATALOG_FILENAME, _CATALOG_FILENAME+'.old')
        fh = open(_CATALOG_FILENAME, 'w')
        fh.write('############################################\n')
        fh.write('#                                          #\n')
        fh.write('# Columns:                                 #\n')
        fh.write('#   1. Name                                #\n')
        fh.write('#   2. RA - HH:MM:SS.SS - J2000            #\n')
        fh.write('#   3. Declication - sDD:MM:SS.S - J2000   #\n')
        fh.write('#   4. Observation Duration - hours        #\n')
        fh.write('#   5. Observing Cadence - days            #\n')
        fh.write('#   6. Last MJD Observed                   #\n')
        fh.write('#                                          #\n')
        fh.write('# Updated:                                 #\n')
        fh.write("#   %s UTC                #\n" % datetime.utcnow().strftime('%Y/%m/%d %H:%M:%S'))
        fh.write('#                                          #\n')
        fh.write('############################################\n')
        for bdy in bdys:
            bdy.to_file(fh)
        fh.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Remove a pulsar from the LK009 observations')
    parser.add_argument('name', type=str, 
                        help='pulsar name')
    parser.add_argument('-f', '--force', action='store_true', 
                        help='do not prompt for confirmation')
    args = parser.parse_args()
    main(args)
