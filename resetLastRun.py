#!/usr/bin/env python3

"""
Script to help manage the LK009 pulsar database.
"""

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
        
    # Load in the backup list, if it exists
    backup_filename = _CATALOG_FILENAME+'.old'
    old_bdys = []
    if os.path.exists(backup_filename):
        ## Read
        fh = open(backup_filename, 'r')
        lines = fh.readlines()
        fh.close()
        ## Parse
        for line in lines:
            line = line.strip().rstrip()
            if line[0] == '#':
                continue
            elif len(line) < 3:
                continue
            bdy = Pulsar.from_line(line)
            old_bdys.append(bdy)
            
    # Find the entry to reset in both catalogs
    idx = None
    for i,bdy in enumerate(bdys):
        if bdy.name == args.name:
            idx = i
            break
    if idx is None:
        raise RuntimeError("Cannot find pulsar '%s' in the database" % args.name)
    old_idx = None
    for i,bdy in enumerate(old_bdys):
        if bdy.name == args.name:
            old_idx = i
            break
            
    # Prompt
    old_mjd = bdys[idx].last_mjd
    if args.mjd is None and old_idx is not None:
        if old_bdys[old_idx].last_mjd < bdys[idx].last_mjd:
            new_mjd = str(old_bdys[old_idx].last_mjd)
        else:
            new_mjd = raw_input('new MJD for last run date for %s? ' % (args.name))
    else:
        new_mjd = raw_input('new MJD for last run date for %s? ' % (args.name))
    new_mjd = int(new_mjd, 10)
    
    yn = raw_input('change last run date for %s from %i to %i? ' % (args.name, old_mjd, new_mjd))
    if yn.lower() in ('y', 'yea', 'yes'):
        # Update
        bdys[idx].last_mjd = new_mjd
        
        # Write out the new version
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
    parser.add_argument('-m', '--mjd', type=int, 
                        help='MJD to reset the last run date to')
    args = parser.parse_args()
    main(args)
