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
        
    # Prompts
    if not args.name:
        args.name = raw_input('name? ')
    if not args.ra:
        args.ra = raw_input('RA [HH:MM:SS.SS; J2000]? ')
    args.ra = ephem.hours(args.ra)
    if not args.dec:
        args.dec = raw_input('Dec [sDD:MM:SS.S; J2000]? ')
    args.dec = ephem.degrees(args.dec)
    if not args.duration:
        args.duration = raw_input('observation duration [hr]? ')
        args.duration = float(args.duration)
    if not args.cadence:
        args.cadence = raw_input('observing cadence [day]? ')
        args.cadence = int(args.cadence, 10)
    new_bdy = Pulsar.from_line("%s %s %s %s %s 0" % (args.name, args.ra, args.dec, args.duration, args.cadence))
    
    # Make sure it isn't already in there
    idx = None
    for i,bdy in enumerate(bdys):
        if bdy.name == new_bdy.name:
            idx = i
            break
    if idx is not None:
        raise RuntimeError("'%s' appears to already be in the database" % args.name)
        
    # One more prompt
    print("Name: %s" % new_bdy.name)
    print("  RA: %s" % new_bdy._ra)
    print("  Dec: %s" % new_bdy._dec)
    print("  Duration: %.3f hr" % (new_bdy.duration/3600.0,))
    print("  Cadence: every %i days" % new_bdy.cadence)
    if args.force:
        yn = 'y'
    else:
        yn = raw_input('add? ')
    if yn.lower() in ('y', 'yea', 'yes'):
        ## Insert
        bdys.append( new_bdy )
        bdys.sort(key=lambda x:x._ra)
        
        ## Write out the new version
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
        
        ## Reminder
        print("")
        print("NOTE: Please remember to make the following directory on the UCF:")
        print("      /data/network/recent_data/kstovall/LK009/%s" % args.name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Add a pulsar to the LK009 observations')
    parser.add_argument('-n', '--name', type=str, 
                        help='pulsar name')
    parser.add_argument('-r', '--ra', type=str, 
                        help='RA [HH:MM:SS.SS, J2000]')
    parser.add_argument('-d', '--dec', type=str, 
                        help='declination [sDD:MM:SS.S; J2000]')
    parser.add_argument('-l', '--duration', type=float, 
                        help='observation duration [hr]')
    parser.add_argument('-c', '--cadence', type=int, 
                        help='observing cadence [days]')
    parser.add_argument('-f', '--force', action='store_true', 
                        help='do not prompt for confirmation')
    args = parser.parse_args()
    main(args)
        
    
