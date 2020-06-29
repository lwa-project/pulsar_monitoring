#!/usr/bin/env python3

"""
Script to help schedule LK009 observations during idle windows.
"""

import os
import sys
import time
import ephem
import numpy
import shutil
import argparse
import subprocess
from datetime import datetime, timedelta
import pytz

from lsl import astro
from lsl.common import busy
from lsl.common.stations import lwa1
from lsl.common import sdf as lslsdf
from lsl.common.mcs import mjdmpm_to_datetime as mjd2dt, datetime_to_mjdmpm as dt2mjd

from lwa_mcs.tp import schedule_sdfs
from lwa_mcs.utils import schedule_at_command


_CATALOG_FILENAME = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'LK009_Pulsars.txt')


UTC = pytz.utc


class Pulsar(ephem.FixedBody):
    """
    Wrapper around the ephem.FixedBody class to allow us to add in custom attributes
    and a few helper methods for determining when things can be observed.
    """
    
    obs = lwa1.get_observer()    # defaults to LWA1
    
    _padding = 10    # total session padding time in seconds
        
    @classmethod
    def from_line(cls, line):
        """
        Return a new Pulsar instance generated from a line in the catalog file.
        """
        
        name, ra, dec, duration, cadence, mjd = line.split(None, 5)
        bdy = cls()
        bdy.name = name
        bdy._ra = ra
        bdy._dec = dec
        bdy._epoch = ephem.J2000
        bdy.duration = float(duration)*3600.0
        bdy.cadence = int(cadence, 10)
        bdy.last_mjd = int(mjd, 10)
        return bdy
        
    def to_line(self):
        """
        Return a text string that is the same format as the catalog file.
        """
        
        return "%-10s  %-11s  %-11s  %-3.1f  %-2i  %-5i" % (self.name, 
                                                            self._ra, 
                                                            self._dec if self._dec < 0 else '+'+str(self._dec),
                                                            self.duration/3600.,
                                                            self.cadence, 
                                                            self.last_mjd)
        
    def to_file(self, fh):
        """
        Similar to to_line(), but writes to an open file handle.
        """
        
        if not isinstance(fh, file):
            raise TypeError("Expected an open filehandle")
        fh.write("%s\n" % self.to_line())
        
    def within_beam(self, other, width_deg=1.5):
        """
        Given an ephem.FixedBody or another Pulsar instance, determin if the other
        target is within the beam of the first target.
        """
        
        if not isinstance(other, (ephem.Body, ephem.FixedBody, Pulsar)):
            raise TypeError("Expected a ephem.Body, ephem.FixedBody, or Pulsar instance")
                
        width = ephem.degrees(width_deg*numpy.pi/180)
        try:
            sep = ephem.separation((self._ra,self._dec), (other.ra,other.dec))
        except (RuntimeError, AttributeError):
            sep = ephem.separation((self._ra,self._dec), (other._ra,other._dec))
        if sep <= width:
            return True
        return False
        
    def set_observer(self, obs):
        """
        Update the ephem.Observer used for this object.
        """
        
        if not isinstance(obs, ephem.Observer):
            raise TypeError("Expected an ephem.Observer instance")
        
        self.obs = obs
        
    def get_start_stop(self, start, stop, padding=True):
        """
        Given a start datetime instance and a stop datetime instance, determine
        when the target should be observed based on its transit time and 
        required observation duration.
        """
        
        if not isinstance(start, datetime):
            raise TypeError("Expected start to be a datetime instance")
        if not isinstance(stop, datetime):
            raise TypeError("Expected stop to be a datetime instance")
        
        self.obs.date = start.strftime('%Y/%m/%d %H:%M:%S')
        self.compute(self.obs)
        
        # Get the transit time
        bdy_transit = self.obs.next_transit(self)
        bdy_transit = mjd2dt(bdy_transit + (astro.DJD_OFFSET - astro.MJD_OFFSET), 0)
        # Round to the nearest second
        if bdy_transit.microsecond >= 5000000:
            bdy_transit += timedelta(seconds=1)
        bdy_transit = bdy_transit.replace(microsecond=0)
        # Compute the start and stop times to center on transit
        bdy_start = bdy_transit - timedelta(seconds=self.duration/2.0)
        bdy_stop  = bdy_transit + timedelta(seconds=self.duration/2.0)
        if padding:
            # Add in the session padding, if needed
            bdy_start -= timedelta(seconds=self._padding/2.0)
            bdy_stop  += timedelta(seconds=self._padding/2.0)
        return bdy_start, bdy_stop
        
    def can_run(self, start, stop):
        """
        Given a start datetime instance and a stop datetime instance, determine
        if the target can be observed within that window.
        """
        try:
            bdy_start, bdy_stop = self.get_start_stop(start, stop, padding=True)
            if bdy_start >= start and bdy_stop <= stop:
                return True
        except ephem.NeverUp:
            pass
        return False
            
    def should_run(self, start, stop):
        """
        Given a start datetime instance and a stop datetime instance, determine
        if the target should be observed given the last time it was observed.
        """
        
        obs_mjd_start, _ = dt2mjd(start)
        if obs_mjd_start >= (self.last_mjd + self.cadence) and self.cadence > 0:
            return True
        return False


_SPACE_CONVERSION_RATE = 19.6e6 / 4096 * 4128 * 2 * 2 * 2    # B/s for two pols, two tunings, and two beams


def get_available_space(user, buffer_factor=0.8, min_free_tb=2.0):
    """
    Given a UCF username, calculate and return how many seconds of dual-beam recording 
    can be stored given the current level of disk usage in /data/network.
    """
    
    df = subprocess.Popen(['ssh', 'mcsdr@lwaucf1', "df -BG /data/network/recent_data/%s" % user], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    space, err = df.communicate()
    try:
        space = space.split('\n')[-2]
        space = space.split(None, 5)[3][:-1]
    except IndexError:
        space = "0"
    space = int(space, 10)*1024**3
    space = min([space*buffer_factor, space-min_free_tb*1024**4])
    space = max([0, space])
    return space / _SPACE_CONVERSION_RATE


def main(args):
    # Get the start and stop times for the window that we are scheduling
    start = datetime.strptime('%s %s' % (args.start_date, args.start_time), '%Y/%m/%d %H:%M:%S')
    stop  = datetime.strptime('%s %s' % (args.stop_date, args.stop_time), '%Y/%m/%d %H:%M:%S')
    print("Scheduling LK009 for %s to %s" % (start.strftime('%Y/%m/%d %H:%M:%S'),  
                                             stop.strftime('%Y/%m/%d %H:%M:%S')))
    print("  Window is %.3f hr long" % ((stop-start).total_seconds()/3600.0,))
    
    # Get the observer and convert to LST
    obs = lwa1.get_observer()
    obs.date = start.strftime('%Y/%m/%d %H:%M:%S')
    lst_start = obs.sidereal_time()
    obs.date = stop.strftime('%Y/%m/%d %H:%M:%S')
    lst_stop = obs.sidereal_time()
    print("  Corresponds to the LST range of %s to %s" % (lst_start, lst_stop))
    
    # Get how much recording we can actually do given the available space
    rec_secs_possible = get_available_space('pulsar')
    print('Recording space available is %.3f hr' % (rec_secs_possible/3600.0,))
    
    # Adjust the start/stop times to allow for an INI
    orig_start, orig_stop = start, stop
    start = start + timedelta(minutes=25)
    
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
    print('Loaded %i targets' % len(bdys))
    
    ## Sort the pulsars by hour angle
    #ha_start = lst_start - max([bdy.duration for bdy in bdys])/86400.0*2*numpy.pi
    #if ha_start < 0:
    #    ha_start += 2*numpy.pi
    #bdys.sort(key=lambda x: ((x._ra - ha_start) if x._ra - ha_start >= 0 else (x._ra - ha_start + 2*numpy.pi)))
    
    # Sort the pulsars by "observability rank"
    mjd_start, _ = dt2mjd(start)
    bdys.sort(key=lambda x: x.last_mjd - mjd_start)
    
    # Come up a list of objects to observe
    bdys_run = []
    rec_secs_allocated = 0.0
    for bdy in bdys:
        ## Is there any space left at this point?
        if bdy.duration >= rec_secs_possible:
            continue
            
        ## Can/should this target be running?
        if not bdy.should_run(start, stop):
            continue
        if not bdy.can_run(start, stop):
            continue
            
        ## Is there a conflict with something that we have previously added to the list?
        css = bdy.get_start_stop(start, stop)
        cta = dt2mjd(start)[0] - (bdy.last_mjd + bdy.cadence)
        conflict = False
        for prv in bdys_run:
            pss = prv.get_start_stop(start, stop)
            pta = dt2mjd(start)[0] - (prv.last_mjd + prv.cadence)
            if (css[0] >= pss[0] and css[0] <= pss[1]) or (css[1] >= pss[0] and css[1] <= pss[1]) \
            or (pss[0] >= css[0] and pss[0] <= css[1]) or (pss[1] >= css[0] and pss[1] <= css[1]):
                conflict = True
                break
        if conflict:
            continue
            
        ## Add it to the list
        bdy.final = bdy.get_start_stop(start, stop, padding=False)
        bdys_run.append(bdy)
        rec_secs_possible -= bdy.duration
        rec_secs_allocated += bdy.duration
    bdys_run.sort(key=lambda x: x.final[0])
    
    # Assign each observation set to beams
    load = {2:0.0, 3:0.0, 4:0.0}
    for bdy in bdys_run:
        beams = [lbp[1] for lbp in sorted([(load[b],b) for b in load])]
        bdy.beams = sorted((beams[0], beams[1]))
        load[beams[0]] += bdy.duration
        load[beams[1]] += bdy.duration
        
    # Report
    print("Identified %i targets (%.3f hr) to observe during this window:" % (len(bdys_run), rec_secs_allocated/3600.0))
    for bdy in bdys_run:
        print("  %s for %.3f hr on beams %i and %i" % (bdy.name, bdy.duration/3600.0, 
                                                       bdy.beams[0], bdy.beams[1]))
        
    # Identify free periods for maintenance scheduling
    ## Times
    freeTime = []
    block_start = start
    block_stop = stop
    tFree = block_start + timedelta(seconds=0)
    while tFree < block_stop:
        isFree = True
        for bdy in bdys_run:
            if tFree < bdy.final[0] - timedelta(minutes=20):
                pass
            elif tFree > bdy.final[1] + timedelta(minutes=20):
                pass
            else:
                isFree = False
                break
                
        if isFree:
            freeTime.append(tFree)
        tFree = tFree + timedelta(minutes=2)
    ## Windows
    freeWindows = []
    if len(freeTime) > 0:
        freeWindows.append( [freeTime[0], freeTime[0]] )
        for tFree in freeTime:
            if tFree-freeWindows[-1][1] <= timedelta(seconds=120):
                freeWindows[-1][1] = tFree
            else:
                freeWindows.append( [tFree, tFree] )
        freeWindows = [freeWindow for freeWindow in freeWindows if freeWindow[1]-freeWindow[0] > timedelta()]
    print("Free Times:")
    for f0,f1 in freeWindows:
        print(" %s to %s, length %s" % (f0.strftime("%m/%d %H:%M:%S"), f1.strftime("%m/%d %H:%M:%S"), f1-f0))
        
    # Indentify busy times
    busyWindows = []
    if len(bdys_run) > 0:
        busyWindows.append( [bdys_run[0].final[0], bdys_run[0].final[1]] )
        for bdy in bdys_run:
            tSDFStart, tSDFStop = bdy.final
            if tSDFStart-busyWindows[-1][1] <= timedelta(minutes=45):
                busyWindows[-1][1] = max([busyWindows[-1][1], tSDFStop])
            else:
                busyWindows.append( [tSDFStart, tSDFStop] )
    print("BusyTimes:")
    for f0,f1 in busyWindows:
        print(" %s to %s, length %s" % (f0.strftime("%m/%d %H:%M:%S"), f1.strftime("%m/%d %H:%M:%S"), f1-f0))
        
    # Load in the next session ID
    try:
        fh = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'state'), 'r')
        session_id = int(fh.read(), 10)
        fh.close()
    except (IOError, ValueError):
        session_id = 1
        
    # Make sure we have a place to put the SDFs
    if not args.dry_run:
        sdf_dir = "/home/op1/MCS/tp/%s/" % start.strftime("%y%m%d")
        if not os.path.exists(sdf_dir):
            print("Creating date directory: %s" % sdf_dir)
            os.mkdir(sdf_dir)
            
    # Build the observations
    fileids, filenames = [], []
    lslobs = lslsdf.Observer('Pratik Kumar', 82)
    for bdy in bdys_run:
        bdy_start, bdy_stop = bdy.final
        bdy_start = UTC.localize(bdy_start)
        bdy_stop = UTC.localize(bdy_stop)
        
        for b,beam in enumerate(bdy.beams):
            targ = lslsdf.DRX(bdy.name, bdy.name, bdy_start, bdy_stop-bdy_start, 
                              bdy._ra, bdy._dec, 
                              35.1e6 if b == 0 else 64.5e6, 
                              49.8e6 if b == 0 else 79.2e6, 
                              7, max_snr=False)                  
            sess = lslsdf.Session('%s, beam %i' % (bdy.name, beam), session_id, [targ,])
            sess.drx_beam = beam
            sess.data_return_method = 'UCF'
            sess.ucf_username = "pulsar/LK009/%s" % bdy.name
            proj = lslsdf.Project(lslobs, 'Continued Regular Monitoring of Pulsars with LWA1', 'LK009', [sess,])
            sdf = proj.render(verbose=False)
            
            if not args.dry_run:
                filename = '%s_%s_%s_%04d_B%i.sdf' % (proj.id, bdy_start.strftime("%y%m%d"), 
                                                      bdy_start.strftime("%H%M"), 
                                                      sess.id, sess.drx_beam)
                filename = os.path.join(sdf_dir, filename)
                
                fh = open(filename, 'w')
                fh.write(sdf)
                fh.close()
                
                fileids.append(session_id)
                filenames.append(filename)
            session_id = session_id + 1
            
    # Submit the SDFs
    print("Submitting SDFs for scheduling")
    if not args.dry_run and filenames:
        bi = busy.BusyIndicator(message="'waiting'")
        bi.start()
        success = schedule_sdfs(filenames)
        bi.stop()
        if not success:
            print("There seems to be an issue with scheduling the SDFs, giving up!")
            print("Script is aborted!")
            sys.exit(1)
            
    print("SDFs successfully scheduled")
    if not args.dry_run:
        # Write out new session id
        fh = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'state'), 'w')
        fh.write("%i\n" % session_id)
        fh.close()
        
        # Update list to get it ready to write back out.  In the process, deal with any 
        # opportunistic targets
        for bdy in bdys_run:
            bdy.last_mjd = dt2mjd(bdy.get_start_stop(start, stop, obs)[0])[0]
            for opt in bdys:
                if opt.cadence <= 0 and bdy.within_beam(opt):
                    print("  Note: observation of %s also contains %s" % (bdy.name, opt.name))
                    opt.last_mjd = bdy.last_mjd
                    
        # Backup the catalog and write out the new version
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
        
    print("Scheduling other commands:")
    atCommands = []
    # Build up the commands
    ## INIdp.sh calls at the start and end of each busy window
    for i,busyWindow in enumerate(busyWindows):
        tINI = busyWindow[0] - timedelta(minutes=20)
        tINI = tINI.replace(second=0, microsecond=0)
        atCommands.append( (tINI, '/home/op1/MCS/sch/INIdp.sh') )
    ## Maintenance commands for the station
    tTBWLast = orig_start - timedelta(hours=12)
    tDRSULast = orig_start - timedelta(hours=12)
    tLWADBLast = orig_start - timedelta(hours=12)
    for i,freeWindow in enumerate(freeWindows):
        length = freeWindow[1] - freeWindow[0]
        
        if length >= timedelta(minutes=45):
            ### Default TBN frequency, TBW health checks, and DRSU scans
            if i == 0 or abs(freeWindow[0] - start) >= timedelta(minutes=3):
                tTBN = freeWindow[0]
                atCommands.append( (tTBN, '/home/op1/MCS/sch/startTBN_split.sh') )
                
            tTBW = freeWindow[0] + timedelta(minutes=4)
            tTBW = tTBW.replace(second=0, microsecond=0)
            tDRSU = freeWindow[0] + timedelta(minutes=4)
            tDRSU = tDRSU.replace(second=0, microsecond=0)
            tLWADB = freeWindow[0] + timedelta(minutes=4)
            tLWADB = tLWADB.replace(second=0, microsecond=0)
            while tTBW <= (freeWindow[1] - timedelta(minutes=6)):
                if tTBW > tTBWLast + timedelta(hours=4):
                    atCommands.append( (tTBW, '/home/op1/MCS/exec/acquireHealthCheckAndProcess.py') )
                    tTBWLast = tTBW
                    
                elif tDRSU > tDRSULast + timedelta(hours=6):
                    atCommands.append( (tDRSU, '/home/op1/MCS/sch/operatorScripts/selectBestDRSU.py --all') )
                    tDRSU = tDRSU + timedelta(minutes=2)
                    atCommands.append( (tDRSU, '/home/op1/MCS/sch/operatorScripts/postDRSUStatus.py') )
                    tDRSULast = tDRSU
                    
                elif tLWADB > tLWADBLast + timedelta(hours=4):
                    if tLWADB.hour >= 8:
                        atCommands.append( (tLWADB, '/home/op1/MCS/sch/operatorScripts/scanDRSUs.sh') )
                        tLWADBLast = tLWADB + timedelta(hours=24)
                        
                tTBW = tTBW + timedelta(minutes=15)
                tDRSU = tDRSU + timedelta(minutes=15)
                tLWADB = tLWADB + timedelta(minutes=15)
                
        elif length >= timedelta(minutes=15):
            ### TBW health checks and DRSU scans
            tTBW = freeWindow[0] + timedelta(minutes=4)
            tTBW = tTBW.replace(second=0, microsecond=0)
            tDRSU = freeWindow[0] + timedelta(minutes=4)
            tDRSU = tDRSU.replace(second=0, microsecond=0)
            if tTBW > tTBWLast + timedelta(hours=4):
                atCommands.append( (tTBW, '/home/op1/MCS/exec/acquireHealthCheckAndProcess.py') )
                tTBWLast = tTBW
                
            elif tDRSU > tDRSULast + timedelta(hours=6):
                atCommands.append( (tDRSU, '/home/op1/MCS/sch/operatorScripts/selectBestDRSU.py --all') )
                tDRSU = tDRSU + timedelta(minutes=2)
                atCommands.append( (tDRSU, '/home/op1/MCS/sch/operatorScripts/postDRSUStatus.py') )
                tDRSULast = tDRSU
                
        elif length >= timedelta(minutes=6):
            ### DRSU scans
            tDRSU = freeWindow[0] + timedelta(minutes=2)
            tDRSU = tDRSU.replace(second=0, microsecond=0)
            if tDRSU > tDRSULast + timedelta(hours=6):
                atCommands.append( (tDRSU, '/home/op1/MCS/sch/operatorScripts/selectBestDRSU.py --all') )
                tDRSU = tDRSU + timedelta(minutes=2)
                atCommands.append( (tDRSU, '/home/op1/MCS/sch/operatorScripts/postDRSUStatus.py') )
                tDRSULast = tDRSU
                
        else:
            ### TBN start if it is the end of the session
            if i == len(freeWindows)-1:
                if freeWindow[0] > max([bdy.final[1] for bdy in bdys_run]):
                    tTBN = freeWindow[0]
                    atCommands.append( (tTBN, '/home/op1/MCS/sch/startTBN_split.sh') )
                    
    ## Implement the commands
    atIDs = []
    for cmd in atCommands:
        if not args.dry_run:
            atID = schedule_at_command(*cmd)
        else:
            atID = -1
        atIDs.append(atID)
        
    print("Done, saving log")
    if not args.dry_run:
        rpt = []
        for cmd,id in zip(atCommands,atIDs):
            if id != -1:
                id = " (#%i)" % id
            else:
                id = ''
            rpt.append( [cmd[0], "%s%s" % (cmd[1], id)] )
        for bdy in bdys_run:
            bdy_start, bdy_stop = bdy.final
            rpt.append( [bdy_start, "%s starts on beams %i and %i" % (bdy.name, bdy.beams[0], bdy.beams[1])] )
            rpt.append( [bdy_stop, "%s stops on beams %i and %i" % (bdy.name, bdy.beams[0], bdy.beams[1])] )
        rpt.sort()
    
        fh = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'runtime.log'), 'a')
        fh.write("Completed Scheduling for UTC %s to %s\n" % (orig_start.strftime('%Y/%m/%d %H:%M:%S'), 
                                                              orig_stop.strftime('%Y/%m/%d %H:%M:%S')))
        fh.write("  Timeline:\n")
        for t,info in rpt:
            fh.write("    %s - %s\n" % (t.strftime("%Y/%m/%d %H:%M:%S"), info))
        fh.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Schedule LK009 observations')
    parser.add_argument('start_date', type=str, 
                        help='scheduling window UTC start date in YYYY/MM/DD format')
    parser.add_argument('start_time', type=str,
                        help='scheduling window UTC start time in HH:MM:SS format')
    parser.add_argument('stop_date', type=str, 
                        help='scheduling window UTC stop date in YYYY/MM/DD format')
    parser.add_argument('stop_time', type=str,
                        help='scheduling window UTC stop time in HH:MM:SS format')
    parser.add_argument('-n', '--dry-run', action='store_true', 
                        help='perform a dry-run only')
    args = parser.parse_args()
    main(args)
    
