import sys
import os
import commands
from commands import getstatusoutput
import datetime
import argparse
import datetime
import math
import ConfigParser
import re
import time
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import math

config = ConfigParser.ConfigParser ()
config.optionxform = str # to preserve the case when reading options
print('reading config file:' + sys.argv[1])
config.read (sys.argv[1])
debugging = 1

################################################################
def getPhantom (config, workingfolder, debugging):#midified from original one, removed unused variables

    # get the precompiled phantom code, untar it and get the name of the pdf libraries
    # used in the compilation
    phantom = config.get ('general', 'package')
    # Check is the package is on http or in local folder
    if "http" in phantom:
        foundIT = execute ('wget ' + phantom, debugging)
        if not foundIT[0] == 0:
            print 'Phantom version: ' + phantom + ' not found, exiting'
            sys.exit (1)
    else:
        # should the phantom code be a local folder
        if not os.path.isfile (phantom):
            print 'Phantom package ' + phantom + ' not found, exiting'
            sys.exit (1)
        execute ('cp ' + phantom + ' ' + workingfolder, debugging)

    execute ('tar xzf ' + phantom.split ('/')[-1], debugging)
    phantomfolder = phantom.split ('/')[-1]
    dummy = '.tar.gz'
    phantomfolder = phantomfolder[0:-len(dummy)] if phantomfolder.endswith(dummy) else phantomfolder
    dummy = '.tgz'
    phantomfolder = phantomfolder[0:-len(dummy)] if phantomfolder.endswith(dummy) else phantomfolder
    
    return [phantomfolder, phantom]
################################################################
def execute (command, verbosity = False) :
    if (verbosity == True):
        print '---> running:'
        print command
    retCode = getstatusoutput (command)
    if (verbosity == True):
        for ri in retCode: print ri
    return retCode
################################################################
# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
def replaceParameterInFile (inputFile, outputFile, substitute):
    f = open (inputFile)
    lines = f.readlines ()
    f.close ()
    f = open (outputFile, 'w')
    for line in lines:
        if line.startswith ('*') or line.startswith (' ') or len (line) < 3 or line.startswith ('.') :
            f.write (line)
        else:
            words = line.split (' ')
            if words[0] in substitute.keys ():
                f.write (words[0] + ' ' + substitute[words[0]])
            else:
                f.write (line)
    f.close ()
# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
################################################################
# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
def addGridsToRin (filename, grids, debug = False):
    if debug : print 'preparing ' + filename + '\n'
    configfile = open (filename, 'read')
    lines = configfile.readlines ()
    configfile.close ()
    configfile = open (filename, 'write')
    for line in lines:
        if line.startswith ('nfiles'):
            configfile.write (line)
            break
        configfile.write (line)
    for gridfile in grids.split ():
        if debug : print 'adding ' + gridfile + '\n'
        configfile.write (gridfile + '\n')
    configfile.close ()
# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
################################################################

#if this works remove this part from the submit.phantom.py script
######variables definition
rootfolder = os.getcwd()
foldername = os.getcwd () + '/' + config.get ('general', 'foldername')

if not os.path.exists (foldername):
    print('gridpack check:' + foldername + 'DOES NOT exist, exiting')
    sys.exit (1)
os.chdir (foldername)
workingfolder = os.getcwd ()


res = getPhantom (config, workingfolder, debugging)
phantomfolder = res[0]
phantom = res[1]

submitfilename = workingfolder + '/'  + config.get ('submission','scheduler') + 'file'

processoutputs = []
submitfile = open (submitfilename, 'read')

if config.get ('submission','scheduler') == 'LSF':
    for line in submitfile.readlines () :
        if 'bsub' in line: processoutputs.append (line.split()[6])
elif config.get ('submission','scheduler') == 'CONDOR':
    for line in submitfile.readlines () :
        if 'condor_submit' in line: 
            path = line.split()[1]
            path = path[:-1] + '/'
            run_out = os.popen('ls '+ path + '*out*').read().split('/')[-1][:-1]
            processoutputs.append (path + run_out)
else: 
    print('Error: submission scheduler not recognized, allowed option are CONDOR/LSF')
    sys.exit(1)
    
submitfile.close ()


finished = True
unfinished = int (0)
for fil in processoutputs:
    if not os.path.exists (fil):
        finished = False
        unfinished += 1
if not finished:
    sys.stdout.write ('Job are not finished yet, EXIT: ' )
    sys.exit(1)

#log file of the generation parameters
logfilename = workingfolder + '/log_GRID.txt'
logfile = open (logfilename, 'write')

# calculate the cross-section
command = 'cd ' + workingfolder + '; grep SIGMA */run.*.out > res ; '
command += workingfolder + '/' + phantomfolder + '/tools/totint.exe > result '
execute (command, debugging)
result = execute ('tail -n 1 ' + workingfolder + '/result', debugging)
Xsection = result[1] + ' pb'

#need to add this function if we want it!
#verifyGridpack (processoutputs, workingfolder, logfile)

logfile.write ('CONFIG FILE\n\n')
for section in config.sections ():
    options = config.options (section)
    for option in options:
        logfile.write ('[' + section + '] ' + option + ' : ' + config.get (section, option) + '\n')
logfile.write ('\nSETUP COMMAND\n\n')
logfile.write (command + '\n')
logfile.write ('\nRESULTING CROSS-SECTION FROM GRID CREATION\n\n')
logfile.write (Xsection + '\n')
logfile.close ()

# NB removing the phantom to save space, it will have to be put back (reasonably in the same place)
#    when generating the events
if not debugging:
    execute ('rm -rf ' + phantomfolder + '*', debugging)
    execute ('rm -rf CMSSW*', debugging)
    execute ('rm -rf LSFJOB*', debugging)

# get all the grids, to be run in workingfolder
#    gridfiles = execute ('for fil in `find -L . -name "phavegas*.dat"` ; do echo `pwd`/$fil ; done', debugging)
gridfiles = execute ('for fil in `find -L . -name "phavegas*.dat"` ; do echo $fil ; done', debugging)

# prepare the r.in file for the event production, starting from the template one
replacement = {'ionesh':'1\n', 'nfiles': str(len (gridfiles[1].split ()))+'\n', 'nunwevts':'EVENTSNUM\n', 'idum':'-RANDOMSEED\n'}
replaceParameterInFile (workingfolder + '/r.in', workingfolder + '/r_GEN.in', replacement)
addGridsToRin (workingfolder + '/r_GEN.in', gridfiles[1], debugging)

execute ('mv r.in r_GRID.in', debugging)
# prepare the script for the event generation
execute ('cp ../' + sys.argv[1] + ' ./', debugging)
execute ('tar cJf ../' + foldername.split ('/')[-1] + '.tar.xz *', debugging)
print('gridpack ' + foldername.split ('/')[-1] + '.tar.xz created')

os.chdir (rootfolder)
sys.exit (0)