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
'''post-processing of the gridpacks: success tests
'''
def verifyGridpack (processoutputs, workingfolder, logfile):

    # verify that all jobs finished successfully
    finished = True
    logfile.write ('\n-->gridpack calculation summary\n')
    plotsFileName = os.path.dirname(logfile.name) + '/' + os.path.basename(logfile.name).split ('.')[0] + '.pdf'
    plotsFile = PdfPages (plotsFileName)
    summary_output = []

    # loop over processes
    for fil in processoutputs:

        if not wordInFile ('SIGMA', fil):
            print 'the following job had issues:\n  ' + fil
            logfile.write ('the following job had issues: ' + fil)
            finished = False

        thisoutput = extractIntegrationResults (fil)
        thisoutput = sorted (thisoutput, key = lambda x: (x[0], x[1], x[2]))

        process = fil.split ('/')[-2]
        NaNcheck =  checkForNaN (thisoutput, logfile, process)
        
        # print a summary of the grid calculation process

        if NaNcheck == False:
            logfile.write ('\nprocess: ' + process + ' \n')
            logfile.write ( 'step'.ljust (21) + '\tchannel\tit\tint             intE            intTOT          intTOTE         chi2\n')
            print 'process:', process
            for elem in thisoutput:
                writing = elem[0] + '\t' + str (elem[1]) + '\t' + str (int (elem[2]))
                for ind in range (3, len(elem)-1):
                   writing = writing + '\t{0:04.2e}'.format(elem[ind])
                writing = writing + '\t' + str (elem[-1])
                logfile.write (writing + '\n')
    
            # prepare drawings of the grid calculation process
    
            # find the unique occurrences in first and second column
            steps = set ([elem[0] for elem in thisoutput])
            # for each combination of those, prepare a plot of evolution with iteration
            for step in steps:
                if 'ALFA' in step : continue
                if 'NORMALIZATION' in step: continue
                chans = set ([elem[1] for elem in thisoutput if elem[0] == step])
                for chan in chans:
                    x_val = [elem[2] for elem in thisoutput if (elem[0] == step and elem[1] == chan)]
                    x_err = [0] * len (x_val)
                    # total integral
                    y_val_tota = [elem[5] for elem in thisoutput if (elem[0] == step and elem[1] == chan)]
                    y_err_tota = [elem[6] for elem in thisoutput if (elem[0] == step and elem[1] == chan)]
                    axes_limits_tota = getAxesLimits (x=x_val, y=y_val_tota, xerr=x_err, yerr=y_err_tota)
                    fig = plt.figure ()
                    plt.errorbar (x=x_val, y=y_val_tota, xerr=x_err, yerr=y_err_tota, label='total')
                    # partial integral
                    y_val_temp = [elem[3] for elem in thisoutput if (elem[0] == step and elem[1] == chan)]
                    y_err_temp = [elem[4] for elem in thisoutput if (elem[0] == step and elem[1] == chan)]
                    axes_limits_temp = getAxesLimits (x=x_val, y=y_val_temp, xerr=x_err, yerr=y_err_temp)
                    plt.errorbar (x=x_val, y=y_val_temp, xerr=x_err, yerr=y_err_temp, color = 'g', label = 'partial')
                    
                    chi2val = [elem[7] for elem in thisoutput if (elem[0] == step and elem[1] == chan and elem[2] == len (y_val_tota))]
                    summary_output.append ((process, step, chan, y_val_temp[-1], y_err_temp[-1], y_val_tota[-1], y_err_tota[-1], chi2val[0]))
                    
                    # drawing
                    leg = plt.legend (loc='best', fancybox=True)
                    leg.get_frame().set_alpha(0.5)
                    plt.title (process + ' ' + step + ' channel ' + str (chan) + '\n')
#                    plt.ticklabel_format (axis='y', style='sci', scilimits=(-2,2))
                    plt.xlabel('iteration')
                    plt.ylabel('total integral')
                    axes_limits = mergeAxesLimits (axes_limits_temp, axes_limits_tota)
                    plt.xlim ((axes_limits[0], axes_limits[1]))
                    fig.savefig (plotsFile, format = 'pdf')
        else:
            finished = False
            
    # end loop over processes

    logfile.write ('\n summary of the results at the last iteration\n\n')
    for elem in summary_output:
        writing = elem[0] + '\t' + str (elem[1]) + '\t' + str (int (elem[2]))
        for ind in range (3, len(elem)-1):
           writing = writing + '\t{0:04.2e}'.format(elem[ind])
        writing = writing + '\t' + str (elem[-1])
        logfile.write (writing + '\n')

    # summary plots
    steps = set ([elem[1] for elem in summary_output])
    for step in steps: 

        # the calculation result

        y_val_tota = [elem[5] for elem in summary_output if (elem[1] == step)]
        y_err_tota = [elem[6] for elem in summary_output if (elem[1] == step)]
        x_label = [elem[0]+' '+str(elem[2]) for elem in summary_output if (elem[1] == step)]
        x_val = range (len (x_label))
        x_err = [0] * len (x_val)
        fig = plt.figure ()
        plt.errorbar (x_val, y_val_tota, xerr=x_err, yerr=y_err_tota)
        plt.xticks (x_val, x_label, rotation='vertical')
        plt.title (step + ' results\n')
        plt.ylabel ('integral')
        plt.yscale ('log')
        plt.xlim((-1,len (x_label)))
        plt.subplots_adjust (bottom=0.5)
        fig.savefig (plotsFile, format = 'pdf')

        # the chisquare

        y_val_tota = [elem[7] for elem in summary_output if (elem[1] == step)]
        x_label = [elem[0]+' '+str(elem[2]) for elem in summary_output if (elem[1] == step)]
        x_val = range (len (x_label))
        fig2 = plt.figure ()
        plt.plot (x_val, y_val_tota, color ='r')
        plt.xticks (x_val, x_label, rotation='vertical')
        plt.title (step + ' chi2\n')
        plt.ylabel ('chi2')
        plt.xlim((-1,len (x_label)))
        plt.subplots_adjust (bottom=0.5)
        fig2.savefig (plotsFile, format = 'pdf')

    # plot the alphas (importance of each integration step)
    alfas = []
    for fil in processoutputs:
        process = fil.split ('/')[-2]
        outfile = open (fil, 'r')
        content = outfile.read ().split ('\n')
        for ind in range(len(content)-1, 0, -1):
            if 'iphs_ind=' in content[ind] and 'alfa=' in content[ind]:
                channel = content[ind].split()[1]
                alfa = content[ind].split()[3]
                alfas.append ([process,channel,alfa])
                break
        outfile.close ()
                                
    y_val = [float (elem[2]) for elem in alfas]
    x_label = [elem[0]+' '+elem[1] for elem in alfas]
    x_val = range (len (x_label))
    fig = plt.figure ()
    plt.plot (x_val, y_val)
    plt.xticks (x_val, x_label, rotation='vertical')
    plt.title ('ALFA VALUE')
    plt.ylabel ('value')
    plt.yscale ('log')
    plt.xlim((-1, len (x_label)))
    plt.ylim(( 0.8 * min (y_val), 2 * max (y_val) ))
    plt.subplots_adjust (bottom=0.5)
    fig.savefig (plotsFile, format = 'pdf')

    # plot the cross-section
    sigmas = []
    for fil in processoutputs:
        process = fil.split ('/')[-2]
        outfile = open (fil, 'r')
        content = outfile.read ().split ('\n')
        for ind in range(len(content)-1, 0, -1):
            if 'SIGMA' in content[ind]:
                sigma_val = content[ind].split()[2]
                sigma_err = content[ind].split()[4]
                sigmas.append ([process,sigma_val,sigma_err])
                break
        outfile.close ()
                                
    y_val = [float (elem[1].replace ('D','E')) for elem in sigmas]
    y_err = [float (elem[2].replace ('D','E')) for elem in sigmas]
    x_label = [elem[0] for elem in sigmas]
    x_val = range (len (x_label))
    x_err = [0] * len (x_val)
    fig = plt.figure ()
    plt.errorbar (x_val, y_val, xerr=x_err, yerr=y_err)
    plt.xticks (x_val, x_label, rotation='vertical')
    plt.title ('SIGMA\n')
    plt.ylabel ('cross-section (pb)')
    plt.yscale ('log')
    plt.xlim((-1,len (x_label)))
    plt.subplots_adjust (bottom=0.5)
    fig.savefig (plotsFile, format = 'pdf')

    # check that the uncertainty on the total XS is small
    # this check stays at the end, to run all the steps before anyhow
    
    summaryfile = open (workingfolder + '/result', 'r')
    totalCrossSection = []
    for line in summaryfile.read ().split ('\n'):
        if 'total cross section=' in line:
            words = line.split ()
            totalCrossSection.append (float (words[3]))
            totalCrossSection.append (float (words[5]))
            break
    print 'TOTAL XS = ' + str (totalCrossSection[0]) + ' +/- ' + str (totalCrossSection[1]) + ' pb'
    logfile.write ('TOTAL XS = ' + str (totalCrossSection[0]) + ' +/- ' + str (totalCrossSection[1]) + ' pb\n') 

    if totalCrossSection[1] / totalCrossSection[0] > 0.02 :
        finished = False
        print 'Uncertainty on the total XS is too large:', str (totalCrossSection[1] / totalCrossSection[0])
        logfile.write ('Uncertainty on the total XS is too large: ' + str (totalCrossSection[1] / totalCrossSection[0]) + '\n') 
    summaryfile.close ()

    plotsFile.close ()
    print 'graphic summary of the calculation written in : ', plotsFileName
    print 'logfile written in : ', logfile.name

    if finished == False :
        print '\n     GRIDPACK PREPARATION FAILED\n'
        logfile.write ('\n     GRIDPACK PREPARATION FAILED\n')
 
    # check the results of each step of the grid integration,
    # for each of the processes that have been generated
  

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