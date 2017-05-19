#!/usr/bin/env python
'''  '''

#####################
# IMPORT OPERATIONS #
#####################

from csv import DictReader
import sys
import os.path

###############
# AUTHOR INFO #
###############

__author__ = 'Michael Gruenstaeudl <m.gruenstaeudl@fu-berlin.de>'
__copyright__ = 'Copyright (C) 2017 Michael Gruenstaeudl'
__info__ = 'DictReplacer'
__version__ = '2017.02.13.2000'

#############
# DEBUGGING #
#############

#import pdb
#pdb.set_trace()


########
# MAIN #
########


# Parse infile names
try:
    inFile = sys.argv[1]
    os.path.isfile(inFile)
except IndexError:
    raise Exception("Could not read inFile.")

try:
    inDict = sys.argv[2]
    os.path.isfile(inDict)
except IndexError:
    
    raise Exception("Could not read inDict.")


# Test if inDict has appropriate title line
raw_inDict = open(inDict, 'rU').readlines()
if 'Target' not in raw_inDict[0] or 'Solution' not in raw_inDict[0]:
    sys.exit("ERROR: Title line does not contain one or both of the keywords `Target`and `Solution`.")


# Open and read in- and outfiles
inf = open(inFile, 'rU')
otf = open(inFile + '__replaced.txt', 'w')

try:
    reader = DictReader(open(inDict, 'rb'), delimiter=',', quotechar='"', skipinitialspace=True)
    repl_dicts = list(reader)
except:
    raise Exception("Could not parse inDict.")


# Do batch replace
targets, solutions = [], []
for dct in repl_dicts:
    if dct['Target'] in targets:
        sys.exit("ERROR: Target `%s` has already been used." %(dct['Target']))
    else:
        targets.append(dct['Target'])
    if dct['Solution'] in solutions:
        sys.exit("ERROR: Solution `%s` has already been used." %(dct['Solution']))
    else:
        solutions.append(dct['Solution'])
if len(targets) != len(solutions):
    sys.exit("ERROR: Number of targets unequal to number of solutions.")


# Do batch replace
hndl = inf.read()
for dct in repl_dicts:
    try:
        hndl = hndl.replace(dct['Target'], dct['Solution'])
    except:
        raise Exception("Could not perform replacement of target `%s` with solution `%s`." %(dct['Target'], dct['Solution']))


# Write and close all files
otf.write(hndl)
otf.close()
