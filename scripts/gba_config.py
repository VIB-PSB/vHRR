
'''
Setup user to run GBA framework
For now: only sets the go.obo file

usage: python3 gba_config.py -obo <go.obo>
'''


import argparse
import configparser

from sys import stderr
from pathlib import Path
from goatools.obo_parser import GODag


#command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("-obo", dest="go_obo", help="set location of go.obo file")
args = parser.parse_args()


# get version info
obodag_version = GODag(args.go_obo).version.split('rel(')[1].split(')')[0]


# init
config = configparser.ConfigParser()

# set values
config['GO'] = {'go_obo': args.go_obo}

# write file
with open('{}/.gbaconfig'.format(Path.home()), 'w') as outf:
    
    outf.write(
"""# GBA config file
# ---------------
#
# go.obo : 
#    - file: {}
#    - date: {}

""".format(args.go_obo, obodag_version)
              )

    config.write(outf)

stderr.write('[INFO] GBA config file writen to ~/.gbaconfig\n')