
 ###########################################################################################
 #                                                                                         #
 # aeropt/cmd.py                                                                           #
 #                                                                                         #
 # author: Ramiro Checa-Garcia                                                             #
 # email:  ramiro.checa-garcia@ecmwf.int                                                   #
 #                                                                                         #
 # history:                                                                                #
 #    - Nov-2022 [Ramiro Checa-Garcia]      1st tested version                             #
 #                                                                                         #
 # info:                                                                                   #
 #        CLASSES                                                                          #
 #        * runinfo         : encapsulate information of current run                       #
 #                                                                                         #
 #        FUNCTIONS                                                                        #
 #        * fill_runinfo    : creates a runinfo object                                     #
 #        * cli             : uses the argparse library to create a cli application        #                
 ###########################################################################################


import os
import sys
import argparse
from datetime import datetime


class runinfo:
    def __init__(self, date, machine, gittag, runpath, version, user, debug, logfile):
        self.date    = date
        self.machine = machine
        self.gittag  = gittag
        self.runpath = runpath
        self.version = version
        self.user    = user
        self.debug   = debug
        self.logfile = logfile
        return

def fill_runinfo(codeversion, debug, logfile):
    today   = datetime.today().strftime('%Y-%m-%d')
    machine = os.uname()[-1]
    user    = os.getlogin()
    version = codeversion
    runpath = os.getcwd()
    gittag  = os.popen("git tag").read()

    return runinfo(today, machine, gittag, runpath, version, user, debug, logfile)


def cli():

    str_prog    = "\033[1m ecaeropt v1.0 \033[0m"
    str_descrip = "\033[95m ECMWF tool to calculate aerosol optical properties (ecaeropt) \033[0m"
    str_epilog  = "Contact: Ramiro Checa-Garcia <ramiro.checa-garcia at ecmwf.int>\n"
    str_setting = "Calculations based on a SETTING toml file (stored results defined in SETTING file)"
    str_config  = "Calculations based on a CONFIG toml file, results stored as OUTNCNAME."
    str_info    = "Information with the list of engines currently implemented"
    str_test    = "Perform tests based on the given setting file"
    str_debug   = 'It is stored debug information in a date-time.log file'

    debug   = None
    setting = None
    config  = None
    info    = None
    test    = None
    parser  = argparse.ArgumentParser(prog       = str_prog,
                                      description= str_descrip,
                                      epilog     = str_epilog)

    parser.add_argument('-d', '--debug',      dest   =debug, 
                        action ='store_true', default=False,
                        help   = str_debug)

    group   = parser.add_mutually_exclusive_group()

    group.add_argument('-s',  '--setting', dest=setting, action='store', help=str_setting,
                       metavar="SETTING",  type=str)

    group.add_argument('-c',  '--config',  dest=config,  action='store', help =str_config,
                       nargs=2,  type=str, metavar=("CONFIG", "OUTNCNAME"))

    group.add_argument('-i',  '--info',    dest = info,    action='store_true',
                       default=False, help=str_info)
    group.add_argument('-t',  '--test',    dest = test,    action='store',
                       default=None,  metavar="SETTING", help=str_test)


    if len(sys.argv) < 2:
        args = parser.parse_args(["-h"])
    else:
        args = parser.parse_args()

    return args
