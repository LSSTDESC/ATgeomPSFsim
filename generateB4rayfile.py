####################################################################################################################
#Generate that particular RAY file used by BeamFour. The convergent beam is generated.
#In this version inclined incident beam is implemented.
#
#The output file is generated using implemented function based on python formatted string.
#
#    author : Sylvie Dagoret-Campagne
#    creation date : August 28th 2020
###################################################################################################################


import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as plt3d
import numpy as np
import os
import matplotlib as mpl
import pandas as pd
import itertools

from astropy import units as u
from astropy.coordinates import Angle

import time
from datetime import datetime,date
import dateutil.parser
import pytz

import argparse

import logging
import coloredlogs
import configparser


if __name__ == "__main__":


    # date
    today = date.today()
    string_date = today.strftime("%Y-%m-%d")


    # time
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")

    tim = time.localtime()
    current_time = time.strftime("%H:%M:%S", tim)



    #timezones
    tz_LA = pytz.timezone('America/Los_Angeles')
    datetime_LA = datetime.now(tz_LA)
    print("LA time:", datetime_LA.strftime("%H:%M:%S"))


    tz_NY = pytz.timezone('America/New_York')
    datetime_NY = datetime.now(tz_NY)
    print("NY time:", datetime_NY.strftime("%H:%M:%S"))

    tz_London = pytz.timezone('Europe/London')
    datetime_London = datetime.now(tz_London)
    print("London time:", datetime_London.strftime("%H:%M:%S"))

    tz_Paris = pytz.timezone('Europe/Paris')
    datetime_Paris = datetime.now(tz_Paris)
    print("Paris time:", datetime_Paris.strftime("%H:%M:%S"))


    # start with logs
    #-----------------
    logging.basicConfig()
    logging.root.setLevel(logging.NOTSET)

    handle = __name__

    logger = logging.getLogger(handle)
    # logging.getLogger().setLevel(logging.INFO)
    logger.setLevel(logging.DEBUG)

    # If you don't want to see log messages from libraries, you can pass a
    # specific logger object to the install() function. In this case only log
    # messages originating from that logger will show up on the terminal.
    coloredlogs.install(level='DEBUG', logger=logger)
    coloredlogs.install(fmt='%(asctime)s,%(msecs)03d %(hostname)s %(name)s[%(process)d] %(levelname)s %(message)s')




    # config file
    #--------------
    configfile=""

    parser = argparse.ArgumentParser()
    parser.add_argument("--config",action="store", dest="configfile",help=f" run generate -config configfilename, with by ex configfilename = default.ini")

    results_args = parser.parse_args()


    print("results_args=",results_args)
    print("configfile=",results_args.configfile)



    msg = f"Start {parser.prog} at date : {string_date} and time :{current_time}"
    logger.info(msg)


