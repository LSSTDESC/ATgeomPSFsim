####################################################################################################################
#Generate that particular RAY file used by BeamFour. The convergent beam is generated.
#In this version inclined incident beam is implemented.
#
#The output file is generated using implemented function based on python formatted string.
#
#    author : Sylvie Dagoret-Campagne
#    creation date : August 28th 2020
#    update        : August 29th 2020
#
# run :
#       python generateB4rayfile.py --config default.ini
#
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

#------------------------------------------------------------------------------------------------------
def GetRayFileString(theta_x,theta_y):
    """
    GetRayFileString(theta_x,theta_y) : convert the angles theta_x, theta_y in armin into a string

    theta_x
    theta_y
    :return: theta_x_str,theta_y_str : output string
    """

    theta_x_num = int(theta_x * 100)
    theta_y_num = int(theta_y * 100)

    if theta_x_num >= 0:
        theta_nstr = '{:0>2}'.format(theta_x_num)
        theta_nstr = theta_nstr.zfill(4)
        theta_x_str = "p" + theta_nstr
    else:
        theta_nstr = '{:0>2}'.format(-theta_x_num)
        theta_nstr = theta_nstr.zfill(4)
        theta_x_str = "m" + theta_nstr

    if theta_y_num >= 0:
        theta_nstr = '{:0>2}'.format(theta_y_num)
        theta_nstr = theta_nstr.zfill(4)
        theta_y_str = "p" + theta_nstr
    else:
        theta_nstr = '{:0>2}'.format(-theta_y_num)
        theta_nstr = theta_nstr.zfill(4)
        theta_y_str = "m" + theta_nstr

    return theta_x_str,theta_y_str
#---------------------------------------------------------------------------------------------------------------------
def set_aspect_equal_3d(ax):
    """Fix equal aspect bug for 3D plots."""

    xlim = ax.get_xlim3d()
    ylim = ax.get_ylim3d()
    zlim = ax.get_zlim3d()

    from numpy import mean
    xmean = mean(xlim)
    ymean = mean(ylim)
    zmean = mean(zlim)

    plot_radius = max([abs(lim - mean_)
                       for lims, mean_ in ((xlim, xmean),
                                           (ylim, ymean),
                                           (zlim, zmean))
                       for lim in lims])

    ax.set_xlim3d([xmean - plot_radius, xmean + plot_radius])
    ax.set_ylim3d([ymean - plot_radius, ymean + plot_radius])
    ax.set_zlim3d([zmean - plot_radius, zmean + plot_radius])

#---------------------------------------------------------------------------------------------------------------------

def EntranceBeamPopulation(nbeamx, nbeamy, f, D, d):
    """
    EntranceBeamPopulation(nbeamx,nbeamy,thetax,thetay,f,D,d)

    input arguments :
      - nbeamx : number of beam along x (should be odd, by example 11)
      - nbeamy : number of beam along y (should be odd, by example 11)
      - thetax : central beam angle along x in radian
      - f : telescope focal length in m
      - D : telescope diameter in m
      - d : distance entrance (beam start) to focal plane. the unit of d is what we want


    output arguments :
      - return X,Y coordinate of beams relative to central beam (position (0,0)) in same unit as d
    """

    numberOfRows = nbeamx * nbeamy

    X = np.linspace(-D * d / 2 / f, D * d / 2 / f, nbeamx)
    Y = np.linspace(-D * d / 2 / f, D * d / 2 / f, nbeamy)

    # need to alloc memory for the dataframe
    df = pd.DataFrame(index=np.arange(0, numberOfRows), columns=['id', 'nx', 'ny', 'X0', 'Y0'])

    idx = 0
    for ny in np.arange(len(Y)):
        for nx in np.arange(len(X)):
            df.iloc[idx] = [idx + 1, nx, ny, X[nx], Y[ny]]
            idx += 1

    df["Z0"] = -d

    return df
#---------------------------------------------------------------------------------------------------------------------

def ComputeAngles(df_in, alphax, alphay, d):
    """
    ComputeAngles(df_in,thetax,thetay,d)

    * input argument :
      - df_in :  input data frame (beam centered at beam entrance)
      - alphax,alphay :  angles in radians of central ray relative to optical beam
      - d : distance entrance -focal point in mm

    * output argument
      -df : modified dataframe


    """
    tan_alphax_0 = np.tan(alphax)
    tan_alphay_0 = np.tan(alphay)
    x = np.asarray(df_in["X0"].values).astype(np.float64)
    y = np.asarray(df_in["Y0"].values).astype(np.float64)

    tan_alphax = np.array(tan_alphax_0 - x / d)
    tan_alphay = np.array(tan_alphay_0 - y / d)

    Uz = 1. / (np.sqrt(1. + tan_alphax ** 2 + tan_alphay ** 2))
    Ux = tan_alphax * Uz
    Uy = tan_alphay * Uz

    df = df_in
    df["U0"] = Ux
    df["V0"] = Uy
    df["W0"] = Uz

    return df




#---------------------------------------------------------------------------------------------------------------------
def ShiftBeamCenter(df_in, alphax, alphay, D_in, D_disp):
    """
    ShiftBeamCenter(df_in,thetax,thetay,f,D_in,D_disp): shift the beam such
    the central ray beam is at center of disperser


    * input arguments:
      - df_in : input data frame (beam centered at beam entrance)
      - alphax,alphay :  angles in radians of central ray relative to optical beam
      - d : distance entrance -focal point in mm

    * output argument:
      -df : modified dataframe

    """

    # shift in X,Y
    # alphax < 0 : move up
    # alphax > 0 : move down
    dx = np.tan(alphax) * np.abs(D_disp - D_in)
    dy = np.tan(alphay) * np.abs(D_disp - D_in)
    df_out = df_in

    df_out["X0"] = df_in.apply(lambda x: x.X0 - dx, axis=1)
    df_out["Y0"] = df_in.apply(lambda x: x.Y0 - dy, axis=1)

    return df_out


#---------------------------------------------------------------------------------------------------------------------
def PlotTransverseBeamViewEntrance(df):

    fig = plt.figure(1,figsize=(12, 6))

    ax = fig.add_subplot(121)
    ax.scatter(df["X0"].values, df["Y0"].values, color="r")
    ax.grid()
    ax.set_xlabel("X0 (mm)")
    ax.set_ylabel("Y0 (mm)")
    ax.set_title("Beams position at Entrance")
    ax.set_aspect('equal')

    ax = fig.add_subplot(122)
    X = np.asarray(df["X0"].values).astype(np.float64)
    Y = np.asarray(df["Y0"].values).astype(np.float64)
    U = np.asarray(df["U0"].values).astype(np.float64)
    V = np.asarray(df["V0"].values).astype(np.float64)

    ax.quiver(X, Y, U, V, color="b")
    ax.grid()
    ax.set_xlabel("X0 (mm)")
    ax.set_ylabel("Y0 (mm)")
    ax.set_title("Beams direction at Entrance")
    ax.set_aspect('equal')
    plt.show()
#---------------------------------------------------------------------------------------------------------------------


def Plot3DBeamView1(df):
    """

    :param df:
    :return:
    """

    fig = plt.figure(2,figsize=(12, 6))

    X0 = np.asarray(df["X0"].values).astype(np.float64)
    Y0 = np.asarray(df["Y0"].values).astype(np.float64)
    Z0 = np.asarray(df["Z0"].values).astype(np.float64)

    U0 = np.asarray(df["U0"].values).astype(np.float64)
    V0 = np.asarray(df["V0"].values).astype(np.float64)
    W0 = np.asarray(df["W0"].values).astype(np.float64)

    X1 = np.asarray(df["X1"].values).astype(np.float64)
    Y1 = np.asarray(df["Y1"].values).astype(np.float64)
    Z1 = np.asarray(df["Z1"].values).astype(np.float64)

    N = len(df)

    all_X = np.concatenate((X0, X1))
    all_Y = np.concatenate((Y0, Y1))
    all_Z = np.concatenate((Z0, Z1))

    ax = fig.add_subplot(121, projection='3d')
    ax.view_init(azim=20)

    for i in np.arange(N):
        xs = (Y0[i], Y1[i])
        ys = (Z0[i], Z1[i])
        zs = (X0[i], X1[i])
        line = plt3d.art3d.Line3D(xs, ys, zs, color="red")
        ax.add_line(line)

    ax.set_xlim3d(all_Y.min(), all_Y.max())
    ax.set_ylim3d(all_Z.min(), all_Z.max())
    ax.set_zlim3d(all_X.min(), all_X.max())
    ax.set_xlabel('Y (mm)')
    ax.set_ylabel('Z (mm)')
    ax.set_zlabel('X (mm)')
    ax.set_title('beams at scale')

    set_aspect_equal_3d(ax)

    ax = fig.add_subplot(122, projection='3d')
    ax.view_init(azim=20)

    for i in np.arange(N):
        xs = (Y0[i], Y1[i])
        ys = (Z0[i], Z1[i])
        zs = (X0[i], X1[i])
        line = plt3d.art3d.Line3D(xs, ys, zs, color="red")
        ax.add_line(line)

    ax.set_xlim3d(all_Y.min(), all_Y.max())
    ax.set_ylim3d(all_Z.min(), all_Z.max())
    ax.set_zlim3d(all_X.min(), all_X.max())
    ax.set_xlabel('Y (mm)')
    ax.set_ylabel('Z (mm)')
    ax.set_zlabel('X (mm)')
    ax.set_title('beams not at scale')

    # set_aspect_equal_3d(ax)

    plt.tight_layout()
    plt.suptitle("Squared beam", Y=1.05, fontsize=25)
    plt.show()
#---------------------------------------------------------------------------------------------------------------------
def Plot3DBeamView2(df):
    fig = plt.figure(3, figsize=(12, 6))

    X0 = np.asarray(df["X0"].values).astype(np.float64)
    Y0 = np.asarray(df["Y0"].values).astype(np.float64)
    Z0 = np.asarray(df["Z0"].values).astype(np.float64)

    U0 = np.asarray(df["U0"].values).astype(np.float64)
    V0 = np.asarray(df["V0"].values).astype(np.float64)
    W0 = np.asarray(df["W0"].values).astype(np.float64)

    X1 = np.asarray(df["X1"].values).astype(np.float64)
    Y1 = np.asarray(df["Y1"].values).astype(np.float64)
    Z1 = np.asarray(df["Z1"].values).astype(np.float64)

    N = len(df)

    all_X = np.concatenate((X0, X1))
    all_Y = np.concatenate((Y0, Y1))
    all_Z = np.concatenate((Z0, Z1))

    ax = fig.add_subplot(121, projection='3d')
    ax.view_init(elev=0, azim=0)

    for i in np.arange(N):
        xs = (Y0[i], Y1[i])
        ys = (Z0[i], Z1[i])
        zs = (X0[i], X1[i])
        line = plt3d.art3d.Line3D(xs, ys, zs, color="red")
        ax.add_line(line)

    ax.set_xlim3d(all_Y.min(), all_Y.max())
    ax.set_ylim3d(all_Z.min(), all_Z.max())
    ax.set_zlim3d(all_X.min(), all_X.max())
    ax.set_xlabel('Y (mm)')
    ax.set_ylabel('Z (mm)')
    ax.set_zlabel('X (mm)')
    ax.set_title('beams at scale')

    set_aspect_equal_3d(ax)

    ax = fig.add_subplot(122, projection='3d')
    ax.view_init(elev=0, azim=0)

    for i in np.arange(N):
        xs = (Y0[i], Y1[i])
        ys = (Z0[i], Z1[i])
        zs = (X0[i], X1[i])
        line = plt3d.art3d.Line3D(xs, ys, zs, color="red")
        ax.add_line(line)

    ax.set_xlim3d(all_Y.min(), all_Y.max())
    ax.set_ylim3d(all_Z.min(), all_Z.max())
    ax.set_zlim3d(all_X.min(), all_X.max())
    ax.set_xlabel('Y (mm)')
    ax.set_ylabel('Z (mm)')
    ax.set_zlabel('X (mm)')
    ax.set_title('beams not at scale')

    # set_aspect_equal_3d(ax)

    plt.tight_layout()
    plt.suptitle("Squared beam (longitudinal view)", Y=1.05, fontsize=25)
    plt.show()







#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
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



    # arguments
    #----------


    parser = argparse.ArgumentParser()
    parser.add_argument("--config",action="store", dest="configfile",help=f" run generate -config configfilename, with by ex configfilename = default.ini")
    results_args = parser.parse_args()


    msg = f"Start {parser.prog} at date : {string_date} and time :{current_time} and with arguments:{results_args}"
    logger.info(msg)



    # config file
    # --------------
    configfile = ""

    config_filename = results_args.configfile
    msg = f"Configuration file : {config_filename}"
    logger.info(msg)

    # 1) CONFIGURATION
    #------------------
    logger.info('1) Configuration')

    config = configparser.ConfigParser()

    if os.path.exists(config_filename):
        config.read(config_filename)
    else:
        msg = f"config file {config_filename} does not exist !"
        logger.error(msg)

    config_section = config.sections()

    if len(config_section) == 0:
        msg = f"empty config file {config_filename} !"
        logger.error(msg)

    if 'GENERAL' in config_section:

        FLAG_DEBUG = bool(config['GENERAL']['FLAG_DEBUG'])
        FLAG_VERBOSE = bool(config['GENERAL']['FLAG_VERBOSE'])
        FLAG_PLOT = bool(config['GENERAL']['FLAG_PLOT'])
        FLAG_PRINT = bool(config['GENERAL']['FLAG_PRINT'])
    else:
        msg = f"empty section GENERAL in config file {config_filename} !"
        logger.error(msg)

    if FLAG_DEBUG:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)


    if 'BEAMFOUR' in config_section:
        output_dir = config['BEAMFOUR']['outputdir']  # output directory
        root_filename = config['BEAMFOUR']['outputrayfile']  # base ray filename
        theta_x = float(config['BEAMFOUR']['Theta_X'])  # target offset in arcmin
        theta_y = float(config['BEAMFOUR']['Theta_Y'])  # target offset in arcmin
        NBEAM_X = int(config['BEAMFOUR']['NBEAM_X'])  # nb of rays should be odd number
        NBEAM_Y = int(config['BEAMFOUR']['NBEAM_Y'])  # nb of rays should be odd number
        WLMIN = float(config['BEAMFOUR']['WLMIN'])  # minimum of wavelength in nm
        WLMAX = float(config['BEAMFOUR']['WLMAX'])  # maximum of wavelength in nm
        WLSTEP = float(config['BEAMFOUR']['WLSTEP'])  # wavelength step in nm
    else:
        msg = f"empty section BEAMFOUR in config file {config_filename} !"
        logger.error(msg)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # number of rays
    NBEAM = NBEAM_X * NBEAM_Y

    Wavelength = np.arange(WLMIN, WLMAX, WLSTEP)
    WL = (Wavelength).astype(int)
    NWL = len(WL)

    theta_x_str, theta_y_str = GetRayFileString(theta_x,theta_y)


    Beam4_Rayfile_RAY = root_filename + "{:d}_nw{:d}_thx{}_thy{}.RAY".format(NBEAM, NWL, theta_x_str, theta_y_str)

    msg = f"Outputfile {Beam4_Rayfile_RAY} "
    logger.info(msg)

    if 'TELESCOPE' in config_section:
        Tel_Focal_Length = float(config['TELESCOPE']['Tel_Focal_Length'])
        Tel_Diameter = float(config['TELESCOPE']['Tel_Diameter'])
    else:
        msg = f"empty section TELESCOPE in config file {config_filename} !"
        logger.error(msg)


    if 'OPTICS' in config_section:
        D_entrance = float(config['OPTICS']['D_entrance']) # mm : distance focal plane - entrance
        D_disperser = float(config['OPTICS']['D_disperser'])# mm : distance focal plane - disperser
    else:
        msg =f"empty section OPTICS in config file {config_filename} !"
        logger.error(msg)

    d= -D_entrance  # positive value

    ################################################
    # 2) Compute Beam Rays location
    #################################################

    logger.info('2) Compute Beam Rays location ')

    df = EntranceBeamPopulation(NBEAM_X, NBEAM_Y, Tel_Focal_Length, Tel_Diameter, -D_entrance)

    if FLAG_VERBOSE or FLAG_DEBUG:
        print(df)

    ################################################
    # 3) Compute Beam Rays orientations
    #################################################

    logger.info('3) Compute Beam Rays orientation ')

    alpha_x = -Angle(theta_x, unit=u.arcmin)
    alpha_y = -Angle(theta_y, unit=u.arcmin)

    df = ComputeAngles(df, alpha_x.radian, alpha_y.radian, -D_entrance)

    if FLAG_VERBOSE or FLAG_DEBUG:
        print(df)

    # compute the norm
    if FLAG_VERBOSE or FLAG_DEBUG:
        print(df.apply(lambda x: np.sqrt(x.U0**2+x.V0**2+x.W0**2), axis=1))

    ################################################
    # 4) Shift Beam
    #################################################

    logger.info('4) Shift Beam ')

    df = ShiftBeamCenter(df, alpha_x.radian, alpha_y.radian, D_entrance, D_disperser)

    if FLAG_VERBOSE or FLAG_DEBUG:
        print(df)

    # plot
    if FLAG_PLOT:
        PlotTransverseBeamViewEntrance(df)


    ################################################
    # 5) Compute Focal Point
    #################################################

    logger.info('5) Compute Focal Point')

    df["X1"] = df.apply(lambda x: (x.X0 + d * x.U0 / x.W0), axis=1)
    df["Y1"] = df.apply(lambda x: (x.Y0 + d * x.V0 / x.W0), axis=1)
    df["Z1"] = 0

    if FLAG_VERBOSE or FLAG_DEBUG:
        print(df)


    if FLAG_PLOT:
        Plot3DBeamView1(df)

        Plot3DBeamView2(df)

    ################################################
    # 6) Write output file
    #################################################

    logger.info(f'6) Write BEAMFOUR Ray output file : {Beam4_Rayfile_RAY}')


    # -----------------------------------------------------------
    def header(N, filename):
        """
        header(N,filename) : provide lines for the header of BeamFour Ray file.
        - input :
            - N : Number of beam rays
            - filename : Ray filename
        - output :
           the three header lines

        The width of the column must be kept fixed. The width = 12 is chosen.
        Note the : is the column separator used by BeamFour.

        """
        wd = ['X0', 'Y0', 'Z0', 'U0', 'V0', 'W0', '@wave', 'X1', 'Y1', 'Z1', 'X2', 'Y2', 'Z2', 'X3', 'Y3', 'Z3',
              'Xgoal', 'Ygoal', 'Xfinal', 'Yfinal', 'Zfinal', 'Notes']
        wd2 = "----------:"

        line1 = "{:d} {}".format(N, filename)
        line1 += os.linesep

        line2 = ""
        Nwd = len(wd)

        for i in np.arange(Nwd):
            line2 += f"{wd[i]:^11}"
        line2 += os.linesep

        line3 = ""
        for i in np.arange(Nwd - 1):
            line3 += f"{wd2:^11}"
        line3 += "------------------"

        line3 += os.linesep

        return line1, line2, line3
    #-----------------------------------------------------------
    def GetLine(x0, y0, z0, u0, v0, w0, wa):
        """
        GetLine(x0,y0,z0,u0,v0,w0,wa) : provide the line format for each ray information.

        - input :
           - (x0,y0,z0) : coordinates of rays at beam entrance in mm unit
           - (u0,v0,w0) : cosinus director of rays at beam entrance
           - wa         : wavelength of the beam ray in mm

        - output :
           - the beam ray line.

        The width of the column must be kept fixed. The width = 12 is chosen.
        Note the : is the column separator used by BeamFour. The column width include the separator :

        """

        wd2 = "----------:"
        wd0 = "          :"

        line = f"{x0: {10.5}}:" + f"{y0: {10.5}}:" + f"{z0: {10}}:" + f"{u0: {10.5}}:" + f"{v0: {10.5}}:" + f"{w0: {10.5}}:" + f"{wa: {10.6}}:"

        for i in np.arange(14):
            line += wd0

        line += "             "
        line += os.linesep

        return line
    #-----------------------------------------------------------------------------

    # open output file
    f = open(os.path.join(output_dir,Beam4_Rayfile_RAY), 'w')

    # header
    line1, line2, line3 = header(NBEAM * NWL, Beam4_Rayfile_RAY)

    f.write(line1)
    f.write(line2)
    f.write(line3)

    # write rows
    # loop on wavelengths
    for iwl in np.arange(NWL):
        wl = WL[iwl] * 1e-6  # wavelength in mm
        # loop on the NBEAM rays
        for idx in np.arange(len(df)):
            line4 = GetLine(df["X0"][idx], df["Y0"][idx], df["Z0"][idx], df["U0"][idx], df["V0"][idx], df["W0"][idx],
                            wl)
            f.write(line4)

    # close file
    f.close()
#----------------------------------------------------------------------------------------------------------------------

