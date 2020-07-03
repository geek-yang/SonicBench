# -*- coding: utf-8 -*-
"""
Copyright Netherlands eScience Center
Function        : Extract Meteorological fields from MERRA2 and array manipulation
Author          : Yang Liu (y.liu@esciencecenter.nl)
First Built     : 2020.07.03
Last Update     : 2020.07.03
Contributor     :
Description     : This module aims to test the speed of input file loading (standard netCDF
                  files) and array manipulation. It conducts quantification of meridional
                  energy transport

                  The sample data is downloaded directly from online data system of NASA,
                  namely MERRA2. MERRA2 is a state-of-the-art atmosphere reanalysis product
                  produced by NASA. It spans from 1980 to 2017. Natively it is generated on
                  a hybrid sigma grid with 72 vertical levels.
Return Values   : netCDF files
Caveat!         : This module is designed to work with a batch of files. Hence, there is
                  pre-requists for the location and arrangement of data. The folder should
                  have the following structure:

                  Please use the default names after downloading from NASA.
                  The files are in netCDF4 format. By default, the latitude
                  ascends.

                  The native grid of MERRA2 is different from normal reanalyses.
                  It is not natively generated on a reduced Gaussian Grid. All the
                  variables are computed on a cubed-sphere grid using GEOS-5 model.
                  As a result, it is not suitable to bring the point data back to
                  spherical harmonics for the calculations of divergence and gradients.
                  More information about its grid set-up is given in the official
                  documentation:
                  https://gmao.gsfc.nasa.gov/pubs/docs/Bosilovich785.pdf
"""

##########################################################################
###########################   Units vacabulory   #########################
# cpT:  [J / kg K] * [K]     = [J / kg]
# Lvq:  [J / kg] * [kg / kg] = [J / kg]
# gz is [m2 / s2] = [ kg m2 / kg s2 ] = [J / kg]

# multiply by v: [J / kg] * [m / s] => [J m / kg s]
# sum over longitudes [J m / kg s] * [ m ] = [J m2 / kg s]

# integrate over pressure: dp: [Pa] = [N m-2] = [kg m2 s-2 m-2] = [kg s-2]
# [J m2 / kg s] * [Pa] = [J m2 / kg s] * [kg / s2] = [J m2 / s3]
# and factor 1/g: [J m2 / s3] * [s2 /m2] = [J / s] = [Wat]
##########################################################################

import sys
import os
import numpy as np
from netCDF4 import Dataset
import logging
import amet
import time

class merra2:
    def __init__(self, path):
        """
        Initialize the extraction of fields from MERRA2.
        The data is on hybrid sigma levels. As the interpolation can introduce
        large errors to the computation of energy transport, we will follow the
        model level. The determination of reference levels is based on the estimation of
        pressure on each level with standard surface pressure 1013.25 hPa
        (see ERA-Interim archive by ECMWF.). For the actual calculation, the varying
        surface pressure should be taken into account.
        param path: the root path of the input fields
        param lat_unit: number of grid boxes meridionally (to calculate the unit width)
        """
        self.path = path
        # 0.75 deg per grid box latitudinally
        self.lat_unit = 360


    @staticmethod
    def defineSigmaLevels():
        """
        Definine sigma levels. For more information, please visit the website of ECMWF.
        Since there are 60 model levels, there are 61 half levels, so it is for A and B values.
        returns: tuple containing arrays with A and B values for the definition of
                 sigma levellist
        rtype: tuple
        """
        # A and B values for the definition of sigma levelist
        # Since there are 72 model levels, there are 73 half levels, so it is for A and B values
        # the unit of A is hPa!!!!!!!!!!!!
        # from surface to TOA
        A = np.array([
            0.000000e+00, 4.804826e-02, 6.593752e+00, 1.313480e+01, 1.961311e+01, 2.609201e+01,
            3.257081e+01, 3.898201e+01, 4.533901e+01, 5.169611e+01, 5.805321e+01, 6.436264e+01,
            7.062198e+01, 7.883422e+01, 8.909992e+01, 9.936521e+01, 1.091817e+02, 1.189586e+02,
            1.286959e+02, 1.429100e+02, 1.562600e+02, 1.696090e+02, 1.816190e+02, 1.930970e+02,
            2.032590e+02, 2.121500e+02, 2.187760e+02, 2.238980e+02, 2.243630e+02, 2.168650e+02,
            2.011920e+02, 1.769300e+02, 1.503930e+02, 1.278370e+02, 1.086630e+02, 9.236572e+01,
            7.851231e+01, 6.660341e+01, 5.638791e+01, 4.764391e+01, 4.017541e+01, 3.381001e+01,
            2.836781e+01, 2.373041e+01, 1.979160e+01, 1.645710e+01, 1.364340e+01, 1.127690e+01,
            9.292942e+00, 7.619842e+00, 6.216801e+00, 5.046801e+00, 4.076571e+00, 3.276431e+00,
            2.620211e+00, 2.084970e+00, 1.650790e+00, 1.300510e+00, 1.019440e+00, 7.951341e-01,
            6.167791e-01, 4.758061e-01, 3.650411e-01, 2.785261e-01, 2.113490e-01, 1.594950e-01,
            1.197030e-01, 8.934502e-02, 6.600001e-02, 4.758501e-02, 3.270000e-02, 2.000000e-02,
            1.000000e-02,],dtype=float)
            # reverse A
        A = A[::-1] * 100 # change unit to Pa
            # the unit of B is 1!!!!!!!!!!!!
            # from surfac eto TOA
        B = np.array([
            1.000000e+00, 9.849520e-01, 9.634060e-01, 9.418650e-01, 9.203870e-01, 8.989080e-01,
            8.774290e-01, 8.560180e-01, 8.346609e-01, 8.133039e-01, 7.919469e-01, 7.706375e-01,
            7.493782e-01, 7.211660e-01, 6.858999e-01, 6.506349e-01, 6.158184e-01, 5.810415e-01,
            5.463042e-01, 4.945902e-01, 4.437402e-01, 3.928911e-01, 3.433811e-01, 2.944031e-01,
            2.467411e-01, 2.003501e-01, 1.562241e-01, 1.136021e-01, 6.372006e-02, 2.801004e-02,
            6.960025e-03, 8.175413e-09, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
            0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
            0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
            0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
            0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
            0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
            0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
            0.000000e+00,],dtype=float)
            # reverse B
        B = B[::-1]

        return (A, B)

    @staticmethod
    def setConstants():
        '''
        Define constants used in the calculations. The constants include:
        const g: gravititional acceleration [m / s2]
        const R: radius of the earth [m]
        const cp: heat capacity of air [J/(Kg*K)]
        const Lv: Latent heat of vaporization [J/Kg]
        const R_dry: gas constant of dry air [J/(kg*K)]
        const R_vap: gas constant for water vapour [J/(kg*K)]
        returns: dictionary with constants for g, R. cp, Lv, R_dry and R_vap
        rtype: dict
        '''
        # define the constant:
        constant = {'g' : 9.80616,      # gravititional acceleration [m / s2]
                    'R' : 6371009,      # radius of the earth [m]
                    'cp': 1004.64,      # heat capacity of air [J/(Kg*K)]
                    'Lv': 2264670,      # Latent heat of vaporization [J/Kg]
                    'R_dry' : 286.9,    # gas constant of dry air [J/(kg*K)]
                    'R_vap' : 461.5,    # gas constant for water vapour [J/(kg*K)]
                   }

        return constant

    def amet_memoryWise(self):
        """
        Quantify Meridional Energy Transport.
        param year_start: the starting time for the calculation
        param year_end: the ending time for the calculation
        param path_uvc: location of the baratropic corrected winds
        param fields: number of fields contained in one file, two options available
        - 1 (default) two seperate files with T,q,u,v on multiple sigma levels and lnsp,z on surface
        - 2 three seperate files, T,q and u,v and lnsp,z
        param example: an example input file for loading dimensions (level)
        return: arrays containing AMET and its components upto differnt pressure levels
        rtype: netCDF4
        """
         # set up logging files to monitor the calculation
        logging.basicConfig(filename = os.path.join(self.path,'history_amet_python.log'),
                            filemode = 'w+', level = logging.DEBUG,
                            format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        # initialize the time span
        # define sigma level
        A, B = self.defineSigmaLevels()
        # use example input file to load the basic dimensions information
        datapath_var = os.path.join(self.path, 'MERRA2_400.inst3_3d_asm_Nv.20160101.nc4.nc')
        var_key = Dataset(datapath_var)
        lat = var_key.variables['lat'][:]
        lon = var_key.variables['lon'][:]
        # calculate the reference levels based on A & B and standard surface pressure
        half_level = A + B * 101325
        level = (half_level[1:] + half_level[:-1]) / 2
        # create space for the output
        # AMET in the entire column
        E = np.zeros((len(lat),len(lon)), dtype=float)
        cpT = np.zeros((len(lat),len(lon)), dtype=float)
        Lvq = np.zeros((len(lat),len(lon)), dtype=float)
        gz = np.zeros((len(lat),len(lon)), dtype=float)
        uv2 = np.zeros((len(lat),len(lon)), dtype=float)
        logging.info("Start retrieving variables T,q,u,v,sp")
        # The shape of each variable is (8,72,361,576)
        T = var_key.variables['T'][:]
        q = var_key.variables['QV'][:]
        sp = var_key.variables['PS'][:] #(8,361,576)
        u = var_key.variables['U'][:]
        v = var_key.variables['V'][:]
        logging.info("Extracting variables successfully!")        
        # compute gz
        z_model = self.calc_gz(var_key)
        # get the basic shape
        tt, hh, yy, xx = q.shape
        AMET = amet.met()
        E, cpT, Lvq, gz, uv2 = AMET.calc_met(T, q, sp, u, v, z_model, A, B,
                              tt, hh, len(lat), len(lon), lat, self.lat_unit)

        return np.mean(E)

    def calc_gz(self, var_key):
        """
        Calculate geopotential on sigma levels.
        The method is given in ECMWF IFS 9220, from section 2.20 - 2.23.
        param T: Absolute Temperature  [K]
        param q: Specific Humidity     [kg/kg]
        param sp: Surface Pressure     [Pa]
        param z: Surface Geopotential  [m2/s2]
        param A: Constant A for Defining Sigma Level [Pa]
        param B: Constant B for Defining Sigma Level
        param t: time dimension of input fields
        param h: level dimension of input fields
        param y: latitude dimension of input fields
        param x: longitude dimension of input fields
        return: An array of geopotential.
        rtype: numpy array
        """
        logging.info("Start the computation of geopotential on model level.")
        # call the function to generate contants
        constant = self.setConstants()
        # define sigma level
        A, B = self.defineSigmaLevels()
        # extract variables
        T = var_key.variables['T'][:]
        q = var_key.variables['QV'][:]
        sp = var_key.variables['PS'][:]
        z = var_key.variables['PHIS'][:]
        # basic dimension
        t, h, y, x = T.shape
        # define the half level pressure matrix
        p_half_plus = np.zeros((t, h, y, x),dtype = float)
        p_half_minus = np.zeros((t, h, y, x),dtype = float)
        # calculate the index of pressure levels
        index_level = np.arange(h)
        # calculate the pressure at each half level
        for i in index_level:
            p_half_plus[:,i,:,:] = A[i+1] + B[i+1] * sp # A is Pa already
            p_half_minus[:,i,:,:] = A[i] + B[i] * sp
        # calculate full pressure level
        #level_full = (p_half_plus + p_half_minus) / 2
        # compute the moist temperature (virtual temperature)
        Tv = T * (1 + (constant['R_vap'] / constant['R_dry'] - 1) * q)
        # initialize the first half level geopotential
        gz_half = np.zeros((t, y, x),dtype =float)
        # initialize the full level geopotential
        gz = np.zeros((t, h, y, x),dtype = float)
        # Calculate the geopotential at each level
        # The integral should be taken from surface level to the TOA
        for i in index_level:
            # reverse the index
            i_inverse = h - 1 - i
            # the ln(p_plus/p_minus) is calculated, alpha is defined
            # an exception lies in the TOA
            # see equation 2.23 in ECMWF IFS 9220
            if i_inverse == 0:
                ln_p = np.log(p_half_plus[:,i_inverse,:,:]/10)
                alpha = np.log(2)
            else:
                ln_p = np.log(p_half_plus[:,i_inverse,:,:]/p_half_minus[:,i_inverse,:,:])
                delta_p = p_half_plus[:,i_inverse,:,:] - p_half_minus[:,i_inverse,:,:]
                alpha = 1 - p_half_minus[:,i_inverse,:,:] / delta_p * ln_p
            # calculate the geopotential of the full level (exclude surface geopotential)
            # see equation 2.22 in ECMWF IFS 9220
            gz_full = gz_half + alpha * constant['R_dry'] * Tv[:,i_inverse,:,:]
            # add surface geopotential to the full level
            # see equation 2.21 in ECMWF IFS 9220
            gz[:,i_inverse,:,:] = z + gz_full
            # renew the half level geopotential for next loop step (from p_half_minus level to p_half_plus level)
            # see equation 2.20 in ECMWF IFS 9220
            gz_half = gz_half + ln_p * constant['R_dry'] * Tv[:,i_inverse,:,:]
        logging.info("Computation of geopotential on model level is finished!")

        return gz

if __name__=="__main__":
    start_time = time.time()
    # sample
    ################################   Input zone  ######################################
    # specify data path
    datapath_MERRA2 = "/home/ESLT0068/NLeSC/Computation_Modeling/SonicBench/data"
    #####################################################################################
    print ('*********************** call functions *************************')
    instance = merra2(datapath_MERRA2)
    E = instance.amet_memoryWise()
    print("Mean total energy transport of a day (TW)",E)
    print("---{} seconds ---".format(time.time() - start_time))

