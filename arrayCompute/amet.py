# -*- coding: utf-8 -*-
"""
Copyright Netherlands eScience Center
Function        : Calculate Meridional Energy Transport in the Atmosphere with Reanalysis
Author          : Yang Liu (y.liu@esciencecenter.nl)
First Built     : 2020.07.03
Last Update     : 2020.07.03
Contributor     :
Description     : This module provides a method to perform the computation
                  of meridional energy transport in the atmosphere.

                  It serves to test the array manipulation speed of the
                  programming language.
Return Values   : netCDF files
Caveat!         : It is highly recommended to compute all the fluxes on native
                  grids.
"""

import numpy as np
from netCDF4 import Dataset
import os
import platform
import sys
import logging

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

class met:
    def __init__(self):
        """
        Quantify the meridional energy transport and its components upto
        certain pressure levels.
        """
        print ("Start quantifying the meridional energy transport in the atmosphere.")

    def calc_met(self, T, q, sp, u, v, gz, A, B, t, h, y, x,
                 lat, lat_unit):
        """
        Calculate the meridional energy transport and its components in the atmosphere.
        All the input files should contain the fields for the entire month.
        Caveat! The integral is taken from TOA to the Surface!
        param q: Specific Humidity [kg/kg]
        param sp: Surface Pressure [Pa]
        param u: Zonal Wind [m/s]
        param v: Meridional Wind [m/s]
        param gz: Geopotential [m2/s2]
        param A: Constant A for Defining Sigma Level
        param B: Constant B for Defining Sigma Level
        param t: time dimension of input fields
        param h: level dimension of input fields
        param y: latitude dimension of input fields
        param x: longitude dimension of input fields
        param lat: latitude
        param lat_unit: number of grid boxes meridionally (to calculate the unit width)
        returns: AMET and its components (internal, latent, geopotential, kinetic energy)
                 upto different heights.
        rtype: numpy array
        """
        # generate the contants
        constant = setConstants()
        # calculate dp
        dp_level = np.zeros((t, h, y, x),dtype = float)
        for i in np.arange(h):
            dp_level[:,i,:,:] =  np.abs((A[i+1] + B[i+1] * sp) - (A[i] + B[i] * sp))
        # calculate each component of total energy and take the vertical integral
        # Internal Energy cpT
        internal_energy = constant['cp'] * (1-q) * T * dp_level / constant['g']
        internal_flux_int = np.mean(np.sum(internal_energy * v,1),0)
        del internal_energy
        logging.info("The calculation of internal energy flux is finished!")
        # Latent heat Lvq
        latent_energy = constant['Lv'] * q * dp_level / constant['g']
        latent_flux_int = np.mean(np.sum(latent_energy * v,1),0)
        del latent_energy
        logging.info("The calculation of latent heat flux is finished!")
        # Geopotential Energy gz
        geopotential_energy = gz * dp_level / constant['g']
        geopotential_flux_int = np.mean(np.sum(geopotential_energy * v,1),0)
        del geopotential_energy
        logging.info("The calculation of geopotential energy flux is finished!")
        # Kinetic Energy u2+v2
        kinetic_energy = 0.5 * (u**2 + v**2) * dp_level / constant['g']
        kinetic_flux_int = np.mean(np.sum(kinetic_energy * v,1),0)
        del kinetic_energy
        logging.info("The calculation of kinetic energy flux is finished!")
        # the earth is taken as a perfect sphere, instead of a ellopsoid
        dx = 2 * np.pi * constant['R'] * np.cos(2 * np.pi * lat / 360) / x
        # plugin the weight of grid box width and apply the correction
        #dx[0] = 0
        # create arrays for each energy transport components after correction
        E_internal = np.zeros((y, x),dtype=float)
        E_latent = np.zeros((y, x),dtype=float)
        E_geopotential = np.zeros((y, x),dtype=float)
        E_kinetic = np.zeros((y, x),dtype=float)
        E_total = np.zeros((y, x),dtype=float)
        # apply correction and weight of grid width
        # also change the unit from Watt to Tera Watt
        for i in np.arange(y):
            E_internal[i,:] = internal_flux_int[i,:] * dx[i]/1e+12
            E_latent[i,:] = latent_flux_int[i,:] * dx[i]/1e+12
            E_geopotential[i,:] = geopotential_flux_int[i,:] * dx[i]/1e+12
            E_kinetic[i,:] = kinetic_flux_int[i,:] * dx[i]/1e+12

        E_total = E_internal + E_latent + E_geopotential + E_kinetic

        logging.info("Weight the energy transport by the width of grid box.")
        logging.info("Computation of meridional energy transport on model level is finished!")

        return E_total, E_internal, E_latent, E_geopotential, E_kinetic
