package com.array;

import java.util.Arrays;

/*
Copyright Netherlands eScience Center
Function      : Extract Meteorological fields from MERRA2 and array manipulation
Author        : Yang Liu (y.liu@esciencecenter.nl)
First Built   : 2020.07.03
Last Update   : 2020.07.03
Contributor   :
Description   : This module aims to test the speed of input file loading (standard netCDF
                files) and array manipulation. It conducts quantification of meridional
                energy transport
                The sample data is downloaded directly from online data system of NASA,
                namely MERRA2. MERRA2 is a state-of-the-art atmosphere reanalysis product
                produced by NASA. It spans from 1980 to 2017. Natively it is generated on
                a hybrid sigma grid with 72 vertical levels.
Return Values : float
Caveat!       : This module is designed to work with a batch of files. Hence, there is
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
 */

// ######################    Units vacabulory    ######################
// cpT:  [J / kg K] * [K]     = [J / kg]
// Lvq:  [J / kg] * [kg / kg] = [J / kg]

class MERRA2 {
    /*

     */
    String path;

    public void path(String path) {
        this.path = path;
        System.out.println(this.path);
    }



}

public class Main {
    public static void main(String[] args) {
        // *****************  input zone  ********************
	    // specify the datapath
        String datapath = "H:\\Creator_Zone\\Script_craft\\SonicBench\\data";
        // *****************  call function  ********************
        MERRA2 instance = new MERRA2();
        instance.path(datapath);
    }
}
