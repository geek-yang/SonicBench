package com.array;

import java.lang.reflect.Array;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

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
    Initialize the extraction of fields from MERRA2.
    The data is on hybrid sigma levels. As the interpolation can introduce
    large errors to the computation of energy transport, we will follow the
    model level. The determination of reference levels is based on the estimation of
    pressure on each level with standard surface pressure 1013.25 hPa
    (see ERA-Interim archive by ECMWF). For the actual calculation, the varying
    surface pressure should be taken into account.
    param path: the root path of the input fields
    param lat_unit: number of grid boxes meridionally (to calculate the unit width)
     */
    // path to the data
    String path;
    // 0.75 deg per grid box latitudinally
    final int lat_unit = 360;

    public void path(String path) {
        this.path = path;
        System.out.println(this.path);
    }

    // Define sigma levels. For more information, please visit the website of ECMWF.
    // Since there are 60 model levels, there are 61 half levels, so it is for A and B values.
    // A and B values for the definition of sigma level list
    // Since there are 72 model levels, there are 73 half levels, so it is for A and B values
    // the unit of A is hPa!!!!!!!!!!!!
    // from surface to TOA
    final double[] A = {
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
            1.000000e-02};

    // reverse A


    final double[] B = {
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
            0.000000e+00};

    // reverse B

    /* Define constants used in the calculations. The constants include:
    const g: gravititional acceleration [m / s2]
    const R: radius of the earth [m]
    const cp: heat capacity of air [J/(Kg*K)]
    const Lv: Latent heat of vaporization [J/Kg]
    const R_dry: gas constant of dry air [J/(kg*K)]
    const R_vap: gas constant for water vapour [J/(kg*K)]
    */
    double g = 9.80616;    // gravititional acceleration [m / s2]
    double R = 6371009;    // radius of the earth [m]
    double cp = 1004.64;   // heat capacity of air [J/(Kg*K)]
    double Lv = 2264670;   // latent heat of vaporization [J/Kg]
    double R_dry = 286.9;  // gas constant of dry air [J/(kg*K)]
    double R_vap = 461.5;  // gas constant for water vapour [J/(kg*K)]

    public double amet_memorywise() {
        /*
        Quantify Meridional Energy Transport.
         */
        // use example input file to load the basic dimensions information
        System.out.println("Extracting input fields");

        double E_mean = 4.0001;

        return E_mean;

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
        // check the A,B matrix
        double[] A = {
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
                1.000000e-02};
        System.out.println(Arrays.toString(A));

        // reverse the array
        for(int i = 0; i < A.length/2; i++)
        {
            double temp = A[i];
            A[i] = A[A.length - i - 1];
            A[A.length - i - 1] = temp;
        }
        System.out.println("After reverse!");

        System.out.println(Arrays.toString(A));
    }
}
