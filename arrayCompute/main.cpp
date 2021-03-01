#include <iostream>
#include <algorithm>
#include <functional>
#include <string>
#include <array>

#include "operators.h"

//using namespace std;

/*
Copyright Netherlands eScience Center
Function         : Extract Meteorological fields from MERRA2 and array manipulation
Author           : Yang Liu (y.liu@esciencecenter.nl)
First Built      : 2020.09.20
Last Update      : 2020.10.11
Contributor      :
Description      : This module aims to test the speed of input file loading (standard netCDF
                   files) and array manipulation. It conducts quantification of meridional
                   energy transport.
                   The sample data is downloaded directly from online data system of NASA,
                   namely MERRA2. MERRA2 is a state-of-the-art atmosphere reanalysis product
                          produced by NASA. It spans from 1980 to 2017. Natively it is generated on
                          a hybrid sigma grid with 72 vertical levels.
Return Values    : float
Caveat!          : This module is designed to work with a batch of files. Hence, there is
                   pre-requists for the location and arrangement of data.

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
// gz is [m2 / s2] = [ kg m2 / kg s2 ] = [J / kg]

// multiply by v: [J / kg] * [m / s] => [J m / kg s]
// sum over longitudes [J m / kg s] * [ m ] = [J m2 / kg s]

// integrate over pressure: dp: [Pa] = [N m-2] = [kg m2 s-2 m-2] = [kg s-2]
// [J m2 / kg s] * [Pa] = [J m2 / kg s] * [kg / s2] = [J m2 / s3]
// and factor 1/g: [J m2 / s3] * [s2 /m2] = [J / s] = [Wat]
// ########################################################

class MERRA2{
public:
    // path to the data
    std::string path;
    // 0.75 deg per grid box latitudinally
    const int lat_unit = 360;
    static const int eta_level = 73; // eta levels of model

    MERRA2()
        : path("H:\\") //variables should be initialized in order.
    {
        std::cout << "Start processing MERRA2 data!" << std::endl;

    }

    void datapath(std::string inpath){
    //this->path = path;
    path = inpath;
    std::cout << "The path to the input files: " << path << std::endl;
   }

   // Define sigma levels. For more information, please visit the website of ECMWF.
   // Since there are 60 model levels, there are 61 half levels, so it is for A and B values.
   // A and B values for the definition of sigma level list
   // Since there are 72 model levels, there are 73 half levels, so it is for A and B values
   // the unit of A is hPa!!!!!!!!!!!!
   // from surface to TOA
    float A[eta_level]={
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
   //reverse(std::begin(A), std::end(A));

   // note that the unit of B is 1!!!!!!!!!!!!
   // from surface to TOA

   float B[eta_level] = {
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
   //reverse(std::begin(B), std::end(B));

    ~MERRA2()
    {
        std::cout << "Finish processing MERRA2 data!" << std::endl;

    }

};

template <typename T>
void reverse_vector(T p_array[], const int numcols){
    // print call for the function
    std::cout << "reverse a 1D vector" << std::endl;
    // an intermediate variable
    T temp;
    for (int i = 0; i < numcols/2; ++i) {
        temp = p_array[numcols-i-1];
        p_array[numcols-i-1] = p_array[i];
        p_array[i] = temp;
    }
}

int main()
{
    std::cout << "Hello world!" << std::endl;

    // data type
    //char singleChar = 'A';
    std::string datapath = "H:\\Creator_Zone\\Script_craft\\SonicBench\\data";
    //int a = 10;
    float b = 2.4;
    double c = 1.24;
    bool check = true;

    std::cout << datapath << std::endl;
    std::cout << datapath.length() << std::endl;
    // indexing from 0
    //cout << datapath[2];

    MERRA2 instance;
    instance.datapath(datapath);

    float AAA[]={
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

    static const int numcol = sizeof(AAA)/sizeof(float);

    std::cout << numcol << std::endl;

    //reverse(std::begin(A), std::end(A));


    //for  (int i = 0; i < end(A)-begin(A); i++) {
    //cout << A[i] << " ";
    //}


    std::cout << "\n";

    // print all the items in AAA
    for (auto item:AAA) {
    std::cout << item << " ";
    }

    std::cout << std::endl;

    //std::transform(A.begin(), A.end(), A.begin(),
    //            std::bind(std::multiplies<T>(), std::placeholders::_1, 100));

    // define an empty pointer
    std::cout << &AAA[0];
    //int* p_AAA = &AAA[0];

    // waiting for any input to continue
    std::cin.get();

    // call the function to reverse a 1D vector
    reverse_vector(AAA, numcol);

     for (auto item:AAA) {
    std::cout << item << " ";
    }


    return 0;
}
