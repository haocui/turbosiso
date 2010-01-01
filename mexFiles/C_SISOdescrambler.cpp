/** \file
 * 
 * \brief Matlab interface for SISO descrambler module
 * 
 * Matlab function declaration: [extrinsic_coded extrinsic_data] = C_SISOdescrambler(intrinsic_coded, apriori_data, scrambler_pattern)
 * 
 * compilation command: mex -litpp -lacml C_SISOdescrambler.cpp
 */

#include "itpp/itmex.h"
#include "SISO.h"

using namespace itpp;

void mexFunction(int n_output, mxArray *output[], int n_input, const mxArray *input[])
{
    // Check the number of inputs and output arguments
    if(n_output!=2) mexErrMsgTxt("Two outputs required!");
    if(n_input!=3) mexErrMsgTxt("Three inputs required!");
    
    // Convert input variables to IT++ format
    vec intrinsic_coded = mxArray2vec(input[0]);
    vec apriori_data = mxArray2vec(input[1]);
    vec scrambler_pattern = mxArray2vec(input[2]);
    
    // ------------------ Start of routine ---------------------------
    tr::SISO siso;
    siso.set_scrambler_pattern(scrambler_pattern);    
    vec extrinsic_data;
    vec extrinsic_coded;
    siso.descrambler(extrinsic_coded, extrinsic_data, intrinsic_coded, apriori_data);
    // ------------------ End of routine -----------------------------
    
    // Convert the IT++ format to Matlab format for output
    output[0] = mxCreateDoubleMatrix(1, extrinsic_coded.length(), mxREAL);
    vec2mxArray(extrinsic_coded, output[0]);
    output[1] = mxCreateDoubleMatrix(1, extrinsic_data.length(), mxREAL);
    vec2mxArray(extrinsic_data, output[1]);
}
