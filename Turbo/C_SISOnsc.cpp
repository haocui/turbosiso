/** \file
 * 
 * \brief Matlab interface for SISO NSC module
 * 
 * Matlab function declaration: [extrinsic_coded extrinsic_data] = C_SISOnsc(intrinsic_coded, apriori_data, bin_gen, scrambler_pattern, tail, map_metric)
 * 
 * compilation command: mex -litpp -lacml C_SISOnsc.cpp
 */

#include "itpp/itmex.h"
#include "/home/bogdan/C++/SISO/SISO.cpp"

using namespace itpp;

void mexFunction(int n_output, mxArray *output[], int n_input, const mxArray *input[])
{
    // Check the number of inputs and output arguments
    if(n_output!=2) mexErrMsgTxt("Two outputs required!");
    if(n_input!=6) mexErrMsgTxt("Six inputs required!");
    
    // Convert input variables to IT++ format
    vec intrinsic_coded = mxArray2vec(input[0]);
    vec apriori_data = mxArray2vec(input[1]);
    bmat bin_gen = mxArray2bmat(input[2]);
    vec scrambler_pattern = mxArray2vec(input[3]);
    bool tail = (mxArray2double(input[4])==1);
    std::string map_metric = mxArray2string(input[5]);
    
    // ------------------ Start of routine ---------------------------
    tr::SISO siso; 
    siso.set_generators(bin_gen);
    siso.set_scrambler_pattern(scrambler_pattern);
    siso.set_tail(tail);
    siso.set_map_metric(map_metric);
    vec extrinsic_data;
    vec extrinsic_coded;
    siso.nsc(extrinsic_coded, extrinsic_data, intrinsic_coded, apriori_data);
    // ------------------ End of routine -----------------------------
    
    // Convert the IT++ format to Matlab format for output
    output[0] = mxCreateDoubleMatrix(1, extrinsic_coded.length(), mxREAL);
    vec2mxArray(extrinsic_coded, output[0]);
    output[1] = mxCreateDoubleMatrix(1, extrinsic_data.length(), mxREAL);
    vec2mxArray(extrinsic_data, output[1]);
}
