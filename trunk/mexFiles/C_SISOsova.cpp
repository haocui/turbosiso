/** \file
 * 
 * \brief Matlab interface for SISO RSC module using SOVA
 * 
 * Matlab function declaration: extrinsic_data =
 * 		C_SISOsova(intrinsic_coded, apriori_data, bin_gen, win_len)
 * 
 * compilation command: mex -litpp -lacml C_SISOsova.cpp
 */

#include "itpp/itmex.h"
#include "SISO.h"

using namespace itpp;

void mexFunction(int n_output, mxArray *output[], int n_input, const mxArray *input[])
{
    // Check the number of inputs and output arguments
    if(n_output!=2) mexErrMsgTxt("Two output required!");
    if(n_input!=8) mexErrMsgTxt("Eight inputs required!");
    
    // Convert input variables to IT++ format
    vec intrinsic_coded = mxArray2vec(input[0]);
    vec apriori_data = mxArray2vec(input[1]);
    bmat bin_gen = mxArray2bmat(input[2]);
    int win_len = mxArray2int(input[3]);
    bool hard_output_flag = (mxArray2double(input[4])==1.0);
    std::string map_metric = mxArray2string(input[5]);
    double matching_scaling_factor = mxArray2double(input[6]);
    double nonmatching_scaling_factor = mxArray2double(input[7]);
    
    // ------------------ Start of routine ---------------------------
    tr::SISO siso; 
    siso.set_generators(bin_gen);
    siso.set_viterbi_win_len(win_len);
    siso.set_viterbi_hard_output_flag(hard_output_flag);
    siso.set_map_metric(map_metric);
    siso.set_viterbi_scaling_factors(matching_scaling_factor, nonmatching_scaling_factor);
    siso.set_sova_scaling_factor(matching_scaling_factor);//last two inputs are used for both Viterbi and SOVA
    siso.set_sova_threshold(nonmatching_scaling_factor);
    vec extrinsic_data;
    vec extrinsic_coded;
    siso.rsc(extrinsic_coded, extrinsic_data, intrinsic_coded, apriori_data);
    // ------------------ End of routine -----------------------------
    
    // Convert the IT++ format to Matlab format for output
    output[0] = mxCreateDoubleMatrix(1, extrinsic_coded.length(), mxREAL);
    vec2mxArray(extrinsic_coded, output[0]);
    output[1] = mxCreateDoubleMatrix(1, extrinsic_data.length(), mxREAL);
    vec2mxArray(extrinsic_data, output[1]);
}
