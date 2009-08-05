/** \file
 * 
 * \brief Matlab interface for SISO equalizer module
 * 
 * Matlab function declaration: extrinsic_data = C_SISOequalizer(rec_sig, apriori_data, ch_imp_response, prec_gen, sigma2)
 * 
 * compilation command: mex -litpp -lacml C_SISOequalizer.cpp
 */

#include "itpp/itmex.h"
#include "SISO.cpp"

using namespace itpp;

void mexFunction(int n_output, mxArray *output[], int n_input, const mxArray *input[])
{
    // Check the number of inputs and output arguments
    if(n_output!=1) mexErrMsgTxt("One output required!");
    if(n_input!=5) mexErrMsgTxt("Five inputs required!");
    
    // Convert input variables to IT++ format
    vec rec_sig = mxArray2vec(input[0]);
    vec apriori_data = mxArray2vec(input[1]);
    vec ch_imp_response = mxArray2vec(input[2]);
    bvec prec_gen = mxArray2bvec(input[3]);
    double sigma2 = mxArray2double(input[4]);
    
    // ------------------ Start of routine ---------------------------
    tr::SISO siso;
    siso.set_precoder_generator(prec_gen);    
    siso.set_impulse_response(ch_imp_response);
    siso.set_noise(sigma2);
    siso.set_map_metric("logMAP");
    siso.set_tail(false);
    vec extrinsic_data;
    siso.equalizer(extrinsic_data, rec_sig, apriori_data);
    // ------------------ End of routine -----------------------------
    
    // Convert the IT++ format to Matlab format for output
    output[0] = mxCreateDoubleMatrix(1, extrinsic_data.length(), mxREAL);
    vec2mxArray(extrinsic_data, output[0]);
}
