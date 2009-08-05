/** \file
 * 
 * \brief Matlab interface for SISO MUD module
 * 
 * Matlab function declaration: extrinsic_data = C_SISOmud(rec_sig, apriori_data, ch_imp_response, prec_gen, sigma2, mud_method)
 * 
 * compilation command: mex -litpp -lacml C_SISOmud.cpp
 */

#include "itpp/itmex.h"
#include "/home/bogdan/C++/SISO/SISO.cpp"

using namespace itpp;

void mexFunction(int n_output, mxArray *output[], int n_input, const mxArray *input[])
{
    // Check the number of inputs and output arguments
    if(n_output!=1) mexErrMsgTxt("One output required!");
    if(n_input!=6) mexErrMsgTxt("Six inputs required!");
    
    // Convert input variables to IT++ format
    vec rec_sig = mxArray2vec(input[0]);
    mat apriori_data = mxArray2mat(input[1]);
    mat ch_imp_response = mxArray2mat(input[2]);
    bvec prec_gen = mxArray2bvec(input[3]);
    double sigma2 = mxArray2double(input[4]);
    std::string mud_method = mxArray2string(input[5]);
    
    // ------------------ Start of routine ---------------------------
    tr::SISO siso;   
    siso.set_impulse_response(ch_imp_response);
    siso.set_precoder_generator(prec_gen);
    siso.set_noise(sigma2);    
    siso.set_mud_method(mud_method);
    mat extrinsic_data;
    siso.mud(extrinsic_data, rec_sig, apriori_data);
    // ------------------ End of routine -----------------------------
    
    // Convert the IT++ format to Matlab format for output
    output[0] = mxCreateDoubleMatrix(extrinsic_data.rows(), extrinsic_data.cols(), mxREAL);
    mat2mxArray(extrinsic_data, output[0]);
}
