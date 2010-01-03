/** \file
 * 
 * \brief Matlab interface for SISO demapper module
 * 
 * Matlab function declaration: extrinsic_data = C_SISOdemapper(rec_sig, apriori_data, ch_imp_response, sigma2, bits_per_symbol, 
 * 														constellation, bin_constellation, demapper_method)
 * 				or
 * 								extrinsic_data = C_SISOdemapper(rec_sig, apriori_data, ch_imp_response, sigma2, bits_per_symbol, 
 * 														constellation, bin_constellation, symbols_per_block, ST_gen1, ST_gen2, demapper_method)
 * 
 * compilation command: mex -litpp -lacml C_SISOdemapper.cpp
 */

#include "itpp/itmex.h"
#include "SISO.h"

using namespace itpp;

void mexFunction(int n_output, mxArray *output[], int n_input, const mxArray *input[])
{
    // Check the number of inputs and output arguments
    if(n_output!=1) mexErrMsgTxt("One output required!");
    if((n_input!=8)&&(n_input!=11)) mexErrMsgTxt("Eigth or eleven inputs required!");
        
	tr::SISO siso;
	vec extrinsic_data;
	if (n_input==8)
	{
		// Convert input variables to IT++ format
    	cvec rec_sig = mxArray2cvec(input[0]);
    	vec apriori_data = mxArray2vec(input[1]);
    	cvec ch_imp_response = mxArray2cvec(input[2]);
    	double sigma2 = mxArray2double(input[3]);    
    	int bits_per_symbol = mxArray2int(input[4]);
    	cvec constellation = mxArray2cvec(input[5]);
    	bmat bin_constellation = mxArray2bmat(input[6]);
    	std::string demapper_method = mxArray2string(input[7]);
    
    	// ------------------ Start of routine ---------------------------
    	siso.set_impulse_response(ch_imp_response);
    	siso.set_noise(sigma2);    
    	siso.set_constellation(bits_per_symbol, constellation, bin_constellation);
    	siso.set_demapper_method(demapper_method);
    	siso.demapper(extrinsic_data, rec_sig, apriori_data);
    	// ------------------ End of routine -----------------------------
	}
	else if (n_input==11)
	{
		// Convert input variables to IT++ format
    	cmat rec_sig = mxArray2cmat(input[0]);
    	vec apriori_data = mxArray2vec(input[1]);
    	cmat ch_imp_response = mxArray2cmat(input[2]);
    	double sigma2 = mxArray2double(input[3]);    
    	int bits_per_symbol = mxArray2int(input[4]);
    	cvec constellation = mxArray2cvec(input[5]);
    	bmat bin_constellation = mxArray2bmat(input[6]);
    	int symbols_per_block = mxArray2int(input[7]);
    	cmat A = mxArray2cmat(input[8]);
    	cmat B = mxArray2cmat(input[9]);
    	std::string demapper_method = mxArray2string(input[10]);
    
    	// ------------------ Start of routine ---------------------------
    	siso.set_impulse_response(ch_imp_response);
    	siso.set_noise(sigma2);    
    	siso.set_constellation(bits_per_symbol, constellation, bin_constellation);
    	siso.set_st_block_code(symbols_per_block, A, B, rec_sig.cols());
    	siso.set_demapper_method(demapper_method);
    	siso.demapper(extrinsic_data, rec_sig, apriori_data);
    	// ------------------ End of routine -----------------------------
	}
    
    // Convert the IT++ format to Matlab format for output
    output[0] = mxCreateDoubleMatrix(1, extrinsic_data.length(), mxREAL);
    vec2mxArray(extrinsic_data, output[0]);
}
