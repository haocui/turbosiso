/*
 * sova_example.cpp
 *
 *  Created on: Oct 13, 2009
 *      Author: bogdan
 */

#include "itpp/itcomm.h"
#include "SISO.h"

using namespace itpp;
using tr::SISO;
using std::cout;
using std::endl;

int main(void)
{
	//Parameters
	ivec gen = "013 015";//octal form
	int constraint_length = 4;
	double EbN0_dB = 1.0;
	double R = 1.0/2.0;
	double Ec = 1;
	int len = 1e4;

    //SISO module
	SISO siso;
	siso.set_generators( gen, constraint_length );
	siso.set_viterbi_win_len( 5*constraint_length );
	siso.set_map_metric("Viterbi");

	//Bits generation & Encoding
    Rec_Syst_Conv_Code cc;
    cc.set_generator_polynomials( gen, constraint_length );
	bvec bits = randb( len );
	bmat parity_bits;
	cc.encode( bits, parity_bits );

	//Multiplexing
	bvec enc( 2*len );
	for ( int n = 0; n < len; n++ )
	{
		enc(2*n) = bits(n);
		enc(2*n+1) = parity_bits(n,0);
	}

	//Modulation & Reception
	BPSK bpsk;
    AWGN_Channel channel;
    double sigma2 = (0.5*Ec/R)*itpp::pow10( -EbN0_dB/10 );
    channel.set_noise( sigma2 );
	vec rec = channel( bpsk.modulate_bits( enc ) );

	//SOVA
	vec apriori_data = zeros( len );
	vec extrinsic_data;
	vec extrinsic_coded;
	tic();
	siso.rsc(extrinsic_coded, extrinsic_data, (-2/sigma2)*rec, apriori_data);
	toc_print();

	//Decision
	BERC berc;
	berc.clear();
	berc.count( bits, bpsk.demodulate_bits( -extrinsic_data ) );
	//berc.count( bits, to_bvec(extrinsic_data) );
	cout << "BER = " << berc.get_errorrate() << endl;
}
