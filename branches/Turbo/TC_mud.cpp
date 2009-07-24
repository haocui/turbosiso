/** \file
 *
 * \brief Transfer Characteristic (TC) for SISO MUD module
 *
 * Reference: L. Liu and L. Ping, ''Iterative detection of chip interleaved CDMA systems in multipath channels,``
 * Electronics letters, vol. 40, pp. 884-886, July 2004
 */

#define TO_FILE
#define USE_PRECODER
#define AWGN_CHANNEL

#include "itpp/itcomm.h"
#include "SISO.h"//SISO class
#include "EXIT.cpp"//EXIT class

using namespace itpp;
using tr::SISO;
using tr::EXIT;
using tr::threshold;
using std::cout;
using std::endl;
using std::string;

int main(void)
{
    //general parameters
    vec sigmaA = "0.01:0.1:7";//standard deviation (sqrt(variance)) of the mutual a priori information
    string mud_method = "maxlogMAP";
    int nb_usr = 8;
    double threshold_value = 50;
#ifndef AWGN_CHANNEL
    int ch_nb_taps = 4;//number of channel multipaths
#else
    int ch_nb_taps = 1;
#endif
    int nb_blocks_lim = 25;
    int perm_len = itpp::pow2i ( 16 );//int(itpp::pow10(5.0));//permutation length
    double EbN0_dB = 8;
    double Ec = 1.0;//chip energy
	double R = 1.0/16.0;//coding rate
#ifdef USE_PRECODER
    ivec prec_gen = "03 02";//octal notation
    int prec_gen_length = 2;
#endif

    //other parameters
    vec sigma2A = sqr(sigmaA);
    int sigma2A_len = sigma2A.length();
    string filename = "TC_mud_"+mud_method+"_"+to_str(nb_usr)+".it";
#ifdef USE_PRECODER
	filename = "prec_"+filename;
#endif
	filename = "Res/" + filename;
    double sigma2 = (0.5*Ec/R)*pow(inv_dB(EbN0_dB), -1.0);//N0/2
    bmat chips(nb_usr,perm_len);
    vec em(perm_len);
    vec rec(perm_len+ch_nb_taps-1);//padding zeros are added
    vec apriori_mutual_info(sigma2A_len);
	vec extrinsic_mutual_info(sigma2A_len);
	register int en,u,nb_blocks;

	//precoder
#ifdef USE_PRECODER
	Rec_Syst_Conv_Code prec;
	prec.set_generator_polynomials(prec_gen, prec_gen_length);
#endif

    //SISO MUD
    mat mud_apriori_data(nb_usr,perm_len);
    mat mud_extrinsic_data;

	//BPSK
	BPSK bpsk;

    //AWGN
    AWGN_Channel awgn;
    awgn.set_noise(sigma2);
    //multipath channel impulse response (Rayleigh fading) with real coefficients
#ifndef AWGN_CHANNEL
    vec single_ch(ch_nb_taps);
    mat ch_imp_response(nb_usr, ch_nb_taps);
    vec ini_state = zeros(ch_nb_taps);
    MA_Filter<double,double,double> multipath_channel(ini_state);
    multipath_channel.set_state(ini_state);//initial state is always 0 due to Zero Padding technique
    vec padding_zeros = zeros(ch_nb_taps-1);
#else
    mat ch_imp_response = ones ( nb_usr, 1 );
#endif

    //SISO blocks
    SISO siso;
    siso.set_mud_method(mud_method);
    siso.set_noise(sigma2);
#ifdef USE_PRECODER
    siso.set_precoder_generator(prec_gen(0), prec_gen_length);
#endif

    //EXIT chart
    EXIT exit;

    //progress timer
    Progress_Timer timer;
    timer.set_max(sigma2A_len);

    //Randomize generators
    RNG_randomize();

    //main loop
    timer.progress(0.0);
    for (en=0;en<sigma2A_len;en++)
    {
    	apriori_mutual_info(en) = exit.apriori_mutual_info(sigma2A(en));//apriori mutual info
        for (nb_blocks=0;nb_blocks<nb_blocks_lim;nb_blocks++)
        {
        	rec.zeros();
        	for (u=0;u<nb_usr;u++)
        	{
	            //chips generation
        		chips.set_row(u, randb(perm_len));

        		//precoder
#ifdef USE_PRECODER
        		bmat parity_bits;
	            prec.encode(chips.get_row(u), parity_bits);
	            em = bpsk.modulate_bits(parity_bits.get_col(0));
#else
	            //BPSK modulation (1->-1,0->+1)
	            em = bpsk.modulate_bits(chips.get_row(u));
#endif

	            //multipath channel
#ifndef AWGN_CHANNEL
	            single_ch = randray(ch_nb_taps);
	            single_ch /= sqrt(sum_sqr(single_ch));//normalized power profile
	            ch_imp_response.set_row(u, single_ch);
	            multipath_channel.set_coeffs(ch_imp_response.get_row(u));
	            rec += multipath_channel(concat(em, padding_zeros));//Zero Padding
#else
	            rec += em;
#endif
        	}
        	rec = awgn(rec);

			//a priori information generation
			for (u=0;u<nb_usr;u++)
			{
				mud_apriori_data.set_row ( u, exit.generate_apriori_info ( chips.get_row ( u ) ) );
			}

            //MUD
            siso.set_impulse_response ( ch_imp_response );
            siso.mud(mud_extrinsic_data, rec, mud_apriori_data);

            //threshold
            mud_extrinsic_data = threshold(mud_extrinsic_data, threshold_value);

            //extrinsic mutual info
            for (u=0;u<nb_usr;u++)
            {
				extrinsic_mutual_info(en) += exit.extrinsic_mutual_info(mud_extrinsic_data.get_row(u), chips.get_row(u));
            }

        }

        //mean extrinsic mutual info over all users and all blocks
        extrinsic_mutual_info(en) /= (nb_usr*nb_blocks_lim);

        //show progress
        timer.progress(1+en);
    }
    timer.toc_print();

    //save results to file
#ifdef TO_FILE
    it_file ff(filename);
    ff << Name("IA") << apriori_mutual_info;
    ff << Name("IE") << extrinsic_mutual_info;
    ff << Name("ch_nb_taps") << ch_nb_taps;
    ff << Name("EbN0_dB") << EbN0_dB;
    ff << Name("perm_len") << perm_len;
    ff << Name("nb_blocks_lim") << nb_blocks_lim;
    ff << Name("R") << R;
    ff.close();
#else
    cout << apriori_mutual_info << endl;
    cout << extrinsic_mutual_info << endl;
#endif

    return 0;
}

