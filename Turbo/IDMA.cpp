/** \file
 *
 * \brief Interleave Division Multiple Access with turbo multiuser detection
 *
 * Reference: L. Liu and L. Ping, ''Iterative detection of chip interleaved CDMA systems in multipath channels,``
 * Electronics letters, vol. 40, pp. 884-886, July 2004
 */

//#define TO_FILE
//#define USE_CC

#include "itpp/itcomm.h"
#include "SISO.h"//SISO class
#include "Progress_Timer.h"

/// Kronecker operator for vectors (both inputs are seen as column or row vectors)
template <class Num_T>
itpp::Vec<Num_T> kron(const itpp::Vec<Num_T> &in, const itpp::Vec<Num_T> &pattern)
{
    int in_len = in.length();
    int pattern_len = pattern.length();
    itpp::Vec<Num_T> out(in_len*pattern_len);
    for (int n=0;n<in_len;n++)
        out.replace_mid(n*pattern_len, in(n)*pattern);
    return out;
}

using namespace itpp;
using tr::SISO;
using tr::threshold;
using std::cout;
using std::endl;
using std::string;

int main(void)
{
    //general parameters
    string mud_method = "maxlogTMAP";
    int nb_usr = 2;
    int spreading_factor = 16;
#ifdef USE_CC
	string map_metric="maxlogMAP";
	ivec gen = "037 021";
    int constraint_length = 5;
    spreading_factor = 8;
#endif
    double threshold_value = 50;
    int ch_nb_taps = 4;//number of channel multipaths
    int nb_errors_lim = 1500;
    int nb_bits_lim = int(1e3);//int(1e6);
    int perm_len = 1024;//38400;//permutation length
    int nb_iter = 15;//number of iterations in the turbo decoder
    vec EbN0_dB = "10";//"0:10:20";
    double Ec = 1.0;//chip energy

#ifdef USE_CC
    int inv_R = spreading_factor*gen.length();
#else
    int inv_R = spreading_factor;
#endif
	double R = 1.0/double(inv_R);//coding rate

    //other parameters
    string filename = "IDMA_"+mud_method+"_"+to_str(nb_usr)+".it";
#ifdef USE_CC
	filename = "cc"+filename;
#endif
	filename = "Res/"+filename;
    int nb_bits = perm_len/inv_R;//number of bits in a block
    vec sigma2 = (0.5*Ec/R)*pow(inv_dB(EbN0_dB), -1.0);//N0/2
    int nb_blocks = 0;//number of blocks
    int nb_errors = 0;//number of errors
    bmat bits(nb_usr,nb_bits);//data bits
#ifdef USE_CC
    bvec coded_bits(nb_bits*gen.length());
    vec mod_bits(nb_bits*gen.length());
#else
    vec mod_bits(nb_bits);
#endif
    vec chips(perm_len);
    imat perm(nb_usr,perm_len);
    imat inv_perm(nb_usr,perm_len);
    vec em(perm_len);
    vec rec(perm_len+ch_nb_taps-1);//padding zeros are added

    //SISO MUD
    mat mud_apriori_data(nb_usr,perm_len);
    mat mud_extrinsic_data;

    //SISO decoder (scrambler or CC)
    vec dec_intrinsic_coded(perm_len);
    vec dec_apriori_data(nb_bits);
    dec_apriori_data.zeros();//always zero
    vec dec_extrinsic_coded;
    vec dec_extrinsic_data;

    //decision
    bvec rec_bits(nb_bits);
    int snr_len = EbN0_dB.length();
    mat ber(nb_iter,snr_len);
    ber.zeros();
    register int en,n,u;

#ifdef USE_CC
    //CC
    Convolutional_Code nsc;
    nsc.set_generator_polynomials(gen, constraint_length);
#endif

	//BPSK
	BPSK bpsk;

    //scrambler pattern
    vec pattern = kron(ones(spreading_factor/2), vec("1.0 -1.0"));

    //AWGN
    AWGN_Channel awgn;
    //multipath channel impulse response (Rayleigh fading) with real coefficients
    vec single_ch(ch_nb_taps);
    mat ch_imp_response(nb_usr, ch_nb_taps);
    vec ini_state = zeros(ch_nb_taps);
    MA_Filter<double,double,double> multipath_channel(ini_state);
    multipath_channel.set_state(ini_state);//initial state is always 0 due to Zero Padding technique
    vec padding_zeros = zeros(ch_nb_taps-1);

    //SISO blocks
    SISO siso;
    siso.set_scrambler_pattern(pattern);
    siso.set_mud_method(mud_method);
#ifdef USE_CC
	siso.set_generators(gen, constraint_length);
	siso.set_map_metric(map_metric);
	siso.set_tail(false);
#endif

	//BER
	BERC berc;

    //progress timer
    tr::Progress_Timer timer;
    timer.set_max(snr_len);

    //Randomize generators
    RNG_randomize();

    //main loop
    timer.progress(0.0);
    for (en=0;en<snr_len;en++)
    {
        awgn.set_noise(sigma2(en));
        siso.set_noise(sigma2(en));
        nb_errors = 0;
        nb_blocks = 0;
        while ((nb_errors<nb_errors_lim) && (nb_blocks*nb_bits<nb_bits_lim))//if at the last iteration the nb. of errors is inferior to lim, then process another block
        {
        	rec.zeros();
        	for (u=0;u<nb_usr;u++)
        	{
	            //permutation
	            perm.set_row(u, sort_index(randu(perm_len)));
	            //inverse permutation
	            inv_perm.set_row(u, sort_index(perm.get_row(u)));

	            //bits generation
	            bits.set_row(u, randb(nb_bits));

#ifdef USE_CC
	            //convolutional code
	            nsc.encode(bits.get_row(u), coded_bits);//no tail

	            //BPSK modulation (1->-1,0->+1)
	            mod_bits = bpsk.modulate_bits(coded_bits);
#else
	            //BPSK modulation (1->-1,0->+1)
	            mod_bits = bpsk.modulate_bits(bits.get_row(u));
#endif

	            //scrambler
	            chips = kron(mod_bits, pattern);

	            //permutation
	            em = chips(perm.get_row(u));

	            //multipath channel
	            single_ch = randray(ch_nb_taps);
	            single_ch /= sqrt(sum_sqr(single_ch));//normalized power profile
	            ch_imp_response.set_row(u, single_ch);
	            multipath_channel.set_coeffs(ch_imp_response.get_row(u));
	            rec += multipath_channel(concat(em, padding_zeros));//Zero Padding
        	}
        	rec = awgn(rec);

            //turbo MUD
            mud_apriori_data.zeros();//a priori LLR of emitted symbols
            siso.set_impulse_response(ch_imp_response);
            for (n=0;n<nb_iter;n++)
            {
            	//MUD
            	siso.mud(mud_extrinsic_data, rec, mud_apriori_data);
            	berc.clear();//mean error rate over all users
            	for (u=0;u<nb_usr;u++)
            	{
	                //deinterleave
	                dec_intrinsic_coded = mud_extrinsic_data.get_row(u)(inv_perm.get_row(u));
#ifdef USE_CC
					//decoder+descrambler
					siso.nsc(dec_extrinsic_coded, dec_extrinsic_data, dec_intrinsic_coded, dec_apriori_data);
#else
	                //descrambler
	                siso.descrambler(dec_extrinsic_coded, dec_extrinsic_data, dec_intrinsic_coded, dec_apriori_data);
#endif
	                //decision
	                rec_bits = bpsk.demodulate_bits(-dec_extrinsic_data);//suppose that a priori info is zero
	                //count errors
	                berc.count(bits.get_row(u), rec_bits);
	                //interleave+threshold
	                mud_apriori_data.set_row(u, threshold(dec_extrinsic_coded(perm.get_row(u)), threshold_value));
            	}
            	ber(n,en) += berc.get_errorrate();
            }//end iterations
            nb_errors += int(berc.get_errors());//get number of errors at the last iteration
            nb_blocks++;
        }//end blocks (while loop)

        //compute BER over all tx blocks
        ber.set_col(en, ber.get_col(en)/nb_blocks);

        //show progress
        timer.progress(1+en);
    }
    timer.toc_print();

    //save results to file
#ifdef TO_FILE
    it_file ff(filename);
    ff << Name("BER") << ber;
    ff << Name("EbN0_dB") << EbN0_dB;
    ff << Name("nb_usr") << nb_usr;
    ff << Name("gen") << spreading_factor;
    ff << Name("nb_iter") << nb_iter;
    ff << Name("total_nb_bits") << nb_bits;
    ff << Name("nb_errors_lim") << nb_errors_lim;
    ff << Name("nb_bits_lim") << nb_bits_lim;
#ifdef USE_CC
	ff << Name("gen") << gen;
#endif
    ff.close();
#else
    //show BER
    cout << ber << endl;
#endif

    return 0;
}
