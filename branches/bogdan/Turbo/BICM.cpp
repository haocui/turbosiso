/** \file
 * 
 * \brief Bit Interleaved Coded Modulation (BICM)
 *
 * Implements BICM using a turbo receiver with a SISO demapper module and a SISO NSC module.
 *
 * Reference: A. Tonello, ''Space-time bit-interleaved coded modulation with an iterative decoding strategy,`` 
 * in Vehicular Technology Conference, vol. 1, pp. 473-478 vol.1, 2000
 */

//#define TO_FILE

#include "itpp/itcomm.h"//for BPSK and BERC classes
#include "SISO.h"//SISO class
#include "Progress_Timer.h"

using namespace itpp;
using tr::SISO;
using tr::threshold;
using std::cout;
using std::endl;
using std::string;

int main(void)
{
    //general parameters
    int const_size = 16;//constellation size
    double threshold_value = 50;
    string map_metric="maxlogMAP";
    ivec gen = "037 021";
    int constraint_length = 5;
    int nb_errors_lim = 1500;
    int nb_bits_lim = int(1e6);
    int perm_len = pow2i(14);//permutation length
    int nb_iter = 10;//number of iterations in the turbo decoder
    vec EbN0_dB = "20";
    double R = 1.0/2.0;//coding rate of FEC
    double Es = 1.0;//mean symbol energy    

    //QAM modulator class
    QAM mod(const_size);

    //other parameters
    string filename = "Res/BICM_"+map_metric+".it";
    int nb_symb = perm_len/mod.bits_per_symbol();
    int nb_bits_tail = perm_len/gen.length();
    int nb_bits = nb_bits_tail-(constraint_length-1);//number of bits in a block (without tail)
    vec sigma2 = (0.5*Es/(R*mod.bits_per_symbol()))*pow(inv_dB(EbN0_dB), -1.0);//N0/2
    int nb_blocks;//number of blocks
    int nb_errors;
    bvec bits(nb_bits);//data bits
    bvec coded_bits(perm_len);//tail is added
    cvec em(perm_len);
    ivec perm(perm_len);
    ivec inv_perm(perm_len);
    cvec rec(perm_len);    
    //SISO demodulator
    vec demod_apriori_data(perm_len);
    vec demod_extrinsic_data;
    //SISO NSC
    vec nsc_intrinsic_coded(perm_len);
    vec nsc_apriori_data(nb_bits_tail);
    nsc_apriori_data.zeros();//always zero
    vec nsc_extrinsic_coded;
    vec nsc_extrinsic_data;
    //decision
    bvec rec_bits(nb_bits_tail);    
    int snr_len = EbN0_dB.length();
    mat ber(nb_iter,snr_len);
    ber.zeros();
    register int en,n;
    
    //CC
    Convolutional_Code nsc;
    nsc.set_generator_polynomials(gen, constraint_length);
    
    //Rayleigh fading
    cvec ch_attenuations(nb_symb);    
    //AWGN
    AWGN_Channel awgn;

    //SISO blocks
    SISO siso;
    siso.set_constellation(mod.bits_per_symbol(), mod.get_symbols(), mod.get_bits2symbols());
    siso.set_generators(gen, constraint_length);
    siso.set_map_metric(map_metric);   
    
    //detector
    BPSK bpsk;
    
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
        awgn.set_noise(2*sigma2(en));
        siso.set_noise(sigma2(en));
        nb_errors = 0;
        nb_blocks = 0;
        while ((nb_errors<nb_errors_lim) && (nb_blocks*nb_bits<nb_bits_lim))//if at the last iteration the nb. of errors is inferior to lim, then process another block
        {
            //permutation
            perm = sort_index(randu(perm_len));
            //inverse permutation
            inv_perm = sort_index(perm);

            //bits generation
            bits = randb(nb_bits);

            //convolutional code
            nsc.encode_tail(bits, coded_bits);//tail is added here to information bits to close the trellis         

            //permutation+QAM modulation
            em = mod.modulate_bits(coded_bits(perm));                       
            
			//flat-fading channel
            ch_attenuations = randn_c(nb_symb);
            rec = awgn(elem_mult(ch_attenuations, em));

            //turbo receiver
            demod_apriori_data.zeros();//a priori information of emitted bits
            siso.set_impulse_response(ch_attenuations);
            for (n=0;n<nb_iter;n++)
            {
                //first decoder
                siso.demapper(demod_extrinsic_data, rec, demod_apriori_data);
               
                //deinterleave+threshold
                nsc_intrinsic_coded = threshold(demod_extrinsic_data(inv_perm), threshold_value);

                //second decoder
                siso.nsc(nsc_extrinsic_coded, nsc_extrinsic_data, nsc_intrinsic_coded, nsc_apriori_data, true);                

                //decision
                rec_bits = bpsk.demodulate_bits(-nsc_extrinsic_data);//suppose that a priori info is zero
                //count errors
                berc.clear();
                berc.count(bits, rec_bits.left(nb_bits));
                ber(n,en) += berc.get_errorrate();

                //interleave
                demod_apriori_data = nsc_extrinsic_coded(perm);
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
    
#ifdef TO_FILE
	//save results to file
    it_file ff(filename);
    ff << Name("BER") << ber;
    ff << Name("EbN0_dB") << EbN0_dB;
    ff << Name("gen") << gen;
    ff << Name("R") << R;
    ff << Name("nb_iter") << nb_iter;
    ff << Name("total_nb_bits") << nb_bits;
    ff << Name("nb_errors_lim") << nb_errors_lim;
    ff << Name("nb_bits_lim") << nb_bits_lim;
    ff.close();
#else
	//show BER
    cout << ber << endl;
#endif

	return 0;
}

