/** \file
 *
 * \brief EXtrinsic Information Transfer (EXIT) chart for Parallel Concatenated Convolutional Codes (PCCCs)
 * of coding rate 1/2
 *
 * Computes the Transfer Characteristics of the two SISO RSC modules used in the turbo decoder for PCCCs
 * at a given Signal to Noise Ratio (SNR).
 * 
 * Reference: S. ten Brink, ''Convergence behavior of iteratively decoded parallel concatenated codes,`` 
 * IEEE Transactions on Communications, vol. 49, pp. 1727-1737, Oct. 2001
 */

#define TO_FILE

#include "itpp/itcomm.h"
#include "SISO.h"//SISO class
#include "EXIT.h"//EXIT class
#include "Progress_Timer.h"

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
    double threshold_value = 50;
    string map_metric="logMAP";
    ivec gen = "023 037";//octal form
    double code_rate = 2.0/3.0;
    bvec puncturing_pattern = "1 0 1 1";//parity bits are punctured
    int constraint_length = 5;
    int nb_blocks_lim = 10;
    int perm_len = int(itpp::pow10(5.0));//total number of bits in a block (with tail)
    double EbN0_dB = 0.8;
    double R = 1.0/2.0;//coding rate of PCCC
    double Ec = 1.0;//coded bit energy
    
    //other parameters
    vec sigma2A = sqr(sigmaA);
    int sigma2A_len = sigma2A.length();
    string filename = "Res/exit_pccc_punctured_"+map_metric+".it";
    int nb_bits = perm_len-(constraint_length-1);//number of bits in a block (without tail)
    double sigma2 = (0.5*Ec/R)*pow(inv_dB(EbN0_dB), -1.0);//N0/2
    double Lc = -2.0/sigma2;//normalisation factor for intrinsic information (take into account the BPSK modulation)
    bvec bits(nb_bits);
    int punctured_bits_len = int(perm_len/code_rate);
    int puncturing_pattern_len = puncturing_pattern.length();
    bvec tail;
    bvec bits_tail(perm_len);
    bmat parity_bits;
    int coded_bits_len = 2*perm_len; 
    bvec coded_bits(coded_bits_len);
    bvec punctured_bits(punctured_bits_len);
    vec mod_bits(punctured_bits_len);
    vec rec_sig(punctured_bits_len);
    vec intrinsic_coded(coded_bits_len);
    vec apriori_data(perm_len);    
    vec extrinsic_coded;
    vec extrinsic_data;
    vec apriori_mutual_info(sigma2A_len);
	vec extrinsic_mutual_info(sigma2A_len);
	extrinsic_mutual_info.zeros();
    register int en,n,k,i,nb_blocks;
    
    //Recursive Systematic Convolutional Code
    Rec_Syst_Conv_Code rsc;
    rsc.set_generator_polynomials(gen, constraint_length);//initial state should be the zero state
    
    //BPSK modulator
    BPSK bpsk;  
    
    //AWGN channel
    AWGN_Channel channel;
    channel.set_noise(sigma2);

    //SISO module
    SISO siso;
    siso.set_generators(gen, constraint_length);
    siso.set_map_metric(map_metric);
    
    //EXIT chart
    EXIT exit;
    
    //Progress timer
    tr::Progress_Timer timer;
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
            //bits generation
            bits = randb(nb_bits);

            //parallel concatenated convolutional code
            rsc.encode_tail(bits, tail, parity_bits);//tail is added
            
            //form coder output
            bits_tail = concat(bits, tail);
            for (n=0;n<perm_len;n++)
            {
            	coded_bits(2*n) = bits_tail(n);//systematic output
            	coded_bits(2*n+1) = parity_bits(n,0);//parity output
            }
            
            //apply puncturing pattern
            i = 0;
            for (n=0;n<coded_bits_len/puncturing_pattern_len;n++)
            {
            	for (k=0;k<puncturing_pattern_len;k++)
            	{
            		if (puncturing_pattern(k))
            		{
            			punctured_bits(i) = coded_bits(n*puncturing_pattern_len+k);
            			i++;
            		}
            	}
            }            	

            //BPSK modulation (1->-1,0->+1)
            mod_bits = bpsk.modulate_bits(punctured_bits);
            
            //AWGN channel
            rec_sig = channel(mod_bits);
            
            //intrinsic info. of coded bits taking into account the puncturing pattern
            i = 0;
            for (n=0;n<coded_bits_len/puncturing_pattern_len;n++)
            {
            	for (k=0;k<puncturing_pattern_len;k++)
            	{
            		if (puncturing_pattern(k))
            		{
            			intrinsic_coded(n*puncturing_pattern_len+k) = Lc*rec_sig(i);
            			i++;
            		}
            		else
            			intrinsic_coded(n*puncturing_pattern_len+k) = 0.0;
            	}
            }
            
            //a priori info. of data bits
            apriori_data = exit.generate_apriori_info(bits_tail);

            //SISO RSC module
            siso.rsc(extrinsic_coded, extrinsic_data, intrinsic_coded, apriori_data, true);
            
            //extract from extrinsic info. of data bits the intrinsic info. corresponding to systematic bits (ten Brink approach is used)
            for (n=0;n<perm_len;n++)
            	extrinsic_data(n) -= intrinsic_coded(2*n);
            
            //threshold
            extrinsic_data = threshold(extrinsic_data, threshold_value);

			//extrinsic mutual info
			extrinsic_mutual_info(en) += exit.extrinsic_mutual_info(extrinsic_data.left(nb_bits), bits);
        }//end blocks (while loop)
        
        //mean extrinsic mutual info over all blocks
        extrinsic_mutual_info(en) /= nb_blocks_lim;

        //show progress
        timer.progress(1+en);
    }
    timer.toc_print();

#ifdef TO_FILE
    //save results to file
    it_file ff(filename);
    ff << Name("IA") << apriori_mutual_info;
    ff << Name("IE") << extrinsic_mutual_info;
    ff << Name("EbN0_dB") << EbN0_dB;
    ff << Name("gen") << gen;
    ff << Name("R") << R;
    ff << Name("perm_len") << perm_len;
    ff << Name("nb_blocks_lim") << nb_blocks_lim;
    ff.close();
#else
    //show BER    
    cout << apriori_mutual_info << endl;
    cout << extrinsic_mutual_info << endl;
#endif    
    
    return 0;
}
