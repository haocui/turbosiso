/** \file
 *
 * \brief Transfer Characteristic (TC) for SISO descrambler module 
 * 
 * Reference: L. Liu and L. Ping, ''Iterative detection of chip interleaved CDMA systems in multipath channels,`` 
 * Electronics letters, vol. 40, pp. 884-886, July 2004
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

/// Kronecker operator for binary vectors
itpp::bvec kron(const itpp::bvec &in, const itpp::bvec &pattern)
{
    int in_len = in.length();
    int pattern_len = pattern.length();
    bvec inv_pattern(pattern_len);
    register int n;
    for (n=0;n<pattern_len;n++)
    	inv_pattern(n) = (itpp::bin(1)^pattern(n));
    itpp::bvec out(in_len*pattern_len);
    for (n=0;n<in_len;n++)
        out.replace_mid(n*pattern_len, (in(n)?inv_pattern:pattern));
    return out;
}

int main(void)
{
    //general parameters
    vec sigmaA = "0.01:0.1:7";//standard deviation (sqrt(variance)) of the mutual a priori information
    int perm_len = int(itpp::pow10(5.0));//permutation length
    int spreading_factor = 16;
    int nb_blocks_lim = 25;   

    //other parameters
    vec sigma2A = sqr(sigmaA);
    int sigma2A_len = sigma2A.length();
    string filename = "Res/TC_descrambler.it";
    int nb_bits = perm_len/spreading_factor;
    bvec bits(nb_bits);    
    bvec pattern = kron(zeros_b(spreading_factor/2), bvec("0 1"));//scrambler pattern
    bvec chips(perm_len);
    vec apriori_mutual_info(sigma2A_len);
	vec extrinsic_mutual_info(sigma2A_len);
	register int en,nb_blocks;
    
    //SISO MUD
    vec dec_apriori_data(nb_bits);
    dec_apriori_data.zeros();//always zero
    vec dec_intrinsic_coded(perm_len);
    vec dec_extrinsic_data;
    vec dec_extrinsic_coded;

    //SISO blocks
    SISO siso;
    siso.set_scrambler_pattern(pattern);
    
    //EXIT chart 
    EXIT exit;

    //progress timer
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
	        
	        //scrambler
	        chips = kron(bits, pattern);	        	            				        	           

			//intrinsic information generation
			dec_intrinsic_coded = exit.generate_apriori_info(chips);

            //SISO descrambler
            siso.descrambler(dec_extrinsic_coded, dec_extrinsic_data, dec_intrinsic_coded, dec_apriori_data);
            
            //extrinsic mutual info
			extrinsic_mutual_info(en) += exit.extrinsic_mutual_info(dec_extrinsic_coded, chips);
        }

        //mean extrinsic mutual info over all users and all blocks
        extrinsic_mutual_info(en) /= nb_blocks_lim;

        //show progress
        timer.progress(1+en);
    }
    timer.toc_print();

    //save results to file
#ifdef TO_FILE
    it_file ff(filename);
    ff << Name("IA") << apriori_mutual_info;
    ff << Name("IE") << extrinsic_mutual_info;
    ff << Name("perm_len") << perm_len;
    ff << Name("nb_blocks_lim") << nb_blocks_lim;
    ff << Name("spreading_factor") << spreading_factor;
    ff.close();
#else
    cout << apriori_mutual_info << endl;
    cout << extrinsic_mutual_info << endl;
#endif

    return 0;
}

