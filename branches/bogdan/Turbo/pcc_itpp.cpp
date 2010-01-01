/*
 * ------------------------------------------------------------------
 * Parallel Concatenated Convolutional (PCC) codes of coding rate 1/3
 * ------------------------------------------------------------------
 * uses Turbo_Codec class provided by IT++ library
 * this version provides the BER at the last iteration only
 */
 
#define TO_FILE

#include "itpp/itbase.h"
#include "itpp/itcomm.h"//for BPSK and BERC classes
#include "Progress_Timer.h"

using namespace itpp;
using std::cout;
using std::endl;
using std::string;

int main(void)
{
    //general parameters
    string metric="LOGMAX";
    ivec gen = "035 021";//generator polynomials
    int constaint_length = 5;
    int nb_errors_lim = 1500;
    int nb_bits_lim = int(1e6);
    int nb_bits = pow2i(14);//total number of bits in a block (permutation length)
    int nb_iter = 10;//number of iterations in the decoder
    vec EbN0_dB = "0:0.1:2";
    double R = double(1)/double(3);//coding rate (without tail bits)
    double Ec = 1.0;//coded bit energy

    //other parameters
    string filename = "pcc_itpp_"+metric+".it";
    int mem_len = constaint_length-1;
    vec N0 = (Ec/R)*pow(inv_dB(EbN0_dB), -1.0);
    double sigma;
    int nb_blocks;
    int nb_errors;
    bvec bits(nb_bits);
    BPSK bpsk;
    int rec_len = int(1/R)*nb_bits+gen.length()*mem_len+gen.length()*mem_len;//The tailbits from the first encoder are placed before the tailbits from the second encoder
    bvec coded_bits(rec_len);
    vec rec(rec_len);
    bvec rec_bits(nb_bits);
    BERC berc;
    int snr_len = EbN0_dB.length();
    vec ber(snr_len);
    ber.zeros();
    register int en;

    //turbo codec
    Turbo_Codec turbo;
    ivec perm(nb_bits);//tail bits are not permuted
    turbo.set_parameters(gen, gen, constaint_length, perm, nb_iter, metric);

    //main loop
    tr::Progress_Timer timer;//progress timer
    timer.set_max(snr_len);
    timer.progress(0.0);
    for (en=0;en<snr_len;en++)
    {
        nb_errors = 0;
        nb_blocks = 0;
        sigma = sqrt(N0(en)/2);
        turbo.set_awgn_channel_parameters(Ec, N0(en));
        while ((nb_errors<nb_errors_lim) && (nb_blocks*nb_bits<nb_bits_lim))//if at the last iteration the nb. of errors is inferior to lim, then process another block
        {
            //bits generation
            bits = randb(nb_bits);
            
            //permutation
            perm = sort_index(randu(nb_bits));

            //parallel concatenated convolutional code
            turbo.set_interleaver(perm);
            turbo.encode(bits, coded_bits);

            //BPSK modulation (1->-1,0->+1) + channel
            rec = sigma*randn(rec_len)+bpsk.modulate_bits(coded_bits);

            //turbo decoder
            turbo.decode(rec, rec_bits);
            
            //count errors
            berc.clear();
            berc.count(bits, rec_bits);
            ber(en) += berc.get_errorrate();
            nb_errors += (int)berc.get_errors();//get number of errors at the last iteration
            nb_blocks++;
        }//end blocks (while loop)

        //compute BER over all tx blocks
        ber(en) /= nb_blocks;

        //show progress
        timer.progress(1+en);
    }
    timer.toc_print();

#ifdef TO_FILE
    //save results
    it_file ff(filename);
    ff << Name("BER") << ber;
    ff << Name("EbN0_dB") << EbN0_dB;
    ff << Name("gen") << gen;
    ff << Name("R") << R;
    ff << Name("nb_iter") << nb_iter;
    ff << Name("nb_bits") << nb_bits;
    ff << Name("nb_errors_lim") << nb_errors_lim;
    ff << Name("nb_bits_lim") << nb_bits_lim;
    ff.close();
#else
    //show BER
    cout << ber << endl;
#endif    
    
    return 0;
}

