/** \file
 * 
 * \brief Soft Input Soft Output (SISO) modules class
 */

#ifndef SISO_CLASS
#define SISO_CLASS

#include <cmath> //use system INFINITY
#include <iostream>
#include "itpp/itbase.h" //IT++ base module

/// turbo receivers workspace
namespace tr
{

/// Soft Input Soft Output (%SISO) modules

/** The following SISO modules are implemented:
 *- decoder for an \f$1/2\f$ Recursive Systematic Convolutional (RSC) code
 *- decoder for an \f$1/r\f$ Non-recursive non-Systematic Convolutional (NSC) code
 *- equalizer (without and with precoding)
 *- descrambler used in Interleave Division Multiple Access (IDMA) systems
 *- Multi User Detectors (MUDs) for IDMA systems
 *- demappers for Bit Interleaved Coded Modulation (BICM) systems
 *- demappers for Space Time (ST) BICM systems
 * 
 * BPSK mapping is realized as follows: 0 -> 1 and 1 -> -1. Thus the xor truth table is preserved when multiplying BPSK symbols.
 */
class SISO
{
public:
    /// %SISO class constructor
    /** Internal variables are initialized with default values:
     * - the trellis used in MAP algorithm is not terminated
     * - maxlogMAP metric is used
     * - simplified GCD is selected
     * - no scrambler is used
     * - GA demapper is selected
     */
    SISO()
    {
        tail = false;
        MAP_metric = "maxlogMAP";
        MUD_method = "sGCD";
        scrambler_pattern = "0";//corresponds to +1 using BPSK mapping
        prec_gen = "1";
        demapper_method = "GA";
    };
    //MAP algorithm setup functions
    /// Sets the metric for MAP algorithm (convolutional codes and multipath channels)
    /** Possible input values are:
     * - logMAP
     * - maxlogMAP
     */
    void set_map_metric(const std::string &in_MAP_metric)
    {
        MAP_metric = in_MAP_metric;
    };
    /// Sets the precoder generator polynomial for turbo equalizer
    /** The precoder is used in turbo equalization.
     * The generator polynomial describes the feedback connections of the precoder.
     */
    void set_precoder_generator(const itpp::bvec &in_prec_gen)//set precoder polynomial
    {
        prec_gen = in_prec_gen;
    };
    void set_precoder_generator(const int &in_prec_gen, const int &constraint_length)//set precoder polynomial
    {
        prec_gen = itpp::dec2bin(constraint_length, in_prec_gen);
    };
    /// Sets convolutional code generator polynomials
    /**
     * The generator polynomials are specified as rows of the input binary matrix.
     */
    void set_generators(const itpp::bmat &in_gen)
    {
        gen = in_gen;
    };
    void set_generators(const itpp::ivec &in_gen, const int &constraint_length)
    {
    	int nb_outputs = in_gen.length();
    	gen.set_size(nb_outputs, constraint_length);
        for (int n=0;n<nb_outputs;n++)
            gen.set_row(n, itpp::dec2bin(constraint_length, in_gen(n)));
    };
    /// Signals whether the trellis used in the MAP algorithm is terminated or not (only for convolutional codes and multipath channels)
    /**
     * If the input value is true, the trellis is terminated. In order to terminate the trellis an ending bit tail is added to the input stream of the encoder.
     */
    void set_tail(const bool &in_tail)
    {
        tail = in_tail;
    };
    //channel setup functions
    /// Sets Additive White Gaussian Noise variance for each dimension
    void set_noise(const double &in_sigma2)
    {
        sigma2 = in_sigma2;
    };
    /// Sets channel impulse response for equalizer
    /** The input is a real channel impulse response.
     */
    void set_impulse_response(const itpp::vec &h)
    {
        impulse_response.set_size(1, h.length());
        impulse_response.set_row(0, h);
    };
    /// Sets channel impulse response for Multi-User Detectors
    /** The input is a real matrix with each row represented by the channel impulse response of one user.
     */
    void set_impulse_response(const itpp::mat &H)
    {
        impulse_response = H;
    };
    /// Sets the channel attenuations for demappers (when only a modulator is used)
    /** The input is the vector of channel attenuations for each received symbol
     */
    void set_impulse_response(const itpp::cvec &h)
    {
        c_impulse_response.set_size(1, h.length());
        c_impulse_response.set_row(0, h);
    };
    /// Sets channel attenuations for demappers (when Space-Time codes are used)
    /** The input is a complex matrix of dimension \f$MN\times tx\_duration\f$, where \f$M\f$ and \f$N\f$ is the number of emission and reception antennas, respectively and \f$tx\_duration\f$ is the transmission duration expressed in symbol durations.
     *
     * This input matrix is formed as follows: 
     * - the starting point is an \f$M\times N\f$ complex matrix (channel matrix) with one row represented by the atenuations seen by each of \f$N\f$ reception antennas during one symbol duration, when the signal is emitted by a given emission antenna. The row number represents the number of emission antenna.
     * - the \f$M\times N\f$ channel matrix is then transformed into a vector of length \f$MN\f$, with the first \f$M\f$ elements the first column of the channel matrix
     * - the vector of length \f$MN\f$ represents one column of the input matrix
     * - in the input matrix, the vector is repeated \f$\tau_c\f$ times and \f$tx\_duration/\tau_c\f$ different vectors are used. Thus, the channel is supposed constant over \f$\tau_c\f$ symbol durations (channel coherence time) and \f$tx\_duration/\tau_c\f$ different channel realisations are used.
     * - in our implementation \f$\tau_c\f$ must be and integer multiple of \f$T\f$, where \f$T\f$ is the ST block code duration expressed in symbol durations. This means that the channel matrix must be constant over at least \f$T\f$ symbol durations.
     */
    void set_impulse_response(const itpp::cmat &cH)
    {
        c_impulse_response = cH;
    };
    /// Sets scrambler pattern
    /** The scrambler pattern must be a sequence of \f$\pm1\f$ and is used by the %SISO NSC module in IDMA systems reception. At emission side the bits are first encoded by an NSC code, BPSK modulated and then scrambled with the given pattern.
     */
    void set_scrambler_pattern(const itpp::vec &phi)
    {
    	int phi_len = phi.length();
    	scrambler_pattern.set_size(phi_len);
        //scrambler_pattern = to_bvec((1-phi)/2);//BPSK mapping: 0->+1 and 1->-1
        register int n;
        for (n=0;n<phi_len;n++)
        	scrambler_pattern(n) = itpp::bin((1-int(phi(n)))/2);//BPSK mapping: 0->+1 and 1->-1
    };
    void set_scrambler_pattern(const itpp::bvec &phi)
    {
        scrambler_pattern = phi;
    };
    /// Sets Multi-User Detector method
    /** Possible input values are:
     * - GCD
     * - sGCD (simplified GCD)
     */
    void set_mud_method(const std::string &method)
    {
        MUD_method = method;
    };
    //demodulator and MIMO demapper setup
    /// Sets symbol constellation
    /** The correspondence between each symbol and its binary representation is given.
     */
    void set_constellation(const int &in_nb_bits_symb, ///< the number of symbols
                           const itpp::cvec &in_constellation, ///< all possible symbols as a complex vector
                           const itpp::bmat &in_bin_constellation ///< binary representations of symbols as a binary matrix (each row corresponds to one symbol)
                          )
    {
        nb_bits_symb = in_nb_bits_symb;
        constellation = in_constellation;
        bin_constellation = in_bin_constellation;
    };
    void set_constellation(const int &in_nb_bits_symb, ///< the number of symbols
                           const itpp::cvec &in_constellation, ///< all possible symbols as a complex vector
                           const itpp::ivec &in_int_constellation ///< integer representations of symbols as a vector
                          )
    {
        nb_bits_symb = in_nb_bits_symb;        
        int nb_symb = in_constellation.length();
        constellation.set_size(nb_symb);
        bin_constellation.set_size(nb_symb, nb_bits_symb);
        for (int n=0;n<nb_symb;n++)
        {
        	constellation(n) = in_constellation(n);
        	bin_constellation.set_row(n, itpp::dec2bin(nb_bits_symb, in_int_constellation(n)));
        }
    };
    /// Sets Space-Time block code parameters
    /** ST block codes are generated using Hassibi's model.
     */
    void set_st_block_code(const int &Q, ///< the number of symbols per block
                           const itpp::cmat &A, ///< generator matrices
                           const itpp::cmat &B, ///< of the ST block code
                           const int &N ///< the number of reception antennas
                          )
    {
        symbols_block = Q;
        nb_em_ant = A.cols();
        nb_rec_ant = N;
        block_duration = A.rows()/Q;
        ST_gen1 = A;
        ST_gen2 = B;
    };
    /// Sets demapper method
    /** Possible input values are:
     * - Hassibi_MAP (maxlogMAP algorithm applied for ST block codes represented using Hassibi's model)
     * - GA
     * - sGA (simplified GA)
     * - mmsePIC
     * - zfPIC (simplified mmsePIC)
     * - Alamouti_MAP (maxlogMAP algorithm applied to Alamouti code using matched-filter reception method)
     */
    void set_demapper_method(const std::string &method)
    {
        demapper_method = method;
    };
    /// %SISO decoder for RSC codes
    void rsc(itpp::vec &extrinsic_parity, ///< extrinsic information of parity bits
             itpp::vec &extrinsic_data, ///< extrinsic information of data bits
             const itpp::vec &intrinsic_coded, ///< intrinsic information of coded bits
             const itpp::vec &apriori_data ///< a priori information of data bits
            );
    /// %SISO decoder for RSC codes (tail is set through input)
    void rsc(itpp::vec &extrinsic_parity, ///< extrinsic information of parity bits
             itpp::vec &extrinsic_data, ///< extrinsic information of data bits
             const itpp::vec &intrinsic_coded, ///< intrinsic information of coded bits
             const itpp::vec &apriori_data, ///< a priori information of data bits
             const bool &tail ///< if true the trellis is terminated
            )
    {
        set_tail(tail);
        rsc(extrinsic_parity, extrinsic_data, intrinsic_coded, apriori_data);
    };
    /// %SISO decoder for NSC codes
    void nsc(itpp::vec &extrinsic_coded, ///< extrinsic information of coded bits
    		 itpp::vec &extrinsic_data, ///< extrinsic information of data bits
    		 const itpp::vec &intrinsic_coded, ///< intrinsic information of coded bits 
             const itpp::vec &apriori_data ///< a priori information of data bits
            );
    /// %SISO decoder for NSC codes (tail is set through input)
    void nsc(itpp::vec &extrinsic_coded, ///< extrinsic information of coded bits
    		 itpp::vec &extrinsic_data, ///< extrinsic information of data bits
    		 const itpp::vec &intrinsic_coded, ///< intrinsic information of coded bits 
             const itpp::vec &apriori_data, ///< a priori information of data bits
             const bool &tail ///< if true the trellis is terminated
            )
    {
    	set_tail(tail);
    	nsc(extrinsic_coded, extrinsic_data, intrinsic_coded, apriori_data);
    };
    /// %SISO equalizer
    /** Channel trellis is generated so that BPSK mapping is assumed: 0->+1 and 1->-1 (xor truth table is preserved)
     */
    void equalizer(itpp::vec &extrinsic_data, ///< extrinsic informations of input symbols
                   const itpp::vec &rec_sig, ///< received signal
                   const itpp::vec &apriori_data ///< a priori informations of input symbols
                  );
    /// %SISO equalizer (tail is set through input)
    /** Channel trellis is generated so that BPSK mapping is assumed: 0->+1 and 1->-1 (xor truth table is preserved)
     */
    void equalizer(itpp::vec &extrinsic_data, ///< extrinsic informations of input symbols
                   const itpp::vec &rec_sig, ///< received signal
                   const itpp::vec &apriori_data,  ///< a priori informations of input symbols
                   const bool &tail ///< if true the trellis is terminated
                  )
    {
    	set_tail(tail);
    	equalizer(extrinsic_data, rec_sig, apriori_data);
    };
    /// %SISO descrambler
    void descrambler(itpp::vec &extrinsic_coded, ///< extrinsic information of scrambled bits
    			   itpp::vec &extrinsic_data, ///< extrinsic information of informational bits
                   const itpp::vec &intrinsic_coded, ///< intrinsic information of scrambled bits    			   
                   const itpp::vec &apriori_data ///< a priori information of informational bits
                  );
    /// %SISO Multi-User Detector
    void mud(itpp::mat &extrinsic_data, ///< extrinsic informations of emitted chips from all users
             const itpp::vec &rec_sig, ///< received signal
             const itpp::mat &apriori_data ///< a priori informations of emitted chips from all users
            );
    /// %SISO demapper (when only a modulator is used)
    void demapper(itpp::vec &extrinsic_data, ///< extrinsic informations of emitted bits
                  const itpp::cvec &rec_sig, ///< received signal
                  const itpp::vec &apriori_data ///< a priori informations of emitted bits
                 );          
    /// %SISO demapper (when Space-Time codes are used)
    void demapper(itpp::vec &extrinsic_data, ///< extrinsic informations of emitted bits
                  const itpp::cmat &rec_sig, ///< received signal
                  const itpp::vec &apriori_data ///< a priori informations of emitted bits
                 );
private:
    /// SISO::rsc using logMAP algorithm
    void rsc_logMAP(itpp::vec &extrinsic_parity, itpp::vec &extrinsic_data, const itpp::vec &intrinsic_coded, const itpp::vec &apriori_data);
    /// SISO::rsc using maxlogMAP algorithm
    void rsc_maxlogMAP(itpp::vec &extrinsic_parity, itpp::vec &extrinsic_data, const itpp::vec &intrinsic_coded, const itpp::vec &apriori_data);
    /// SISO::nsc using logMAP algorithm
    void nsc_logMAP(itpp::vec &extrinsic_coded, itpp::vec &extrinsic_data, const itpp::vec &intrinsic_coded, const itpp::vec &apriori_data);
    /// SISO::nsc using maxlogMAP algorithm
    void nsc_maxlogMAP(itpp::vec &extrinsic_coded, itpp::vec &extrinsic_data, const itpp::vec &intrinsic_coded, const itpp::vec &apriori_data);
    /// SISO::equalizer using logMAP algorithm
    void equalizer_logMAP(itpp::vec &extrinsic_data, const itpp::vec &rec_sig, const itpp::vec &apriori_data);
    /// SISO::equalizer using maxlogMAP algorithm
    void equalizer_maxlogMAP(itpp::vec &extrinsic_data, const itpp::vec &rec_sig, const itpp::vec &apriori_data);
    /// SISO::mud using maxlogMAP algorithm
    void mud_maxlogMAP(itpp::mat &extrinsic_data, const itpp::vec &rec_sig, const itpp::mat &apriori_data);
    /// SISO::mud using Gaussian Chip Detector (GCD)
    void GCD(itpp::mat &extrinsic_data, const itpp::vec &rec_sig, const itpp::mat &apriori_data);
    /// SISO::mud using simplified Gaussian Chip Detector (sGCD)
    void sGCD(itpp::mat &extrinsic_data, const itpp::vec &rec_sig, const itpp::mat &apriori_data);
    /// SISO::demapper using maxlogMAP algorithm for ST block codes described using Hassibi's model
    void Hassibi_maxlogMAP(itpp::vec &extrinsic_data, const itpp::cmat &rec_sig, const itpp::vec &apriori_data);
    /// SISO::demapper using Gaussian Approximation (GA) algorithm
    void GA(itpp::vec &extrinsic_data, const itpp::cmat &rec_sig, const itpp::vec &apriori_data);
    /// SISO::demapper using simplified Gaussian Approximation (sGA) algorithm
    void sGA(itpp::vec &extrinsic_data, const itpp::cmat &rec_sig, const itpp::vec &apriori_data);
    /// SISO::demapper using MMSE Parallel Interference Canceller (PIC)
    void mmsePIC(itpp::vec &extrinsic_data, const itpp::cmat &rec_sig, const itpp::vec &apriori_data);
    /// SISO::demapper using ZF Parallel Interference Canceller (PIC)
    void zfPIC(itpp::vec &extrinsic_data, const itpp::cmat &rec_sig, const itpp::vec &apriori_data);
    /// SISO::demapper using maxlogMAP algorithm and matched filter receiver for Alamouti ST code
    void Alamouti_maxlogMAP(itpp::vec &extrinsic_data, const itpp::cmat &rec_sig, const itpp::vec &apriori_data);
    /// SISO::demapper using logMAP algorithm for complex modulators
    void demodulator_logMAP(itpp::vec &extrinsic_data, const itpp::cvec &rec_sig, const itpp::vec &apriori_data);
    /// SISO::demapper using maxlogMAP algorithm for complex modulators
    void demodulator_maxlogMAP(itpp::vec &extrinsic_data, const itpp::cvec &rec_sig, const itpp::vec &apriori_data);
    /// Prints an error message to standard output
    /** If the %SISO class is used in a mex file, this function ensures that the proper function is used for displaying the error message
     */
    void print_err_msg(const std::string &msg) const
    {
#ifdef mex_h
		mexErrMsgTxt(msg.c_str());
#else
		std::cout << msg << std::endl;
#endif    	
    };

    // MAP algorithm variables
    /// MAP algorithm metric
    std::string MAP_metric;
    /// Generator polynomials for convolutional codes (CC)
    itpp::bmat gen;
    /// Precoder generator polynomial
    itpp::bvec prec_gen;
    /// True if trellis of CC is terminated
    bool tail;
    //channel variables
    /// AWGN noise variance
    double sigma2;
    /// Real channel impulse response
    itpp::mat impulse_response;
    /// Complex channel impulse response
    itpp::cmat c_impulse_response;
    /// Scrambler pattern
    itpp::bvec scrambler_pattern;
    /// MUD method
    std::string MUD_method;
    //constellation variables
    /// Number of bits/symbol
    int nb_bits_symb;
    /// Complex constellation
    itpp::cvec constellation;
    /// Binary constellation
    itpp::bmat bin_constellation;
    //Space Time block code variables
    /// Number of symbols/block
    int symbols_block;
    /// Number of emission antennas
    int nb_em_ant;
    /// Number of reception antennas
    int nb_rec_ant;
    /// ST code block duration
    int block_duration;
    /// ST generator matrix 1
    itpp::cmat ST_gen1;
    /// ST generator matrix 2
    itpp::cmat ST_gen2;
    /// Demapper method
    std::string demapper_method;

    //internal variables and functions
    /// FIR filter for a zero padded signals (\f$L\f$ zeros are added at the end of the signal, where \f$L\f$ is the order of the filter)
    void zpFIRfilter(itpp::vec& filt, ///< filtered signal
                     const itpp::vec &h, ///< filter impulse response
                     const itpp::vec &sig ///< signal to filter
                    );
    /// Generates (precoded) channel trellis
    void gen_chtrellis(void);
    /// Generates (precoded) hyper channel trellis
    void gen_hyperTrellis(void);
    /// (Hyper) Channel trellis
    struct
    {
    	int numInputSymbols;///< number of input symbols
        int stateNb;///< number of states
        int* prevState;///< previous states
        int* nextState;///< next states
        double* output;///< output
        int* input;///< input
    } chtrellis;
    /// Generates Recursive and Systematic Convolutional (RSC) code trellis
    void gen_rsctrellis(void);
    /// RSC code trellis
    struct
    {
        int numStates;///< number of states
        int* prevStates;///< previous states
        int* nextStates;///< next states
        double* PARout;///< parity output bit
        itpp::bin* fm;/// feedback memory
    } rsctrellis;
    /// Generates Non recursive and non Systematic Convolutional (NSC) code trellis
    void gen_nsctrellis(void);
    /// NSC code trellis
    struct
    {
        int stateNb;///< number of states
        int* prevState;///< previous states
        int* nextState;///< next states
        double* output;///< output
        int* input;///< input
    } nsctrellis;
    /// Finds half constellations
    void find_half_const(int &select_half, itpp::vec &re_part, itpp::bmat &re_bin_part, itpp::vec &im_part, itpp::bmat &im_bin_part);
    /// Finds equivalent received signal with real coefficients
    void EquivRecSig(itpp::vec &x_eq, const itpp::cmat &rec_sig);
    /// Finds equivalent channel with real coefficients
    void EquivCh(itpp::mat &H_eq, const itpp::cvec &H);
};

void SISO::gen_rsctrellis(void)
//generates 1/2 RSC trellis structure for binary symbols
//the states are numbered from 0
{
    int mem_len = gen.cols()-1;
    register int n,k,j;
    itpp::bin feedback,out;
    int buffer;

    rsctrellis.numStates = (1<<mem_len);
    rsctrellis.prevStates = new int[2*rsctrellis.numStates];
    rsctrellis.nextStates = new int[2*rsctrellis.numStates];
    rsctrellis.PARout = new double[2*rsctrellis.numStates];
    rsctrellis.fm = new itpp::bin[rsctrellis.numStates];

    itpp::bvec cases(mem_len);
    for (n=0;n<2;n++)
    {
        for (k=0;k<rsctrellis.numStates;k++)
        {
            cases = itpp::dec2bin(mem_len, k);
            //feedback
            feedback = (itpp::bin)n;
            for (j=1;j<(mem_len+1);j++)
                feedback ^= (gen(0,j)*cases[j-1]);
            //out
            out = feedback*gen(1,0);
            for (j=1;j<(mem_len+1);j++)
                out ^= (gen(1,j)*cases[j-1]);
            rsctrellis.PARout[k+n*rsctrellis.numStates] = (out?1.0:0.0);//parity bit
            rsctrellis.fm[k] = n^out;
            //shift
            for (j=mem_len-1;j>0;j--)
                cases[j] = cases[j-1];
            cases[0] = feedback;
            //next and previous state
            buffer = itpp::bin2dec(cases, true);
            rsctrellis.nextStates[k+n*rsctrellis.numStates] = buffer;//next state
            rsctrellis.prevStates[buffer+n*rsctrellis.numStates] = k;//previous state
        }
    }
}

inline void SISO::rsc(itpp::vec &extrinsic_parity, itpp::vec &extrinsic_data, const itpp::vec &intrinsic_coded, const itpp::vec &apriori_data)
{
	if (gen.size()==0)
	{
		print_err_msg("SISO::rsc: generator polynomials not initialized");
		return;
	} 
	
	if (MAP_metric=="logMAP")
		rsc_logMAP(extrinsic_parity, extrinsic_data, intrinsic_coded, apriori_data);
	else if (MAP_metric=="maxlogMAP")
		rsc_maxlogMAP(extrinsic_parity, extrinsic_data, intrinsic_coded, apriori_data);
	else
		print_err_msg("SISO::rsc: unknown MAP metric. The MAP metric should be either logMAP or maxlogMAP");	
}

void SISO::rsc_logMAP(itpp::vec &extrinsic_parity, itpp::vec &extrinsic_data, const itpp::vec &intrinsic_coded, const itpp::vec &apriori_data)
/*
 logMAP (SISO) decoder for RSC of rate 1/2
 extrinsic_parity - extrinsic information of parity bits
 extrinsic_data - extrinsic information of data (informational) bits
 intrinsic_coded - intrinsic information of coded (systematic and parity) bits
 apriori_data - a priori information of data (informational) bits
 Reference: Steven S. Pietrobon and Adrian S. Barbulescu, "A simplification of the modified Bahl decoding algorithm for systematic convolutional codes", Proc. ISITA, 1994
 */
{
    //get parameters
    int N = apriori_data.length();
    //other parameters
    register int n,k;
    double buffer;
    int index;
    double A_min, A_max;
    double sum0, sum1;

    //trellis generation
    gen_rsctrellis();

    //parameter initialization
    double* Lc1I = new double[N];
    double* Lc2I = new double[N];
    for (n=0;n<N;n++)
    {
        Lc1I[n] = intrinsic_coded[2*n];
        Lc2I[n] = intrinsic_coded[2*n+1];
    }
    double* A0 = new double[rsctrellis.numStates*N];
    double* A1 = new double[rsctrellis.numStates*N];
    double* A_mid = new double[N];
    double* B0 = new double[rsctrellis.numStates*N];
    double* B1 = new double[rsctrellis.numStates*N];
    buffer = (tail?-INFINITY:0);//log(buffer)
    for (n=0;n<N;n++)
    {
        for (k=0;k<rsctrellis.numStates;k++)
        {
            A0[k+n*rsctrellis.numStates] = -INFINITY;
            A1[k+n*rsctrellis.numStates] = -INFINITY;
            B0[k+n*rsctrellis.numStates] = buffer;
            B1[k+n*rsctrellis.numStates] = buffer;
        }
        A_mid[n] = 0;
    }

    //A
    A0[0] = Lc2I[0]*rsctrellis.PARout[0];//i=0
    A1[0] = Lc1I[0]+apriori_data[0]+Lc2I[0]*rsctrellis.PARout[rsctrellis.numStates];//i=1
    for (n=1;n<N;n++)
    {
        A_min = INFINITY;
        A_max = 0;
        for (k=0;k<rsctrellis.numStates;k++)
        {
            buffer = itpp::log_add(A0[rsctrellis.prevStates[k]+(n-1)*rsctrellis.numStates], 
                     A1[rsctrellis.prevStates[k+rsctrellis.numStates]+(n-1)*rsctrellis.numStates]);//log(alpha0+alpha1)
            A0[k+rsctrellis.numStates*n] = Lc2I[n]*rsctrellis.PARout[k]+buffer;
            A1[k+rsctrellis.numStates*n] = Lc1I[n]+apriori_data[n]+Lc2I[n]*rsctrellis.PARout[k+rsctrellis.numStates]+buffer;
            //find min A(:,n)
            A_min = std::min(A_min, A0[k+rsctrellis.numStates*n]);
            //find max A(:,n)
            A_max = std::max(A_max, A0[k+rsctrellis.numStates*n]);
        }
        //normalization
        A_mid[n] = (A_min+A_max)/2;
        if (isinf(A_mid[n]))
            continue;
        for (k=0;k<rsctrellis.numStates;k++)
        {
            A0[k+rsctrellis.numStates*n] -= A_mid[n];
            A1[k+rsctrellis.numStates*n] -= A_mid[n];
        }
    }
    //B
    B0[rsctrellis.prevStates[0]+(N-1)*rsctrellis.numStates] = 0;
    B1[rsctrellis.prevStates[rsctrellis.numStates]+(N-1)*rsctrellis.numStates] = 0;
    for (n=N-2;n>=0;n--)
    {
        for (k=0;k<rsctrellis.numStates;k++)
        {
            index = rsctrellis.nextStates[k];
            B0[k+rsctrellis.numStates*n] = itpp::log_add(B0[index+(n+1)*rsctrellis.numStates]+Lc2I[n+1]*rsctrellis.PARout[index], 
                                              B1[index+(n+1)*rsctrellis.numStates]+Lc1I[n+1]+apriori_data[n+1]+Lc2I[n+1]*rsctrellis.PARout[index+rsctrellis.numStates]);
            index = rsctrellis.nextStates[k+rsctrellis.numStates];
            B1[k+rsctrellis.numStates*n] = itpp::log_add(B0[index+(n+1)*rsctrellis.numStates]+Lc2I[n+1]*rsctrellis.PARout[index], 
                                              B1[index+(n+1)*rsctrellis.numStates]+Lc1I[n+1]+apriori_data[n+1]+Lc2I[n+1]*rsctrellis.PARout[index+rsctrellis.numStates]);

        }
        if (isinf(A_mid[n+1]))
            continue;
        for (k=0;k<rsctrellis.numStates;k++)
        {
            B0[k+rsctrellis.numStates*n] -= A_mid[n+1];
            B1[k+rsctrellis.numStates*n] -= A_mid[n+1];
        }
    }

    //updated LLR for information bits
    extrinsic_data.set_size(N);
    for (n=0;n<N;n++)
    {
        sum0 = 0;
        sum1 = 0;
        for (k=0;k<rsctrellis.numStates;k++)
        {
            sum1 += exp(A1[k+n*rsctrellis.numStates]+B1[k+n*rsctrellis.numStates]);
            sum0 += exp(A0[k+n*rsctrellis.numStates]+B0[k+n*rsctrellis.numStates]);
        }
        extrinsic_data[n] = log(sum1/sum0)-apriori_data[n];//updated information must be independent of input LLR
    }

    //updated LLR for coded (parity) bits
    extrinsic_parity.set_size(N);
    for (n=0;n<N;n++)
    {
        sum0 = 0;
        sum1 = 0;
        for (k=0;k<rsctrellis.numStates;k++)
        {
            if (rsctrellis.fm[k])
            {
                sum1 += exp(A1[k+n*rsctrellis.numStates]+B1[k+n*rsctrellis.numStates]);
                sum0 += exp(A0[k+n*rsctrellis.numStates]+B0[k+n*rsctrellis.numStates]);
            }
            else
            {
                sum0 += exp(A1[k+n*rsctrellis.numStates]+B1[k+n*rsctrellis.numStates]);
                sum1 += exp(A0[k+n*rsctrellis.numStates]+B0[k+n*rsctrellis.numStates]);
            }
        }
        extrinsic_parity[n] = log(sum0/sum1)-Lc2I[n];//updated information must be independent of input LLR
    }
    //destroy trellis
    delete[] rsctrellis.prevStates;
    delete[] rsctrellis.nextStates;
    delete[] rsctrellis.PARout;
    delete[] rsctrellis.fm;
    //destroy MAP parameters
    delete[] Lc1I;
    delete[] Lc2I;
    delete[] A0;
    delete[] A1;
    delete[] A_mid;
    delete[] B0;
    delete[] B1;
}

void SISO::rsc_maxlogMAP(itpp::vec &extrinsic_parity, itpp::vec &extrinsic_data, const itpp::vec &intrinsic_coded, const itpp::vec &apriori_data)
/*
 maxlogMAP (SISO) decoder for RSC of rate 1/2
 extrinsic_parity - extrinsic information of parity bits
 extrinsic_data - extrinsic information of data (informational) bits
 intrinsic_coded - intrinsic information of coded (systematic and parity) bits
 apriori_data - a priori information of data (informational) bits
 Reference: Steven S. Pietrobon and Adrian S. Barbulescu, "A simplification of the modified Bahl decoding algorithm for systematic convolutional codes", Proc. ISITA, 1994*/
{
    //get parameters
    int N = apriori_data.length();
    //other parameters
    register int n,k;
    double buffer;
    int index;
    double A_min, A_max;
    double sum0, sum1;

    //trellis generation
    gen_rsctrellis();

    //parameter initialization
    double* Lc1I = new double[N];
    double* Lc2I = new double[N];
    for (n=0;n<N;n++)
    {
        Lc1I[n] = intrinsic_coded[2*n];
        Lc2I[n] = intrinsic_coded[2*n+1];
    }
    double* A0 = new double[rsctrellis.numStates*N];
    double* A1 = new double[rsctrellis.numStates*N];
    double* A_mid = new double[N];
    double* B0 = new double[rsctrellis.numStates*N];
    double* B1 = new double[rsctrellis.numStates*N];
    buffer = (tail?-INFINITY:0);//log(buffer)
    for (n=0;n<N;n++)
    {
        for (k=0;k<rsctrellis.numStates;k++)
        {
            A0[k+n*rsctrellis.numStates] = -INFINITY;
            A1[k+n*rsctrellis.numStates] = -INFINITY;
            B0[k+n*rsctrellis.numStates] = buffer;
            B1[k+n*rsctrellis.numStates] = buffer;
        }
        A_mid[n] = 0;
    }

    //A
    A0[0] = Lc2I[0]*rsctrellis.PARout[0];//i=0
    A1[0] = Lc1I[0]+apriori_data[0]+Lc2I[0]*rsctrellis.PARout[rsctrellis.numStates];//i=1
    for (n=1;n<N;n++)
    {
        A_min = INFINITY;
        A_max = 0;
        for (k=0;k<rsctrellis.numStates;k++)
        {
            buffer = std::max(A0[rsctrellis.prevStates[k]+(n-1)*rsctrellis.numStates], 
                     A1[rsctrellis.prevStates[k+rsctrellis.numStates]+(n-1)*rsctrellis.numStates]);//log(alpha0+alpha1)
            A0[k+rsctrellis.numStates*n] = Lc2I[n]*rsctrellis.PARout[k]+buffer;
            A1[k+rsctrellis.numStates*n] = Lc1I[n]+apriori_data[n]+Lc2I[n]*rsctrellis.PARout[k+rsctrellis.numStates]+buffer;
            //find min A(:,n)
            A_min = std::min(A_min, A0[k+rsctrellis.numStates*n]);
            //find max A(:,n)
            A_max = std::max(A_max, A0[k+rsctrellis.numStates*n]);
        }
        //normalization
        A_mid[n] = (A_min+A_max)/2;
        if (isinf(A_mid[n]))
            continue;
        for (k=0;k<rsctrellis.numStates;k++)
        {
            A0[k+rsctrellis.numStates*n] -= A_mid[n];
            A1[k+rsctrellis.numStates*n] -= A_mid[n];
        }
    }
    //B
    B0[rsctrellis.prevStates[0]+(N-1)*rsctrellis.numStates] = 0;
    B1[rsctrellis.prevStates[rsctrellis.numStates]+(N-1)*rsctrellis.numStates] = 0;
    for (n=N-2;n>=0;n--)
    {
        for (k=0;k<rsctrellis.numStates;k++)
        {
            index = rsctrellis.nextStates[k];
            B0[k+rsctrellis.numStates*n] = std::max(B0[index+(n+1)*rsctrellis.numStates]+Lc2I[n+1]*rsctrellis.PARout[index], 
                                              B1[index+(n+1)*rsctrellis.numStates]+Lc1I[n+1]+apriori_data[n+1]+Lc2I[n+1]*rsctrellis.PARout[index+rsctrellis.numStates]);
            index = rsctrellis.nextStates[k+rsctrellis.numStates];
            B1[k+rsctrellis.numStates*n] = std::max(B0[index+(n+1)*rsctrellis.numStates]+Lc2I[n+1]*rsctrellis.PARout[index], 
                                              B1[index+(n+1)*rsctrellis.numStates]+Lc1I[n+1]+apriori_data[n+1]+Lc2I[n+1]*rsctrellis.PARout[index+rsctrellis.numStates]);

        }
        if (isinf(A_mid[n+1]))
            continue;
        for (k=0;k<rsctrellis.numStates;k++)
        {
            B0[k+rsctrellis.numStates*n] -= A_mid[n+1];
            B1[k+rsctrellis.numStates*n] -= A_mid[n+1];
        }
    }

    //updated LLR for information bits
    extrinsic_data.set_size(N);
    for (n=0;n<N;n++)
    {
        sum0 = -INFINITY;
        sum1 = -INFINITY;
        for (k=0;k<rsctrellis.numStates;k++)
        {
            sum1 = std::max(sum1, A1[k+n*rsctrellis.numStates]+B1[k+n*rsctrellis.numStates]);
            sum0 = std::max(sum0, A0[k+n*rsctrellis.numStates]+B0[k+n*rsctrellis.numStates]);
        }
        extrinsic_data[n] = (sum1-sum0)-apriori_data[n];//updated information must be independent of input LLR
    }

    //updated LLR for coded (parity) bits
    extrinsic_parity.set_size(N);
    for (n=0;n<N;n++)
    {
        sum0 = -INFINITY;
        sum1 = -INFINITY;
        for (k=0;k<rsctrellis.numStates;k++)
        {
            if (rsctrellis.fm[k])
            {
                sum1 = std::max(sum1, A1[k+n*rsctrellis.numStates]+B1[k+n*rsctrellis.numStates]);
                sum0 = std::max(sum0, A0[k+n*rsctrellis.numStates]+B0[k+n*rsctrellis.numStates]);
            }
            else
            {
                sum0 = std::max(sum0, A1[k+n*rsctrellis.numStates]+B1[k+n*rsctrellis.numStates]);
                sum1 = std::max(sum1, A0[k+n*rsctrellis.numStates]+B0[k+n*rsctrellis.numStates]);
            }
        }
        extrinsic_parity[n] = (sum0-sum1)-Lc2I[n];//updated information must be independent of input LLR
    }
    //destroy trellis
    delete[] rsctrellis.prevStates;
    delete[] rsctrellis.nextStates;
    delete[] rsctrellis.PARout;
    delete[] rsctrellis.fm;
    //destroy MAP parameters
    delete[] Lc1I;
    delete[] Lc2I;
    delete[] A0;
    delete[] A1;
    delete[] A_mid;
    delete[] B0;
    delete[] B1;
}

void SISO::gen_nsctrellis(void)
/*
 generate trellis for non systematic non recursive convolutional codes
 r - number of outputs of the CC
 mem_len - memory length of the CC
 */
{
    //get parameters
    int r = gen.rows();
    int mem_len = gen.cols()-1;
    //other parameters
    register int n,k,j,p;
    itpp::bin inputs[] = {0,1};
    int index;

    nsctrellis.stateNb = (1<<mem_len);
    nsctrellis.output = new double[nsctrellis.stateNb*2*r];
    nsctrellis.nextState = new int[nsctrellis.stateNb*2];
    nsctrellis.prevState = new int[nsctrellis.stateNb*2];
    nsctrellis.input = new int[nsctrellis.stateNb*2];

    itpp::bvec enc_mem(mem_len);
    itpp::bin out;
    for (n=0;n<nsctrellis.stateNb;n++) //initial state
    {
        enc_mem = itpp::dec2bin(mem_len, n);
        //output
        for (k=0;k<2;k++)
            for (j=0;j<r;j++)
            {
                out = inputs[k]*gen(j,0);
                for (p=1;p<=mem_len;p++)
                    out ^= (enc_mem[p-1]*gen(j,p));//0 or 1
                nsctrellis.output[n+k*nsctrellis.stateNb+j*nsctrellis.stateNb*2] = double(out);
            }
        //next state
        for (j=(mem_len-1);j>0;j--)
            enc_mem[j] = enc_mem[j-1];
        for (k=0;k<2;k++)
        {
            enc_mem[0] = inputs[k];
            nsctrellis.nextState[n+k*nsctrellis.stateNb] = itpp::bin2dec(enc_mem, true);
        }
    };

    for (j=0;j<nsctrellis.stateNb;j++)
    {
        index = 0;
        for (n=0;n<nsctrellis.stateNb;n++)
            for (k=0;k<2;k++)
                if (nsctrellis.nextState[n+k*nsctrellis.stateNb]==j)
                {
                    nsctrellis.prevState[j+index*nsctrellis.stateNb] = n;
                    nsctrellis.input[j+index*nsctrellis.stateNb]  = k;//0 or 1
                    index++;
                }
    }
}

inline void SISO::nsc(itpp::vec &extrinsic_coded, itpp::vec &extrinsic_data, const itpp::vec &intrinsic_coded, const itpp::vec &apriori_data)
{
	if (gen.size()==0)
	{
		print_err_msg("SISO::nsc: generator polynomials not initialized");
		return;
	}
	
    if (MAP_metric=="logMAP")
        nsc_logMAP(extrinsic_coded, extrinsic_data, intrinsic_coded, apriori_data);
    else if (MAP_metric=="maxlogMAP")
        nsc_maxlogMAP(extrinsic_coded, extrinsic_data, intrinsic_coded, apriori_data);
    else
		print_err_msg("SISO::nsc: unknown MAP metric. The MAP metric should be either logMAP or maxlogMAP");
}

void SISO::nsc_logMAP(itpp::vec &extrinsic_coded, itpp::vec &extrinsic_data, const itpp::vec &intrinsic_coded, const itpp::vec &apriori_data)
/*
 * generalized decoder for NSC codes (after the NSC code a scrambler of pattern phi is used) using log MAP alg.
 * extrinsic_coded - extrinsic information of coded bits (output if the shaping filter)
 * extrinsic_data - extrinsic information of data bits
 * intrinsic_coded - intrinsic information of coded bits
 * apriori_data - a priori information of data bits
 */
{
    //get parameters
    int N = apriori_data.length();
    int Nc = scrambler_pattern.length();
    int r = gen.rows();//number of outputs of the CC
    //other parameters
    register int n,k,m,mp,j,i;
    int pstates[2];
    int nstates[2];
    int inputs[2];
    double C[2];//log(gamma)
    double sum;
    double sumbis;
    int index;

    //initialize trellis
    gen_nsctrellis();
    //initialize log(alpha) and log(beta)
    double* A = new double[nsctrellis.stateNb*(N+1)];
    double* B = new double[nsctrellis.stateNb*(N+1)];
    A[0] = 0;
    for (n=1;n<nsctrellis.stateNb;n++)
        A[n] = -INFINITY;
    B[N*nsctrellis.stateNb] = 0;
    sum = (tail?-INFINITY:0);
    for (n=1;n<nsctrellis.stateNb;n++)
        B[n+N*nsctrellis.stateNb] = sum;//if tail==false the final state is not known

    //forward recursion
    for (n=1;n<(N+1);n++)
    {
    	sum = 0;//normalization factor
        for (m=0;m<nsctrellis.stateNb;m++) //final state
        {
            for (k=0;k<2;k++)
            {
                pstates[k] = nsctrellis.prevState[m+nsctrellis.stateNb*k];//determine previous states
                inputs[k] = nsctrellis.input[m+nsctrellis.stateNb*k];//determine input
                C[k] = (inputs[k]?apriori_data[n-1]:0);//compute log of gamma
            }
            for (i=0;i<2;i++)//for each C[i]
            {
                for (k=0;k<r;k++)
                    for (j=0;j<Nc;j++)
                        C[i] += nsctrellis.output[pstates[i]+inputs[i]*nsctrellis.stateNb+k*nsctrellis.stateNb*2]*(1-2*double(scrambler_pattern[j]))*intrinsic_coded[j+k*Nc+(n-1)*Nc*r];
            }
            A[m+n*nsctrellis.stateNb] = itpp::log_add(A[pstates[0]+(n-1)*nsctrellis.stateNb]+C[0], A[pstates[1]+(n-1)*nsctrellis.stateNb]+C[1]);
            sum += exp(A[m+n*nsctrellis.stateNb]);
        }
        //normalization
        sum = log(sum);
        if (isinf(sum))
        	continue;
        for (m=0;m<nsctrellis.stateNb;m++)
            A[m+n*nsctrellis.stateNb] -= sum;
    }

    //backward recursion
    for (n=N-1;n>=0;n--)
    {
    	sum = 0;//normalisation factor
        for (m=0;m<nsctrellis.stateNb;m++) //initial state
        {
            for (k=0;k<2;k++)
            {
                nstates[k] = nsctrellis.nextState[m+k*nsctrellis.stateNb];//determine next states
                C[k] = (k?apriori_data[n]:0);//compute log of gamma
            }
            for (i=0;i<2;i++)
            {
                for (k=0;k<r;k++)
                    for (j=0;j<Nc;j++)
                        C[i] += nsctrellis.output[m+i*nsctrellis.stateNb+k*nsctrellis.stateNb*2]*(1-2*double(scrambler_pattern[j]))*intrinsic_coded[j+k*Nc+n*Nc*r];
            }
            B[m+n*nsctrellis.stateNb] = itpp::log_add(B[nstates[0]+(n+1)*nsctrellis.stateNb]+C[0], B[nstates[1]+(n+1)*nsctrellis.stateNb]+C[1]);
            sum += exp(B[m+n*nsctrellis.stateNb]);
        }
        //normalisation            
        sum = log(sum);
        if (isinf(sum))
        	continue;
        for (m=0;m<nsctrellis.stateNb;m++)
            B[m+n*nsctrellis.stateNb] -= sum;
    }

    //compute extrinsic_data
    extrinsic_data.set_size(N);
    for (n=1;n<(N+1);n++)
    {
        sum = 0;
        sumbis = 0;
        for (k=0;k<(nsctrellis.stateNb/2);k++)
        {
            sum += exp(A[k+(nsctrellis.stateNb/2)+n*nsctrellis.stateNb]+B[k+(nsctrellis.stateNb/2)+n*nsctrellis.stateNb]);//nominator
            sumbis += exp(A[k+n*nsctrellis.stateNb]+B[k+n*nsctrellis.stateNb]);//denominator
        }
        extrinsic_data[n-1] = log(sum/sumbis)-apriori_data[n-1];
    }

    //compute extrinsic_coded
    double *sum0 = new double[r];
    double *sum1 = new double[r];
    extrinsic_coded.set_size(N*Nc*r);
    for (n=0;n<N;n++)
    {
        for (k=0;k<r;k++)
        {
            sum0[k] = 0;
            sum1[k] = 0;
        }
        for (mp=0;mp<nsctrellis.stateNb;mp++) //initial state
            for (i=0;i<2;i++)
            {
                m = nsctrellis.nextState[mp+i*nsctrellis.stateNb];//final state
                //compute log of sigma
                index = (m>=(nsctrellis.stateNb/2));//0 if input is 0, 1 if input is 1
                C[0] = (index?apriori_data[n]:0);
                for (k=0;k<r;k++)
                    for (j=0;j<Nc;j++)
                        C[0] += nsctrellis.output[mp+index*nsctrellis.stateNb+k*nsctrellis.stateNb*2]*(1-2*double(scrambler_pattern[j]))*intrinsic_coded[j+k*Nc+n*Nc*r];
                C[1] = A[mp+n*nsctrellis.stateNb]+C[0]+B[m+(n+1)*nsctrellis.stateNb];//this is only a buffer
                //compute sums
                for (k=0;k<r;k++)
                    if (nsctrellis.output[mp+index*nsctrellis.stateNb+k*nsctrellis.stateNb*2])
                        sum1[k] += exp(C[1]);
                    else
                        sum0[k] += exp(C[1]);
            }
        for (k=0;k<r;k++)
            for (j=0;j<Nc;j++)
            {
            	index = j+k*Nc+n*r*Nc;
                extrinsic_coded[index] = (1-2*double(scrambler_pattern[j]))*log(sum1[k]/sum0[k])-intrinsic_coded[index];
            }
    }

    //free memory
    delete[] nsctrellis.output;
    delete[] nsctrellis.nextState;
    delete[] nsctrellis.prevState;
    delete[] nsctrellis.input;
    delete[] A;
    delete[] B;
    delete[] sum0;
    delete[] sum1;
}

void SISO::nsc_maxlogMAP(itpp::vec &extrinsic_coded, itpp::vec &extrinsic_data, const itpp::vec &intrinsic_coded, const itpp::vec &apriori_data)
/*
 * generalized decoder for NSC codes (after the NSC code a scrambler of pattern phi is used) using max log MAP alg.
 * extrinsic_coded - extrinsic information of coded bits (output if the shaping filter)
 * extrinsic_data - extrinsic information of data bits
 * intrinsic_coded - intrinsic information of coded bits
 * apriori_data - a priori information of data bits
 */
{
    //get parameters
    int N = apriori_data.length();
    int Nc = scrambler_pattern.length();
    int r = gen.rows();//number of outputs of the CC
    //other parameters
    register int n,k,m,mp,j,i;
    int pstates[2];
    int nstates[2];
    int inputs[2];
    double C[2];//log(gamma)
    double sum;
    double sumbis;
    int index;

    //initialize trellis
    gen_nsctrellis();
    //initialize log(alpha) and log(beta)
    double* A = new double[nsctrellis.stateNb*(N+1)];
    double* B = new double[nsctrellis.stateNb*(N+1)];
    A[0] = 0;
    for (n=1;n<nsctrellis.stateNb;n++)
        A[n] = -INFINITY;
    B[N*nsctrellis.stateNb] = 0;
    sum = (tail?-INFINITY:0);
    for (n=1;n<nsctrellis.stateNb;n++)
        B[n+N*nsctrellis.stateNb] = sum;//if tail==false the final state is not known

    //forward recursion
    for (n=1;n<(N+1);n++)
    {
    	sum = -INFINITY;//normalisation factor
        for (m=0;m<nsctrellis.stateNb;m++) //final state
        {
            for (k=0;k<2;k++)
            {
                pstates[k] = nsctrellis.prevState[m+nsctrellis.stateNb*k];//determine previous states
                inputs[k] = nsctrellis.input[m+nsctrellis.stateNb*k];//determine input
                C[k] = (inputs[k]?apriori_data[n-1]:0);//compute log of gamma
            }
            for (i=0;i<2;i++)//for each C[i]
            {
                for (k=0;k<r;k++)
                    for (j=0;j<Nc;j++)
                        C[i] += nsctrellis.output[pstates[i]+inputs[i]*nsctrellis.stateNb+k*nsctrellis.stateNb*2]*(1-2*double(scrambler_pattern[j]))*intrinsic_coded[j+k*Nc+(n-1)*Nc*r];
            }
            A[m+n*nsctrellis.stateNb] = std::max(A[pstates[0]+(n-1)*nsctrellis.stateNb]+C[0], A[pstates[1]+(n-1)*nsctrellis.stateNb]+C[1]);
            sum = std::max(sum, A[m+n*nsctrellis.stateNb]);
        }
        //normalization
        if (isinf(sum))
        	continue;
        for (m=0;m<nsctrellis.stateNb;m++)
            A[m+n*nsctrellis.stateNb] -= sum;
    }

    //backward recursion
    for (n=N-1;n>=0;n--)
    {
    	sum = -INFINITY;//normalisation factor
        for (m=0;m<nsctrellis.stateNb;m++) //initial state
        {
            for (k=0;k<2;k++)
            {
                nstates[k] = nsctrellis.nextState[m+k*nsctrellis.stateNb];//determine next states
                C[k] = (k?apriori_data[n]:0);//compute log of gamma
            }
            for (i=0;i<2;i++)
            {
                for (k=0;k<r;k++)
                    for (j=0;j<Nc;j++)
                        C[i] += nsctrellis.output[m+i*nsctrellis.stateNb+k*nsctrellis.stateNb*2]*(1-2*double(scrambler_pattern[j]))*intrinsic_coded[j+k*Nc+n*Nc*r];
            }
            B[m+n*nsctrellis.stateNb] = std::max(B[nstates[0]+(n+1)*nsctrellis.stateNb]+C[0], B[nstates[1]+(n+1)*nsctrellis.stateNb]+C[1]);
            sum = std::max(sum, B[m+n*nsctrellis.stateNb]);
        }
        //normalisation   
        if (isinf(sum))
        	continue;           
        for (m=0;m<nsctrellis.stateNb;m++)
            B[m+n*nsctrellis.stateNb] -= sum;
    }

    //compute extrinsic_data
    extrinsic_data.set_size(N);
    for (n=1;n<(N+1);n++)
    {
        sum = -INFINITY;
        sumbis = -INFINITY;
        for (k=0;k<(nsctrellis.stateNb/2);k++)
        {
            sum = std::max(sum, A[k+(nsctrellis.stateNb/2)+n*nsctrellis.stateNb]+B[k+(nsctrellis.stateNb/2)+n*nsctrellis.stateNb]);//nominator
            sumbis = std::max(sumbis, A[k+n*nsctrellis.stateNb]+B[k+n*nsctrellis.stateNb]);//denominator
        }
        extrinsic_data[n-1] = (sum-sumbis)-apriori_data[n-1];
    }

    //compute extrinsic_coded
    double *sum0 = new double[r];
    double *sum1 = new double[r];
    extrinsic_coded.set_size(N*Nc*r);
    for (n=0;n<N;n++)
    {
        for (k=0;k<r;k++)
        {
            sum0[k] = -INFINITY;
            sum1[k] = -INFINITY;
        }
        for (mp=0;mp<nsctrellis.stateNb;mp++) //initial state
            for (i=0;i<2;i++)
            {
                m = nsctrellis.nextState[mp+i*nsctrellis.stateNb];//final state
                //compute log of sigma
                index = (m>=(nsctrellis.stateNb/2));//0 if input is 0, 1 if input is 1
                C[0] = (index?apriori_data[n]:0);
                for (k=0;k<r;k++)
                    for (j=0;j<Nc;j++)
                        C[0] += nsctrellis.output[mp+index*nsctrellis.stateNb+k*nsctrellis.stateNb*2]*(1-2*double(scrambler_pattern[j]))*intrinsic_coded[j+k*Nc+n*Nc*r];
                C[1] = A[mp+n*nsctrellis.stateNb]+C[0]+B[m+(n+1)*nsctrellis.stateNb];//this is only a buffer
                //compute sums
                for (k=0;k<r;k++)
                    if (nsctrellis.output[mp+index*nsctrellis.stateNb+k*nsctrellis.stateNb*2])
                        sum1[k] = std::max(sum1[k], C[1]);
                    else
                        sum0[k] = std::max(sum0[k], C[1]);
            }
        for (k=0;k<r;k++)
            for (j=0;j<Nc;j++)
            {
            	index = j+k*Nc+n*r*Nc;
                extrinsic_coded[index] = (1-2*double(scrambler_pattern[j]))*(sum1[k]-sum0[k])-intrinsic_coded[index];
            }
    }

    //free memory
    delete[] nsctrellis.output;
    delete[] nsctrellis.nextState;
    delete[] nsctrellis.prevState;
    delete[] nsctrellis.input;
    delete[] A;
    delete[] B;
    delete[] sum0;
    delete[] sum1;
}

void SISO::gen_chtrellis(void)
// generate trellis for precoded FIR channels with real coefficients
{
    //get parameters
    int mem_len = impulse_response.cols()-1;//memory length of the channel
    int p_order = prec_gen.length()-1;//order of the precoder polynomial

    //other variables
    register int n,k,j;
    double inputs[] = {1.0,-1.0};//1->-1, 0->+1
    int index;
    double feedback[2];

    //create channel trellis
    int equiv_ch_mem_len = std::max(mem_len, p_order);
    chtrellis.stateNb = (1<<equiv_ch_mem_len);
    chtrellis.output = new double[2*chtrellis.stateNb];
    chtrellis.nextState = new int[2*chtrellis.stateNb];
    chtrellis.prevState = new int[2*chtrellis.stateNb];
    chtrellis.input = new int[2*chtrellis.stateNb];

    //initialize trellis
    itpp::ivec enc_mem(equiv_ch_mem_len);
    for (n=0;n<chtrellis.stateNb;n++) //initial state
    {
        enc_mem = 1-2*itpp::to_ivec(itpp::dec2bin(equiv_ch_mem_len, n));//1->-1, 0->+1
        //output
        for (k=0;k<2;k++)
        {
            feedback[k] = inputs[k];
            for (j=1;j<=p_order;j++)
                if (prec_gen[j])
                    feedback[k] *= enc_mem[j-1];//xor truth table must remain the same
            chtrellis.output[n+k*chtrellis.stateNb] = feedback[k]*impulse_response(0,0);
            for (j=1;j<=mem_len;j++)
                chtrellis.output[n+k*chtrellis.stateNb] += (enc_mem[j-1]*impulse_response(0,j));
        }
        //next state
        for (j=(equiv_ch_mem_len-1);j>0;j--)
            enc_mem[j] = enc_mem[j-1];
        for (k=0;k<2;k++)
        {
            enc_mem[0] = int(feedback[k]);
            chtrellis.nextState[n+k*chtrellis.stateNb] = itpp::bin2dec(itpp::to_bvec((1-enc_mem)/2), true);//-1->1, +1->0
        }
    }

    for (j=0;j<chtrellis.stateNb;j++)
    {
        index = 0;
        for (n=0;n<chtrellis.stateNb;n++)
        {
            for (k=0;k<2;k++)
            {
                if (chtrellis.nextState[n+k*chtrellis.stateNb]==j)
                {
                    chtrellis.prevState[j+index*chtrellis.stateNb] = n;
                    chtrellis.input[j+index*chtrellis.stateNb] = k;//this is an index to the channel input
                    index++;
                }
            }
        }
    }
}

inline void SISO::equalizer(itpp::vec &extrinsic_data, const itpp::vec &rec_sig, const itpp::vec &apriori_data)
{
	if (impulse_response.size()==0)
	{
		print_err_msg("SISO::equalizer: channel impulse response not initialized");
		return;
	}
	if ((impulse_response.size()==1)&&(prec_gen.length()==1))
	{
		print_err_msg("SISO::equalizer: flat fading channel and no precoder. Use the soft output of the channel (no need for a priori information)");
		return;
	}
	
    if (MAP_metric=="logMAP")
        equalizer_logMAP(extrinsic_data, rec_sig, apriori_data);
    else if (MAP_metric=="maxlogMAP")
        equalizer_maxlogMAP(extrinsic_data, rec_sig, apriori_data);
    else
		print_err_msg("SISO::equalizer: unknown MAP metric. The MAP metric should be either logMAP or maxlogMAP");
}

void SISO::equalizer_logMAP(itpp::vec &extrinsic_data, const itpp::vec &rec_sig, const itpp::vec &apriori_data)
/*
  extrinsic_data - extrinsic information of channel input symbols
  rec_sig - received symbols
  apriori_data - a priori information of channel input symbols
 */
{
    //get parameters
    int N = rec_sig.length();//length of the received frame
    //other parameters
    register int n,k,m;
    int pstates[2];
    int nstates[2];
    int inputs[2];
    double C[2];//log(gamma)
    double sum;
    double sumbis;
    double buffer;

    //initialize trellis
    gen_chtrellis();
    //initialize log(alpha) and log(beta)
    double* A = new double[chtrellis.stateNb*(N+1)];
    double* B = new double[chtrellis.stateNb*(N+1)];
    A[0] = 0;
    B[N*chtrellis.stateNb] = 0;
    sum = (tail?-INFINITY:0);
    for (n=1;n<chtrellis.stateNb;n++)
    {
        A[n] = -INFINITY;
        B[n+N*chtrellis.stateNb] = sum;//if tail==false the final state is not known
    }

    //forward recursion
    for (n=1;n<=N;n++)
    {
    	sum = 0;//normalisation factor
        for (m=0;m<chtrellis.stateNb;m++) //final state
        {
            for (k=0;k<2;k++)
            {
                pstates[k] = chtrellis.prevState[m+chtrellis.stateNb*k];//determine previous states
                inputs[k] = chtrellis.input[m+chtrellis.stateNb*k];//determine input
                C[k] = (inputs[k])*apriori_data[n-1]-itpp::sqr(rec_sig[n-1]-chtrellis.output[pstates[k]+chtrellis.stateNb*inputs[k]])/(2*sigma2);//compute log of gamma
            }
            A[m+n*chtrellis.stateNb] = itpp::log_add(A[pstates[0]+(n-1)*chtrellis.stateNb]+C[0], A[pstates[1]+(n-1)*chtrellis.stateNb]+C[1]);
            sum += exp(A[m+n*chtrellis.stateNb]);
        }
        //normalisation
        sum = log(sum);
        for (m=0;m<chtrellis.stateNb;m++)
        {
            A[m+n*chtrellis.stateNb] -= sum;
        }
    }

    //backward recursion
    for (n=N-1;n>=0;n--)
    {
    	sum = 0;//normalisation factor
        for (m=0;m<chtrellis.stateNb;m++) //initial state
        {
            for (k=0;k<2;k++)
            {
                nstates[k] = chtrellis.nextState[m+k*chtrellis.stateNb];//determine next states
                C[k] = (k)*apriori_data[n]-itpp::sqr(rec_sig[n]-chtrellis.output[m+k*chtrellis.stateNb])/(2*sigma2);//compute log of gamma
            }
            B[m+n*chtrellis.stateNb] = itpp::log_add(B[nstates[0]+(n+1)*chtrellis.stateNb]+C[0], B[nstates[1]+(n+1)*chtrellis.stateNb]+C[1]);
            sum += exp(B[m+n*chtrellis.stateNb]);
        }
        //normalisation
        sum = log(sum);
        for (m=0;m<chtrellis.stateNb;m++)
        {
            B[m+n*chtrellis.stateNb] -= sum;
        }
    }

    //compute extrinsic_data
    extrinsic_data.set_size(N);
    for (n=1;n<=N;n++)
    {
        sum = 0;//could be replaced by a vector
        sumbis = 0;
        for (m=0;m<chtrellis.stateNb;m++) //initial state
        {
            for (k=0;k<2;k++) //input index
            {
                nstates[k] = chtrellis.nextState[m+k*chtrellis.stateNb];//determine next states
                C[k] = (k)*apriori_data[n-1]-itpp::sqr(rec_sig[n-1]-chtrellis.output[m+k*chtrellis.stateNb])/(2*sigma2);//compute log of gamma
                buffer = exp(A[m+(n-1)*chtrellis.stateNb]+C[k]+B[nstates[k]+n*chtrellis.stateNb]);
                if (k)
                    sum += buffer;//1
                else
                    sumbis += buffer;//0
            }
        }
        extrinsic_data[n-1] = log(sum/sumbis)-apriori_data[n-1];
    }

    //free memory
    delete[] chtrellis.output;
    delete[] chtrellis.nextState;
    delete[] chtrellis.prevState;
    delete[] chtrellis.input;
    delete[] A;
    delete[] B;
}

void SISO::equalizer_maxlogMAP(itpp::vec &extrinsic_data, const itpp::vec &rec_sig, const itpp::vec &apriori_data)
/*
  extrinsic_data - extrinsic information for channel input symbols
  rec_sig - received symbols
  apriori_data - a priori information for channel input symbols
 */
{
    //get parameters
    int N = rec_sig.length();//length of the received frame
    //other parameters
    register int n,k,m;
    int pstates[2];
    int nstates[2];
    int inputs[2];
    double C[2];//log(gamma)
    double sum;
    double sumbis;
    double buffer;

    //initialize trellis
    gen_chtrellis();
    //initialize log(alpha) and log(beta)
    double* A = new double[chtrellis.stateNb*(N+1)];
    double* B = new double[chtrellis.stateNb*(N+1)];
    A[0] = 0;
    for (n=1;n<chtrellis.stateNb;n++)
        A[n] = -INFINITY;
    B[N*chtrellis.stateNb] = 0;
    sum = (tail?-INFINITY:0);
    for (n=1;n<chtrellis.stateNb;n++)
        B[n+N*chtrellis.stateNb] = sum;//if tail==false the final state is not known

    //forward recursion
    for (n=1;n<=N;n++)
    {
    	sum = -INFINITY;//normalization factor
        for (m=0;m<chtrellis.stateNb;m++) //final state
        {
            for (k=0;k<2;k++)
            {
                pstates[k] = chtrellis.prevState[m+chtrellis.stateNb*k];//determine previous states
                inputs[k] = chtrellis.input[m+chtrellis.stateNb*k];//determine input
                C[k] = (inputs[k])*apriori_data[n-1]-itpp::sqr(rec_sig[n-1]-chtrellis.output[pstates[k]+chtrellis.stateNb*inputs[k]])/(2*sigma2);//compute log of gamma
            }
            A[m+n*chtrellis.stateNb] = std::max(A[pstates[0]+(n-1)*chtrellis.stateNb]+C[0], A[pstates[1]+(n-1)*chtrellis.stateNb]+C[1]);
            sum = std::max(sum, A[m+n*chtrellis.stateNb]);
        }
        for (m=0;m<chtrellis.stateNb;m++)
            A[m+n*chtrellis.stateNb] -= sum;
    }

    //backward recursion
    for (n=N-1;n>=0;n--)
    {
    	sum = -INFINITY;//normalisation factor
        for (m=0;m<chtrellis.stateNb;m++) //initial state
        {
            for (k=0;k<2;k++)
            {
                nstates[k] = chtrellis.nextState[m+k*chtrellis.stateNb];//determine next states
                C[k] = (k)*apriori_data[n]-itpp::sqr(rec_sig[n]-chtrellis.output[m+k*chtrellis.stateNb])/(2*sigma2);//compute log of gamma
            }
            B[m+n*chtrellis.stateNb] = std::max(B[nstates[0]+(n+1)*chtrellis.stateNb]+C[0], B[nstates[1]+(n+1)*chtrellis.stateNb]+C[1]);
            sum = std::max(sum, B[m+n*chtrellis.stateNb]);
        }
        for (m=0;m<chtrellis.stateNb;m++)
            B[m+n*chtrellis.stateNb] -= sum;
    }

    //compute extrinsic_data
    extrinsic_data.set_size(N);
    for (n=1;n<=N;n++)
    {
        sum = -INFINITY;
        sumbis = -INFINITY;
        for (m=0;m<chtrellis.stateNb;m++) //initial state
        {
            for (k=0;k<2;k++) //input index
            {
                nstates[k] = chtrellis.nextState[m+k*chtrellis.stateNb];//determine next states
                C[k] = (k)*apriori_data[n-1]-itpp::sqr(rec_sig[n-1]-chtrellis.output[m+k*chtrellis.stateNb])/(2*sigma2);//compute log of gamma
                buffer = A[m+(n-1)*chtrellis.stateNb]+C[k]+B[nstates[k]+n*chtrellis.stateNb];
                if (k)
                    sum = std::max(sum, buffer);//1
                else
                    sumbis = std::max(sumbis, buffer);//0
            }
        }
        extrinsic_data[n-1] = (sum-sumbis)-apriori_data[n-1];
    }

    //free memory
    delete[] chtrellis.output;
    delete[] chtrellis.nextState;
    delete[] chtrellis.prevState;
    delete[] chtrellis.input;
    delete[] A;
    delete[] B;
}

void SISO::descrambler(itpp::vec &extrinsic_coded, itpp::vec &extrinsic_data, const itpp::vec &intrinsic_coded, const itpp::vec &apriori_data)
/*
  inputs:
  intrinsic_coded - intrinsic information of coded bits (repetition code output)  
  apriori_data - a priori information of informational bits (repetition code input)
  outputs:
  extrinsic_coded - extrinsic information of coded bits
  extrinsic_data - Logarithm of Likelihood Ratio of informational bits  
*/
{
    //get parameters
    int nb_bits = apriori_data.length();
    int Nc = scrambler_pattern.length();
    //implementation
    extrinsic_data.set_size(nb_bits);
    extrinsic_coded.set_size(nb_bits*Nc);
    register int n,k;    
    for (k=0;k<nb_bits;k++)
    {
        extrinsic_data[k] = 0;//apriori_data[k];//add a priori info
        for (n=0;n<Nc;n++)
            extrinsic_data[k] += (1-2*double(scrambler_pattern[n]))*intrinsic_coded[n+k*Nc];
        for (n=0;n<Nc;n++)
            extrinsic_coded[n+k*Nc] = (1-2*double(scrambler_pattern[n]))*extrinsic_data[k]-intrinsic_coded[n+k*Nc];
    }
}

inline void SISO::zpFIRfilter(itpp::vec &filt, const itpp::vec &h, const itpp::vec &sig)
//FIR filter for a zero padded signal (L zeros are added at the end of the signal)
{
    //get parameters
    int L = h.length()-1;
    int N = sig.length();
    //implementation
    register int n,l;
    for (n=0;n<(N+L);n++)
    {
        filt[n] = 0;
        for (l=0;l<=L;l++)
        {
            if ((n-l)<0)
                break;//channel has state 0 at the beginning
            if ((n-l)>=N)
                continue;//channel has state 0 in the end
            filt[n] += (h[l]*sig[n-l]);
        }
    }
}

void SISO::gen_hyperTrellis(void)
/* generates channel hyper trellis for binary symbols
 * the channel is a MISO system
 * BPSK mapping: 0->+1, 1->-1
 */
{
	//get parameters
	int nb_usr = impulse_response.rows();
	int ch_order = impulse_response.cols()-1;
	int p_order = prec_gen.length()-1;
	int max_order = std::max(ch_order, p_order);

	//initialize hypertrellis
	chtrellis.numInputSymbols = itpp::pow2i(nb_usr);
	int mem_len = nb_usr*max_order;
	if (mem_len>=(int)(8*sizeof(int)))
	{
		std::string msg = "SISO::gen_hyperTrellis: memory length of the hyperchannel should be at most ";
		msg += itpp::to_str(8*sizeof(int)-1);
		msg += ". Try to lower the number of users, channel order or the order of the precoding polynomial (if any)";
		print_err_msg(msg);
		return;
	}
	chtrellis.stateNb = itpp::pow2i(mem_len);
	try
	{	
		long long int len =  static_cast<long long int>(chtrellis.stateNb)*static_cast<long long int>(chtrellis.numInputSymbols);
		chtrellis.nextState = new int[len];	
		chtrellis.prevState = new int[len];
		chtrellis.output = new double[len];
		chtrellis.input = new int[len];
	} catch (std::bad_alloc)
	{
		std::string msg = "SISO::gen_hyperTrellis: not enough memory for the channel trellis variables. The number of states is ";
		msg += itpp::to_str(chtrellis.stateNb);
		msg += " and the number of input symbols ";
		msg += itpp::to_str(chtrellis.numInputSymbols);
		print_err_msg(msg);
		return;
	}
	itpp::ivec index(chtrellis.stateNb);	    
    index.zeros();
	itpp::bvec hyper_ch_mem(mem_len);
	itpp::bvec hyper_ch_in(nb_usr);
	itpp::bvec hyper_states(mem_len);
	itpp::bin feedback;
    
    //create hypertrellis
    register int n,k,p,r;
    int buffer;
	double hyper_ch_out;
	for(k=0;k<chtrellis.stateNb;k++)
	{
		hyper_ch_mem = itpp::dec2bin(mem_len, k);//initial state
		for(n=0;n<chtrellis.numInputSymbols;n++)
		{
			hyper_ch_in = itpp::dec2bin(nb_usr, n);//MISO channel input
            hyper_ch_out = 0;
            for(r=0;r<nb_usr;r++)
            {
            	buffer = r*max_order;
            	//precoder
            	feedback = hyper_ch_in[r];
            	for (p=1;p<=p_order;p++)
            		if (prec_gen(p))
            			feedback ^= hyper_ch_mem[p-1+buffer];//xor
            	//FIR channel output                
                hyper_ch_out += (1-2*double(feedback))*impulse_response(r,0);
                for(p=0;p<ch_order;p++)
                    hyper_ch_out += (1-2*double(hyper_ch_mem[p+buffer]))*impulse_response(r,p+1);//hyper channel output for user r
                //(equivalent) channel next state
                hyper_states[buffer] = feedback;
                for(p=0;p<(max_order-1);p++)
                	hyper_states[p+buffer+1] = hyper_ch_mem[p+buffer];//next hyper state for user r
            }
			chtrellis.output[k+n*chtrellis.stateNb] = hyper_ch_out;            
			buffer = itpp::bin2dec(hyper_states);//next state from an initial state and a given input 
			chtrellis.nextState[k+n*chtrellis.stateNb] = buffer;
			chtrellis.prevState[buffer+index[buffer]*chtrellis.stateNb] = k;
			chtrellis.input[buffer+index[buffer]*chtrellis.stateNb] = n;
			index[buffer]++;
		}
	}
}

inline void SISO::mud(itpp::mat &extrinsic_data, const itpp::vec &rec_sig, const itpp::mat &apriori_data)
{
	if (impulse_response.size()==0)
	{
		print_err_msg("SISO::mud: channel impulse response not initialized");
		return;
	}
	
	if (MUD_method=="maxlogMAP")
		mud_maxlogMAP(extrinsic_data, rec_sig, apriori_data);
    else if (MUD_method=="GCD")
        GCD(extrinsic_data, rec_sig, apriori_data);
    else if (MUD_method=="sGCD")
        sGCD(extrinsic_data, rec_sig, apriori_data);
    else
		print_err_msg("SISO::mud: unknown MUD method. The MUD method should be either GCD or sGCD");
}

/// Maximum A Posteriori algorithm for Multi-User Detection in IDMA systems
/** uses max log MAP algorithm
 * use with care for large number of users and/or FIR channel order
 */
void SISO::mud_maxlogMAP(itpp::mat &extrinsic_data, const itpp::vec &rec_sig, const itpp::mat &apriori_data)
/* output:
 * extrinsic_data - extrinsic information for the chips (usr_nb x block_len)
 * inputs:
 * rec_sig - received signal (1 x block_len)
 * apriori_data - a priori information for the chips (usr_nb x block_len)
 * impulse_response - channel matrix (usr_nb x (L+1))
 * prec_gen - precoder generator polynomial
 * sigma2 - noise variance
 * tail - if true the final state is known
 */
{
    //get parameters
    int nb_usr = apriori_data.rows();
    int block_len = apriori_data.cols();
    
    //init trellis
    gen_hyperTrellis();

    //initial conditions for A = log(alpha) and B = log(beta)
    double *A = NULL,*B = NULL;
    try
    {
    	A = new double[chtrellis.stateNb*(block_len+1)];
    	B = new double[chtrellis.stateNb*(block_len+1)];
    } catch (std::bad_alloc)
    {
    	std::string msg = "SISO::mud_maxlogMAP: Not enough memory for alphas and betas. The number of states is ";
    	msg += itpp::to_str(chtrellis.stateNb);
    	msg += " and the block length ";
    	msg += itpp::to_str(block_len);
    	print_err_msg(msg);
    }
    register int n;    
    A[0] = 0;
    B[block_len*chtrellis.stateNb] = 0;
    double buffer = (tail?-INFINITY:0);
    for (n=1;n<chtrellis.stateNb;n++)
    {
        A[n] = -INFINITY;
        B[n+block_len*chtrellis.stateNb] = buffer;//if tail==false the final state is not known
    }
    
    //compute log(alpha) (forward recursion)
    register int s,k;
    int sp,i;
    itpp::bvec in_chips(nb_usr);
    for(n=1;n<=block_len;n++)
    {
        buffer = -INFINITY;//normalization factor
        for(s=0;s<chtrellis.stateNb;s++)
        {
        	A[s+n*chtrellis.stateNb] = -INFINITY;
            for(k=0;k<chtrellis.numInputSymbols;k++)
            {
                sp = chtrellis.prevState[s+k*chtrellis.stateNb];
                i = chtrellis.input[s+k*chtrellis.stateNb];
                in_chips = itpp::dec2bin(nb_usr, i);
                A[s+n*chtrellis.stateNb] = std::max(A[s+n*chtrellis.stateNb], \
                								A[sp+(n-1)*chtrellis.stateNb]-itpp::sqr(rec_sig[n-1]-chtrellis.output[sp+i*chtrellis.stateNb])/(2*sigma2)+\
                                                itpp::to_vec(in_chips)*apriori_data.get_col(n-1));
            }
            buffer = std::max(buffer, A[s+n*chtrellis.stateNb]);
        }
        //normalization
        for(s=0;s<chtrellis.stateNb;s++)
            A[s+n*chtrellis.stateNb] -= buffer;
    }
    
    //compute log(beta) (backward recursion)
    for(n=block_len-1;n>=0;n--)
    {
        buffer = -INFINITY;//normalization factor
        for(s=0;s<chtrellis.stateNb;s++)
        {
        	B[s+n*chtrellis.stateNb] = -INFINITY;
            for(k=0;k<chtrellis.numInputSymbols;k++)
            {
                sp = chtrellis.nextState[s+k*chtrellis.stateNb];
                in_chips = itpp::dec2bin(nb_usr, k);
                B[s+n*chtrellis.stateNb] = std::max(B[s+n*chtrellis.stateNb], \
                								B[sp+(n+1)*chtrellis.stateNb]-itpp::sqr(rec_sig[n]-chtrellis.output[s+k*chtrellis.stateNb])/(2*sigma2)+\
                                                itpp::to_vec(in_chips)*apriori_data.get_col(n));
            }
            buffer = std::max(buffer, B[s+n*chtrellis.stateNb]);
        }
        //normalization
        for(s=0;s<chtrellis.stateNb;s++)
            B[s+n*chtrellis.stateNb] -= buffer;
    }
    
    //compute extrinsic information
    double nom, denom;
    extrinsic_data.set_size(nb_usr,block_len);
    register int u;
    for(u=0;u<nb_usr;u++)
    {
        for(n=1;n<=block_len;n++)
        {
            nom = -INFINITY;
            denom = -INFINITY;
            for(s=0;s<chtrellis.stateNb;s++)
            {
                for(i=0;i<chtrellis.numInputSymbols;i++)
                {
                    in_chips = itpp::dec2bin(nb_usr, i);
                    buffer = A[s+(n-1)*chtrellis.stateNb]+B[chtrellis.nextState[s+i*chtrellis.stateNb]+n*chtrellis.stateNb]-\
                             itpp::sqr(rec_sig[n-1]-chtrellis.output[s+i*chtrellis.stateNb])/(2*sigma2)+
                             itpp::to_vec(in_chips)*apriori_data.get_col(n-1);
                    if(in_chips[u])
                        nom = std::max(nom, buffer);
                    else
                        denom = std::max(denom, buffer);
                }
            }
            extrinsic_data(u,n-1) = (nom-denom)-apriori_data(u,n-1);
        }
    }
    //free memory
    delete[] chtrellis.nextState;
    delete[] chtrellis.prevState;
    delete[] chtrellis.output;
    delete[] chtrellis.input;
    delete[] A;
    delete[] B;
}

/// Gaussian Chip Detector for IDMA systems
/** Use with care for large size of interleavers.
 */
void SISO::GCD(itpp::mat &extrinsic_data, const itpp::vec &rec_sig, const itpp::mat &apriori_data)
/* Gaussian Chip Detector
 * output:
 * extrinsic_data - extrinsic information of emitted chips
 * inputs:
 * rec_sig - received signal
 * apriori_data - a priori information of emitted chips
 */
{
    //get parameters
    int N = apriori_data.cols();//emitted frames of non-zero samples
    int K = apriori_data.rows();//number of users
    int L = impulse_response.cols()-1;//channel order
    //other parameters
    register int n,k;

    //mean and variance of each chip (NxK)
    itpp::mat Ex = -itpp::tanh(apriori_data/2.0);//take into account BPSK mapping
    itpp::mat Vx = 1.0-sqr(Ex);

    //expectation and variance of the received signal
    itpp::vec Er(N+L);
    Er.zeros();
    itpp::mat Covr;
    try
    {
    	Covr.set_size(N+L,N+L);
    } catch (std::bad_alloc)
    {
    	std::string msg = "SISO::GCD: not enough memory when allocating space for the covariance matrix. The interleaver length is ";
    	msg += itpp::to_str(N);
    	msg += ". Use sGCD instead or try to reduce the interleaver length";
    	print_err_msg(msg);
    	return;
    }    
    //create eye function
    Covr.zeros();
    for (n=0;n<(N+L);n++)
    	Covr(n,n) = 1;
    itpp::vec col(N+L);
    col.zeros();
    itpp::vec row(N);
    row.zeros();
    itpp::mat h_eq(N+L,N);
    for (k=0;k<K;k++)
    {
        col.replace_mid(0, impulse_response.get_row(k));
        row(0) = impulse_response(k,0);
        h_eq = itpp::toeplitz(col, row);
        Er += h_eq*Ex.get_row(k);
        Covr += (h_eq*itpp::diag(Vx.get_row(k)))*h_eq.T();
    }

    //extrinsic information
    double nom,denom;
    Er = rec_sig-Er;
    itpp::mat inv_Covr(N+L,N+L);
    inv_Covr = itpp::inv(Covr);
    itpp::mat inv_cov_zeta(N+L,N+L);
    itpp::vec rec_chip(N+L);
    extrinsic_data.set_size(K,N);
    for (k=0;k<K;k++)
        for (n=0;n<N;n++)
        {
            col.replace_mid(n, impulse_response.get_row(k));
            inv_cov_zeta = inv_Covr+itpp::outer_product(inv_Covr*col, inv_Covr.T()*(col*Vx(k,0)))/(1-(col*Vx(k,0))*(inv_Covr*col));
            rec_chip = Er-col*(+1-Ex(k,n));
            nom = -0.5*rec_chip*(inv_cov_zeta*rec_chip);
            rec_chip = Er-col*(-1-Ex(k,n));
            denom = -0.5*rec_chip*(inv_cov_zeta*rec_chip);
            extrinsic_data(k,n) = denom-nom;//take into account BPSK mapping
            col(n) = 0;
        }
}

/// simplified Gaussian Chip Detector for IDMA systems
/** This algorithm is simplified and uses much less memory than its counterpart, the GCD.
 * Recommended to use in most cases.
 */
void SISO::sGCD(itpp::mat &extrinsic_data, const itpp::vec &rec_sig, const itpp::mat &apriori_data)
/* simplified GCD
 * output:
 * extrinsic_data - extrinsic information of emitted chips
 * inputs:
 * rec_sig - received signal
 * apriori_data - a priori information of emitted chips
 */
{
    //get parameters
    int N = apriori_data.cols();//emitted frames of non-zero samples
    int K = apriori_data.rows();//number of users
    int L = impulse_response.cols()-1;//channel order
    //other parameters
    register int n,k;

    //mean and variance of each chip (NxK)
    itpp::mat Ex = -itpp::tanh(apriori_data/2.0);//take into account BPSK mapping
    itpp::mat Vx = 1.0-itpp::sqr(Ex);

    //mean and variance for the samples of the received signal
    itpp::vec Er(N+L);
    Er.zeros();
    itpp::vec Vr = sigma2*itpp::ones(N+L);
    itpp::vec buffer(N+L);
    for (k=0;k<K;k++)
    {
        zpFIRfilter(buffer, impulse_response.get_row(k), Ex.get_row(k));
        Er += buffer;
        zpFIRfilter(buffer, itpp::sqr(impulse_response.get_row(k)), Vx.get_row(k));
        Vr += buffer;
    }

    //extrinsic information for the samples of the received signal
    Er = rec_sig-Er;
    itpp::vec ch(L+1);
    extrinsic_data.set_size(K,N);
    for (k=0;k<K;k++)
    {
        ch = impulse_response.get_row(k);
        for (n=0;n<N;n++)
            extrinsic_data(k,n) = -2*itpp::elem_div(ch, Vr.mid(n,L+1)-sqr(ch)*Vx(k,n))*(Er.mid(n,L+1)+ch*Ex(k,n));//take into account BPSK mapping
    }
}

void SISO::find_half_const(int &select_half, itpp::vec &re_part, itpp::bmat &re_bin_part, itpp::vec &im_part, itpp::bmat &im_bin_part)
/* finds real in imaginary parts of the constellation and its corresponding bits
 * this approach is used for equivalent channel according to Hassibi's model 
 * the constellation must be quadratic and the number of bits per symbol must be a multiple of two
 */
{
    //values needed for initializations
    int const_size = itpp::pow2i(nb_bits_symb);//constellation size
    int half_nb_bits_symb = nb_bits_symb/2;
    int half_len = itpp::pow2i(half_nb_bits_symb);//number of values of real(imaginary) part    
    //initialize output variables
    select_half = 0;
    re_part.set_size(half_len);    
    re_bin_part.set_size(half_len, half_nb_bits_symb);
	re_part.zeros();
	re_part(0) = constellation(0).real();    
    im_part.set_size(half_len);        
    im_bin_part.set_size(half_len, half_nb_bits_symb);
	im_part.zeros();
	im_part(0) = constellation(0).imag();    
    //select half for real (imaginary) to binary correspondence
    if (nb_bits_symb%2)
	{
		print_err_msg("SISO::find_half_const: number of bits per symbol must be a multiple of two");
		return;
	}
    const double min_diff = 1e-3;
	itpp::ivec idx = itpp::find(itpp::abs(itpp::real(constellation)-re_part(0))<min_diff);
	if (idx.length()!=half_len)
	{
		print_err_msg("SISO::find_half_const: the constellation must be quadratic");
		return;
	}
    itpp::bvec temp(nb_bits_symb);
    register int n;
    for (n=0;n<2;n++)
    {
    	temp = bin_constellation.get_row(idx(n));        
    	re_bin_part.set_row(n,temp(0,half_nb_bits_symb-1));
    }
    select_half = (re_bin_part.get_row(0)==re_bin_part.get_row(1))?0:1;
    //algorithm    
    double buffer;
	temp = bin_constellation.get_row(0);
    re_bin_part.set_row(0,temp(select_half*half_nb_bits_symb,(1+select_half)*half_nb_bits_symb-1));
    im_bin_part.set_row(0,temp((1-select_half)*half_nb_bits_symb,(2-select_half)*half_nb_bits_symb-1));   
    int re_idx = 0;
    int im_idx = 0;
    for (n=1;n<const_size;n++)
    {
        temp = bin_constellation.get_row(n);
        buffer = constellation(n).real();
        if (itpp::prod(re_part-buffer)>min_diff)
        {
            re_idx++;
            re_part(re_idx) = buffer;
            re_bin_part.set_row(re_idx, temp(select_half*half_nb_bits_symb,(1+select_half)*half_nb_bits_symb-1));
        }
        buffer = constellation(n).imag();
        if (itpp::prod(im_part-buffer)>min_diff)
        {
            im_idx++;
            im_part(im_idx) = buffer;
            im_bin_part.set_row(im_idx, temp((1-select_half)*half_nb_bits_symb,(2-select_half)*half_nb_bits_symb-1));
        }
    }
}

inline void SISO::EquivRecSig(itpp::vec &x_eq, const itpp::cmat &rec_sig)
//finds equivalent received signal with real coefficients
//the equivalent received signal follows the model of Hassibi's paper
//ouput:
//x_eq - equivalent received signal with real coefficients
//inputs:
//rec_sig - received signal
{
    for (int k=0;k<nb_rec_ant;k++)
    {
        x_eq.set_subvector(k*2*block_duration, itpp::real(rec_sig.get_col(k)));
        x_eq.set_subvector(k*2*block_duration+block_duration, itpp::imag(rec_sig.get_col(k)));
    }
}

void SISO::EquivCh(itpp::mat &H_eq, const itpp::cvec &H)
//finds equivalent channel with real coefficients following the model of Hassibi's paper
//output:
//H_eq - equivalent channel
//input:
//H - channel matrix
{
    itpp::mat Aq(2*block_duration,2*nb_em_ant);
    itpp::mat Bq(2*block_duration,2*nb_em_ant);
    itpp::cmat temp(block_duration,nb_em_ant);
    itpp::vec h(2*nb_em_ant);
    itpp::mat AhBh(2*block_duration,2);
    register int n,k;
    for (k=0;k<symbols_block;k++)
    {
        temp = ST_gen1.get(k*block_duration,k*block_duration+block_duration-1,0,nb_em_ant-1);
        Aq.set_submatrix(0, 0, itpp::real(temp));
        Aq.set_submatrix(0, nb_em_ant, -itpp::imag(temp));
        Aq.set_submatrix(block_duration, 0, itpp::imag(temp));
        Aq.set_submatrix(block_duration, nb_em_ant, itpp::real(temp));
        temp = ST_gen2.get(k*block_duration,k*block_duration+block_duration-1,0,nb_em_ant-1);
        Bq.set_submatrix(0, 0, -itpp::imag(temp));
        Bq.set_submatrix(0, nb_em_ant, -itpp::real(temp));
        Bq.set_submatrix(block_duration, 0, itpp::real(temp));
        Bq.set_submatrix(block_duration, nb_em_ant, -itpp::imag(temp));
        for (n=0;n<nb_rec_ant;n++)
        {
            h.set_subvector(0, real(H.mid(n*nb_em_ant,nb_em_ant)));
            h.set_subvector(nb_em_ant, imag(H.mid(n*nb_em_ant,nb_em_ant)));
            AhBh.set_col(0, Aq*h);
            AhBh.set_col(1, Bq*h);
            H_eq.set_submatrix(2*block_duration*n, 2*k, AhBh);
        }
    }
}

inline void SISO::demapper(itpp::vec &extrinsic_data, const itpp::cvec &rec_sig, const itpp::vec &apriori_data)
{
	if (c_impulse_response.size()==0)
	{
		print_err_msg("SISO::demapper: channel impulse response not initialized");
		return;
	}
	if ((constellation.size()==0) || (bin_constellation.size()==0))
	{
		print_err_msg("SISO::demapper: constellation not initialized");		
		return;
	}
	if (MAP_metric=="logMAP")
		demodulator_logMAP(extrinsic_data, rec_sig, apriori_data);
	else if (MAP_metric=="maxlogMAP")
		demodulator_maxlogMAP(extrinsic_data, rec_sig, apriori_data);
	else
		print_err_msg("SISO::demapper: unknown MAP metric. The MAP metric should be either logMAP or maxlogMAP");
}

inline void SISO::demapper(itpp::vec &extrinsic_data, const itpp::cmat &rec_sig, const itpp::vec &apriori_data)
{
	if (c_impulse_response.size()==0)
	{
		print_err_msg("SISO::demapper: channel impulse response not initialized");
		return;
	}
	if ((ST_gen1.size()==0) || (ST_gen2.size()==0))
	{
		print_err_msg("SISO::demapper: Space-Time generator polynomials not initialized");
		return;
	}
	if ((constellation.size()==0) || (bin_constellation.size()==0))
	{
		print_err_msg("SISO::demapper: constellation not initialized");
		return;
	}
	
    if (demapper_method=="Hassibi_maxlogMAP")
        Hassibi_maxlogMAP(extrinsic_data, rec_sig, apriori_data);
    else if (demapper_method=="GA")
        GA(extrinsic_data, rec_sig, apriori_data);
    else if (demapper_method=="sGA")
        sGA(extrinsic_data, rec_sig, apriori_data);
    else if (demapper_method=="mmsePIC")
        mmsePIC(extrinsic_data, rec_sig, apriori_data);
    else if (demapper_method=="zfPIC")
        zfPIC(extrinsic_data, rec_sig, apriori_data);
    else if (demapper_method=="Alamouti_maxlogMAP")
        Alamouti_maxlogMAP(extrinsic_data, rec_sig, apriori_data);
    else
		print_err_msg("SISO::demapper: unknown demapper method. The demapper method should be either Hassibi_maxlogMAP, GA, sGA, mmsePIC, zfPIC or Alamouti_maxlogMAP");
}

void SISO::Hassibi_maxlogMAP(itpp::vec &extrinsic_data, const itpp::cmat &rec_sig, const itpp::vec &apriori_data)
//maxlogMAP algorithm for ST block codes using Hassibi's model
{
    //general parameters
    int nb_subblocks = rec_sig.rows()/block_duration;//number of subblocks of ST matrices/int period
    double N0 = 2*sigma2;//noise DSP
    int nb_all_symb = itpp::pow2i(nb_bits_symb*symbols_block);//nb. of all possible input symbols as a binary vector
    double nom,denom;//nominator and denominator of extrinsic information
    itpp::bvec bin_frame(nb_bits_symb*symbols_block);//binary frame at channel input
    itpp::bmat mat_bin_frame(nb_bits_symb, symbols_block);
    itpp::vec symb_frame_eq(2*symbols_block);//frame of symbols at equivalent channel input
    double temp;
    itpp::mat H_eq(2*nb_rec_ant*block_duration,2*symbols_block);//equivalent channel matrix
    itpp::vec x_eq(2*block_duration*nb_rec_ant);//equivalent received signal
    register int ns,q,nb,n,k;
	int index;
	extrinsic_data.set_size(nb_bits_symb*nb_subblocks*symbols_block);
    //main loop
    for (ns=0;ns<nb_subblocks;ns++)//for each subblock
    {
        //find equivalent channel matrix
        EquivCh(H_eq, c_impulse_response.get_col(ns));
        //find equivalent received signal
        EquivRecSig(x_eq, rec_sig(ns*block_duration,(ns+1)*block_duration-1,0,nb_rec_ant-1));
        //compute the LLR of each bit in a frame of symbols_block symbols
        for (q=0;q<symbols_block;q++)//for each symbol in a subblock
        {
            for (nb=0;nb<nb_bits_symb;nb++)//for a given bit try all possible sollutions for the input symbol vector
            {
                nom = -INFINITY;
                denom = -INFINITY;
                for (n=0;n<nb_all_symb;n++)//all possible symbols
                {
                    bin_frame = itpp::dec2bin(nb_bits_symb*symbols_block, n);
                    mat_bin_frame = itpp::reshape(bin_frame, nb_bits_symb, symbols_block);
                    for (k=0;k<symbols_block;k++)
                    {
                        symb_frame_eq(2*k) = constellation(itpp::bin2dec(mat_bin_frame.get_col(k))).real();
                        symb_frame_eq(1+2*k) = constellation(itpp::bin2dec(mat_bin_frame.get_col(k))).imag();
                    }
                    temp = -itpp::sum_sqr(x_eq-H_eq*symb_frame_eq)/N0+\
                           itpp::to_vec(bin_frame)*apriori_data.mid(ns*nb_bits_symb*symbols_block,nb_bits_symb*symbols_block);
                    if (bin_frame(nb+q*nb_bits_symb))
                        nom = std::max(nom, temp);
                    else
                        denom = std::max(denom, temp);
                }
                index = nb+q*nb_bits_symb+ns*nb_bits_symb*symbols_block;
                extrinsic_data(index) = (nom-denom)-apriori_data(index);
            }//bits/symbol
        }//symbols/subblock
    }//subblocks
}

void SISO::GA(itpp::vec &extrinsic_data, const itpp::cmat &rec_sig, const itpp::vec &apriori_data)
// Gaussian Approximation algorithm for ST codes using Hassibi's model
{
    //general parameters
    int nb_subblocks = rec_sig.rows()/block_duration;//number of subblocks
    int half_nb_bits_symb = nb_bits_symb/2;
    int half_len = itpp::pow2i(half_nb_bits_symb);//number of values of real(imaginary) part

    //correspondence between real and imaginary part of symbols and their binary representations
    int select_half;
    itpp::vec re_part;
    itpp::bmat re_bin_part;
    itpp::vec im_part;
    itpp::bmat im_bin_part;
    find_half_const(select_half, re_part, re_bin_part, im_part, im_bin_part);

    //equivalent channel
    itpp::mat H_eq(2*nb_rec_ant*block_duration,2*symbols_block);
    itpp::vec E_re_s(symbols_block);
    itpp::vec E_im_s(symbols_block);
    itpp::vec Var_re_s(symbols_block);
    itpp::vec Var_im_s(symbols_block);
    itpp::vec Ey(2*block_duration*nb_rec_ant);
    itpp::mat Cy(2*block_duration*nb_rec_ant,2*block_duration*nb_rec_ant);
    itpp::mat Cy_inv(2*block_duration*nb_rec_ant,2*block_duration*nb_rec_ant);
    itpp::vec x_eq(2*block_duration*nb_rec_ant);
    itpp::vec EZeta(2*block_duration*nb_rec_ant);
    itpp::mat CZeta_inv(2*block_duration*nb_rec_ant,2*block_duration*nb_rec_ant);
    double nom,denom;
    double temp;
    register int ns,q,k,p,cs;
    int index;
    extrinsic_data.set_size(nb_bits_symb*nb_subblocks*symbols_block);
    for (ns=0;ns<nb_subblocks;ns++)//subblock by subblock
    {
        //mean and variance of real and imaginary parts of emitted symbols
        E_re_s.zeros();
        E_im_s.zeros();
        Var_re_s.zeros();
        Var_im_s.zeros();
        for (q=0;q<symbols_block;q++)
        {
        	index = q*nb_bits_symb+ns*symbols_block*nb_bits_symb;
            for (k=0;k<half_len;k++)
            {
                E_re_s(q) += re_part(k)*itpp::prod(itpp::elem_div(exp(itpp::elem_mult(itpp::to_vec(re_bin_part.get_row(k)),\
                                                      apriori_data.mid(select_half*half_nb_bits_symb+index,half_nb_bits_symb))),\
                                                      1+exp(apriori_data.mid(select_half*half_nb_bits_symb+index,half_nb_bits_symb))));
                E_im_s(q) += im_part(k)*itpp::prod(itpp::elem_div(exp(itpp::elem_mult(itpp::to_vec(im_bin_part.get_row(k)),\
                                                      apriori_data.mid((1-select_half)*half_nb_bits_symb+index,half_nb_bits_symb))),\
                                                      1+exp(apriori_data.mid((1-select_half)*half_nb_bits_symb+index,half_nb_bits_symb))));
                Var_re_s(q) += itpp::sqr(re_part(k))*itpp::prod(itpp::elem_div(exp(itpp::elem_mult(itpp::to_vec(re_bin_part.get_row(k)),\
                                                    apriori_data.mid(select_half*half_nb_bits_symb+index,half_nb_bits_symb))),\
                                                    1+exp(apriori_data.mid(select_half*half_nb_bits_symb+index,half_nb_bits_symb))));
                Var_im_s(q) += itpp::sqr(im_part(k))*itpp::prod(itpp::elem_div(exp(itpp::elem_mult(itpp::to_vec(im_bin_part.get_row(k)),\
                                                    apriori_data.mid((1-select_half)*half_nb_bits_symb+index,half_nb_bits_symb))),\
                                                    1+exp(apriori_data.mid((1-select_half)*half_nb_bits_symb+index,half_nb_bits_symb))));
            }
            Var_re_s(q) -= itpp::sqr(E_re_s(q));
            Var_im_s(q) -= itpp::sqr(E_im_s(q));
        }

        //find equivalent channel
        EquivCh(H_eq, c_impulse_response.get_col(ns));

        //compute E[y] and Cov[y]
        Ey.zeros();
        Cy = sigma2*itpp::eye(2*block_duration*nb_rec_ant);
        for (q=0;q<symbols_block;q++)
        {
            //real & imaginary
            Ey += (H_eq.get_col(2*q)*E_re_s(q)+H_eq.get_col(1+2*q)*E_im_s(q));
            Cy += (itpp::outer_product(H_eq.get_col(2*q), H_eq.get_col(2*q)*Var_re_s(q))+\
            	itpp::outer_product(H_eq.get_col(1+2*q), H_eq.get_col(1+2*q)*Var_im_s(q)));
        }

        //inverse of Cov[y]
        Cy_inv = itpp::inv(Cy);

        //find equivalent received signal
        EquivRecSig(x_eq, rec_sig(ns*block_duration,(ns+1)*block_duration-1,0,nb_rec_ant-1));

        //compute extrinsic information of coded bits
        for (q=0;q<symbols_block;q++)
        {
            //real part
            EZeta = Ey-H_eq.get_col(2*q)*E_re_s(q);
            CZeta_inv = Cy_inv+itpp::outer_product(Cy_inv*\
                                ((Var_re_s(q)/(1-(((H_eq.get_col(2*q)).transpose()*Cy_inv)*(H_eq.get_col(2*q)*Var_re_s(q)))(0)))*\
                                 H_eq.get_col(2*q)), Cy_inv.transpose()*H_eq.get_col(2*q));
			index = select_half*half_nb_bits_symb+q*nb_bits_symb+ns*symbols_block*nb_bits_symb;			                                 
            for (p=0;p<half_nb_bits_symb;p++)
            {
                nom = -INFINITY;
                denom = -INFINITY;                
                for (cs=0;cs<half_len;cs++)
                {
                    temp = -0.5*((x_eq-H_eq.get_col(2*q)*re_part(cs)-EZeta).transpose()*CZeta_inv*(x_eq-H_eq.get_col(2*q)*re_part(cs)-EZeta))(0)+\
                           itpp::to_vec(re_bin_part.get_row(cs))*apriori_data.mid(index,half_nb_bits_symb);
                    if (re_bin_part(cs,p))
                        nom = std::max(nom, temp);
                    else
                        denom = std::max(denom, temp);
                }
                extrinsic_data(index+p) = (nom-denom)-apriori_data(index+p);
            }
            //imaginary part
            EZeta = Ey-H_eq.get_col(1+2*q)*E_im_s(q);
            CZeta_inv = Cy_inv+itpp::outer_product(Cy_inv*\
                                ((Var_im_s(q)/(1-(((H_eq.get_col(1+2*q)).transpose()*Cy_inv)*(H_eq.get_col(1+2*q)*Var_im_s(q)))(0)))*\
                                 H_eq.get_col(1+2*q)), Cy_inv.transpose()*H_eq.get_col(1+2*q));
			index = (1-select_half)*half_nb_bits_symb+q*nb_bits_symb+ns*symbols_block*nb_bits_symb;                                 
            for (p=0;p<half_nb_bits_symb;p++)
            {
                nom = -INFINITY;
                denom = -INFINITY;                
                for (cs=0;cs<half_len;cs++)
                {
                    temp = -0.5*((x_eq-H_eq.get_col(1+2*q)*im_part(cs)-EZeta).transpose()*CZeta_inv*(x_eq-H_eq.get_col(1+2*q)*im_part(cs)-EZeta))(0)+\
                           itpp::to_vec(im_bin_part.get_row(cs))*apriori_data.mid(index,half_nb_bits_symb);
                    if (im_bin_part(cs,p))
                        nom = std::max(nom, temp);
                    else
                        denom = std::max(denom, temp);
                }
                extrinsic_data(index+p) = (nom-denom)-apriori_data(index+p);
            }
        }
    }//subblock by subblock
}

void SISO::sGA(itpp::vec &extrinsic_data, const itpp::cmat &rec_sig, const itpp::vec &apriori_data)
//simplified Gaussian Approximation algorithm for ST codes using Hassibi's model
{
    //general parameters
    int nb_subblocks = (int)(rec_sig.rows()/block_duration);//number of subblocks
    int half_nb_bits_symb = (int)(nb_bits_symb/2);
    int half_len = itpp::pow2i(half_nb_bits_symb);//number of values of real(imaginary) part

    //correspondence between real and imaginary part of symbols and their binary representations
    int select_half;
    itpp::vec re_part;
    itpp::bmat re_bin_part;
    itpp::vec im_part;
    itpp::bmat im_bin_part;
    find_half_const(select_half, re_part, re_bin_part, im_part, im_bin_part);

    //equivalent channel
    itpp::mat H_eq(2*nb_rec_ant*block_duration,2*symbols_block);

    itpp::vec E_re_s(symbols_block);
    itpp::vec E_im_s(symbols_block);
    itpp::vec Var_re_s(symbols_block);
    itpp::vec Var_im_s(symbols_block);
    itpp::vec Ey(2*block_duration*nb_rec_ant);
    itpp::mat Cy(2*block_duration*nb_rec_ant,2*block_duration*nb_rec_ant);
    itpp::vec x_eq(2*block_duration*nb_rec_ant);
    itpp::vec EZeta(2*block_duration*nb_rec_ant);
    itpp::vec CZeta(2*block_duration*nb_rec_ant);
    double nom,denom;
    double temp;
    register int ns,q,k,p,cs;
    int index;
    extrinsic_data.set_size(nb_bits_symb*nb_subblocks*symbols_block);
    for (ns=0;ns<nb_subblocks;ns++)//subblock by subblock
    {
        //mean and variance of real and imaginary parts of emitted symbols
        E_re_s.zeros();
        E_im_s.zeros();
        Var_re_s.zeros();
        Var_im_s.zeros();
        for (q=0;q<symbols_block;q++)
        {
        	index = q*nb_bits_symb+ns*symbols_block*nb_bits_symb;
            for (k=0;k<half_len;k++)
            {
                E_re_s(q) += re_part(k)*itpp::prod(itpp::elem_div(exp(itpp::elem_mult(itpp::to_vec(re_bin_part.get_row(k)),\
                                                      apriori_data.mid(select_half*half_nb_bits_symb+index,half_nb_bits_symb))),\
                                                      1+exp(apriori_data.mid(select_half*half_nb_bits_symb+index,half_nb_bits_symb))));
                E_im_s(q) += im_part(k)*itpp::prod(itpp::elem_div(exp(itpp::elem_mult(itpp::to_vec(im_bin_part.get_row(k)),\
                                                      apriori_data.mid((1-select_half)*half_nb_bits_symb+index,half_nb_bits_symb))),\
                                                      1+exp(apriori_data.mid((1-select_half)*half_nb_bits_symb+index,half_nb_bits_symb))));
                Var_re_s(q) += itpp::sqr(re_part(k))*itpp::prod(itpp::elem_div(exp(itpp::elem_mult(itpp::to_vec(re_bin_part.get_row(k)),\
                                                    apriori_data.mid(select_half*half_nb_bits_symb+index,half_nb_bits_symb))),\
                                                    1+exp(apriori_data.mid(select_half*half_nb_bits_symb+index,half_nb_bits_symb))));
                Var_im_s(q) += itpp::sqr(im_part(k))*itpp::prod(itpp::elem_div(exp(itpp::elem_mult(itpp::to_vec(im_bin_part.get_row(k)),\
                                                    apriori_data.mid((1-select_half)*half_nb_bits_symb+index,half_nb_bits_symb))),\
                                                    1+exp(apriori_data.mid((1-select_half)*half_nb_bits_symb+index,half_nb_bits_symb))));
            }
            Var_re_s(q) -= itpp::sqr(E_re_s(q));
            Var_im_s(q) -= itpp::sqr(E_im_s(q));
        }

        //find equivalent channel
        EquivCh(H_eq, c_impulse_response.get_col(ns));

        //compute E[y] and Cov[y]
        Ey.zeros();
        Cy = sigma2*itpp::eye(2*block_duration*nb_rec_ant);
        for (q=0;q<symbols_block;q++)
        {
            //real & imaginary
            Ey += (H_eq.get_col(2*q)*E_re_s(q)+H_eq.get_col(1+2*q)*E_im_s(q));
            Cy += (itpp::outer_product(H_eq.get_col(2*q), H_eq.get_col(2*q)*Var_re_s(q))+\
            	 itpp::outer_product(H_eq.get_col(1+2*q), H_eq.get_col(1+2*q)*Var_im_s(q)));
        }

        //find equivalent received signal
        EquivRecSig(x_eq, rec_sig(ns*block_duration,(ns+1)*block_duration-1,0,nb_rec_ant-1));

        //compute extrinsic INFINITYormation of coded bits
        for (q=0;q<symbols_block;q++)
        {
            //real part
            EZeta = Ey-H_eq.get_col(2*q)*E_re_s(q);
            CZeta = diag(Cy-itpp::outer_product(H_eq.get_col(2*q), H_eq.get_col(2*q)*Var_re_s(q)));
            index = select_half*half_nb_bits_symb+q*nb_bits_symb+ns*symbols_block*nb_bits_symb;
            for (p=0;p<half_nb_bits_symb;p++)
            {
                nom = -INFINITY;
                denom = -INFINITY;
                for (cs=0;cs<half_len;cs++)
                {
                    temp = -0.5*itpp::sum(itpp::elem_div(sqr(x_eq-H_eq.get_col(2*q)*re_part(cs)-EZeta), CZeta))+\
                           itpp::to_vec(re_bin_part.get_row(cs))*apriori_data.mid(index,half_nb_bits_symb);
                    if (re_bin_part(cs,p))
                        nom = std::max(nom, temp);
                    else
                        denom = std::max(denom, temp);
                }
                extrinsic_data(index+p) = (nom-denom)-apriori_data(index+p);
            }
            //imaginary part
            EZeta = Ey-H_eq.get_col(1+2*q)*E_im_s(q);
            CZeta = itpp::diag(Cy-itpp::outer_product(H_eq.get_col(1+2*q), H_eq.get_col(1+2*q)*Var_im_s(q)));
            index = (1-select_half)*half_nb_bits_symb+q*nb_bits_symb+ns*symbols_block*nb_bits_symb;
            for (p=0;p<half_nb_bits_symb;p++)
            {
                nom = -INFINITY;
                denom = -INFINITY;
                for (cs=0;cs<half_len;cs++)
                {
                    temp = -0.5*itpp::sum(itpp::elem_div(sqr(x_eq-H_eq.get_col(1+2*q)*im_part(cs)-EZeta), CZeta))+\
                           itpp::to_vec(im_bin_part.get_row(cs))*apriori_data.mid(index,half_nb_bits_symb);
                    if (im_bin_part(cs,p))
                        nom = std::max(nom, temp);
                    else
                        denom = std::max(denom, temp);
                }
                extrinsic_data(index+p) = (nom-denom)-apriori_data(index+p);
            }
        }
    }//subblock by subblock
}

void SISO::mmsePIC(itpp::vec &extrinsic_data, const itpp::cmat &rec_sig, const itpp::vec &apriori_data)
//MMSE Parallel Interference Canceller
{
    //general parameters
    int nb_subblocks = rec_sig.rows()/block_duration;//number of subblocks
    int half_nb_bits_symb = nb_bits_symb/2;
    int half_const_len = itpp::pow2i(half_nb_bits_symb);
    int nb_bits_subblock = nb_bits_symb*symbols_block;//number of coded bits in an ST block
    itpp::vec Es(2*symbols_block);
    itpp::vec Vs(2*symbols_block);
    itpp::mat H_eq(2*nb_rec_ant*block_duration,2*symbols_block);
    itpp::mat K(2*nb_rec_ant*block_duration,2*nb_rec_ant*block_duration);
    itpp::mat K_inv(2*nb_rec_ant*block_duration,2*nb_rec_ant*block_duration);
    itpp::vec x_eq(2*nb_rec_ant*block_duration);
    itpp::vec interf(2*symbols_block);
    itpp::vec temp(2*nb_rec_ant*block_duration);
    itpp::vec w(2*nb_rec_ant*block_duration);//filter impulse response
    double s_tilde;
    double mu_res;
    double sigma2_res;
    double nom,denom;
    double tmp;
    register int ns,q,k,s;
    int index;

    //correspondence between real and imaginary part of symbols and their binary representations
    int select_half;
    itpp::vec re_part;
    itpp::bmat re_bin_part;
    itpp::vec im_part;
    itpp::bmat im_bin_part;
    find_half_const(select_half, re_part, re_bin_part, im_part, im_bin_part);
    double part_var = 1/(double)(2*nb_em_ant);//real and imaginary part variance

	extrinsic_data.set_size(nb_bits_symb*nb_subblocks*symbols_block);
    for (ns=0;ns<nb_subblocks;ns++)//compute block by block
    {
        //mean and variance of real and imaginary parts of emitted symbols
        Es.zeros();
        Vs.zeros();
        for (q=0;q<symbols_block;q++)
        {
        	index = q*nb_bits_symb+ns*symbols_block*nb_bits_symb;
            for (k=0;k<half_const_len;k++)
            {
                Es(2*q) += re_part(k)*itpp::prod(itpp::elem_div(exp(itpp::elem_mult(itpp::to_vec(re_bin_part.get_row(k)), \
                                                        apriori_data.mid(select_half*half_nb_bits_symb+index, half_nb_bits_symb))), \
                                                    (1+exp(apriori_data.mid(select_half*half_nb_bits_symb+index, half_nb_bits_symb)))));
                Es(1+2*q) += im_part(k)*itpp::prod(itpp::elem_div(exp(itpp::elem_mult(itpp::to_vec(im_bin_part.get_row(k)), \
                                                      apriori_data.mid((1-select_half)*half_nb_bits_symb+index, half_nb_bits_symb))), \
                                                      (1+exp(apriori_data.mid((1-select_half)*half_nb_bits_symb+index, half_nb_bits_symb)))));
                Vs(2*q) += itpp::sqr(re_part(k))*itpp::prod(itpp::elem_div(exp(itpp::elem_mult(itpp::to_vec(re_bin_part.get_row(k)), \
                                                apriori_data.mid(select_half*half_nb_bits_symb+index, half_nb_bits_symb))), \
                                                (1+exp(apriori_data.mid(select_half*half_nb_bits_symb+index, half_nb_bits_symb)))));
                Vs(1+2*q) += itpp::sqr(im_part(k))*itpp::prod(itpp::elem_div(exp(itpp::elem_mult(itpp::to_vec(im_bin_part.get_row(k)), \
                                                  apriori_data.mid((1-select_half)*half_nb_bits_symb+index, half_nb_bits_symb))), \
                                                  (1+exp(apriori_data.mid((1-select_half)*half_nb_bits_symb+index, half_nb_bits_symb)))));
            }
            Vs(2*q) -= itpp::sqr(Es(2*q));
            Vs(1+2*q) -= itpp::sqr(Es(1+2*q));
        }

        //find equivalent channel matrix
        EquivCh(H_eq, c_impulse_response.get_col(ns));
        //compute invariant inverse
        K = H_eq*diag(Vs)*H_eq.transpose()+sigma2*itpp::eye(2*block_duration*nb_rec_ant);
        K_inv = itpp::inv(K);
        //find equivalent received signal
        EquivRecSig(x_eq, rec_sig(ns*block_duration,(ns+1)*block_duration-1,0,nb_rec_ant-1));
        for (q=0;q<symbols_block;q++)//symbols/block
        {
            //compute the extrinsic information of coded bits
            //real part
            //IC + filtering (real and imaginary parts of one symbol)
            interf = Es;
            interf(2*q) = 0;//this is the symbol to recover
            temp = H_eq.get_col(2*q);
			w = (part_var*temp.transpose())*(K_inv-itpp::outer_product(K_inv*(((part_var-Vs(2*q))/ \
                                                    (1+((temp.transpose()*K_inv)*(temp*(part_var-Vs(2*q))))(0)))*temp), K_inv.transpose()*temp));                                                    
            s_tilde = w*(x_eq-H_eq*interf);
            mu_res = w*temp;//mean of the filtered signal
            index = select_half*half_nb_bits_symb+nb_bits_symb*q+ns*nb_bits_subblock;
            for (k=0;k<half_nb_bits_symb;k++)
            {
                nom = -INFINITY;
                denom = -INFINITY;
                for (s=0;s<half_const_len;s++)
                {
                    sigma2_res = ((w.transpose()*(K+itpp::outer_product(temp, temp*(itpp::sqr(re_part(s))-Vs(2*q)))))*w)(0)-itpp::sqr(re_part(s)*mu_res);//variance of the filtered signal
                    tmp = -itpp::sqr(s_tilde-mu_res*re_part(s))/(2*sigma2_res)+ \
                          itpp::to_vec(re_bin_part.get_row(s))*apriori_data.mid(index, half_nb_bits_symb);
                    if (re_bin_part(s,k))
                        nom = std::max(nom, tmp);
                    else
                        denom = std::max(denom, tmp);
                }
                extrinsic_data(index+k) = (nom-denom)-apriori_data(index+k);
            }
            //end real part
            //imaginary part
            //IC + filtering (real and imaginary parts of one symbol)
            interf = Es;
            interf(2*q+1) = 0;//this is the symbol to recover
            temp = H_eq.get_col(2*q+1);
			w = (part_var*temp.transpose())*(K_inv-itpp::outer_product(K_inv*(((part_var-Vs(1+2*q))/ \
                                                    (1+((temp.transpose()*K_inv)*(temp*(part_var-Vs(1+2*q))))(0)))*temp), K_inv.transpose()*temp));                                                    
            s_tilde = w*(x_eq-H_eq*interf);
            mu_res = w*temp;//mean of the filtered signal
            index = (1-select_half)*half_nb_bits_symb+nb_bits_symb*q+ns*nb_bits_subblock;
            for (k=0;k<half_nb_bits_symb;k++)
            {
                nom = -INFINITY;
                denom = -INFINITY;
                for (s=0;s<half_const_len;s++)
                {
                    sigma2_res = ((w.transpose()*(K+itpp::outer_product(temp, temp*(itpp::sqr(im_part(s))-Vs(1+2*q)))))*w)(0)-itpp::sqr(im_part(s)*mu_res);//variance of the filtered signal
                    tmp = -itpp::sqr(s_tilde-mu_res*im_part(s))/(2*sigma2_res)+ \
                          itpp::to_vec(im_bin_part.get_row(s))*apriori_data.mid(index, half_nb_bits_symb);
                    if (im_bin_part(s,k))
                        nom = std::max(nom, tmp);
                    else
                        denom = std::max(denom, tmp);
                }
                extrinsic_data(index+k) = (nom-denom)-apriori_data(index+k);
            }
            //end imaginary part
        }//symbols/block
    }//block by block
}

void SISO::zfPIC(itpp::vec &extrinsic_data, const itpp::cmat &rec_sig, const itpp::vec &apriori_data)
//ZF Parallel Interference Canceller
{
    //general parameters
    int nb_subblocks = rec_sig.rows()/block_duration;//number of subblocks
    int half_nb_bits_symb = nb_bits_symb/2;
    int half_const_len = itpp::pow2i(half_nb_bits_symb);
    int nb_bits_subblock = nb_bits_symb*symbols_block;//number of coded bits in an ST block
    itpp::vec Es(2*symbols_block);
    itpp::vec Vs(2*symbols_block);
    itpp::mat H_eq(2*nb_rec_ant*block_duration,2*symbols_block);
    itpp::mat K(2*nb_rec_ant*block_duration,2*nb_rec_ant*block_duration);
    itpp::vec x_eq(2*nb_rec_ant*block_duration);
    itpp::vec interf(2*symbols_block);
    itpp::vec temp(2*nb_rec_ant*block_duration);
    itpp::vec w(2*nb_rec_ant*block_duration);//filter impulse response
    double s_tilde;
    double mu_res;
    double sigma2_res;
    double nom,denom;
    double tmp;
    register int ns,q,k,s;
    int index;

    //correspondence between real and imaginary part of symbols and their binary representations
    int select_half;
    itpp::vec re_part;
    itpp::bmat re_bin_part;
    itpp::vec im_part;
    itpp::bmat im_bin_part;
    find_half_const(select_half, re_part, re_bin_part, im_part, im_bin_part);

	extrinsic_data.set_size(nb_bits_symb*nb_subblocks*symbols_block);
    for (ns=0;ns<nb_subblocks;ns++)//compute block by block
    {
        //mean and variance of real and imaginary parts of emitted symbols
        Es.zeros();
        Vs.zeros();
        for (q=0;q<symbols_block;q++)
        {
        	index = q*nb_bits_symb+ns*symbols_block*nb_bits_symb;
            for (k=0;k<half_const_len;k++)
            {
                Es(2*q) += re_part(k)*itpp::prod(itpp::elem_div(exp(itpp::elem_mult(itpp::to_vec(re_bin_part.get_row(k)), \
                                                        apriori_data.mid(select_half*half_nb_bits_symb+index, half_nb_bits_symb))), \
                                                    (1+exp(apriori_data.mid(select_half*half_nb_bits_symb+index, half_nb_bits_symb)))));
                Es(1+2*q) += im_part(k)*itpp::prod(itpp::elem_div(exp(itpp::elem_mult(itpp::to_vec(im_bin_part.get_row(k)), \
                                                      apriori_data.mid((1-select_half)*half_nb_bits_symb+index, half_nb_bits_symb))), \
                                                      (1+exp(apriori_data.mid((1-select_half)*half_nb_bits_symb+index, half_nb_bits_symb)))));
                Vs(2*q) += itpp::sqr(re_part(k))*itpp::prod(itpp::elem_div(exp(itpp::elem_mult(itpp::to_vec(re_bin_part.get_row(k)), \
                                                apriori_data.mid(select_half*half_nb_bits_symb+index, half_nb_bits_symb))), \
                                                (1+exp(apriori_data.mid(select_half*half_nb_bits_symb+index, half_nb_bits_symb)))));
                Vs(1+2*q) += itpp::sqr(im_part(k))*itpp::prod(itpp::elem_div(exp(itpp::elem_mult(to_vec(im_bin_part.get_row(k)), \
                                                  apriori_data.mid((1-select_half)*half_nb_bits_symb+index, half_nb_bits_symb))), \
                                                  (1+exp(apriori_data.mid((1-select_half)*half_nb_bits_symb+index, half_nb_bits_symb)))));
            }
            Vs(2*q) -= itpp::sqr(Es(2*q));
            Vs(1+2*q) -= itpp::sqr(Es(1+2*q));
        }

        //find equivalent channel matrix
        EquivCh(H_eq, c_impulse_response.get_col(ns));
        //compute invariant inverse
        K = H_eq*itpp::diag(Vs)*H_eq.transpose()+sigma2*itpp::eye(2*block_duration*nb_rec_ant);
        //find equivalent received signal
        EquivRecSig(x_eq, rec_sig(ns*block_duration,(ns+1)*block_duration-1,0,nb_rec_ant-1));
        for (q=0;q<symbols_block;q++)//symbols/block
        {
            //compute the extrinsic information of coded bits
            //real part
            //IC + filtering (real and imaginary parts of one symbol)
            interf = Es;
            interf(2*q) = 0;//this is the symbol to recover
            temp = H_eq.get_col(2*q);
            w = temp/(temp*temp);//filter impulse response
            s_tilde = w*(x_eq-H_eq*interf);
            mu_res = w*temp;//mean of the filtered signal
            index = select_half*half_nb_bits_symb+nb_bits_symb*q+ns*nb_bits_subblock;
            for (k=0;k<half_nb_bits_symb;k++)
            {
                nom = -INFINITY;
                denom = -INFINITY;
                for (s=0;s<half_const_len;s++)
                {
                    sigma2_res = ((w.transpose()*(K+itpp::outer_product(temp, temp*(itpp::sqr(re_part(s))-Vs(2*q)))))*w)(0)-itpp::sqr(re_part(s)*mu_res);
                    tmp = -itpp::sqr(s_tilde-mu_res*re_part(s))/(2*sigma2_res)+ \
                          itpp::to_vec(re_bin_part.get_row(s))*apriori_data.mid(index, half_nb_bits_symb);
                    if (re_bin_part(s,k))
                        nom = std::max(nom, tmp);
                    else
                        denom = std::max(denom, tmp);
                }
                extrinsic_data(index+k) = (nom-denom)-apriori_data(index+k);
            }
            //end real part
            //imaginary part
            //IC + filtering (real and imaginary parts of one symbol)
            interf = Es;
            interf(2*q+1) = 0;//this is the symbol to recover
            temp = H_eq.get_col(2*q+1);
            w = temp/(temp*temp);//filter impulse response
            s_tilde = w*(x_eq-H_eq*interf);
            mu_res = w*temp;//mean of the filtered signal
            index = (1-select_half)*half_nb_bits_symb+nb_bits_symb*q+ns*nb_bits_subblock;
            for (k=0;k<half_nb_bits_symb;k++)
            {
                nom = -INFINITY;
                denom = -INFINITY;
                for (s=0;s<half_const_len;s++)
                {
                    sigma2_res = ((w.transpose()*(K+itpp::outer_product(temp, temp*(itpp::sqr(im_part(s))-Vs(1+2*q)))))*w)(0)-itpp::sqr(im_part(s)*mu_res);
                    tmp = -itpp::sqr(s_tilde-mu_res*im_part(s))/(2*sigma2_res)+ \
                          itpp::to_vec(im_bin_part.get_row(s))*apriori_data.mid(index, half_nb_bits_symb);
                    if (im_bin_part(s,k))
                        nom = std::max(nom, tmp);
                    else
                        denom = std::max(denom, tmp);
                }
                extrinsic_data(index+k) = (nom-denom)-apriori_data(index+k);
            }
            //end imaginary part
        }//symbols/block
    }//block by block
}

void SISO::Alamouti_maxlogMAP(itpp::vec &extrinsic_data, const itpp::cmat &rec_sig, const itpp::vec &apriori_data)
//maxlogMAP algorithm for Alamouti ST code
{
    //matched filter
    int int_len = apriori_data.length();//interleaver length
    int nb_symb = (int)(int_len/nb_bits_symb);//number of symbols/block
    itpp::cvec comb_sig(nb_symb);
    comb_sig.zeros();
    itpp::cmat conj_H = itpp::conj(c_impulse_response);
    itpp::cmat conj_X = itpp::conj(rec_sig);
    register int nr,n,cs;
    for (nr=0;nr<nb_rec_ant;nr++)
    {
        for (n=0;n<(nb_symb/2);n++)
        {
            comb_sig(2*n) += (conj_H(2*nr,n)*rec_sig(2*n,nr)+c_impulse_response(1+2*nr,n)*conj_X(1+2*n,nr));
            comb_sig(1+2*n) += (conj_H(1+2*nr,n)*rec_sig(2*n,nr)-c_impulse_response(2*nr,n)*conj_X(1+2*n,nr));
        }
    }

    //extrinsic information of coded bits
    int const_size = itpp::pow2i(nb_bits_symb);//constellation size
    double buffer;
    double nom,denom;
    double temp;
    int index;
    extrinsic_data.set_size(nb_bits_symb*nb_symb);
    for (n=0;n<nb_symb;n++)
    {
        buffer = itpp::sum_sqr(itpp::abs(c_impulse_response.get_col(n/2)));
        for (nr=0;nr<nb_bits_symb;nr++)
        {
            nom = -INFINITY;
            denom = -INFINITY;
            for (cs=0;cs<const_size;cs++)
            {
                temp = -itpp::sqr(comb_sig(n)-buffer*constellation(cs))/(2*buffer*sigma2)+ \
                       itpp::to_vec(bin_constellation.get_row(cs))*apriori_data.mid(n*nb_bits_symb, nb_bits_symb);
                if (bin_constellation(cs,nr))
                    nom = std::max(nom, temp);
                else
                    denom = std::max(denom, temp);
            }
            index = n*nb_bits_symb+nr;
            extrinsic_data(index) = (nom-denom)-apriori_data(index);//extrinsic information
        }
    }
}

void SISO::demodulator_logMAP(itpp::vec &extrinsic_data, const itpp::cvec &rec_sig, const itpp::vec &apriori_data)
/// logMAP demodulator
{
	int nb_symb = rec_sig.length();
	int const_size = itpp::pow2i(nb_bits_symb);
	double nom,denom,temp;
	register int k,i,cs;
	int index;
	extrinsic_data.set_size(nb_bits_symb*nb_symb);
    for(k=0;k<nb_symb;k++)
    {
        for(i=0;i<nb_bits_symb;i++)
        {
            nom = 0;
            denom = 0;
            for(cs=0;cs<const_size;cs++)
            {
                temp = -itpp::sqr(rec_sig(k)-c_impulse_response(0,k)*constellation(cs))/(2*sigma2)+\
                    itpp::to_vec(bin_constellation.get_row(cs))*apriori_data.mid(k*nb_bits_symb, nb_bits_symb);
                if(bin_constellation(cs,i))
                    nom += exp(temp);
                else
                    denom += exp(temp);
            }
            index = k*nb_bits_symb+i;
            extrinsic_data(index) = log(nom/denom)-apriori_data(index);//extrinsic information
        }
    }
}

void SISO::demodulator_maxlogMAP(itpp::vec &extrinsic_data, const itpp::cvec &rec_sig, const itpp::vec &apriori_data)
/// maxlogMAP demodulator
{
	int nb_symb = rec_sig.length();
	int const_size = itpp::pow2i(nb_bits_symb);
	double nom,denom,temp;
	register int k,i,cs;
	int index;
	extrinsic_data.set_size(nb_bits_symb*nb_symb);
    for(k=0;k<nb_symb;k++)
    {
        for(i=0;i<nb_bits_symb;i++)
        {
            nom = -INFINITY;
            denom = -INFINITY;
            for(cs=0;cs<const_size;cs++)
            {
                temp = -itpp::sqr(rec_sig(k)-c_impulse_response(0,k)*constellation(cs))/(2*sigma2)+\
                    itpp::to_vec(bin_constellation.get_row(cs))*apriori_data.mid(k*nb_bits_symb, nb_bits_symb);
                if(bin_constellation(cs,i))
                    nom = std::max(nom, temp);
                else
                    denom = std::max(denom, temp);
            }
            index = k*nb_bits_symb+i;
            extrinsic_data(index) = (nom-denom)-apriori_data(index);//extrinsic information
        }
    }
}

// threshold functions
inline double threshold(const double &x, const double &value)
{
    if ((x>value)||(x<-value))
        return (x>0?value:-value);
    return x;
}

itpp::vec threshold(const itpp::vec &in, const double &value)
{
    itpp::vec out(in.length());
    register int n;
    for (n=0;n<in.length();n++)
        out(n) = threshold(in(n), value);
    return out;
}

itpp::mat threshold(const itpp::mat &in, const double &value)
{
    itpp::mat out(in.rows(),in.cols());
    register int n;
    for (n=0;n<in.rows();n++)
        out.set_row(n, threshold(in.get_row(n), value));
    return out;
}

}//end turbo receivers workspace

#endif
