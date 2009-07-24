/** \file
 * 
 * \brief EXtrinsic Information Transfer (EXIT) chart class
 */

#ifndef EXIT_CLASS
#define EXIT_CLASS

#include "itpp/comm/modulator.h" //BPSK class for a priori information generation
#include "itpp/stat/histogram.h" //histogram class for mutual information computation

namespace tr 
{

/// EXtrinsic Information Transfer (%EXIT) chart
/** Computes the A priori Mutual Information assuming a Gaussian distribution of the a priori information and
 * the Extrinsic Mutual Information between the emitted bits and their extrinsic information
 *
 * Description:
 * - the a priori mutual information is computed using relation (14)
 * - the extrinsic mutual information is computed by estimating first the conditional Probability Density Functions (PDF),
 * given the emitted bits, and then numerically integrating according to relation (19)
 *
 * Reference: 
 * Stephan ten Brink, ''Convergence behavior of iteratively decoded parallel concatenated codes,`` 
 * IEEE Transactions on Communications, oct. 2001
 */ 
class EXIT
{
public:
	/// Computes the a priori mutual information
	/** It is assumed that the a priori information has a Gaussian distribution
	 */
	double apriori_mutual_info(const double &in_sigma2A, ///< variance of the a priori information
							   const double &lim=100 ///< [-lim,+lim] is the integration interval (theoretically it should be [-inf,+inf])
							   )
	{
    	sigma2A = in_sigma2A;
    	return double(1)-itpp::quad(&gaussian_fct, -lim, lim);
	};
	/// Generates a priori information assuming a Gaussian distribution of the a priori information
	/** The variance of the a priori information must be already initialized through EXIT::apriori_mutual_info function.
	 * The information generated in this way is used sometimes as intrinsic information at the SISO module input.
	 */
	itpp::vec generate_apriori_info(const itpp::bvec &bits)
	{
		itpp::BPSK bpsk;
		return (-sigma2A/2)*bpsk.modulate_bits(bits)+sqrt(sigma2A)*itpp::randn(bits.length());
	};
	/// Computes the extrinsic mutual information
	/** The conditional Probability Density Function (PDF) of the extrinsic information is estimated using the histogram of the
	 * extrinsic information and the knowledge of the emitted bits corresponding to the extrinsic information.
	 */
	double extrinsic_mutual_info(const itpp::vec &obs, ///< extrinsic information obtained from the SISO module output
								 const itpp::bvec &cond, ///< emitted bits corresponding to the extrinsic information
								 const int &N=100 ///< number of subintervals used to compute the histogram
								 );
private:	
	static double sigma2A;
	friend double itpp::quad(double (*f)(double), double a, double b, double tol);
	static double gaussian_fct(const double x)
	{
    	return (1.0/sqrt(sigma2A*itpp::m_2pi))*exp(-itpp::sqr(x-(sigma2A/2.0))/(2.0*sigma2A))*log2(1+exp(-x));
	};
};

double EXIT::sigma2A;//allocate memory for static member variable

double EXIT::extrinsic_mutual_info(const itpp::vec &obs, const itpp::bvec &cond, const int &N)
{
    //initialize histogram
    itpp::Histogram<double> hist(itpp::min(obs), itpp::max(obs), N);//common definition interval for both PDFs

    //conditional PDF knowing that a bit of 0 was emitted
    itpp::ivec idx = itpp::find(cond==itpp::bin(0));
    itpp::vec cond_obs = obs(idx);
    hist.reset();//start counting
    hist.update(cond_obs);
    itpp::vec left_pdf = hist.get_pdf();
    itpp::ivec left_int = itpp::find(left_pdf!=0);//integration interval for the left PDF

    //conditional PDF knowing that a bit of 1 was emitted
    idx = itpp::find(cond==itpp::bin(1));
    cond_obs = obs(idx);
    hist.reset();//restart counting
    hist.update(cond_obs);
    itpp::vec right_pdf = hist.get_pdf();
    itpp::ivec right_int = itpp::find(right_pdf!=0);//integration interval for the right PDF

    //mutual extrinsic information
    itpp::vec left_half = itpp::elem_mult(left_pdf(left_int), itpp::log2(itpp::elem_div(2.0*left_pdf(left_int), left_pdf(left_int)+right_pdf(left_int))));
    double IE = itpp::sum(left_half)-0.5*(left_half(0)+left_half(left_half.length()-1));//numerical integration
    itpp::vec right_half = itpp::elem_mult(right_pdf(right_int), itpp::log2(itpp::elem_div(2.0*right_pdf(right_int), left_pdf(right_int)+right_pdf(right_int))));
    IE += itpp::sum(right_half)-0.5*(right_half(0)+right_half(right_half.length()-1));//numerical integration
    IE *= 0.5;

    return IE;
}

}
#endif
