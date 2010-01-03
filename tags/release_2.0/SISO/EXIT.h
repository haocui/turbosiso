/** \file
 *
 * \brief Implementation of EXtrinsic Information Transfer (EXIT) chart class header
 */

#ifndef EXIT_H_
#define EXIT_H_

#include "itpp/itbase.h"
#include "itpp/comm/modulator.h" //BPSK class for a priori information generation
#include "itpp/stat/histogram.h" //histogram class for mutual information computation

namespace tr
{

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

}
#endif /* EXIT_H_ */
