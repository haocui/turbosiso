//element wise multiplication of the two vectors followed by xor sum of the resulting vector
inline bin xor_sum(const bvec &in1, const bvec &in2)
{
    bin out = 0;
    for(int n=0;n<in1.length();n++)
    {
        out ^= (in1(n)*in2(n));
    };
    return out;
}

//non systematic non recursive convolutional code
//class CC with members: set_generator, set_initial_state, etc.
//member function nsc should have two versions: nsc(in) and nsc(in, gen)
bvec nsc(const bvec &in, const bmat &gen)
{
    int r = gen.rows();//inverse of coding rate
    int reg_len = gen.cols()-1;//register length
    int nb_bits = in.length();

	bvec out(nb_bits*r);
    bvec reg(reg_len);
    bvec in_reg(reg_len+1);
    register int k,n;
    for(k=0;k<r;k++)
    {
        reg.zeros();
        for(n=0;n<nb_bits;n++)
        {
            in_reg(0) = in(n);
            in_reg.replace_mid(1, reg);
            out(k+r*n) = xor_sum(in_reg, gen.get_row(k));
            reg.shift_right(in(n));
        }
    }
    return out;
}

//precoder (recursive CC of rate 1) (operates on BPSK symbols)
vec precoder(vec in, bvec precoder_pol)
{
	int in_len = in.length();
	vec out(in_len);
	int reg_len = precoder_pol.length()-1;
	vec reg(reg_len);//precoder register
	reg.ones();//intial state
	ivec reg_idx = find(precoder_pol.right(reg_len));//feedback index (precoder_pol(0) is omitted since it represents the register input)
	register int n;
	for(n=0;n<in_len;n++)
	{
		out(n) = in(n)*prod(reg(reg_idx));
		reg.shift_right(out(n), 1);
	}
	return out;
}

//Kronecker operator for vectors (make a patch)
template <class Num_T>
Vec<Num_T> kron(const Vec<Num_T> &in, const Vec<Num_T> &pattern)
{
    int in_len = in.length();
    int pattern_len = pattern.length();
    Vec<Num_T> out(in_len*pattern_len);
    for (int n=0;n<in_len;n++)
        out.replace_mid(n*pattern_len, in(n)*pattern);
    return out;
}

/*
 *-------------------------------------------------------------------------------------------------------------
 * Computes the A priori Mutual Information supposing a Gaussian distribution of the a priori information and
 * the Extrinsic Mutual Information between the emitted BPSK symbols and their extrinsic information
 *-------------------------------------------------------------------------------------------------------------

 *Description:
 * - the a priori mutual information is computed using relation (14) [1]
 * - the extrinsic mutual information is computed by estimating first the conditional Probability Density Functions (PDF), 
 * given the emitted symbol -1 and +1, then numerically integrating according to relation (19) [1]

 *Reference: [1] Stephan ten Brink, "Convergence behavior of iteratively decoded parallel concatenated codes", 
 * IEEE Transactions on Communications, oct. 2001
 */

#include "itpp/stat/histogram.h"

double s2A;
double gaussian_fct(const double x)
{
	return (1/sqrt(s2A*m_2pi))*exp(-sqr(x-(s2A/2))/(2*s2A))*log2(1+exp(-x));
}

/*
 * A templated version of quad needed to accept a pointer to a member function
 
class EXIT_Chart
{
	public:
	double apriori_mutual_info(const double &sigma2A, const double lim=100);
	double extrinsic_mutual_info(const vec &obs, const vec &cond, const int N=100);
};
*/

inline double apriori_mutual_info(const double &sigma2A, const double lim=100)
{
	s2A = sigma2A;
	return 1-quad(&gaussian_fct, -lim, lim);
}

inline double extrinsic_mutual_info(const vec &obs, const vec &cond, const int N=100)
/*
 *inputs
 * obs - extrinsic information of the emitted symbols obtained at a SISO module output
 * cond - emitted BPSK symbols
 * N - number of intervals for the histogram used to approximate the conditional probability density function of the extrinsic information for a given emitted symbol (-1 or +1)
 * output
 *  mutual extrinsic information
 */
{
	//initialize histogram
	Histogram<double> hist(min(obs), max(obs), N);//common definition interval for both PDFs
	
	//conditional PDF knowing that a -1 was emitted
	ivec idx = find(cond==-1.0);
	vec cond_obs = obs(idx);
	hist.reset();//start counting
	hist.update(cond_obs);
	vec left_pdf = hist.get_pdf();
	ivec left_int = find(left_pdf!=0);//integration interval for left PDF

	//conditional PDF knowing that a +1 was emitted
	idx = find(cond==+1.0);
	cond_obs = obs(idx);
	hist.reset();//restart counting
	hist.update(cond_obs);
	vec right_pdf = hist.get_pdf();
	ivec right_int = find(right_pdf!=0);//integration interval for right PDF
	
	//mutual extrinsic information
	vec left_half = elem_mult(left_pdf(left_int), log2(elem_div(2*left_pdf(left_int), left_pdf(left_int)+right_pdf(left_int))));
	double IE = sum(left_half)-0.5*(left_half(0)+left_half(left_half.length()-1));//numerical integration
	vec right_half = elem_mult(right_pdf(right_int), log2(elem_div(2*right_pdf(right_int), left_pdf(right_int)+right_pdf(right_int))));
	IE += sum(right_half)-0.5*(right_half(0)+right_half(right_half.length()-1));//numerical integration
	IE *= 0.5;
	
	return IE;
}

