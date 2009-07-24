/** \file
 *
 * \brief Implementation of SISO modules for RSC codes
 */

#include "SISO.h"

namespace tr
{
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
#pragma omp parallel for private(k, j, cases, feedback, out, buffer)
		for (k=0;k<rsctrellis.numStates;k++)
		{
			cases = itpp::dec2bin(mem_len, k);
			//feedback
			feedback = (itpp::bin)n;
			for (j=1;j<(mem_len+1);j++)
			{
				feedback ^= (gen(0,j)*cases[j-1]);
			}
			//out
			out = feedback*gen(1,0);
			for (j=1;j<(mem_len+1);j++)
			{
				out ^= (gen(1,j)*cases[j-1]);
			}
			rsctrellis.PARout[k+n*rsctrellis.numStates] = (out?1.0:0.0);//parity bit
			rsctrellis.fm[k] = n^out;
			//shift
			for (j=mem_len-1;j>0;j--)
			{
				cases[j] = cases[j-1];
			}
			cases[0] = feedback;
			//next and previous state
			buffer = itpp::bin2dec(cases, true);
			rsctrellis.nextStates[k+n*rsctrellis.numStates] = buffer;//next state
			rsctrellis.prevStates[buffer+n*rsctrellis.numStates] = k;//previous state
		}
	}
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
#pragma omp parallel for private(n)
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
#pragma omp parallel for private(n,k)
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
		{
			continue;
		}
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
		{
			continue;
		}
		for (k=0;k<rsctrellis.numStates;k++)
		{
			B0[k+rsctrellis.numStates*n] -= A_mid[n+1];
			B1[k+rsctrellis.numStates*n] -= A_mid[n+1];
		}
	}

	//updated LLR for information bits
	extrinsic_data.set_size(N);
#pragma omp parallel for private(n, k, sum0, sum1)
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
#pragma omp parallel for private(n, k, sum0, sum1)
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
#pragma omp parallel for private(n)
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
#pragma omp parallel for private(n,k)
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
#pragma omp parallel for private(n, k, sum0, sum1)
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
#pragma omp parallel for private(n, k, sum0, sum1)
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
}//namespace tr
