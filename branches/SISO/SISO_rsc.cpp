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

void SISO::rsc_sova(itpp::vec &extrinsic_data, const itpp::vec &intrinsic_coded,
		const itpp::vec &apriori_data, const int &win_len)
/* Soft Output Viterbi Algorithm (SOVA) optimized for 1/N encoders
 * Output: extrinsic_data - extrinsic information of data bits
 * Inputs: intrinsic_coded - intrinsic information of data bits
 * 		   apriori_data - a priori information of data bits
 *		   win_len - window length used to represent the code trellis
 *
 * The original version has been written by Adrian Bohdanowicz (2003).
 * It is assumed that the BPSK mapping is: 0 -> +1, 1 -> -1.
 * Changes have been made to adapt the code for RSC codes of rate 1/2
 * and for soft input informations.
 */
{
	//number of code outputs
	int nb_outputs = gen.rows();

	//setup internal variables based on RSC trellis
	register int i,j,s;
	gen_rsctrellis();//trellis generation for 1/2 RSC codes
	int nb_states = rsctrellis.numStates;
	itpp::Array<itpp::mat> bin_out(2);//contains code output for each initial state and code input
	itpp::imat next_state(nb_states,2);//next state in the trellis
	for(i=0;i<2;i++)
	{
		bin_out(i).set_size(nb_states, nb_outputs);
		for(j=0;j<nb_states;j++)
		{
			bin_out(i)(j,0) = double(i);//systematic bit
			bin_out(i)(j,1) = rsctrellis.PARout[j+i*nb_states];//parity bit
			next_state(j,i) = rsctrellis.nextStates[j+i*nb_states];
		}
	}
	itpp::vec bin_inp("0 1");//binary code inputs

	int len = apriori_data.length();//number of decoding steps (total)

	//allocate memory for the trellis window
	itpp::mat metr(nb_states,win_len+1);//path metric buffer
	metr.zeros();
	metr += -INFINITY;
	metr(0,0) = 0;//initial state => (0,0)
	itpp::mat surv(nb_states,win_len+1);//survivor state buffer
	surv.zeros();
	itpp::mat inpt(nb_states,win_len+1);//survivor input buffer (dec. output)
	inpt.zeros();
	itpp::mat diff(nb_states,win_len+1);//path metric difference
	diff.zeros();
	itpp::mat comp(nb_states,win_len+1);//competitor state buffer
	comp.zeros();
	itpp::mat inpc(nb_states,win_len+1);//competitor input buffer
	inpc.zeros();
	//soft output (sign with reliability)
	itpp::vec sft(len);
	sft.zeros();
	sft += INFINITY;

	//decode all the bits
	int Cur,Nxt,nxt,sur,b,tmp,idx;
	itpp::vec buf(nb_outputs);
	double llb,mtr,dif,cmp,inc,srv,inp;
	itpp::vec bin(nb_outputs);
	itpp::ivec cyclic_buffer(win_len);
	extrinsic_data.set_size(len);
	for(i = 0; i < len; i++)
	{
	    //indices + precalculations
		Cur = i%(win_len+1);//curr trellis (cycl. buf) position
		Nxt = (i+1)%(win_len+1);//next trellis (cycl. buf) position
		buf = intrinsic_coded(i*nb_outputs,(i+1)*nb_outputs-1);//intrinsic_info portion to be processed
	    llb = apriori_data(i);//SOVA: apriori_info portion to be processed
	    metr.set_col(Nxt, -INFINITY*itpp::ones(nb_states));

	    //forward recursion
	    for(s = 0; s<nb_states; s++)
	    {
	        for(j = 0; j<2; j++)
	        {
	            nxt = next_state(s,j);//state after transition
	            bin = bin_out(j).get_row(s);//transition output (encoder)
	            mtr = bin*buf+metr(s,Cur);//transition metric
	            mtr += bin_inp(j)*llb;//SOVA

	            if(metr(nxt,Nxt) < mtr)
	            {
	                diff(nxt,Nxt) = mtr-metr(nxt,Nxt);//SOVA
	                comp(nxt,Nxt) = surv(nxt,Nxt);//SOVA
	                inpc(nxt,Nxt) = inpt(nxt,Nxt);//SOVA

	                metr(nxt,Nxt) = mtr;//store the metric
	                surv(nxt,Nxt) = s;//store the survival state
	                inpt(nxt,Nxt) = j;//store the survival input
	            }
	            else
	            {
	                dif = metr(nxt,Nxt)-mtr;
	                if(dif <= diff(nxt,Nxt))
	                {
	                    diff(nxt,Nxt) = dif;//SOVA
	                    comp(nxt,Nxt) = s;//SOVA
	                    inpc(nxt,Nxt) = j;//SOVA
	                }
	            }
	        }
	    }

	    //trace backwards
	    if(i < (win_len-1))
	    {
	        continue;
	    }//proceed if the buffer has been filled
	    mtr = itpp::max(metr.get_col(Nxt), sur);//find the initial state (max metric)
	    b = i;//temporary bit index
        for(j=0; j<win_len; j++)//indices in a 'cyclic buffer' operation
        {
        	cyclic_buffer(j) = (Nxt-j)%(win_len+1);
            cyclic_buffer(j) = (cyclic_buffer(j)<0)?(cyclic_buffer(j)+win_len+1):cyclic_buffer(j);
        }

	    for(j=0; j<win_len; j++)//for all the bits in the buffer
	    {
	        inp = inpt(sur,cyclic_buffer(j));//current bit-decoder output (encoder input)
	        extrinsic_data(b) = inp;//store the hard output

	        tmp = cyclic_buffer(j);
	        cmp = comp(sur, tmp);//SOVA: competitor state (previous)
	        inc = inpc(sur, tmp);//SOVA: competitor bit output
	        dif = diff(sur, tmp);//SOVA: corresp. path metric difference
	        srv = surv(sur, tmp);//SOVA: temporary survivor path state

	        for(s=j+1; s<=win_len; s++)//check all buffer bits srv and cmp paths
	        {
	        	if(inp != inc)
	        	{
	        		idx = b-((s-1)-j);//calculate index: [enc.k*(b-(k-1)+j-1)+1:enc.k*(b-(k-1)+j)]
	        		sft(idx) = std::min(sft(idx), dif);//update LLRs for bits that are different
	        	}
	        	if(srv == cmp)
	        	{
	        		break;
	        	}//stop if surv and comp merge (no need to continue)
	            if(s == win_len)
	            {
	                break;
	            }//stop if the end (otherwise: error)
	            tmp = cyclic_buffer(s);
	            inp = inpt(srv, tmp);//previous surv bit
	            inc = inpt(cmp, tmp);//previous comp bit
	            srv = surv(srv, tmp);//previous surv state
	            cmp = surv(cmp, tmp);//previous comp state
	        }
	        sur = surv(sur, cyclic_buffer(j));//state for the previous surv bit
	        b--;//update bit index
	    }
	}

	// provide soft output with +/- sign:
	extrinsic_data = itpp::elem_mult((2.0*extrinsic_data-1.0), sft)-apriori_data;

	//free trellis memory
	delete[] rsctrellis.prevStates;
	delete[] rsctrellis.nextStates;
	delete[] rsctrellis.PARout;
	delete[] rsctrellis.fm;
}

}//namespace tr
