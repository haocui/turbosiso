/** \file
 * 
 * \brief Implementation of SISO modules for equalizers
 */

#include "SISO.h"

namespace tr
{
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
}//end namespace tr
