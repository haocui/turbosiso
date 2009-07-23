/** \file
 * 
 * \brief Space Time block Codes (STC) class
 * 
 * Implements Space Time block Codes using Hassibi's model
 * 
 * Reference: B. Hassibi and B. M. Hochwald, ''High-rate codes that are linear in space and time,`` 
 * IEEE Transactions on Information Theory, vol. 48, pp. 1804-1824, July 2002
 */

#ifndef STC_CLASS
#define STC_CLASS

#include <iostream>
#include <complex>
#include "itpp/itbase.h" //IT++ base module

namespace tr
{

/// Space Time block Codes (%STC)
/** ST block codes are implemented following Hassibi's approach.
 * 
 * Reference:
 * B. Hassibi and B. M. Hochwald, ''High-rate codes that are linear in space and time,`` IEEE Transactions on Information Theory, 
 * vol. 48, pp. 1804-1824, July 2002
 */
class STC
{
public:
	/// Setup ST block codes (Hassibi's method is used)
	void setup(const int &in_em_antennas, const int &in_channel_uses, const std::string &in_code_name, const int &in_const_size)
	{
		em_antennas = in_em_antennas;
		channel_uses = in_channel_uses;
		code_name = in_code_name;
		const_size = in_const_size;
		Hassibi_block_code();
	};
	/// Encodes input symbols according to specified ST code
	itpp::cmat encode(const itpp::cvec &symb)
	{
		return Hassibi_encode(symb);
	};
	/// Gets the number of symbols per ST code block
	const int get_nb_symbols_per_block(void) const
	{
		return symb_block;
	};
	/// Gets the first generator matrix of the ST code following Hassibi's approach
	const itpp::cmat get_1st_gen_matrix(void) const
	{
		return A;
	};
	/// Gets the second generator matrix of the ST code following Hassibi's approach
	const itpp::cmat get_2nd_gen_matrix(void) const
	{
		return B;
	};
	/// Gets the number of emission antennas
	const int get_nb_em_antennas(void) const
	{
		return em_antennas;
	};
	/// Gets the number of channel uses (ST block code duration [symbols])
	const int get_channel_uses(void) const
	{
		return channel_uses;
	};
private:
	void Hassibi_block_code(void);
	itpp::cmat Hassibi_encode(const itpp::cvec &symb);
	itpp::cmat diag_pow(const itpp::cmat &in_mat, const double &in_exp);
	itpp::mat mat_pow(const itpp::mat &in_mat, const int &in_exp);
	int symb_block;
	int const_size;
	itpp::cmat A;
	itpp::cmat B;
	int em_antennas;
	int channel_uses;
	std::string code_name;	
};

void STC::Hassibi_block_code(void)
/* this function implements the A and B matrices needed for Space-Time block codes generation following the Hassibi's approach:
 * S = sum_{q=1}^symb_block (A_q alpha_q + jB_q beta_q),
 * where s_q = alpha_q+jbeta_q is the symbol after modulation
 * each A_q and B_q matrix has dimension TxM
 * different A_q and B_q matrices are stacked one below the other, e.g. [A_1;A_2;...;A_Q]

 * input: code_name - code name whose generator matrices are to be generated
 *        const_size - constellation size (used in Damen code)
 * ouputs: symb_block - number of symbols per block
 *         A, B - generator matrices
 * inputs/ouputs: for some codes these are inputs for others they are
 * predefined, so they are outputs only
 *               em_antennas - number of emission antenna
 *               channel_uses - channel uses


 * author: Bogdan Cristea
 * revision date 19/02/07: Golden code was changed to have TxM generator matrices
 * revision date 18/01/08: translation to C++
 */
{ 
	if (code_name=="V-BLAST_MxN")//classical V-BLAST
	{
    	symb_block = channel_uses*em_antennas;//number of symbols/block
    	std::cout << "STC::LDcode: Warning! For " << code_name << " the following parameter is predefined:" << std::endl;
    	std::cout << "symb_block = channel_uses*em_antennas = " << symb_block << std::endl;
    	A.set_size(symb_block*channel_uses, em_antennas);
    	A.zeros();
    	itpp::mat temp(channel_uses, em_antennas);
    	temp.zeros();
    	register int tau,m;
    	for(tau=0;tau<channel_uses;tau++)
    	{
        	for(m=0;m<em_antennas;m++)
        	{
            	temp(tau,m) = 1;
            	A.set_submatrix((em_antennas*(tau-1)+m-1), 0, to_cmat(temp));
            	temp(tau,m) = 0;
        	}
    	}
    	B = A;
	} 
	else if (code_name=="imp_V-BLAST_MxN")//improved V-BLAST (code (31) in Hassibi's paper)
	{
    	if (channel_uses!=em_antennas)
    	{
        	std::cout << "STC::LDcode: Warning! For " << code_name << " channel_uses and em_antennas must be equal. Choosing channel_uses=em_antennas" << std::endl;
        	channel_uses = em_antennas;
    	}
    	symb_block = channel_uses*em_antennas;//number of symbols/block
    	std::cout << "STC::LDcode: Warning! For " << code_name << " the following parameter is predefined:" << std::endl;
    	std::cout << "symb_block = " << symb_block << std::endl;
    	std::complex<double> j(0,1);
    	itpp::cmat D = itpp::diag(exp(j*(2*itpp::pi/em_antennas)*itpp::linspace(0, em_antennas-1, em_antennas)));
    	itpp::mat P = itpp::diag(itpp::ones(em_antennas-1), -1);
    	P(0,em_antennas-1) = 1;
    	A.set_size(symb_block*channel_uses, em_antennas);
    	A.zeros();
    	register int k,l;
    	for (k=0;k<channel_uses;k++)    	
        	for (l=0;l<em_antennas;l++)        	
        		A.set_submatrix((em_antennas*k+l)*channel_uses, 0, diag_pow(D, k)*itpp::to_cmat(mat_pow(P, l))/sqrt(em_antennas));            	        	    	
    	B = A;
	}
	else if (code_name=="Alamouti_2xN")//Alamouti's orthogonal code
	{    
    	em_antennas = 2;//emission antenna
    	channel_uses = 2;//channel uses
    	symb_block = 2;//number of symbols/block
    	std::cout << "STC::LDcode: Warning! For " << code_name << " the following parameters are predefined:" << std::endl;
    	std::cout << "em_antennas = " << em_antennas << ", channel_uses = " << channel_uses << ", symb_block = " << symb_block << std::endl;
    	A = "1  0;\
         	 0  1;\
         	 0  1;\
        	-1  0";//A_1; A_2
    	B = "1  0;\
         	 0 -1;\
         	 0  1;\
         	 1  0";//B_1; B_2
	}
	else if (code_name=="Switched_Alamouti_4xN")
	{
    	em_antennas = 4;//emission antenna
    	channel_uses = 4;//channel uses
    	symb_block = 4;//number of symbols/block
    	std::cout << "STC::LDcode: Warning! For " << code_name << " the following parameters are predefined:" << std::endl;
    	std::cout << "em_antennas = " << em_antennas << ", channel_uses = " << channel_uses << ", symb_block = " << symb_block << std::endl;
    	A = "1  0  0  0;\
         	 0  1  0  0;\
	         0  0  0  0;\
         	 0  0  0  0;\
	         0  1  0  0;\
	         -1 0  0  0;\
	         0  0  0  0;\
	         0  0  0  0;\
	         0  0  0  0;\
	         0  0  0  0;\
	         0  0  1  0;\
	         0  0  0  1;\
	         0  0  0  0;\
	         0  0  0  0;\
	         0  0  0  1;\
	         0  0  -1 0";//A_1; A_2; A_3; A_4
        A *= sqrt(2);//normalization
	    B = "1  0  0  0;\
	         0  -1 0  0;\
	         0  0  0  0;\
	         0  0  0  0;\
	         0  1  0  0;\
	         1  0  0  0;\
	         0  0  0  0;\
	         0  0  0  0;\
	         0  0  0  0;\
	         0  0  0  0;\
	         0  0  1  0;\
	         0  0  0  -1;\
	         0  0  0  0;\
	         0  0  0  0;\
	         0  0  0  1;\
           	 0  0  1  0";//B_1; B_2; B_3; B_4
		 B *= sqrt(2);         
	}
	else if (code_name=="Double_Alamouti_4xN")
	{
    	em_antennas = 4;//emission antenna
    	channel_uses = 2;//channel uses
    	symb_block = 4;//number of symbols/block
    	std::cout << "STC::LDcode: Warning! For " << code_name << " the following parameters are predefined:" << std::endl;
    	std::cout << "em_antennas = " << em_antennas << ", channel_uses = " << channel_uses << ", symb_block = " << symb_block << std::endl;
    	A = "1  0  0  0;\
	         0  1  0  0;\
	         0  0  1  0;\
	         0  0  0  1;\
	         0  1  0  0;\
	         -1 0  0  0;\
	         0  0  0  1;\
	         0  0  -1 0";//A_1; A_2; A_3; A_4
    	B = "1  0  0  0;\
	         0  -1 0  0;\
	         0  0  1  0;\
	         0  0  0 -1;\
	         0  1  0  0;\
	         1  0  0  0;\
	         0  0  0  1;\
         	 0  0  1  0";//B_1; B_2; B_3; B_4
	}
	else if (code_name=="Jafarkhani_4xN")//Jafarkhani's quasi-orthogonal code
	{
    	em_antennas = 4;//emission antenna
    	channel_uses = 4;//channel uses
    	symb_block = 4;//number of symbols/block
    	std::cout << "STC::LDcode: Warning! For " << code_name << " the following parameters are predefined:" << std::endl;
	    std::cout << "em_antennas = " << em_antennas << ", channel_uses = " << channel_uses << ", symb_block = " << symb_block << std::endl;
	    A = "1  0  0  0;\
	         0  1  0  0;\
	         0  0  1  0;\
	         0  0  0  1;\
	         0  1  0  0;\
	         -1 0  0  0;\
	         0  0  0  1;\
	         0  0  -1 0;\
	         0  0  1  0;\
	         0  0  0  1;\
	         -1 0  0  0;\
	         0  -1 0  0;\
	         0  0  0  1;\
	         0  0  -1 0;\
	         0  -1 0  0;\
	         1  0  0  0";//A_1; A_2; A_3; A_4
	    B = "1  0  0  0;\
	         0  -1  0  0;\
	         0  0  -1  0;\
	         0  0  0  1;\
	         0  1  0  0;\
	         1  0  0  0;\
	         0  0  0  -1;\
	         0  0  -1 0;\
	         0  0  1  0;\
	         0  0  0  -1;\
	         1  0  0  0;\
	         0  -1 0  0;\
	         0  0  0  1;\
	         0  0  1  0;\
	         0  1  0  0;\
	         1  0  0  0";//B_1; B_2; B_3; B_4
	}
	else if (code_name=="Golden_2x2")//Golden code as proposed by Belfiore
	{
    	em_antennas = 2;//emission antenna
    	channel_uses = 2;//channel uses
    	symb_block = 4;//number of symbols/block
	    std::cout << "STC::LDcode: Warning! For " << code_name << " the following parameters are predefined:" << std::endl;
	    std::cout << "em_antennas = " << em_antennas << ", channel_uses = " << channel_uses << ", symb_block = " << symb_block << std::endl;
    	std::complex<double> theta((1+sqrt(5))/2,0);
    	std::complex<double> theta_b((1-sqrt(5))/2,0);
    	std::complex<double> j(0,1);
    	std::complex<double> one(1,0);
    	std::complex<double> alpha = one+j*(one-theta);
    	std::complex<double> alpha_b = one+j*(one-theta_b);
    	std::complex<double> gamma = j;
    	A.set_size(8,2);
    	A(0,0) = alpha/sqrt(5); 	  A(0,1) = 0;
    	A(1,0) = 0;			    	  A(1,1) =  alpha_b/sqrt(5);//A_1
    	A(2,0) = alpha*theta/sqrt(5); A(2,1) = 0;
    	A(3,0) = 0;					  A(3,1) = alpha_b*theta_b/sqrt(5);//A_2
    	A(4,0) = 0;					  A(4,1) = gamma*alpha_b/sqrt(5);
    	A(5,0) = alpha/sqrt(5);		  A(5,1) = 0;//A_3
    	A(6,0) = 0;					  A(6,1) = gamma*alpha_b*theta_b/sqrt(5);
    	A(7,0) = alpha*theta/sqrt(5); A(7,1) = 0;//A_4
     	B = A;
	}
	else if (code_name=="Damen_2x2")//ST code based on number theory as proposed by Damen
	{
    	em_antennas = 2;//emission antenna
    	channel_uses = 2;//channel uses
    	symb_block = 4;//number of symbols/block
    	std::cout << "STC::LDcode: Warning! For " << code_name << " the following parameters are predefined:" << std::endl;
    	std::cout << "em_antennas = " << em_antennas << ", channel_uses = " << channel_uses << ", symb_block = " << symb_block << std::endl;
    	double lambda;
    	if (const_size==4)
        	lambda = 0.5;
    	else if (const_size==16)
	        lambda = 0.521;
    	else if (const_size>=256)
	        lambda = itpp::pi/4;
    	else
    	{
        	lambda = itpp::pi/4;
        	std::cout << "STC::LDcode: Warning! For " << code_name << " and const. size " << const_size << ", lambda has the value " << lambda << std::endl;
    	}
    	std::complex<double> j(0,1);
    	std::complex<double> phi = exp(j*lambda);
    	std::complex<double> theta = exp(j*(lambda/2));
    	A.set_size(8, 2);
    	A(0,0) = 1/sqrt(2);     	A(0,1) = 0;
    	A(1,0) = 0;			    	A(1,1) = 1/sqrt(2);//A_1
    	A(2,0) = phi/sqrt(2);   	A(2,1) = 0;
    	A(3,0) = 0;			    	A(3,1) = -phi/sqrt(2);//A_2
    	A(4,0) = 0;			    	A(4,1) = theta/sqrt(2);
    	A(5,0) = theta/sqrt(2); 	A(5,1) = 0;//A_3
    	A(6,0) = 0;					A(6,1) = -theta*phi/sqrt(2);
    	A(7,0) = theta*phi/sqrt(2); A(7,1) = 0;//A_4
    	B = A;
	}
	else if (code_name=="34ortho_3xN")//rate 3/4 orthogonal code (mutual information 5.13 bits/channel use at rho=20 dB)
	{
    	em_antennas = 3;//emission antenna
    	channel_uses = 4;//channel uses
    	symb_block = 3;//number of symbols/block
    	std::cout << "STC::LDcode: Warning! For " << code_name << " the following parameters are predefined:" << std::endl;
    	std::cout << "em_antennas = " << em_antennas << ", channel_uses = " << channel_uses << ", symb_block = " << symb_block << std::endl;
    	A = "1 0 0;\
         	 0 1 0;\
	         0 0 1;\
	         0 0 0;\
	         0 1 0;\
	         -1 0 0;\
	         0 0 0;\
	         0 0 1;\
	         0 0 1;\
	         0 0 0;\
	         -1 0 0;\
	         0 -1 0";//A_1; A_2; A_3
        A /= sqrt(double(4)/double(3));
     	B = "1 0 0;\
	         0 -1 0;\
	         0 0 -1;\
	         0 0 0;\
	         0 1 0;\
	         1 0 0;\
	         0 0 0;\
	         0 0 -1;\
	         0 0 1;\
	         0 0 0;\
	         1 0 0;\
          	 0 1 0";//B_1; B_2; B_3
         B /= sqrt(double(4)/double(3));
	}
	else if (code_name=="36LD_3xN")//(36) LD code with mutual info. 6.25bits/channel use at rho=20dB
	{
    	em_antennas = 3;//emission antenna
    	channel_uses = 4;//channel uses
    	symb_block = 4;//number of symbols/block
    	std::cout << "STC::LDcode: Warning! For " << code_name << " the following parameters are predefined:" << std::endl;
    	std::cout << "em_antennas = " << em_antennas << ", channel_uses = " << channel_uses << ", symb_block = " << symb_block << std::endl;
    	A.set_size(16, 3);
    	A(0,0) = 1; 		  A(0,1) = 0;	 	    A(0,2) = 0;
    	A(1,0) = 1; 		  A(1,1) = 1; 		    A(1,2) = 0;
    	A(2,0) = 0; 		  A(2,1) = 0; 		    A(2,2) = 1;
    	A(3,0) = 0; 		  A(3,1) = 0; 		    A(3,2) = 0;//A_1
    	A(4,0) = 0; 		  A(4,1) = 1/sqrt(2);   A(4,2) = 0;
    	A(5,0) = -1/sqrt(2);  A(5,1) = 0;	 	    A(5,2) = -1/sqrt(2);
    	A(6,0) = 0; 		  A(6,1) = 1/sqrt(2);   A(6,2) = 0;
    	A(7,0) = 1/sqrt(2);   A(7,1) = 0;	 	    A(7,2) = -1/sqrt(2);//A_2
		A(8,0) = 1; 		  A(8,1) = 0;	 	    A(8,2) = 0;
		A(9,0) = 0; 		  A(9,1) = 0;	 	    A(9,2) = 0;
		A(10,0) = 0; 		  A(10,1) = 0;	 	    A(10,2) = -1;
		A(11,0) = 0; 		  A(11,1) = -1;	 	    A(11,2) = 0;//A_3
		A(12,0) = 0; 		  A(12,1) = -1/sqrt(2); A(12,2) = 0;
		A(13,0) = 1/sqrt(2);  A(13,1) = 0;	 	    A(13,2) = -1/sqrt(2);
		A(14,0) = 0; 		  A(14,1) = 1/sqrt(2);  A(14,2) = 0;
		A(15,0) = -1/sqrt(2); A(15,1) = 0;	 	    A(15,2) = -1/sqrt(2);//A_4
		B.set_size(16, 3);
		B(0,0) = 0; 		  		    B(0,1) = -1/sqrt(2);   		   B(0,2) = 0;
		B(1,0) = -1/sqrt(2);  		    B(1,1) = 0;  			 	   B(1,2) = 1/sqrt(2);
		B(2,0) = 0; 		  		    B(2,1) = 1/sqrt(2);    		   B(2,2) = 0;
		B(3,0) = 1/sqrt(2);   		    B(3,1) = 0;  			 	   B(3,2) = 1/sqrt(2);//B_1
		B(4,0) = 1/sqrt(2);             B(4,1) = double(-1)/double(2); B(4,2) = 0;
		B(5,0) = double(-1)/double(2);  B(5,1) = -1/sqrt(2);  		   B(5,2) = double(-1)/double(2);
		B(6,0) = 0; 		  		    B(6,1) = double(-1)/double(2); B(6,2) = 1/sqrt(2);
		B(7,0) = double(1)/double(2);   B(7,1) = 0;   		  		   B(7,2) = double(-1)/double(2);//B_2
		B(8,0) = 1/sqrt(2); 		    B(8,1) = double(1)/double(2);  B(8,2) = 0;
		B(9,0) = double(1)/double(2);   B(9,1) = -1/sqrt(2);   		   B(9,2) = double(1)/double(2);
		B(10,0) = 0; 		  		    B(10,1) = double(1)/double(2); B(10,2) = 1/sqrt(2);
		B(11,0) = double(-1)/double(2); B(11,1) = 0;   		  		   B(11,2) = double(1)/double(2);//B_3
		B(12,0) = 1; 		  		    B(12,1) = 0;   		   		   B(12,2) = 0;
		B(13,0) = 0; 		  		    B(13,1) = 0;   		  		   B(13,2) = 0;
		B(14,0) = 0; 		  		    B(14,1) = 0;   		   		   B(14,2) = -1;
		B(15,0) = 0; 		  		    B(15,1) = 1;   		   		   B(15,2) = 0;//B_4
	}
	else if (code_name=="37LD_3xN")//(37) LD code 3-antenna LD code obtained from the symetrical concatenation of 3 2-antenna orthogonal design
	{
    	em_antennas = 3;//emission antenna
    	channel_uses = 6;//channel uses
    	symb_block = 6;//number of symbols/block
    	std::cout << "STC::LDcode: Warning! For " << code_name << " the following parameters are predefined:" << std::endl;
    	std::cout << "em_antennas = " << em_antennas << ", channel_uses = " << channel_uses << ", symb_block = " << symb_block << std::endl;
    	A = "1  0  0;\
         	 0  1  0;\
         	 0  0  0;\
         	 0  0  0;\
         	 0  0  0;\
         	 0  0  0;\
         	 0  1  0;\
         	 -1 0  0;\
          	 0  0  0;\
         	 0  0  0;\
         	 0  0  0;\
         	 0  0  0;\
         	 0  0  0;\
         	 0  0  0;\
         	 0  1  0;\
         	 0  0  1;\
         	 0  0  0;\
         	 0  0  0;\
         	 0  0  0;\
         	 0  0  0;\
         	 0  0  1;\
         	 0 -1  0;\
         	 0  0  0;\
         	 0  0  0;\
         	 0  0  0;\
         	 0  0  0;\
         	 0  0  0;\
         	 0  0  0;\
         	 1  0  0;\
          	 0  0  1;\
         	 0  0  0;\
         	 0  0  0;\
         	 0  0  0;\
         	 0  0  0;\
         	 0  0  1;\
         	 -1 0  0";//A_1; A_2; A_3; A_4; A_5; A_6
        A *= sqrt(double(3)/double(2));
    	B = "1  0  0;\
         	 0  -1  0;\
         	 0  0  0;\
         	 0  0  0;\
         	 0  0  0;\
         	 0  0  0;\
         	 0  1  0;\
         	 1 0  0;\
         	 0  0  0;\
         	 0  0  0;\
         	 0  0  0;\
         	 0  0  0;\
         	 0  0  0;\
         	 0  0  0;\
         	 0  1  0;\
         	 0  0  -1;\
         	 0  0  0;\
         	 0  0  0;\
         	 0  0  0;\
         	 0  0  0;\
         	 0  0  1;\
         	 0  1  0;\
         	 0  0  0;\
         	 0  0  0;\
         	 0  0  0;\
         	 0  0  0;\
         	 0  0  0;\
         	 0  0  0;\
         	 1  0  0;\
         	 0  0  -1;\
         	 0  0  0;\
         	 0  0  0;\
         	 0  0  0;\
         	 0  0  0;\
         	 0  0  1;\
         	 1  0  0";//B_1; B_2; B_3; B_4; B_5; B_6
        B *= sqrt(double(3)/double(2));
	}
	else if (code_name=="39LD_3xN")
	{
    	em_antennas = 3;//emission antenna
    	channel_uses = 6;//channel uses
    	symb_block = 6;//number of symbols/block
    	std::cout << "STC::LDcode: Warning! For " << code_name << " the following parameters are predefined:" << std::endl;
    	std::cout << "em_antennas = " << em_antennas << ", channel_uses = " << channel_uses << ", symb_block = " << symb_block << std::endl;
    	A.set_size(36, 3);
    	A(0,0) = 1/sqrt(2);   	  		A(0,1) = 0; 		   		    A(0,2) = 0;
    	A(1,0) = 0; 		  	  		A(1,1) = 1/sqrt(2);  		    A(1,2) = 0;
    	A(2,0) = 0; 		  	  		A(2,1) = 1/sqrt(2);  		    A(2,2) = 0;
    	A(3,0) = 0; 		  	 		A(3,1) = 0; 		   		    A(3,2) = 1/sqrt(2);
    	A(4,0) = 1/sqrt(2);   	  		A(4,1) = 0; 		   		    A(4,2) = 0;
    	A(5,0) = 0; 		  	  		A(5,1) = 0; 		   		    A(5,2) = 1/sqrt(2);//A_1
    	A(6,0) = 0; 		  	  		A(6,1) = 1/sqrt(2);  		    A(6,2) = 0;
    	A(7,0) = -1/sqrt(2);  	  		A(7,1) = 0; 		   		    A(7,2) = 0;
    	A(8,0) = 0; 		  	  		A(8,1) = 0; 		   		    A(8,2) = 1/sqrt(2);
    	A(9,0) = 0; 		  	  		A(9,1) = -1/sqrt(2); 		    A(9,2) = 0;
    	A(10,0) = 0; 		  	  		A(10,1) = 0; 		   		    A(10,2) = 1/sqrt(2);
    	A(11,0) = -1/sqrt(2); 	  		A(11,1) = 0; 		   		    A(11,2) = 0;//A_2
    	A(12,0) = 1/sqrt(2);  	  		A(12,1) = 0; 		   		    A(12,2) = 0;
    	A(13,0) = 0; 		  	  		A(13,1) = 1/sqrt(2);      	    A(13,2) = 0;
    	A(14,0) = 0; 		  	  		A(14,1) = -1/(2*sqrt(2)); 	    A(14,2) = -sqrt(3)/(2*sqrt(2));
    	A(15,0) = 0; 		  	  		A(15,1) = sqrt(3)/(2*sqrt(2));  A(15,2) = -1/(2*sqrt(2));
    	A(16,0) = -1/(2*sqrt(2)); 		A(16,1) = 0; 	 			    A(16,2) = sqrt(3)/(2*sqrt(2));
		A(17,0) = -sqrt(3)/(2*sqrt(2)); A(17,1) = 0; 	 			    A(17,2) = -1/(2*sqrt(2));//A_3
		A(18,0) = 0; 		  	  		A(18,1) = 1/sqrt(2);  		    A(18,2) = 0;
		A(19,0) = -1/sqrt(2);   	  	A(19,1) = 0; 		   		    A(19,2) = 0;
		A(20,0) = 0; 		  	  		A(20,1) = sqrt(3)/(2*sqrt(2));  A(20,2) = -1/(2*sqrt(2));
		A(21,0) = 0; 		  	  		A(21,1) = 1/(2*sqrt(2)); 	    A(21,2) = sqrt(3)/(2*sqrt(2));
		A(22,0) = -sqrt(3)/(2*sqrt(2)); A(22,1) = 0; 	 			    A(22,2) = -1/(2*sqrt(2));
		A(23,0) = 1/(2*sqrt(2)); 		A(23,1) = 0; 	 			    A(23,2) = -sqrt(3)/(2*sqrt(2));//A_4
    	A(24,0) = 1/sqrt(2);   	  		A(24,1) = 0; 		   		    A(24,2) = 0;
    	A(25,0) = 0; 		  	  		A(25,1) = 1/sqrt(2);  		    A(25,2) = 0;
    	A(26,0) = 0; 		  	  		A(26,1) = -1/(2*sqrt(2)); 	    A(26,2) = sqrt(3)/(2*sqrt(2));
    	A(27,0) = 0; 		  	  		A(27,1) = -sqrt(3)/(2*sqrt(2)); A(27,2) = -1/(2*sqrt(2));
    	A(28,0) = -1/(2*sqrt(2)); 		A(28,1) = 0; 	 			    A(28,2) = -sqrt(3)/(2*sqrt(2));
    	A(29,0) = sqrt(3)/(2*sqrt(2));  A(29,1) = 0; 	 			    A(29,2) = -1/(2*sqrt(2));//A_5
		A(30,0) = 0; 		  	  		A(30,1) = 1/sqrt(2);  		    A(30,2) = 0;
		A(31,0) = -1/sqrt(2);  	  		A(31,1) = 0; 		   		    A(31,2) = 0;
		A(32,0) = 0; 		  	  		A(32,1) = -sqrt(3)/(2*sqrt(2)); A(32,2) = -1/(2*sqrt(2));
		A(33,0) = 0; 		  	  		A(33,1) = 1/(2*sqrt(2)); 	    A(33,2) = -sqrt(3)/(2*sqrt(2));
		A(34,0) = sqrt(3)/(2*sqrt(2));  A(34,1) = 0; 	 			    A(34,2) = -1/(2*sqrt(2));
		A(35,0) = 1/(2*sqrt(2)); 		A(35,1) = 0; 	 			    A(35,2) = sqrt(3)/(2*sqrt(2));//A_6
		B.set_size(36, 3);
		B(0,0) = 1/sqrt(2);   	  		B(0,1) = 0; 		   		    B(0,2) = 0;
		B(1,0) = 0; 		  	  		B(1,1) = -1/sqrt(2);  		    B(1,2) = 0;
		B(2,0) = 0; 		  	  		B(2,1) = 1/sqrt(2);  		    B(2,2) = 0;
		B(3,0) = 0; 		  	 		B(3,1) = 0; 		   		    B(3,2) = -1/sqrt(2);
		B(4,0) = 1/sqrt(2);   	  		B(4,1) = 0; 		   		    B(4,2) = 0;
		B(5,0) = 0; 		  	 		B(5,1) = 0; 		   		    B(5,2) = -1/sqrt(2);//B_1
		B(6,0) = 0; 		  	  		B(6,1) = 1/sqrt(2);  		    B(6,2) = 0;
		B(7,0) = 1/sqrt(2);   	  		B(7,1) = 0; 		   		    B(7,2) = 0;
		B(8,0) = 0; 		  	 		B(8,1) = 0; 		   		    B(8,2) = 1/sqrt(2);
		B(9,0) = 0; 		  	  		B(9,1) = 1/sqrt(2);  		    B(9,2) = 0;
		B(10,0) = 0; 		  	 		B(10,1) = 0; 		   		    B(10,2) = 1/sqrt(2);
		B(11,0) = 1/sqrt(2);   	  		B(11,1) = 0; 		   		    B(11,2) = 0;//B_2
		B(12,0) = 1/sqrt(2);   	  		B(12,1) = 0; 		   		    B(12,2) = 0;
		B(13,0) = 0; 		  	  		B(13,1) = -1/sqrt(2);  		    B(13,2) = 0;
		B(14,0) = 0; 		  	  		B(14,1) = -1/(2*sqrt(2)); 	    B(14,2) = -sqrt(3)/(2*sqrt(2));
		B(15,0) = 0; 		  	  		B(15,1) = -sqrt(3)/(2*sqrt(2)); B(15,2) = 1/(2*sqrt(2));
		B(16,0) = -1/(2*sqrt(2)); 		B(16,1) = 0; 	 			    B(16,2) = sqrt(3)/(2*sqrt(2));
		B(17,0) = sqrt(3)/(2*sqrt(2));  B(17,1) = 0; 	 			    B(17,2) = 1/(2*sqrt(2));//B_3
		B(18,0) = 0; 		  	  		B(18,1) = 1/sqrt(2);  		    B(18,2) = 0;
		B(19,0) = 1/sqrt(2);   	  		B(19,1) = 0; 		   		    B(19,2) = 0;
		B(20,0) = 0; 		  	  		B(20,1) = sqrt(3)/(2*sqrt(2));  B(20,2) = -1/(2*sqrt(2));
		B(21,0) = 0; 		  	  		B(21,1) = -1/(2*sqrt(2)); 	    B(21,2) = -sqrt(3)/(2*sqrt(2));
		B(22,0) = -sqrt(3)/(2*sqrt(2)); B(22,1) = 0; 	 			    B(22,2) = -1/(2*sqrt(2));
		B(23,0) = -1/(2*sqrt(2)); 		B(23,1) = 0; 	 			    B(23,2) = sqrt(3)/(2*sqrt(2));//B_4
		B(24,0) = 1/sqrt(2);   	  		B(24,1) = 0; 		   		    B(24,2) = 0;
		B(25,0) = 0; 		  	  		B(25,1) = -1/sqrt(2);  		    B(25,2) = 0;
		B(26,0) = 0; 		  	  		B(26,1) = -1/(2*sqrt(2)); 	    B(26,2) = sqrt(3)/(2*sqrt(2));
		B(27,0) = 0; 		  	  		B(27,1) = sqrt(3)/(2*sqrt(2));  B(27,2) = 1/(2*sqrt(2));
		B(28,0) = -1/(2*sqrt(2)); 		B(28,1) = 0; 	 			    B(28,2) = -sqrt(3)/(2*sqrt(2));
		B(29,0) = -sqrt(3)/(2*sqrt(2)); B(29,1) = 0; 	 			    B(29,2) = 1/(2*sqrt(2));//B_5
		B(30,0) = 0; 		  	  		B(30,1) = 1/sqrt(2);  		    B(30,2) = 0;
		B(31,0) = 1/sqrt(2);   	  		B(31,1) = 0; 		   		    B(31,2) = 0;
		B(32,0) = 0; 		  	  		B(32,1) = -sqrt(3)/(2*sqrt(2)); B(32,2) = -1/(2*sqrt(2));
		B(33,0) = 0; 		  	  		B(33,1) = -1/(2*sqrt(2)); 	    B(33,2) = sqrt(3)/(2*sqrt(2));
		B(34,0) = sqrt(3)/(2*sqrt(2));  B(34,1) = 0; 	 			    B(34,2) = -1/(2*sqrt(2));
		B(35,0) = -1/(2*sqrt(2)); 		B(35,1) = 0; 	 			    B(35,2) = -sqrt(3)/(2*sqrt(2));//B_6
	}
	else
		std::cout << "STC::LDcode: unknown code name. Available codes are: V-BLAST_MxN, imp_V-BLAST_MxN, Alamouti_2xN, \
		Switched_Alamouti_4xN, Double_Alamouti_4xN, Jafarkhani_4xN, Golden_2x2, Damen_2x2, 34ortho_3xN, 36LD_3xN, 37LD_3xN, 39LD_3xN" << std::endl;
}

itpp::cmat STC::Hassibi_encode(const itpp::cvec &symb)
//LD code generation (symb_block symbols go to an channel_uses x em_antennas matrix) following Hassibi's approach
{	
	int nb_subblocks = symb.length()/symb_block;
	int tx_duration = channel_uses*nb_subblocks;
    itpp::cmat S(tx_duration,em_antennas);
    itpp::cmat temp(channel_uses,em_antennas);
    std::complex<double> j(0,1);
    register int ns,k;    
    for (ns=0;ns<nb_subblocks;ns++)//encode block by block (symb_block symbols)
    {
    	temp.zeros();
        for (k=0;k<symb_block;k++)//sum over all symb_block matrices
            temp += (A(k*channel_uses,(k+1)*channel_uses-1,0,em_antennas-1)*static_cast< std::complex<double> >(symb(k+ns*symb_block).real())\
            	    +j*B(k*channel_uses,(k+1)*channel_uses-1,0,em_antennas-1)*static_cast< std::complex<double> >(symb(k+ns*symb_block).imag()));
        S.set_submatrix(ns*channel_uses, 0, temp);
    }
    return S;
}

inline itpp::cmat STC::diag_pow(const itpp::cmat &in_mat, const double &in_exp)
//first input should be a diagonal square matrix with complex elements
{
	register int n;
	int dim = in_mat.rows();
	itpp::cmat out_mat(dim,dim);
	for (n=0;n<dim;n++)
		out_mat(n,n) = std::pow(in_mat(n,n), in_exp);
	return out_mat;
}

inline itpp::mat STC::mat_pow(const itpp::mat &in_mat, const int &in_exp)
//square matrix power of integer exponent
{
	if (in_exp==0)
		return itpp::eye(in_mat.rows());
	itpp::mat out = in_mat;
	int abs_in_exp = std::abs(in_exp);
	register int n;
	for (n=1;n<abs_in_exp;n++)
		out *= in_mat;
	return (in_exp>0)?out:itpp::inv(out);
}

}
#endif
