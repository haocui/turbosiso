/** \file
 *
 * \brief Space Time Codes (STC) class header
 */

#ifndef STC_H_
#define STC_H_

#include <iostream>
#include <complex>
#include "itpp/itbase.h" //IT++ base module

namespace tr
{

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

}
#endif /* STC_H_ */
