/*
 * Convolutional Codes Class
 * Non recursive non Systematic Convolutional (NSC) code
 * Recursive Systematic Convolutional (RSC) code
 * Precoder (recursive non systematic convolutional code with coding rate 1)
 */

#ifndef CC_CLASS
#define CC_CLASS

#include <iostream>
#include <complex>
#include "itpp/itbase.h" //IT++ base module

namespace tr
{

class CC
{
public:
    CC()
    {
        tail = false;
    };
    void set_generators(const itpp::bmat &gen)//set CC generators
    {
        generators = gen;
        nb_outputs = gen.rows();
        reg_len = gen.cols()-1;
    };
    void set_generators(const itpp::ivec &gen, const int &constraint_length)
    {
    	nb_outputs = gen.length();
    	reg_len = constraint_length-1;
    	generators.set_size(nb_outputs, constraint_length);
        for (int n=0;n<nb_outputs;n++)
        	generators.set_row(n, itpp::dec2bin(constraint_length, gen(n)));
    };
    void set_generator(const itpp::bvec &gen)//set precoder generator
    {
        generators = gen;
        nb_outputs = 1;
        reg_len = gen.length()-1;
    };
    void set_generator(const int &gen, const int &constraint_length)
    {
    	nb_outputs = 1;
    	reg_len = constraint_length-1;
    	generators.set_size(1, constraint_length);
        generators.set_row(0, itpp::dec2bin(constraint_length, gen));
    };
    void set_initial_state(const itpp::bvec &reg)
    {
        ini_state = reg;
    };
    void set_tail(const bool &t)
    {
        tail = t;
    };
    itpp::bvec nsc(const itpp::bvec &in);
    itpp::bvec nsc(const itpp::bvec &in, const bool &tail)
    {
        set_tail(tail);
        return nsc(in);
    };
    itpp::bvec nsc(const itpp::bvec &in, const itpp::bmat &gen, const bool &tail)
    {
        set_generators(gen);
        set_tail(tail);
        return nsc(in);
    };
    itpp::bvec rsc(const itpp::bvec &in);
    itpp::bvec rsc(const itpp::bvec &in, const bool &tail)
    {
        set_tail(tail);
        return rsc(in);
    };
    itpp::bvec rsc(const itpp::bvec &in, const itpp::bmat &gen, const bool &tail)
    {
        set_generators(gen);
        set_tail(tail);
        return rsc(in);
    };
    itpp::bvec precoder(const itpp::bvec &in);
    itpp::bvec precoder(const itpp::bvec &in, const itpp::bvec &precoder_pol)
    {
        set_generators(precoder_pol);
        return precoder(in);
    };
    int outputs(void) const
    {
        return nb_outputs;
    };
    int memory_len(void) const
    {
        return reg_len;
    };
    itpp::bvec get_tail_input_bits(void) const
    {
        return tail_input_bits;
    };
private:
    itpp::bmat generators;//when appliable, first polynomial defines the feedback
    bool tail;
    int nb_outputs;
    itpp::bvec ini_state;
    int reg_len;
    itpp::bvec tail_input_bits;
    itpp::bin xor_sum(const itpp::bvec &in1, const itpp::bvec &in2);
};

//element wise multiplication of the two vectors followed by xor sum of the resulting vector
inline itpp::bin CC::xor_sum(const itpp::bvec &in1, const itpp::bvec &in2)
{
    itpp::bin out = 0;
    for (int n=0;n<in1.length();n++)
        out ^= (in1(n)*in2(n));
    return out;
}

//non recursive non systematic convolutional code
itpp::bvec CC::nsc(const itpp::bvec &in)
{
    int nb_bits = in.length();
    int len = nb_bits*nb_outputs;
    itpp::bvec out(len);
    itpp::bvec reg = ini_state;
    itpp::bvec in_reg(reg_len+1);
    register int k,n;
    for (n=0;n<nb_bits;n++)
    {
        in_reg(0) = in(n);
        in_reg.replace_mid(1, reg);
        for (k=0;k<nb_outputs;k++)
            out(k+nb_outputs*n) = xor_sum(in_reg, generators.get_row(k));
        reg.shift_right(in(n));
    }
    if (tail)//compute tail in order to finish in zero state
    {
        out.set_length(len+reg_len*nb_outputs, true);
        tail_input_bits.set_length(reg_len);
        for (n=0;n<reg_len;n++)
        {
            in_reg(0) = 0;
            in_reg.replace_mid(1, reg);
            for (k=0;k<nb_outputs;k++)
                out(len+k+nb_outputs*n) = xor_sum(in_reg, generators.get_row(k));
            reg.shift_right((itpp::bin)0);
            tail_input_bits(n) = 0;
        }
    }
    return out;
}

//recursive systematic convolutional code
itpp::bvec CC::rsc(const itpp::bvec &in)
{
    int nb_bits = in.length();
    int len = nb_bits*nb_outputs;
    itpp::bvec out(len);
    itpp::bvec reg = ini_state;
    itpp::bvec in_reg(reg_len+1);
    itpp::bin feedback;
    register int k,n;
    for (n=0;n<nb_bits;n++)
    {
        //systematic output
        out(nb_outputs*n) = in(n);
        //feedback
        in_reg(0) = in(n);
        in_reg.replace_mid(1, reg);
        feedback = xor_sum(in_reg, generators.get_row(0));
        //parity outputs
        in_reg(0) = feedback;
        for (k=1;k<nb_outputs;k++)
            out(k+nb_outputs*n) = xor_sum(in_reg, generators.get_row(k));
        reg.shift_right(feedback);
    }
    if (tail)//compute tail in order to finish in zero state
    {
        out.set_length(len+reg_len*nb_outputs, true);
        tail_input_bits.set_length(reg_len);
        for (n=0;n<reg_len;n++)
        {
            in_reg(0) = 0;
            feedback = xor_sum(in_reg, generators.get_row(0));
            out(len+nb_outputs*n) = feedback;//systematic output
            for (k=1;k<nb_outputs;k++)
                out(len+k+nb_outputs*n) = xor_sum(in_reg, generators.get_row(k));//parity output
            reg.shift_right((itpp::bin)0);
            tail_input_bits(n) = feedback;
        }
    }
    return out;
}

//precoder (recursive non systematic convolutional code with coding rate 1)
itpp::bvec CC::precoder(const itpp::bvec &in)
{
    int nb_bits = in.length();
    itpp::bvec reg = ini_state;
    itpp::bvec in_reg(reg_len+1);
    itpp::bvec out(nb_bits);
    register int n;
    //parity output
    for (n=0;n<nb_bits;n++)
    {
        in_reg(0) = in(n);
        in_reg.replace_mid(1, reg);
        out(n) = xor_sum(in_reg, generators.get_row(0));
        reg.shift_right(out(n));
    }
    return out;
}

/* Utilities functions
 */
//Kronecker operator for vectors (both seen as column or row vectors)
template <class Num_T>
itpp::Vec<Num_T> kron(const itpp::Vec<Num_T> &in, const itpp::Vec<Num_T> &pattern)
{
    int in_len = in.length();
    int pattern_len = pattern.length();
    itpp::Vec<Num_T> out(in_len*pattern_len);
    for (int n=0;n<in_len;n++)
        out.replace_mid(n*pattern_len, in(n)*pattern);
    return out;
}

}

#endif
