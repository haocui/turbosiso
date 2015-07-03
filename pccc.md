# Parallel Concatenated Convolutional Codes of coding rate 1/3 #

The Parallel Concatenated Convolutional Codes (PCCCs) are Recursive and Systematic Convolutional (RSC) codes with the same generator polynomials ![http://turbosiso.googlecode.com/svn/wiki/form_26.png](http://turbosiso.googlecode.com/svn/wiki/form_26.png) and separated by an interleaver. The first RSC code has the trellis terminated (tail bits are added) and the second has unterminated trellis (no tail bits). The separating interleaver is random and has length ![http://turbosiso.googlecode.com/svn/wiki/form_27.png](http://turbosiso.googlecode.com/svn/wiki/form_27.png). The output is formed by using the systematic bits and the parity bits from both RSC codes. The symbols are BPSK modulated before being send into an Additive White Gaussian Noise (AWGN) channel.

![http://turbosiso.googlecode.com/svn/wiki/PCCC.jpg](http://turbosiso.googlecode.com/svn/wiki/PCCC.jpg)

**PCCC of coding rate 1/3**

> The turbo decoder for PCCCs uses two SISO RSC modules. The input of the first SISO RSC module is represented by the intrinsic information of data and parity bits. The input of the second RSC module is represented only by the intrinsic information of parity bits since the intrinsic information of data bits has been used in the first SISO RSC module. The extrinsic information of data bits is exchanged between the two SISO RSC modules.

![http://turbosiso.googlecode.com/svn/wiki/PCCC_tdec.jpg](http://turbosiso.googlecode.com/svn/wiki/PCCC_tdec.jpg)

**Turbo decoder for PCCC**

> The performance of PCCCs is shown below when using the max log MAP and the log MAP algorithms. It can be seen that the use of the max log MAP algorithm has worser performance at low Signal to Noise Ratio (SNR) than the log MAP algorithm. This is due to the fact that the approximation ![http://turbosiso.googlecode.com/svn/wiki/form_28.png](http://turbosiso.googlecode.com/svn/wiki/form_28.png), used in the max log MAP algorithm, does not hold at low SNR.

![http://turbosiso.googlecode.com/svn/wiki/pccc_maxlogMAP.jpg](http://turbosiso.googlecode.com/svn/wiki/pccc_maxlogMAP.jpg)

**Performance of the turbo decoder for PCCC using max log MAP algorithm**

![http://turbosiso.googlecode.com/svn/wiki/pccc_logMAP.jpg](http://turbosiso.googlecode.com/svn/wiki/pccc_logMAP.jpg)

**Performance of the turbo decoder for PCCC using log MAP algorithm**

> Simulation times of programs written in C++ and Matlab are shown in the table below. Both programs use the SISO class. Note that the time difference between simulation times (C++ and Matlab) is almost constant for both MAP algorithms since the most complex routine is written in C++.

Simulation times for PCCC in C++ and Matlab

|  | **C++** | **Matlab** |
|:-|:--------|:-----------|
| **max log MAP** | 1 hr, 0 min, 22 sec | 1 hr, 35 min, 40 sec |
| **log MAP** | 3 hr, 3 min, 57 sec | 3 hr, 32 min, 55 sec |

The EXIT diagram of the turbo decoder is shown below at a SNR ![http://turbosiso.googlecode.com/svn/wiki/form_29.png](http://turbosiso.googlecode.com/svn/wiki/form_29.png). Both SISO RSC modules use the log MAP algorihm. Note that in ten Brink's approach the extrinsic information of data bits is computed differently by substracting from the LLR both a priori and intrinsic information, not only the a priori information as in our approach.

![http://turbosiso.googlecode.com/svn/wiki/EXIT_pccc_logMAP.jpg](http://turbosiso.googlecode.com/svn/wiki/EXIT_pccc_logMAP.jpg)

**EXIT diagram of the turbo decoder**

> Programs used to obtain the figures: [pccc.cpp](http://turbosiso.googlecode.com/svn/trunk/pccc.cpp), [pccc.m](http://turbosiso.googlecode.com/svn/trunk/pccc.m) and [EXIT\_pccc.cpp](http://turbosiso.googlecode.com/svn/trunk/EXIT_pccc.cpp).

**References:**

S. Benedetto, D. Divsalar, G. Motorsi and F. Pollara, "A Soft-Input Soft-Output Maximum A posteriori (MAP) Module to Decode Parallel and Serial Concatenated Codes", TDA Progress Report, nov. 1996

S. ten Brink, ''Convergence behavior of iteratively decoded parallel concatenated codes,`` IEEE Transactions on Communications, vol. 49, pp. 1727-1737, Oct. 2001