# Serial Concatenated Convolutional Codes of coding rate 1/4 #

The Serial Concatenated Convolutional Codes (SCCCs) are realized with a Non-recursive non-Systematic Convolutional (NSC) code followed, after interleaving, by a Recursive and Systematic Convolutional (RSC) code. This configuration ensures the best performance for SCCCs. Both convolutional codes have the same generator polynomials ![http://turbosiso.googlecode.com/svn/wiki/form_26.png](http://turbosiso.googlecode.com/svn/wiki/form_26.png). The interleaver length is ![http://turbosiso.googlecode.com/svn/wiki/form_27.png](http://turbosiso.googlecode.com/svn/wiki/form_27.png). The symbols are BPSK modulated before being send into an Additive White Gaussian Noise (AWGN) channel.

![http://turbosiso.googlecode.com/svn/wiki/SCCC.jpg](http://turbosiso.googlecode.com/svn/wiki/SCCC.jpg)

**SCCC of coding rate 1/4**

> The turbo decoder for SCCCs uses a SISO RSC module and a SISO NSC module. The input of the SISO RSC module is represented by the intrinsic information of coded bits (the output of the RSC code). Its output, the extrinsic information of data bits of the RSC code, is used as intrinsic information for the SISO NSC module. This is the only information used by the SISO NSC module to compute the extrinsic information of coded bits (used in the next iteration as a priori information for the SISO RSC module) and the extrinsic information of data bits (used at the last iteration to recover the data bits).

![http://turbosiso.googlecode.com/svn/wiki/SCCC_tdec.jpg](http://turbosiso.googlecode.com/svn/wiki/SCCC_tdec.jpg)

**Turbo decoder for SCCC**

> The performance for SCCCs is shown below using max log MAP and log MAP algorithms. Worser performance of the max log MAP algorithm are due to the fact that the approximation ![http://turbosiso.googlecode.com/svn/wiki/form_28.png](http://turbosiso.googlecode.com/svn/wiki/form_28.png), used in the max log MAP algorithm, does not hold at low Signal to Noise Ratio (SNR).

![http://turbosiso.googlecode.com/svn/wiki/sccc_maxlogMAP.jpg](http://turbosiso.googlecode.com/svn/wiki/sccc_maxlogMAP.jpg)

**Performance of the turbo decoder for SCCC using max log MAP algorithm**

![http://turbosiso.googlecode.com/svn/wiki/sccc_logMAP.jpg](http://turbosiso.googlecode.com/svn/wiki/sccc_logMAP.jpg)

**Performance of the turbo decoder for SCCC using log MAP algorithm**

> Program used to obtain the figures: [sccc.cpp](http://turbosiso.googlecode.com/svn/trunk/sccc.cpp).

**Reference:**

S. Benedetto, D. Divsalar, G. Motorsi and F. Pollara, "A Soft-Input Soft-Output Maximum A posteriori (MAP) Module to Decode Parallel and Serial Concatenated Codes", TDA Progress Report, nov. 1996