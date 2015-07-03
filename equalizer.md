# Turbo equalizer #

The coded transmission is realized using a Non-recursive and non-Systematic Convolutional (NSC) code whose bits, after interleaving (of length ![http://turbosiso.googlecode.com/svn/wiki/form_27.png](http://turbosiso.googlecode.com/svn/wiki/form_27.png)), are BPSK modulated and send into a multipath channel. The multipath channel has real coefficients and must have at least ![http://turbosiso.googlecode.com/svn/wiki/form_23.png](http://turbosiso.googlecode.com/svn/wiki/form_23.png) multipaths in order to use an iterative algorithm at the receiver side. The multipath channel can be seen as a coder of coding rate ![http://turbosiso.googlecode.com/svn/wiki/form_21.png](http://turbosiso.googlecode.com/svn/wiki/form_21.png), thus the coded transmission being a type of Serial Concatenated Code (SCC). In our case, the multipath channel has Rayleigh distributed coefficients and has ![http://turbosiso.googlecode.com/svn/wiki/form_30.png](http://turbosiso.googlecode.com/svn/wiki/form_30.png) multipaths. Optionally a precoder can be used at channel input in order to obtain an equivalent recursive channel. The precoder polynomial is ![http://turbosiso.googlecode.com/svn/wiki/form_31.png](http://turbosiso.googlecode.com/svn/wiki/form_31.png). At the channel output, the received signal is affected by an Additive White Gaussian Noise (AWGN).

![http://turbosiso.googlecode.com/svn/wiki/coded_tx.jpg](http://turbosiso.googlecode.com/svn/wiki/coded_tx.jpg)

**Coded transmission**

> The turbo equalizer uses a SISO equalizer module and a SISO NSC module. The SISO equalizer module uses at its input the received signal and the a priori information of coded bits. It delivers at its output the extrinsic information of channel input symbols (coded bits of the NSC code). After deinterleaving, the extrinsic information is used as intrinsic information for the SISO NSC module. The SISO NSC module computes the extrinsic information of coded bits (used in the next iteration as a priori information at the SISO equalizer module input) and the extrinsic information of data bits (used at the last iteration to recover the data bits).

![http://turbosiso.googlecode.com/svn/wiki/turbo_eq.jpg](http://turbosiso.googlecode.com/svn/wiki/turbo_eq.jpg)

**Turbo equalizer**

> The performance of the turbo equalizer without precoder is shown below. As for the SCCCs, the max log MAP algorithm has worser performance than the log MAP algorithm at low SNR.

![http://turbosiso.googlecode.com/svn/wiki/equalizer_maxlogMAP.jpg](http://turbosiso.googlecode.com/svn/wiki/equalizer_maxlogMAP.jpg)

**Performances of the turbo equalizer without precoder using max log MAP algorithm**

![http://turbosiso.googlecode.com/svn/wiki/equalizer_logMAP.jpg](http://turbosiso.googlecode.com/svn/wiki/equalizer_logMAP.jpg)

**Performances of the turbo equalizer without precoder using log MAP algorithm**

> The performance of the turbo equalizer with precoder is shown below. It can be seen that, if the SNR is high enough, the performance is better than in the case when no precoder is used.

![http://turbosiso.googlecode.com/svn/wiki/prec_equalizer_maxlogMAP.jpg](http://turbosiso.googlecode.com/svn/wiki/prec_equalizer_maxlogMAP.jpg)

**Performance of the turbo equalizer with precoder using max log MAP algorithm**

![http://turbosiso.googlecode.com/svn/wiki/prec_equalizer_logMAP.jpg](http://turbosiso.googlecode.com/svn/wiki/prec_equalizer_logMAP.jpg)

**Performance of the turbo equalizer with precoder using log MAP algorithm**

> Program used to obtain the figures: [equalizer.cpp](http://turbosiso.googlecode.com/svn/trunk/equalizer.cpp).

**Reference:**

R. Koetter, A. C. Singer, and M. Tuchler, ''Turbo equalization: an iterative equalization and decoding technique for coded data transmision,`` IEEE Signal Processing Magazine, pp. 67-80, Jan. 2004