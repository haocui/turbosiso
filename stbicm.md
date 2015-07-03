# Space Time Bit Interleaved Coded Modulation #

The Space Time Bit Interleaved Coded Modulation (ST-BICM) system is a generalisation of the BICM system (here the generator polynomials are _0133_ _0171_), using, after the modulator, an ST block code. Thus, transmit diversity can be achieved. Further, capacity increase of communication system is achieved by transmitting different symbols through different emission antennas. We have considered two ST codes: Golden code (among the best known ST block codes) and Alamouti code (orthogonal ST block code). The number of emission antennas is ![http://turbosiso.googlecode.com/svn/wiki/form_23.png](http://turbosiso.googlecode.com/svn/wiki/form_23.png) (equal to the number of reception antennas). The order of the QAM modulation is chosen so that the spectral efficiency is constant ![http://turbosiso.googlecode.com/svn/wiki/form_38.png](http://turbosiso.googlecode.com/svn/wiki/form_38.png): for Golden code QAM modulation is used and for Alamouti code ![http://turbosiso.googlecode.com/svn/wiki/form_36.png](http://turbosiso.googlecode.com/svn/wiki/form_36.png) QAM modulation is used. The Multiple-Input Multiple-Output (MIMO) channel is considered flat fading with complex attenuations.

![http://turbosiso.googlecode.com/svn/wiki/stbicm.jpg](http://turbosiso.googlecode.com/svn/wiki/stbicm.jpg)

**ST-BICM system model**

The turbo receiver uses a SISO demapper module and a SISO NSC module. The SISO demapper module decodes both the ST block code and the MIMO channel. It uses at its input the received complex signal and the a priori information of coded bits. Its output is represented by the extrinsic information of coded bits which, after deinterleaving, becomes intrinsic information for the SISO NSC. The SISO NSC module outputs an extrinsic information of coded bits (used as a priori information in the next iteration at the SISO demapper input) and an extrinsic information of data bits (used to recover the data bits at the last iteration).

![http://turbosiso.googlecode.com/svn/wiki/bicm_rec.jpg](http://turbosiso.googlecode.com/svn/wiki/bicm_rec.jpg)

**Turbo receiver for ST-BICM systems**

Performance of the ST-BICM using different algorithms in the SISO demapper and different ST block codes are presented below. Note that the peformances of GA and MMSE PIC algorithms are very close, but the algorithm complexity is a function of system parameters.

![http://turbosiso.googlecode.com/svn/wiki/Golden_maxlogMAP_Hassibi.jpg](http://turbosiso.googlecode.com/svn/wiki/Golden_maxlogMAP_Hassibi.jpg)

**Performance of the ST-BICM using Golden code and max log MAP algorithm**

![http://turbosiso.googlecode.com/svn/wiki/Golden_GA.jpg](http://turbosiso.googlecode.com/svn/wiki/Golden_GA.jpg)

**Performance of the ST-BICM using Golden code and GA algorithm**

![http://turbosiso.googlecode.com/svn/wiki/Golden_sGA.jpg](http://turbosiso.googlecode.com/svn/wiki/Golden_sGA.jpg)

**Performance of the ST-BICM using Golden code and simplified GA algorithm**

![http://turbosiso.googlecode.com/svn/wiki/Golden_mmsePIC.jpg](http://turbosiso.googlecode.com/svn/wiki/Golden_mmsePIC.jpg)

**Performance of the ST-BICM using Golden code and MMSE PIC algorithm**

![http://turbosiso.googlecode.com/svn/wiki/Golden_zfPIC.jpg](http://turbosiso.googlecode.com/svn/wiki/Golden_zfPIC.jpg)

**Performance of the ST-BICM using Golden code and ZF PIC algorithm**

![http://turbosiso.googlecode.com/svn/wiki/Alamouti_maxlogMAP.jpg](http://turbosiso.googlecode.com/svn/wiki/Alamouti_maxlogMAP.jpg)

**Performance of the ST-BICM using Alamouti code and max log MAP algorithm**

Program used to obtain the figures: [STBICM.cpp](http://turbosiso.googlecode.com/svn/trunk/STBICM.cpp).

**Reference:**

B. Cristea, ''Turbo receivers for Space-Time BICM``, to be published in IEEE Transactions on Wireless Communications