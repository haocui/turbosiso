# Interleave Division Multiple Access with turbo multiuser detection #

The Interleave Division Multiple Access (IDMA) system shown below uses encoders (ENC) at emitter side, followed by interleavers (![http://turbosiso.googlecode.com/svn/wiki/form_32.png](http://turbosiso.googlecode.com/svn/wiki/form_32.png)), different for each user. The encoder can be a scrambler of length ![http://turbosiso.googlecode.com/svn/wiki/form_33.png](http://turbosiso.googlecode.com/svn/wiki/form_33.png) or a NSC code (generator polynomials ![http://turbosiso.googlecode.com/svn/wiki/form_26.png](http://turbosiso.googlecode.com/svn/wiki/form_26.png)) followed by a scrambler (![http://turbosiso.googlecode.com/svn/wiki/form_34.png](http://turbosiso.googlecode.com/svn/wiki/form_34.png)). After interleaving, a Zero Padding technique is used in order to eliminate the interference between adjacent blocks of chips (permuted coded bits). The block length is defined by the interleaver length (![http://turbosiso.googlecode.com/svn/wiki/form_27.png](http://turbosiso.googlecode.com/svn/wiki/form_27.png)). The chips are BPSK modulated and send into multipath channels (![http://turbosiso.googlecode.com/svn/wiki/form_30.png](http://turbosiso.googlecode.com/svn/wiki/form_30.png) multipaths). The channels are different for each user and have real, Rayleigh distributed, attenuations. The sum of all emitted signals is received and the received signal is affected by an Additive White Gaussian Noise (AWGN).

![http://turbosiso.googlecode.com/svn/wiki/IDMA_model.jpg](http://turbosiso.googlecode.com/svn/wiki/IDMA_model.jpg)

**IDMA system model**

> The turbo multi-user receiver uses a SISO MUD module and SISO descrambler modules, corresponding to each user. When the NSC code is followed by a scrambler, a SISO NSC module is used. This turbo multi-user receiver is a generalisation of turbo equalization principles. When there is a single user in the system, we have a turbo equalizer.

![http://turbosiso.googlecode.com/svn/wiki/IDMA_receiver.jpg](http://turbosiso.googlecode.com/svn/wiki/IDMA_receiver.jpg)

**Turbo multi-user receiver for IDMA systems**

Below the performance of the IDMA systems is shown in two different configurations: using only a scrambler at each user side and using a NSC code followed by a scrambler. In both cases the performance is limited by the performance of a mono-user transmission in an AWGN channel.

![http://turbosiso.googlecode.com/svn/wiki/IDMA_sGCD_8.jpg](http://turbosiso.googlecode.com/svn/wiki/IDMA_sGCD_8.jpg)

**Performance of the IDMA system with scrambler and using the simplified GCD for 8 users**

![http://turbosiso.googlecode.com/svn/wiki/ccIDMA_sGCD_5.jpg](http://turbosiso.googlecode.com/svn/wiki/ccIDMA_sGCD_5.jpg)

**Performance of the IDMA system with NSC code followed by a scrambler and using the simplified GCD for 5 users**

> It is also shown the EXIT diagram of the turbo multi-user receiver based on the simplified GCD and the scrambler at a SNR ![http://turbosiso.googlecode.com/svn/wiki/form_35.png](http://turbosiso.googlecode.com/svn/wiki/form_35.png). Thus it can be seen the iterative exchange of informations between the two SISO modules.

![http://turbosiso.googlecode.com/svn/wiki/EXIT_sGCD_descr.jpg](http://turbosiso.googlecode.com/svn/wiki/EXIT_sGCD_descr.jpg)

**EXIT diagram of the turbo multi-user receiver**

> Programs used to obtain the figures: [IDMA.cpp](http://turbosiso.googlecode.com/svn/trunk/IDMA.cpp), [TC\_mud.cpp](http://turbosiso.googlecode.com/svn/trunk/TC_mud.cpp) and [TC\_descrambler.cpp](http://turbosiso.googlecode.com/svn/trunk/TC_descrambler.cpp)

**Reference:**

L. Liu and L. Ping, ''Iterative detection of chip interleaved CDMA systems in multipath channels,`` Electronics letters, vol. 40, pp. 884-886, July 2004.