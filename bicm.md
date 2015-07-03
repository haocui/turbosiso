# Bit Interleaved Coded Modulation #

The Bit Interleaved Coded Modulation (BICM) system uses a Non-recursive non-Systematic Convolutional (NSC) code whose bits are interleaved (interleaver length is ![http://turbosiso.googlecode.com/svn/wiki/form_27.png](http://turbosiso.googlecode.com/svn/wiki/form_27.png)) and then mapped to complex symbols (![http://turbosiso.googlecode.com/svn/wiki/form_36.png](http://turbosiso.googlecode.com/svn/wiki/form_36.png) QAM modulation is used). Channel capacity is increased by using higher order modulations. Gray and binary mapping (following Matlab's terminology) are used. The channel is assumed flat fading, defined only by an attenuation. The attenuations are generated from a complex Gaussian distribution with variance ![http://turbosiso.googlecode.com/svn/wiki/form_14.png](http://turbosiso.googlecode.com/svn/wiki/form_14.png) on each dimension. The received signal is affected by an Additive White Gaussian Noise (AWGN). Note that the BICM system is a type of Serial Concatenated Code (SCC).

![http://turbosiso.googlecode.com/svn/wiki/bicm.jpg](http://turbosiso.googlecode.com/svn/wiki/bicm.jpg)

**BICM system model**

The turbo receiver uses a SISO demapper module and a SISO NSC module. The SISO demapper module is a soft demodulator using at its input the received complex signal and the a priori information of coded bits. Its output is the extrinsic information of coded bits. The extrinsic information becomes, after deinterleaving, intrinsic information at the SISO NSC module input. The SISO NSC module outputs an extrinsic information of coded bits (used as a priori information in the next iteration at the SISO demapper input) and an extrinsic information of data bits (used to recover the data bits at the last iteration).

![http://turbosiso.googlecode.com/svn/wiki/bicm_rec.jpg](http://turbosiso.googlecode.com/svn/wiki/bicm_rec.jpg)

**Turbo receiver for BICM systems**

Performance of BICM are shown below when using Gray or binary mapping with log MAP and max log MAP algorithms in both SISO modules. Note that Gray mapping has worse performance than binary mapping when a turbo receiver is used.

![http://turbosiso.googlecode.com/svn/wiki/BICM_Gray_logMAP.jpg](http://turbosiso.googlecode.com/svn/wiki/BICM_Gray_logMAP.jpg)

**Performances of the BICM using Gray mapping and logMAP algorithm**

![http://turbosiso.googlecode.com/svn/wiki/BICM_Gray_maxlogMAP.jpg](http://turbosiso.googlecode.com/svn/wiki/BICM_Gray_maxlogMAP.jpg)

**Performances of the BICM using Gray mapping and maxlogMAP algorithm**

![http://turbosiso.googlecode.com/svn/wiki/BICM_binary_logMAP.jpg](http://turbosiso.googlecode.com/svn/wiki/BICM_binary_logMAP.jpg)

**Performances of the BICM using binary mapping and logMAP algorithm**

![http://turbosiso.googlecode.com/svn/wiki/BICM_binary_maxlogMAP.jpg](http://turbosiso.googlecode.com/svn/wiki/BICM_binary_maxlogMAP.jpg)

**Performances of the BICM using binary mapping and maxlogMAP algorithm**

Simulation times for programs written in C++ and Matlab are shown below. We have used a Matlab program in order to modulate with binary mapping, since in IT++ only Gray mapping is provided. The simulation times in both cases (with log MAP and max log MAP algorithms) are very close due to the fact that the Matlab program uses several ''mex`` functions (besides the functions implementing the SISO modules).

Simulation times for BICM in C++ and Matlab

|  | **C++** | **Matlab** |
|:-|:--------|:-----------|
| **max log MAP** | 43 min, 51 sec | 52 min, 20 sec |
| **log MAP** | 1 hr, 10 min, 53 sec | 1 hr, 16 min, 33 sec |

Programs used to obtain the figures: [BICM.cpp](http://turbosiso.googlecode.com/svn/trunk/BICM.cpp) and [BICM.m](http://turbosiso.googlecode.com/svn/trunk/BICM.m).

**Reference:**

A. Tonello, ''Space-time bit-interleaved coded modulation with an iterative decoding strategy,`` in Vehicular Technology Conference, vol. 1, pp. 473-478 vol.1, 2000