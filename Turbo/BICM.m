%
% Implements BICM using a turbo receiver with a SISO demapper module and a SISO NSC module.
%
% Reference: A. Tonello, ''Space-time bit-interleaved coded modulation with an iterative decoding strategy,`` 
% in Vehicular Technology Conference, vol. 1, pp. 473-478 vol.1, 2000

%general parameters
const_size = 16;%constellation size
select_mapping = 'binary';
threshold_value = 50;
demapper_method = 'maxlogMAP';
map_metric = 'maxlogMAP';
gen = [37 21];
constraint_length = 5;
nb_errors_lim = 1500;
nb_bits_lim = 1e6;
perm_len = pow2(14);%permutation length
nb_iter = 10;%number of iterations in the turbo decoder
EbN0_dB = 0:20;
R = 1/2;%coding rate of FEC
Es = 1;%mean symbol energy    

%QAM modulator class
bits_per_symbol = log2(const_size);
h_QAMmod = modem.qammod('M', const_size, 'PhaseOffset', 0, 'SymbolOrder', select_mapping, 'InputType', 'bit');
constellation = h_QAMmod.Constellation;%complex constellation
bin_constellation = de2bi(h_QAMmod.SymbolMapping, bits_per_symbol, 'left-msb');%binary constellation
scale = modnorm(constellation, 'avpow', 1);%normalisation factor    

%other parameters
filename = ['Res/BICM_' map_metric '_' select_mapping];
nb_symb = perm_len/bits_per_symbol;
nb_bits_tail = perm_len/length(gen);
nb_bits = nb_bits_tail-(constraint_length-1);%number of bits in a block (without tail)
sigma2 = (0.5*Es/(R*bits_per_symbol))*(10.^(-EbN0_dB/10));%N0/2
%SISO NSC
nsc_apriori_data = zeros(1, nb_bits_tail);
%decision
snr_len = length(EbN0_dB);
BER = zeros(nb_iter,snr_len);
    
%CC
bin_gen = de2bi(base2dec(int2str(gen.'), 8), 'left-msb');
constraint_len = size(bin_gen, 2);
trellis = poly2trellis(constraint_len, gen);
    
%main loop
txtprogressbar;
for en=1:snr_len
        nb_errors = 0;
        nb_blocks = 0;
        while ((nb_errors<nb_errors_lim) && (nb_blocks*nb_bits<nb_bits_lim))%if at the last iteration the nb. of errors is inferior to lim, then process another block
            %permutation
            perm = randperm(perm_len);
            %inverse permutation
            [ignore, inv_perm] = sort(perm);

            %bits generation
            bits = randint(1, nb_bits, 2);

            %convolutional code
            [coded_bits, final_state] = convenc(bits, trellis, 0);%initial state is zero
            tail = convenc(zeros(1, constraint_length-1), trellis, final_state);%tail is added to information bits to close the trellis         

            %permutation+QAM modulation
            temp = [coded_bits tail];
            em = scale*modulate(h_QAMmod, temp(perm)').';            
            
			%flat-fading channel
            ch_attenuations = sqrt(0.5)*(randn(1, nb_symb)+j*randn(1, nb_symb));
            rec = ch_attenuations.*em+sqrt(sigma2(en))*(randn(1, nb_symb)+j*randn(1, nb_symb));

            %turbo receiver
            demod_apriori_data = zeros(1, perm_len);%a priori information of emitted bits
            for n=1:nb_iter
                %first decoder
                demod_extrinsic_data = C_SISOdemapper(rec, demod_apriori_data, ch_attenuations, sigma2(en), bits_per_symbol, ...
                    scale*constellation, bin_constellation, demapper_method);
               
                %deinterleave
                nsc_intrinsic_coded = demod_extrinsic_data(inv_perm);
                
                %threshold
                nsc_intrinsic_coded(nsc_intrinsic_coded>threshold_value) = threshold_value;
                nsc_intrinsic_coded(nsc_intrinsic_coded<-threshold_value) = -threshold_value;

                %second decoder
                [nsc_extrinsic_coded nsc_extrinsic_data] = C_SISOnsc(nsc_intrinsic_coded, nsc_apriori_data, bin_gen, ...
                    1, 1, map_metric);               

                %decision
                rec_bits = double(nsc_extrinsic_data>0);%suppose that a priori info is zero
                %count errors
                BER(n,en) = BER(n,en)+nnz(bits-rec_bits(1:nb_bits))/nb_bits;

                %interleave
                demod_apriori_data = nsc_extrinsic_coded(perm);
            end %end iterations
            nb_errors = nb_errors+nnz(bits-rec_bits(1:nb_bits));%get number of errors at the last iteration
            nb_blocks = nb_blocks+1;
        end %end blocks (while loop)

        %compute BER over all tx blocks
        BER(:,en) = BER(:,en)/nb_blocks;

        %show progress
        txtprogressbar(en/snr_len);
end

save(filename, 'BER', 'EbN0_dB', 'gen', 'R', 'nb_iter', 'nb_errors_lim', 'nb_bits_lim', 'perm_len', 'const_size')
