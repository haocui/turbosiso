% -------------------------------------------------------------------
% Parallel Concatenated Convolutional Codes (PCCC) of coding rate 1/3
% -------------------------------------------------------------------
% Reference: S. Benedetto, D. Divsalar, G. Motorsi and F. Pollara, "A Soft-Input Soft-Output Maximum A posteriori (MAP) Module
% to Decode Parallel and Serial Concatenated Codes", TDA Progress Report, nov. 1996

clear all
clc

%general parameters
threshold = 50;
map_metric = 'logMAP';
gen = [37 21];%first polynomial is the feedback polynomial
constraint_len = 5;
nb_errors_lim = 1500;
nb_bits_lim = 1e6;
perm_len = pow2(14);%total number of bits in a block (with tail)
nb_iter = 10;%number of iterations in the turbo decoder
EbN0_dB = 0:0.1:5;
R = 1/3;%coding rate (non punctured PCCC)
Ec = 1;%coded bit energy

%other parameters
bin_gen = de2bi(base2dec(int2str(gen.'), 8), 'left-msb');
filename = ['Res/pccc_' map_metric];
mem_len = constraint_len-1;
nb_bits = perm_len-mem_len;%number of bits in a block (without tail)
sigma2 = (0.5*Ec/R)*10.^(-EbN0_dB/10);%N0/2
perm = zeros(1, perm_len);
inv_perm = zeros(1, perm_len);
tail = zeros(1, mem_len);
cod_bits_len = perm_len*length(gen);
rec_len = perm_len/R;
coded_bits = zeros(1, rec_len);
rec = zeros(1, rec_len);
dec1_intrinsic_coded = zeros(1, cod_bits_len);
dec2_intrinsic_coded = zeros(1, cod_bits_len);
extrinsic_coded = zeros(1, perm_len);
extrinsic_data = zeros(1, perm_len);
rec_bits = zeros(1, perm_len);
snr_len = length(EbN0_dB);
BER = zeros(nb_iter,snr_len);

%Recursive Systematic Convolutional Code trellis
trellis = poly2trellis(constraint_len, gen, gen(1));

%main loop
txtprogressbar;
for en=1:snr_len
    Lc = -2/sigma2(en);%normalisation factor for intrinsic information (take into account the BPSK mapping)
    nb_errors = 0;
    nb_blocks = 0;
    while ((nb_errors<nb_errors_lim) && (nb_blocks*nb_bits<nb_bits_lim))%if at the last iteration the nb. of errors is inferior to lim, then process another block
        %permutation
        perm = randperm(perm_len);
        %inverse permutation
        [ignore, inv_perm] = sort(perm);

        %bits generation
        bits = randint(1, nb_bits, 2);

        %parallel concatenated convolutional code
        [cod1_bits final_state] = convenc(bits, trellis, 0);%initial state is zero
        %tail is added here to information bits to close the trellis
        cod_mem = de2bi(final_state, mem_len, 'left-msb');
        for n=1:mem_len
            tail(n) = mod(cod_mem*bin_gen(1,2:end)', 2);%xor sum
            cod1_bits = [cod1_bits tail(n)];%systematic bit
            cod1_bits = [cod1_bits mod(cod_mem*bin_gen(2,2:end)', 2)];%parity output
            cod_mem(2:end) = cod_mem(1:end-1);%shift memory
            cod_mem(1) = 0;
        end
        cod2_input = [bits tail];
        cod2_bits = convenc(cod2_input(perm), trellis, 0);%initial state is zero, no tail
        for n=1:perm_len %output with no puncturing
            coded_bits(3*n-2) = cod2_input(n);%systematic output
            coded_bits(3*n-1) = cod1_bits(2*n);%first parity output
            coded_bits(3*n) = cod2_bits(2*n);%second parity output
        end
        clear cod1_bits

        %BPSK modulation (1->-1,0->+1) + AWGN channel
        rec = (1-2*coded_bits)+sqrt(sigma2(en))*randn(1, rec_len);

        %form input for SISO blocks
        for n=1:perm_len
            dec1_intrinsic_coded(2*n-1) = Lc*rec(3*n-2);
            dec1_intrinsic_coded(2*n) = Lc*rec(3*n-1);
            dec2_intrinsic_coded(2*n-1) = 0;%systematic output of the CC is already used in decoder1
            dec2_intrinsic_coded(2*n) = Lc*rec(3*n);
        end
        %turbo decoder
        apriori_data = zeros(1, perm_len);%a priori LLR for information bits
        for n=1:nb_iter
            %first decoder
            [extrinsic_coded extrinsic_data] = C_SISOrsc(dec1_intrinsic_coded, apriori_data, bin_gen, 1, map_metric);%with tail
            %interleave
            apriori_data = extrinsic_data(perm);
            %threshold
            apriori_data(apriori_data>threshold) = threshold;
            apriori_data(apriori_data<-threshold) = -threshold;
            %second decoder
            [extrinsic_coded extrinsic_data] = C_SISOrsc(dec2_intrinsic_coded, apriori_data, bin_gen, 0, map_metric);%no tail

            %decision
            apriori_data = apriori_data+extrinsic_data;%a posteriori information
            rec_bits = double(apriori_data(inv_perm)>0);%take into account the BPSK mapping
            %count errors
            BER(n,en) = BER(n,en)+nnz(bits-rec_bits(1:nb_bits))/nb_bits;

            %deinterleave for the next iteration
            apriori_data = extrinsic_data(inv_perm);
        end%end iterations
        nb_errors = nb_errors+nnz(bits-rec_bits(1:nb_bits));%get number of errors at the last iteration
        nb_blocks = nb_blocks+1;
    end%end blocks (while loop)

    %compute BER over all tx blocks
    BER(:,en) = BER(:,en)/nb_blocks;

    %show progress
    txtprogressbar(en/snr_len);
end

% figure
% semilogy(EbN0_dB, BER, 'o-')
% grid on
% xlabel('E_b/N_0 [dB]')
% ylabel('BER')

save(filename, 'BER', 'EbN0_dB', 'gen', 'R', 'nb_iter', 'perm_len', 'nb_errors_lim', 'nb_bits_lim')


