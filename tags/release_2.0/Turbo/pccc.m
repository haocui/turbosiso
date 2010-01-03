% -------------------------------------------------------------------
% Parallel Concatenated Convolutional Codes (PCCC) of coding rate 1/3
% -------------------------------------------------------------------
% Reference: S. Benedetto, D. Divsalar, G. Motorsi and F. Pollara, 
% "A Soft-Input Soft-Output Maximum A posteriori (MAP) Module
% to Decode Parallel and Serial Concatenated Codes", TDA Progress Report, Nov. 1996

clear all
clc

%general parameters
threshold = 100;
map_metric = 'maxlogMAP';
%gen = [13 15];%first polynomial is the feedback polynomial
gen = [37 21];
%gen = [7 5];
nb_errors_lim = 3000;
perm_len_lim = 1e6;
perm_len = pow2(16);%total number of bits in a block (with tail)
nb_iter = 10;%number of iterations in the turbo decoder
EbN0_dB = 0:0.1:5;
R = 1/3;%coding rate (non punctured PCCC)
Ec = 1;%coded bit energy

%other parameters
bin_gen = de2bi(base2dec(int2str(gen.'), 8), 'left-msb');
constraint_len = size ( bin_gen, 2 );
sova_win_len = 5*constraint_len;
filename = ['Res/pccc_' map_metric  '_' num2str(gen(1)) '_' num2str(gen(2)) '.mat'];
sigma2 = (0.5*Ec/R)*10.^(-EbN0_dB/10);%N0/2
perm = zeros(1, perm_len);
inv_perm = zeros(1, perm_len);
cod_bits_len = perm_len*length(gen);
rec_len = perm_len/R;
coded_bits = zeros(1, rec_len);
rec = zeros(1, rec_len);
dec1_intrinsic_coded = zeros(1, cod_bits_len);
dec2_intrinsic_coded = zeros(1, cod_bits_len);
extrinsic_coded = zeros(1, perm_len);
extrinsic_data = zeros(1, perm_len);
snr_len = length(EbN0_dB);
BER = zeros(nb_iter, snr_len);

%Recursive Systematic Convolutional Code trellis
trellis = poly2trellis(constraint_len, gen, gen(1));

%main loop
txtprogressbar;
for en=1:snr_len
    nb_errors = 0;
    nb_blocks = 0;
    while ((nb_errors<nb_errors_lim) && (nb_blocks*perm_len<perm_len_lim))
        %permutation
        perm = randperm(perm_len);
        %inverse permutation
        [ignore, inv_perm] = sort(perm);

        %bits generation
        bits = randint(1, perm_len, 2);

        %parallel concatenated convolutional code
        cod1_bits = convenc(bits, trellis, 0);%initial state is zero
        cod2_bits = convenc(bits(perm), trellis, 0);
        
        %multiplexer
        coded_bits(1:3:end) = bits;
        coded_bits(2:3:end) = cod1_bits(2:2:end);
        coded_bits(3:3:end) = cod2_bits(2:2:end);

        %BPSK modulation (1->-1,0->+1) + AWGN channel
        rec = (1-2*coded_bits)+sqrt(sigma2(en))*randn(1, rec_len);

        %form input for SISO blocks
        dec1_intrinsic_coded(1:2:end) = -2/sigma2(en)*rec(1:3:end);
        dec1_intrinsic_coded(2:2:end) = -2/sigma2(en)*rec(2:3:end);
        dec2_intrinsic_coded(2:2:end) = -2/sigma2(en)*rec(3:3:end);
        
        %turbo decoder
        apriori_data = zeros(1, perm_len);%a priori LLR for information bits
        for n=1:nb_iter
            
            %first decoder
            [extrinsic_coded extrinsic_data] = ...
                C_SISOrsc(dec1_intrinsic_coded, apriori_data, bin_gen, 0, map_metric);%no tail
            %extrinsic_data = C_SISOsova(dec1_intrinsic_coded, apriori_data, bin_gen, sova_win_len);
            
            %interleave
            apriori_data = extrinsic_data(perm);
            
            %second decoder
            [extrinsic_coded extrinsic_data] = ...
                C_SISOrsc(dec2_intrinsic_coded, apriori_data, bin_gen, 0, map_metric);%no tail
            %extrinsic_data = C_SISOsova(dec2_intrinsic_coded, apriori_data, bin_gen, sova_win_len);

            %decision on a posteriori information
            apriori_data = apriori_data+extrinsic_data;
            rec_bits = double(apriori_data(inv_perm)>0);
            
            %deinterleave for the next iteration
            apriori_data = extrinsic_data(inv_perm);

            %threshold
            apriori_data(apriori_data>threshold) = threshold;
            apriori_data(apriori_data<-threshold) = -threshold;
            
            %count errors
            BER(n,en) = BER(n,en)+nnz(bits-rec_bits)/perm_len;

        end%end iterations
        nb_errors = nb_errors+nnz(bits-rec_bits(1:perm_len));%get number of errors at the last iteration
        nb_blocks = nb_blocks+1;
    end%end blocks (while loop)

    %compute BER over all tx blocks
    BER(:,en) = BER(:,en)/nb_blocks;

    %show progress
    txtprogressbar(en/snr_len);
end

figure
semilogy(EbN0_dB, BER, 'o-')
grid on
xlabel('E_b/N_0 [dB]')
ylabel('BER')

save(filename, 'BER', 'EbN0_dB', 'gen', 'R', 'nb_iter', ...
    'perm_len', 'nb_errors_lim', 'perm_len_lim')
