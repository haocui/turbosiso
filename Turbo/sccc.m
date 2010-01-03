% -------------------------------------------------------------------
% Serial Concatenated Convolutional Codes (SCCC) of coding rate 1/3
% -------------------------------------------------------------------

clear all
clc

%general parameters
threshold = 100;
map_metric = 'logMAP';
gen = [7 5];%recursive polynomial first
outer_code_puncturing_pattern = [1 1 1 1]; 
nb_errors_lim = 3000;
nb_bits_lim = 1e6;
perm_len = 3*pow2(15);%total number of bits in a block
nb_iter = 10;%number of iterations in the turbo decoder
EbN0_dB = 0:0.1:5;
Ec = 1;%coded bit energy

%other parameters
bin_gen = de2bi(base2dec(int2str(gen.'), 8), 'left-msb');
code_outputs = length(gen);
outer_code_rate = (length(outer_code_puncturing_pattern)/code_outputs)/...
    length(find(outer_code_puncturing_pattern~=0));
inner_code_rate = 1/length(gen);
R = outer_code_rate*inner_code_rate;%coding rate of the SCCC
constraint_len = size ( bin_gen, 2 );
sova_win_len = 5*constraint_len;
filename = ['Res/sccc5_' map_metric  '_' num2str(gen(1)) '_' num2str(gen(2)) '.mat'];
sigma2 = (0.5*Ec/R)*10.^(-EbN0_dB/10);%N0/2
perm = zeros(1, perm_len);
inv_perm = zeros(1, perm_len);
nb_bits = perm_len*outer_code_rate;
rec_len = perm_len/inner_code_rate;
if (int32(nb_bits)~=nb_bits) || (int32(rec_len)~=rec_len)
    error('sccc', 'specified permutation length, %d, does not accept these code rates, %d, %d', ...
        outer_code_rate, inner_code_rate);
end
rec = zeros(1, rec_len);
snr_len = length(EbN0_dB);
BER = zeros(nb_iter, snr_len);
expanded_puncturing_pattern = logical(kron(ones(1, code_outputs*nb_bits/length(outer_code_puncturing_pattern)), ...
    outer_code_puncturing_pattern));
inner_code_extrinsic_data = zeros(1, perm_len);
outer_code_apriori_data = zeros(1, nb_bits);
outer_code_intrinsic_coded = zeros(1, nb_bits*code_outputs);

%Recursive Systematic Convolutional Code trellis
trellis = poly2trellis(constraint_len, gen, gen(1));

%main loop
txtprogressbar;
for en=1:snr_len
    nb_errors = 0;
    nb_blocks = 0;
    while ((nb_errors<nb_errors_lim) && (nb_blocks*perm_len<nb_bits_lim))
        %permutation
        perm = randperm(perm_len);
        %inverse permutation
        [ignore, inv_perm] = sort(perm);

        %bits generation
        bits = randint(1, nb_bits, 2);

        %serial concatenated convolutional code
        temp = convenc(bits, trellis, 0);%initial state is zero
        outer_code_bits = temp(expanded_puncturing_pattern);%puncturing
        inner_code_bits = convenc(outer_code_bits(perm), trellis, 0);

        %BPSK modulation (1->-1,0->+1) + AWGN channel
        rec = (1-2*inner_code_bits)+sqrt(sigma2(en))*randn(1, rec_len);
        
        %turbo decoder
        inner_code_apriori_data = zeros(1, perm_len);%a priori LLR for information bits
        inner_code_intrinsic_coded = (-2/sigma2(en))*rec;
        for n=1:nb_iter
            
            %inner decoder
            [inner_code_extrinsic_parity inner_code_extrinsic_data] = ...
                C_SISOrsc(inner_code_intrinsic_coded, inner_code_apriori_data, bin_gen, 0, map_metric);%no tail            
            
            %deinterleave
            outer_code_intrinsic_coded(expanded_puncturing_pattern) = ...
                inner_code_extrinsic_data(inv_perm);
            
            %outer decoder
            [outer_code_extrinsic_parity outer_code_extrinsic_data] = ...
                C_SISOrsc(outer_code_intrinsic_coded, outer_code_apriori_data, bin_gen, 0, map_metric);%no tail
            
            %deinterleave for the next iteration
            temp(1:2:end) = outer_code_extrinsic_data;
            temp(2:2:end) = outer_code_extrinsic_parity;
            temp2 = temp(expanded_puncturing_pattern);
            inner_code_apriori_data = temp2(perm);
            
            %threshold
            inner_code_apriori_data(inner_code_apriori_data>threshold) = threshold;
            inner_code_apriori_data(inner_code_apriori_data<-threshold) = -threshold;
            
            %count errors
            rec_bits = double(outer_code_extrinsic_data>0);
            BER(n,en) = BER(n,en)+nnz(bits-rec_bits)/nb_bits;

        end%end iterations
        nb_errors = nb_errors+nnz(bits-rec_bits);%get number of errors at the last iteration
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
    'perm_len', 'nb_errors_lim', 'nb_bits_lim', 'outer_code_puncturing_pattern')
