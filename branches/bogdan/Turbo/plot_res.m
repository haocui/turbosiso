%simple script for plotting the results (BER vs EbN0_dB) from an *.it file

clear all

filename = input('File to load: ', 's');
itload(filename);
figure
semilogy(EbN0_dB,  ber, 'o-')
grid on
xlabel('E_b/N_0 [dB]')
ylabel('BER')
