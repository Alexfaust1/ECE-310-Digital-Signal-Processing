%% Alexander Faust
%
% ECE-310 Problem Set 4 - Quantization
%
% November 26, 2023
clc; clear; close all;

%% Question 4 
ellip_order = 8;                        % Elliptic filter order
Rp_ripple = 1.5;                        % Passband ripple [dB]
Rs_ripple = 30;                         % Stopband ripple [dB]
Wn = [0.3 0.6];                        % Passband length

[z, p, k] = ellip((ellip_order/2), Rp_ripple, Rs_ripple, Wn);

% Part a) - Convert to transfer function form
N = 2048;                           % Number of points for filter
[NUM, DEN] = zp2tf(z, p, k);
[tf, W] = freqz(NUM, DEN, N);

% List frequency to plot against:
fs = 10;              
frequencies = linspace(0, fs/2, 1000);

% Create figure of transfer function
figure;
plot(W/pi, 20*log10(abs(tf)));
title("Magnitude Response of Transfer Function");
xlabel("Normalized digital radian frequency");
ylabel("Magnitude [dB]");

% Part b) - Obtain SOS form
[SOS_up, G_up] = zp2sos(z, p, k, 'UP','inf');
[SOS_down, G_down] = zp2sos(z, p, k, 'DOWN', 'inf');

% Part c) - Check magnitude of poles and zeros from SOS
[SOS_up_poles, SOS_up_zeros, SOS_up_k] = sos2zp(SOS_up);
[SOS_down_poles, SOS_down_zeros, SOS_down_k] = sos2zp(SOS_down);

MAG_poles_up = abs(SOS_up_poles);
MAG_poles_down = abs(SOS_down_poles);
MAG_zeros_up = abs(SOS_up_zeros);
MAG_zeros_down = abs(SOS_down_zeros);

% Part d) - Create H1(z), H2(z), ... H_L(z):
NUM1_up = SOS_up(1, 1:3); %B's
NUM2_up = SOS_up(2, 1:3);
NUM3_up = SOS_up(3, 1:3);
NUM4_up = SOS_up(4, 1:3);

DEN1_up = SOS_up(1, 4:end); % A's
DEN2_up = SOS_up(2, 4:end);
DEN3_up = SOS_up(3, 4:end);
DEN4_up = SOS_up(4, 4:end);

% Compute convolution terms for numerator terms i.e. B's
conv_NUM1_NUM2_up = conv(NUM1_up, NUM2_up);
conv_NUM1_NUM2_NUM3_up = conv(conv_NUM1_NUM2_up, NUM3_up);
conv_NUM1_NUM2_NUM3_NUM4_up = conv(conv_NUM1_NUM2_NUM3_up, NUM4_up);

% Compute convolution terms for denominator terms i.e. A's
conv_DEN1_DEN2_up = conv(DEN1_up, DEN2_up);
conv_DEN1_DEN2_DEN3_up = conv(conv_DEN1_DEN2_up, DEN3_up);
conv_DEN1_DEN2_DEN3_DEN4_up = conv(conv_DEN1_DEN2_DEN3_up, DEN4_up);

% Construct transfer functions for each stage:
H0_up = freqz(G_up, DEN1_up, N);
H1_up = freqz(G_up * NUM1_up, conv_DEN1_DEN2_up, N);
H2_up = freqz(G_up * conv_NUM1_NUM2_up, conv_DEN1_DEN2_DEN3_up, N);
H3_up = freqz(G_up * conv_NUM1_NUM2_NUM3_up, conv_DEN1_DEN2_DEN3_DEN4_up, N);
H4_up = freqz(G_up * conv_NUM1_NUM2_NUM3_NUM4_up, conv_DEN1_DEN2_DEN3_DEN4_up, N);

figure;
hold on;
plot(W/pi, 20*log10(abs(H0_up)));
plot(W/pi, 20*log10(abs(H1_up)));
plot(W/pi, 20*log10(abs(H2_up)));
plot(W/pi, 20*log10(abs(H3_up)));
plot(W/pi, 20*log10(abs(H4_up)));
hold off;
title("Plot of H1(z) H2(z) ... for SOS UP");
xlabel("Normalized digital radian frequency");
ylabel("Magnitude [dB]");

% Do the same as above but for SOS_down: %
NUM1_down = SOS_down(1, 1:3);
NUM2_down = SOS_down(2, 1:3);
NUM3_down = SOS_down(3, 1:3);
NUM4_down = SOS_down(4, 1:3);

DEN1_down = SOS_down(1, 4:end);
DEN2_down = SOS_down(2, 4:end);
DEN3_down = SOS_down(3, 4:end);
DEN4_down = SOS_down(4, 4:end);

% Compute convolution terms for numerator terms
conv_NUM1_NUM2_down = conv(NUM1_down, NUM2_down);
conv_NUM1_NUM2_NUM3_down = conv(conv_NUM1_NUM2_down, NUM3_down);
conv_NUM1_NUM2_NUM3_NUM4_down = conv(conv_NUM1_NUM2_NUM3_down, NUM4_down);

% Compute convolution terms for denominator terms
conv_DEN1_DEN2_down = conv(DEN1_down, DEN2_down);
conv_DEN1_DEN2_DEN3_down = conv(conv_DEN1_DEN2_down, DEN3_down);
conv_DEN1_DEN2_DEN3_DEN4_down = conv(conv_DEN1_DEN2_DEN3_down, DEN4_down);

% Construct transfer functions for each stage:
H0_down = freqz(G_down, DEN1_down, N);
H1_down = freqz(G_down * NUM1_down, conv_DEN1_DEN2_down, N);
H2_down = freqz(G_down * conv_NUM1_NUM2_down, conv_DEN1_DEN2_DEN3_down, N);
H3_down = freqz(G_down * conv_NUM1_NUM2_NUM3_down, conv_DEN1_DEN2_DEN3_DEN4_down, N);
H4_down = freqz(G_down * conv_NUM1_NUM2_NUM3_NUM4_down, conv_DEN1_DEN2_DEN3_DEN4_down, N);

figure;
hold on;
plot(W/pi, 20*log10(abs(H0_down)));
plot(W/pi, 20*log10(abs(H1_down)));
plot(W/pi, 20*log10(abs(H2_down)));
plot(W/pi, 20*log10(abs(H3_down)));
plot(W/pi, 20*log10(abs(H4_down)));
hold off;
title("Plot of H1(z) H2(z) ... for SOS DOWN");
xlabel("Normalized digital radian frequency");
ylabel("Magnitude [dB]");

% Part e) - Check denominators are the same
check1 = DEN1_up == DEN4_down;
check3 = NUM1_up./NUM4_down;


%% Question 5 - Sensitivity Properties of Parallel Allpass Realizations
% Filter parameters:
ellip_lp_ord = 3;                     % Third-order elliptic lowpass filter
Pb_edge = 0.3*pi;                     % Passband edge
Pb_ripple = 0.92;                     % Passband ripple [dB]
Sb_ripple = 20;                       % Stopband ripple [dB]

% Original transfer function coefficients:
original_num = [0.1336, 0.0568, 0.0563, 0.1336];
original_den = [1, -1.5055, 1.2630, -0.3778];

% Coefficients for sum of two all pass realization:
H1_num = [-0.4954, 1];  
H1_den = [1, -0.4954]; 

H2_num = [0.7626, -1.0101, 1];
H2_den = [1, -1.0101, 0.7626]; 

% Quantize coefficients to compare resulting transfer functions...
%       *** Can interpret 0.5 as a shift for quantization ***

DEN_Q5 = zeros(4, 1);
NUM_Q5 = zeros(4, 1);   
DEN_Q5(1) = 1;
for i = 2:4
    DEN_y = fi(original_den(i), 1, 5, 3);
    DEN_Q5(i) = DEN_y.data;
end

for i = 1:4
    NUM_y = fi(original_num(i), 1, 5, 3);
    NUM_Q5(i) = NUM_y.data;
end

% Quantize the H! and H2 num and den's:
H1_den_quant = fi(H1_den(2), 1, 5, 3);
H1_den_quant = H1_den_quant.data;

H1_num_quant = fi(H1_num(2), 1, 5, 3);
H1_num_quant = H1_num_quant.data;

H2_den_quant = fi(H2_den(2), 1, 5, 3);
H2_den_quant = H2_den_quant.data;

H2_num_quant = fi(H2_num(1), 1, 5, 3);
H2_num_quant = H2_num_quant.data;

H2_den_quant_pos3 = fi(H2_den(3), 1, 5, 3);
H2_den_quant_pos3 = H2_den_quant_pos3.data;

H2_num_quant_pos3 = fi(H2_num(2), 1, 5, 3); 
H2_num_quant_pos3 = H2_num_quant_pos3.data;

% part b) - Check ideal gain 1 at omega = 0 and 0 at omega = pi:
% Effie Bluestone helped me understand this section!
gainAtZeroOriginal = sum(original_num)/sum(original_den);
t1 = [1, -1, 1, -1];
gainAtPiOriginal = sum(original_num.*t1) / sum(original_den.*t1);

gainAtZeroTwoAllP = (1/2) * (sum(H1_num) / sum(H1_den) + sum(H2_num) / sum(H2_den));
t2 = [1, -1];
i = [1, -1, 1];
gainAtPiTwoAllp = (1/2) * (sum(H1_num.*t2) / sum(H1_den.*t2) + sum(H2_num.*i) / sum(H2_den.*i) );

gainAtZeroQuant = sum(NUM_Q5) / sum(DEN_Q5);
gainAtPiQuant = sum(NUM_Q5'.*t1) / sum(DEN_Q5'.*t1);

H1_num_quant1 = [H1_num_quant, 1];
H2_num_quant2 = [H2_num_quant_pos3, H2_num_quant, 1];

H1_den_quant1 = [1, H1_den_quant];
H2_den_quant2 = [1, H2_den_quant, H2_den_quant_pos3];

gainAtZero_quant_two_allp = (1/2) * ...
    (sum(H1_num_quant1) / sum(H1_den_quant1) + sum(H2_num_quant2) / sum(H2_den_quant2));

gainAtPi_quant_two_allp = (1/2) * ...
    (sum(H1_num_quant1.*t2) / sum(H1_den_quant1.*t2) + sum(H2_num_quant2.*i) / sum(H2_den_quant2.*i));

% Plot gain
errorZero = 20*log10(gainAtZeroQuant / gainAtZeroOriginal);

errorPi = 20*log10(gainAtPiQuant / gainAtPiOriginal);

errorZeroTwoAllpass = 20*log10(gainAtZero_quant_two_allp / gainAtZeroTwoAllP);

errorPiTwoAllpass = 20*log10(gainAtPi_quant_two_allp / gainAtPiTwoAllp);

% Part c)
[H, omega] = freqz(original_num, original_den, 1000);
HQZero = freqz(NUM_Q5, DEN_Q5, 1000);
Omega = omega / pi;
HQa = (1/2) * (freqz(H1_num_quant1, H1_den_quant1, 1000) + freqz(H2_num_quant2, H2_den_quant2, 1000));

% Find max:
maxQZero = max(abs(H-HQZero));
maxQA = max(abs(H-HQa));

figure;
hold on
plot(Omega, 20*log10(abs(H)));
plot(Omega, 20*log10(abs(HQZero)));
plot(Omega, 20*log10(abs(HQa)));
xlim([0 1]);
ylim([-40 10]);
xlabel("Normalized frequency");
ylabel("Magnitude [dB]");
title("Superposition of each filter");
legend("H(\omega)", "H_{Q0}(\omega)", "H_{QA}(\omega)");
hold off


% Part d)
% 1.
maxDevPBGainZero = max(abs(20*log10(abs(H(1:0.2*1000))) - ...
    20*log10(abs(HQZero(1:0.2*1000)))));

maxDevPBGainA = max(abs(20*log10(abs(H(1:0.2*1000))) - ...
    20*log10(abs(HQa(1:0.2*1000)))));

% 2.
%% Comment
%
% Similar to the first filter in terms of the equiripple
%
%%
%3. 
maxSBGain = max(20*log10(abs(H(0.4*1000:end))));

maxSBGainZero = max(20*log10(abs(HQZero(0.4*1000:end))));

maxSBGainA = max(20*log10(abs(HQa(0.4*1000:end))));

% Part e)
[group_delay, Wp] = grpdelay(original_num, original_den, 1000);

group_delay_zero = grpdelay(original_num, original_den, 1000);

grou_delay_A = grpdelay(H1_num_quant1, H1_den_quant1, 1000) + ...
    grpdelay(H2_num_quant2, H2_den_quant2, 1000);

omega_freqs = Wp(1:0.2*1000) / pi;

% Create figure
figure;
hold on
plot(omega_freqs, group_delay(1:0.2*1000));
plot(omega_freqs, group_delay_zero(1:0.2*1000));
plot(omega_freqs, grou_delay_A(1:0.2*1000));

xlabel("Normalized \omega");
ylabel("Group Delay");
title("Superimposed graphs of group delay in PB");
legend("H(\omega)", "H_{Q0}(\omega)", "H_{QA}(\omega)");
hold off;
% We can see that the all pass is more sensitive

% Part f) - Compute poles and zeros

% Poles and zeros for original filter
figure;
zplane(original_num, original_den);
title("Original filter poles and zeros");

figure;
zplane(NUM_Q5', DEN_Q5');
title("Quantized filter poles and zeros");

figure;
zplane(conv(H2_num_quant2, H1_den_quant1) + conv(H1_num_quant1, H2_den_quant2), ...
    conv(H1_den_quant1, H1_den_quant1));
title("Reduced filter poles and zeros");
%% Comment
%
% The zero has moved further inside the unit circle towards the imaginary
% axis for the reduced filter. For the quantized filter it seems to have
% stayed at z = -1








