%% Alexander Faust
%
% ECE-310 Problem Set 3 - Multirate
%
% November 18, 2023
clc; clear; close all;

%% Question 5 - Wavelets and Filter Banks
% Produce filter coeff for Daubechies wavelets order N:
N= 1;                       % substitute a value between 1 and 45
wname= ['db',int2str(N)];
[h0,h1,f0,f1]= wfilters(wname);

% Frequency vector:
freq = linspace(0, pi, 10^4);

% Create formula for the impulse response for H0 and H1:
[H0, w] = freqz(h0, 1, freq);
[H1, ~] = freqz(h1, 1, freq);

figure;
hold on
plot(w, abs(H0));
plot(w, abs(H1));
title("Magnitude of Haar Wavelet Filters");
xlabel("Frequency (\pi radians/sample)");
ylabel("Magnitude [linear]");
legend("|H_{0}(\omega)|", "|H_{1}(\omega)|");
xlim([0 pi]);
hold off

% Construct E(z):
E = [h0(1:2:end), h1(1:2:end) ; h0(2:2:end) , h1(2:2:end)];
E_conj = conj(E');
disp("Verifying paraunitary property: ");
disp(E*E_conj);
disp("Furthermore, determinant of E: ");
disp(det(E));


% part b)
% Compute power and some other values:
% * Use filter parameters for power NOT freqz *
Power_omega = (abs(h0)).^2 + (abs(h1)).^2;
difference = max(Power_omega) - min(Power_omega);
mean_power = mean(Power_omega);

% Check variation small:
disp("Variation in power: " + difference);
disp("Mean Power: " + mean_power);

% part c) - generate db5 filter
[h0_db5, h1_db5, f0_db5, f1_db5] = wfilters('db5');
disp("db5 H0 coefficients: ");
disp(h0_db5);
disp("db5 H1 coefficients: ");
disp(h1_db5);

% 1. Break down db5 filters h1 and h0 into type I polyphase components:
h0_db5_even = h0_db5(1:2:end);              % Obtain even indexed coeff
h0_db5_odd = h0_db5(2:2:end);               % Obtain odd indexed coeff

h1_db5_even = h1_db5(1:2:end);              % Obtain even indexed coeff
h1_db5_odd = h1_db5(2:2:end);               % Obtain odd indexed coeff

% 2. Create 3-D array with each entry containing the polyphase coefficients
E_db5 = cat(3, h0_db5_even, h1_db5_even , h0_db5_odd, h1_db5_odd);
E_db5(:, :, 1);                 % displays h0_db5_even for example


% 3. Compute ~E*E
E_tild = fliplr(E_db5);
% Preallocate ~E*E matrix as a 3D array
filler_matrix = zeros(1, 9);
E_tild_E = cat(3, filler_matrix, filler_matrix, filler_matrix, filler_matrix);
E_tild_E(:, :, 1) = conv(E_tild(:,:,1), E_db5(:,:,1)) + conv(E_tild(:,:,2), E_db5(:,:,2));
E_tild_E(:, :, 2) = conv(E_tild(:,:,3), E_db5(:,:,3)) + conv(E_tild(:,:,4), E_db5(:,:,4));
E_tild_E(:, :, 3) = conv(E_tild(:,:,1), E_db5(:,:,3)) + conv(E_tild(:,:,2), E_db5(:,:,4));
E_tild_E(:, :, 4) = conv(E_tild(:,:,3), E_db5(:,:,1)) + conv(E_tild(:,:,4), E_db5(:,:,2));

% 4. Find end-to-end delay for analysis/synthesis filter banks
delay = (length(h0_db5)) / 2;
disp("End-to-end delay for filter bank: " + delay + newline);

% 5. Superimpose plots of |H0| and |H1| for db5:
[H0_db5, w_db5] = freqz(h0_db5, 1, 1024, 'whole');
[H1_db5, ~] = freqz(h1_db5, 1, 1024, 'whole');
W_db5 = w_db5/pi;

figure;
hold on
plot(W_db5, abs(H0_db5)');
plot(W_db5, abs(H1_db5)');
xlabel("Normalized frequency [rad]");
ylabel("Magnitude [linear]");
title("Magnitude plots for H_{1}(\omega) and H_{0}(\omega) for db5");
hold off;
legend("H_{0}(\omega)", "H_{1}(\omega)");

% 6.
%
% Paraunitary property for db5 wavelet filter confirmed in 3. above
%

% 7. Create the new one stage filter bank
frequencies = linspace(0, pi, 1000);
H1Z = freqz(h1_db5, 1, frequencies);
H1Z_decimate2 = freqz(h1_db5, 1, 2*frequencies);
H0Z = freqz(h0_db5, 1, frequencies);
H0Z_decimate2 = freqz(h0_db5, 1, 2*frequencies);
H1Z_decimate4 = freqz(h1_db5, 1, 4*frequencies);
H0Z_decimate4 = freqz(h0_db5, 1, 4*frequencies);

% Plot superimposed magnitude response for each filter bank:
figure;
hold on
plot(frequencies, abs(H1Z), 'DisplayName', "H_{1}(z)");
plot(frequencies, abs(H0Z.*H1Z_decimate2), "DisplayName", "H_{0}(z)H_{1}(z^{2})");
plot(frequencies, abs(H0Z.*H0Z_decimate2.*H1Z_decimate4), "DisplayName", "H_{0}(z)H_{0}(z^{2})H_{1}(z^{4})");
plot(frequencies, abs(H0Z.*H0Z_decimate2.*H0Z_decimate4), "DisplayName", "H_{0}(z)H_{0}(z^{2})H_{0}(z^{4})");
legend;
title("Superimposed magnitude response of each filter bank");
xlabel("\omega / \pi");
ylabel("Magnitude [linear]");
xlim([0 pi]);
hold off

% 8. Compute power associated with each filter bank channel
Power_k1 = (1/2) * abs(H1Z).^2;
Power_k2 = (1/4) * abs(H0Z.*H1Z_decimate2).^2;
Power_k3_channel1 = (1/8) * abs(H0Z.*H0Z_decimate2.*H1Z_decimate4).^2;
Power_k3_channel2 = (1/8) * abs(H0Z.*H0Z_decimate2.*H0Z_decimate4).^2;

% Total power is sum of channel power
total_power_G = Power_k1 + Power_k2 + Power_k3_channel1 + Power_k3_channel2;

power_difference_G = max(total_power_G) - min(total_power_G);
disp("Power difference between max and min is very small: " + power_difference_G + newline);

% Compute nominal power:
nominal_power = mean(total_power_G);

% Finally, compute sum of 1/Mk terms:
M = [2, 4, 8, 8];
inv_M_sum = sum(1 ./ (M));
disp("Sum of inverse of M terms: " + inv_M_sum);








