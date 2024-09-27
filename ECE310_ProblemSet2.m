%% Alexander Faust
%
% ECE 310 - Problem Set 2
% October 5, 2023
%
clc; clear; close all;

%% Question 1
% List a and b values:
a = [-2 5 -2];
b = [6 13 6 0];
% Compute zpk from num and den:
[z,p,k] = tf2zp(a, b);

figure;
zplane(z, p);

text(real(z)+0.1,imag(z),"Zero")
text(real(p)+0.1,imag(p),"Pole")

%% Question 2
% List sampling and filter parameters:

fs = 40e6;          % Sampling rate [Hz]
Ny = fs/2;          % Nyquist [Hz]
fpass1 = 9e6;       % Lower passband edge [Hz]
fpass2 = 12.5e6;    % Upper passband edge [Hz]
fstop1 = 9.5e6;     % Lower stopband edge [Hz]
fstop2 = 12e6;      % Upper stopband edge [Hz]
Ap = 1.5;           % Passband variation [dB]
Ast = 40;           % Stopband attenuation [dB]

% Question 2 - part a)
% Obtain the order of the butterworth filter for each case:
[n_butter_dig, Wn_butter_digital] = buttord([fpass1, fpass2]/(fs/2), [fstop1, fstop2]/(fs/2), Ap, Ast);
[n_butter_analog, Wn_butter_analog] = buttord([fpass1, fpass2], [fstop1, fstop2], Ap, Ast, 's');

% Obtain the order of the Chebyshev I filter for each case:
[n_cheby1_dig, Wn_cheby1_digital] = cheb1ord([fpass1, fpass2]/(fs/2), [fstop1, fstop2]/(fs/2),Ap, Ast); 
[n_cheby1_analog, Wn_cheby1_analog] = cheb1ord([fpass1, fpass2], [fstop1, fstop2], Ap, Ast, 's');

% Obtain the order of the Chebyshev II filter for each case:
[n_cheby2_dig, Wn_cheby2_digital] = cheb2ord([fpass1, fpass2]/(fs/2), [fstop1, fstop2]/(fs/2), Ap, Ast);
[n_cheby2_analog, Wn_cheby2_analog] = cheb2ord([fpass1, fpass2], [fstop1, fstop2], Ap, Ast, 's');

% Obtain the order of the Elliptic filter for each case:
[n_ellip_dig, Wn_ellip_digital] = ellipord([fpass1, fpass2]/(fs/2), [fstop1, fstop2]/(fs/2), Ap, Ast);
[n_ellip_analog, Wn_ellip_analog] = ellipord([fpass1, fpass2], [fstop1, fstop2], Ap, Ast, 's');

% Question 2 - part b) & part c)
% Frequency range to plot over:
w = linspace(0, 20e6, 1000);

%% Create the Butterworth filter given filter order in part a):
[B_BW_dig, A_BW_dig] = butter(n_butter_dig, Wn_butter_digital, 'stop');
[B_BW_analog, A_BW_analog] = butter(n_butter_analog, Wn_butter_analog, 'stop', 's');
[Z_BW_dig, P_BW_dig, ~] = butter(n_butter_dig, Wn_butter_digital, 'stop');
[Z_BW_analog, P_BW_analog, ~] = butter(n_butter_analog, Wn_butter_analog, 'stop', 's');

% Compute frequency response at frequencies w for both cases:
H_BW_dig = freqz(B_BW_dig, A_BW_dig, w, fs);
H_BW_analog = freqs(B_BW_analog, A_BW_analog, w);
H_BW_dig_dB = 20*log10(abs(H_BW_dig));
H_BW_analog_dB = 20*log10(abs(H_BW_analog));

% Plot the poles and zeros plot / magnitude and phase response:
poles_zeros(Z_BW_dig, P_BW_dig, Z_BW_analog, P_BW_analog, "Butterworth");
plotFreq(H_BW_dig, H_BW_analog, w, "Butterworth");

%% Create the Chebyshev I filter given filter order in part a):
[B_cheby1_dig, A_cheby1_dig] = cheby1(n_cheby1_dig, Ap, Wn_cheby1_digital, 'stop');
[B_cheby1_analog, A_cheby1_analog] = cheby1(n_cheby1_analog, Ap, Wn_cheby1_analog, 'stop', 's');
[Z_cheby1_dig, P_cheby1_dig, ~] = cheby1(n_cheby1_dig, Ap, Wn_cheby1_digital, 'stop');
[Z_cheby1_analog, P_cheby1_analog, ~] = cheby1(n_cheby1_analog, Ap, Wn_cheby1_analog, 'stop', 's');

% Compute frequency response at frequencies w for both cases:
H_cheby1_dig = freqz(B_cheby1_dig, A_cheby1_dig, w, fs);
H_cheby1_analog = freqs(B_cheby1_analog, A_cheby1_analog, w);
H_cheby1_dig_dB = 20*log10(abs(H_cheby1_dig));
H_cheby1_analog_dB = 20*log10(abs(H_cheby1_analog));

% Plot the poles and zeros plot / magnitude and phase response:
poles_zeros(Z_cheby1_dig, P_cheby1_dig, Z_cheby1_analog, P_cheby1_analog, "Chebyshev I");
plotFreq(H_cheby1_dig, H_cheby1_analog, w, "Chebyshev I");

%% Create the Chebyshev II filter given filter order in part a):
[B_cheby2_dig, A_cheby2_dig] = cheby2(n_cheby2_dig, Ast, Wn_cheby2_digital, 'stop');
[B_cheby2_analog, A_cheby2_analog] = cheby2(n_cheby2_analog, Ast, Wn_cheby2_analog, 'stop', 's');
[Z_cheby2_dig, P_cheby2_dig, ~] = cheby2(n_cheby2_dig, Ast, Wn_cheby2_digital, 'stop');
[Z_cheby2_analog, P_cheby2_analog, ~] = cheby2(n_cheby2_analog, Ast, Wn_cheby2_analog, 'stop', 's');

% Compute frequency response at frequencies w for both cases:
H_cheby2_dig = freqz(B_cheby2_dig, A_cheby2_dig, w, fs);
H_cheby2_analog = freqs(B_cheby2_analog, A_cheby2_analog, w);
H_cheby2_dig_dB = 20*log10(abs(H_cheby2_dig));
H_cheby2_analog_dB = 20*log10(abs(H_cheby2_analog));

% Plot the poles and zeros plot / magnitude and phase response:
poles_zeros(Z_cheby2_dig, P_cheby2_dig, Z_cheby2_analog, P_cheby2_analog, "Chebyshev II");
plotFreq(H_cheby2_dig, H_cheby2_analog, w, "Chebyshev II");

% List Filter order:
Cases = ["Analog" ; "Digital"];
Butterworth_Order = [2*n_butter_analog ; 2*n_butter_dig];
ChebyshevI_Order = [2*n_cheby1_analog ; 2*n_cheby1_dig];
ChebyshevII_Order = [2*n_cheby2_analog ; 2*n_cheby2_dig];
Elliptic_order = [2*n_ellip_analog ; 2*n_ellip_dig];

Filter_Order_Table = table(Cases, Butterworth_Order, ChebyshevI_Order, ChebyshevII_Order, Elliptic_order);
disp("Displaying filter order for each type:" + newline);
disp(Filter_Order_Table);
%% Create the Elliptic filter given the filter order in part a):
[B_ellip_dig, A_ellip_dig] = ellip(n_ellip_dig, Ap, Ast, Wn_ellip_digital, 'stop');
[B_ellip_analog, A_ellip_analog] = ellip(n_ellip_analog, Ap, Ast, Wn_ellip_analog, 'stop', 's');
[Z_ellip_dig, P_ellip_dig, ~] = ellip(n_ellip_dig, Ap, Ast, Wn_ellip_digital, 'stop');
[Z_ellip_analog, P_ellip_analog, ~] = ellip(n_ellip_analog, Ap, Ast, Wn_ellip_analog, 'stop', 's');

% Compute frequency response at frequencies w for both cases:
H_ellip_dig = freqz(B_ellip_dig, A_ellip_dig, w, fs);
H_ellip_analog = freqs(B_ellip_analog, A_ellip_analog, w);
H_ellip_dig_dB = 20*log10(abs(H_ellip_dig));
H_ellip_analog_dB = 20*log10(abs(H_ellip_analog));

% Plot the poles and zeros plot / magnitude and phase response:
poles_zeros(Z_ellip_dig, P_ellip_dig, Z_ellip_analog, P_ellip_analog, "Elliptic");
plotFreq(H_ellip_dig, H_ellip_analog, w, "Elliptic");

% Question 2 - part d)
% Find index of edges:
indstop1 = find(w <= fstop1);
indstop2 = find(w <= fstop2);
indpass1 = find(w <= fpass1);
indpass2 = find(w <= fpass2);


Stopband_edges = [9.5 ; 12 ; "Passband Edges" ; 9 ; 12.5];           % Stop band frequencies in Hz
BW_dig = [H_BW_dig_dB(indstop1(end)) ; H_BW_dig_dB(indstop2(end)) ; "=====" ; H_BW_dig_dB(indpass1(end)) ; H_BW_dig_dB(indpass2(end))];
BW_analog = [H_BW_analog_dB(indstop1(end)) ; H_BW_analog_dB(indstop2(end)) ; "=====" ; H_BW_analog_dB(indpass1(end)) ; H_BW_analog_dB(indpass2(end))];
Cheby1_dig = [H_cheby1_dig_dB(indstop1(end)) ; H_cheby1_dig_dB(indstop2(end)) ; "=====" ; H_cheby1_dig_dB(indpass1(end)) ; H_cheby1_dig_dB(indpass2(end))];
Cheby1_analog = [H_cheby1_analog_dB(indstop1(end)) ; H_cheby1_analog_dB(indstop2(end)) ; "=====" ; H_cheby1_analog_dB(indpass1(end)) ; H_cheby1_analog_dB(indpass2(end))];
Cheby2_dig = [H_cheby2_dig_dB(indstop1(end)) ; H_cheby2_dig_dB(indstop2(end)) ; "=====" ; H_cheby2_dig_dB(indpass1(end)) ; H_cheby2_dig_dB(indpass2(end))];
Cheby2_analog = [H_cheby2_analog_dB(indstop1(end)) ; H_cheby2_analog_dB(indstop2(end)) ; "=====" ; H_cheby2_analog_dB(indpass1(end)) ; H_cheby2_analog_dB(indpass2(end))];
Elliptic_dig = [H_ellip_dig_dB(indstop1(end)) ; H_ellip_dig_dB(indstop2(end)) ; "=====" ; H_ellip_dig_dB(indpass1(end)) ; H_ellip_dig_dB(indpass2(end))];
Elliptic_analog = [H_ellip_analog_dB(indstop1(end)) ; H_ellip_analog_dB(indstop2(end)) ; "=====" ; H_ellip_analog_dB(indpass1(end)) ; H_ellip_analog_dB(indpass2(end))];

Gain_table = table(Stopband_edges, BW_dig, BW_analog, Cheby1_dig, Cheby1_analog, Cheby2_dig, Cheby2_analog, Elliptic_dig, Elliptic_analog);
disp("Displaying table with gain values in [dB] at the stop band or passband edges:" + newline);
disp(Gain_table);
%% Problem 3
% Create the Kaiser Window:
% List the delta pass and stop parameters on a LINEAR scale
del_pass = 10 ^ (-Ap / 20);
del_stop = 10 ^ (-Ast / 20);

% Create stopband and passband edge frequeny vector [Hz]:
F = [fpass1, fstop1, fstop2, fpass2];
% Specify the passband and stop band via the A vector:
A = [1 0 1];
% Specify maximum ripples (ON A LINEAR SCALE)
DEV = [del_pass, del_stop, del_pass];
% Determine kaiser order based on the specs:
[kais_ord, Wn, bta, filtype] = kaiserord(F, A, DEV, fs);
kais_window = kaiser(kais_ord + 1, bta);
% Determine Equiripple FIR design based on Parks-McClellan:
[pm_ord, Fp, As, W] = firpmord(F, A, DEV, fs);

% Determine the filter coefficients of the Kaiser filter:
B_kaiser = fir1(kais_ord, Wn, filtype, kais_window, 'noscale');

% Determine pm FIR filter:
B_pm = firpm(pm_ord, Fp, As, W);


% Question 3 - part b)
% Indicate length of each filter:
kaiser_length = length(B_kaiser);
pm_length = length(B_pm);

% Create stem plot of each filter:
figure;
subplot(2, 1, 1);
stem(0 : (kaiser_length - 1), B_kaiser);
title("Filter Coefficients for Kaiser");
xlabel("Index");
ylabel("Value");

subplot(2, 1, 2);
stem(0 : (pm_length - 1), B_pm);
title("Filter Coefficients for Equiripple");
xlabel("Index");
ylabel("Value");


% Create pole and zero plot for each filter:
figure;
subplot(2, 1, 1);
zplane(B_kaiser, 1);
title("Pole Zero Plot for Kaiser");

subplot(2, 1, 2);
zplane(B_pm);
title("Pole Zero Plot for Equiripple");

% Plot the magnitude response for each filter:
[H_kaiser, Wn_kaiser] = freqz(B_kaiser, 1, w, fs);
[H_pm, Wn_pm] = freqz(B_pm, 1, w, fs);

figure;
subplot(2, 1, 1);
plot(Wn_kaiser / 1e6, 20*log10(abs(H_kaiser)));
title("Frequency Response of Kaiser");
ylim([-50 2]);
xlabel("Frequency [MHz]");
ylabel("Magnitude [dB]");

% Plot phase degree
phase_kaiser_radians_dig = angle(H_kaiser);
phase_kaiser_deg_dig = unwrap(phase_kaiser_radians_dig) * 180/pi;
subplot(2, 1, 2);
plot(w/1e6, phase_kaiser_deg_dig);
title("Phase Response of Digital Kaiser FIR Filter");
xlabel("Frequency [MHz]");
ylabel("Phase [degrees]");

figure;
subplot(2, 1, 1);
plot(Wn_pm / 1e6, 20*log10(abs(H_pm)));
title("Frequency Response of Equiripple");
ylim([-50 2]);
xlabel("Frequency [MHz]");
ylabel("Magnitude [dB]");

phase_PM_radians_dig = angle(H_pm);
phase_PM_deg_dig = unwrap(phase_PM_radians_dig) * 180/pi;
subplot(2, 1, 2);
plot(w/1e6, phase_PM_deg_dig);
title("Phase Response of digital Equiripple FIR Filter");
xlabel("Frequency [MHz]");
ylabel("Phase [degrees]");

% Question 3 - part c)
W_ratio = W(1) / W(2);
del_ratio = del_stop / del_pass;

disp("Checking the ratio of weights to the ratio of tolerance to look for equivalence..." + newline);
disp("W_pass / W_stop = " + W_ratio + ", and del_stop / del_pass = " + del_ratio + newline + ...
     "Clearly, the values are equal." + newline);

% Question 3 - part d)
% Examine the passband range:
passband_indices = find(w < 9e6);

% Determine max and mingain in passband:
maxgain_PB_pm = max(20*log10(abs(H_pm(1 : passband_indices(end)))));
mingain_PB_pm = min(20*log10(abs(H_pm(1 : passband_indices(end)))));
disp("Max gain in passband for Parks McClellan: " + maxgain_PB_pm + ", Min gain: " + mingain_PB_pm);

maxgain_PB_kaiser = max(20*log10(abs(H_kaiser(1 : passband_indices(end)))));
mingain_PB_kaiser = min(20*log10(abs(H_kaiser(1 : passband_indices(end)))));
disp("Max gain in passband for Kaiser: " + maxgain_PB_kaiser + ", Min gain: " + mingain_PB_kaiser + newline);

% Determine if the max gain in the stop band is always less than or equal
% to 30dB:
ind1 = find(w < 12e6);
ind2 = find(w < 9.5e6);
stopband_indices = ind1(ind2(end) + 1 : end);

% For the PM filter:
for i = 1:length(stopband_indices)
    % Create a way to tell at the end if the process completed
    % successfully:
    err = 0;

    if 20*log10(abs(H_pm(stopband_indices(i)))) > -30
        disp("Equiripple filter design in the stop band is above 30dB at frequency: " + w(stopband_indices(i)));
        err = err + 1;
    end
    
    if i == length(stopband_indices)
        if err == 0
            disp("Equiripple filter design in the stop band is always below 30dB!");
        else 
            disp("Equiripple filter design failed to meet specs at more than one frequency...");
        end
    end
end

% For the Kaiser Filter:
for i = 1:length(stopband_indices)
    % Create a way to tell at the end if the process completed
    % successfully:
    err = 0;

    if 20*log10(abs(H_kaiser(stopband_indices(i)))) > -30
        disp("Kaiser filter design in the stop band is below 30dB at frequency: " + w(stopband_indices(i)));
        err = err + 1;
    end
    
    if i == length(stopband_indices)
        if err == 0
            disp("Kaiser filter design in the stop band is always below 30dB!");
        else
            disp("Kaiser filter design failed to meet specs at more than one frequency");
        end
    end
end

% Question 3 - part e)
disp(newline);
disp("According to the simulations in the passband and the difference between the max and min" + newline + ...
    "the filter design for the Kaiser filter, it meets the specs. However, the equiripple has a little more than a 1.5dB" + newline + ...
    "difference in the max and min gain. This means that it did not meet the specs exactly in the passband. I would fix this" + newline + ...
    "by perhaps lowering the ripple factor when passing variables through the filter.");

%% Functions Created:
% Create function to simply make freqs/freqz plots of each filter created
% the same way for different A and B coefficients:
function plotFreq(H1, H2, w, filter_name)
    % Create digital magnitude response with the given w frequency range 
    % and the phase response
    figure;
    subplot(2, 1, 1);
    plot(w/1e6, 20*log10(abs(H1)));
    % freqz(B_digital, A_digital);
    title("Magnitude response of digital " + filter_name + " IIR Filter");
    xlabel("Frequency [MHz]");
    ylabel("Magnitude [dB]");
    ylim([-50 2]);
    
    % Create plot of phase response:
    phase_radians_dig = angle(H1);
    phase_deg_dig = unwrap(phase_radians_dig) * 180/pi;
    subplot(2, 1, 2);
    plot(w/1e6, phase_deg_dig);
    title("Phase Response of digital " + filter_name + " IIR Filter");
    xlabel("Frequency [MHz]");
    ylabel("Phase [degrees]");
    
    % Create analog magnitude response:
    figure;
    subplot(2, 1, 1);
    plot(w/1e6, 20*log10(abs(H2)));
    title("Magnitude Response of analog " + filter_name + " IIR Filter");
    xlabel("Frequency [MHz]");
    ylabel("Magnitude [dB]");
    ylim([-50 2]);

    % Create subplot for analog phas response:
    phase_radians_analog = angle(H2);
    phase_deg_analog = unwrap(phase_radians_analog) * 180/pi;
    subplot(2, 1, 2);
    plot(w/1e6, phase_deg_analog);
    title("Phase Response of analog " + filter_name + " IIR Filter");
    xlabel("Frequency [MHz]");
    ylabel("Phase [degrees]");
end

% Function to create pole and zero plots for each filter:
function poles_zeros(Z_dig, P_dig, Z_analog, P_analog, filter_name)
    figure;
    subplot(2, 1, 1);
    zplane(Z_dig, P_dig);
    title("Poles and Zero Plot for Digital " + filter_name + " Filter");

    subplot(2, 1, 2);
    zplane(Z_analog, P_analog);
    title("Poles and Zero Plot for Analog " + filter_name + " Filter");
end

