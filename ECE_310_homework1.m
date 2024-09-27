%% Alexander Faust
%% ECE-310 Homework 1 
clear; clc; close all;

%% Problem 2
% List variables:
f_s = 20e6;                 % Sample frequency in Hz (20 MHz)
n = 500;                    % Number of samples
t = 0:1/f_s:(n-1)/f_s;      % Time vector
x = sin(2*pi*(6e6)*t);      % Real valued sinusoid at 6 MHz
N = 512;                    % Number of points for DFT

% Question 2a:
bin_spacing = f_s/N;    % Bin spacing in [Hz]
disp("Question 2a:" + newline);
disp("Bin spacing in Hertz: " + bin_spacing + newline);

% Question 2b & 2c:
k = 0:N-1;              % Index for the 512-point DFT
freqs = k*f_s/N;        % Bin frequency vector 

% Compute the DFT of this sequence using a rectangular window:
Rect_win = rectwin(n);
Rect_win_scaled = sqrt(1/n) * Rect_win; 
x_Rwin = x .* Rect_win_scaled';

% Take a 512-point DFT                                
X1 = fft(x_Rwin, N);      % DFT of signal x using N points 

% Find the values of k where peaks occur i.e. bin frequencies
[PKS, LOCS] = findpeaks(abs(X1),"MinPeakHeight",0.2);
bin_freqs(:) = LOCS*(f_s/N)/1e6;    % In [MHz]
disp("Question 2b:" + newline);
disp("The indices k = " + LOCS(1) + " and k = " + LOCS(2) + " are the indices where a significant");
disp("peak is observed. Furthermore, this would correlate to the bin frequencies"); 
disp("w = " + bin_freqs(1) + "MHz, and w = " + bin_freqs(2) + "MHz.");

% Create a plot of the computed DFT of the signal using the rectangular window
figure;
subplot(2,1,1);
plot(freqs, abs(X1));      
title("512-point DFT of 6MHz signal using a Rectangular window");  
xlabel("Frequency [Hz]");
ylabel("Magnitude");

% Compute DFT using Chebyshev window with peak sidelobe level 30dB
Cheb_win = chebwin(n, 30);
energy_Cwin = sum(Cheb_win .^ 2);
Cheb_win_scaled = Cheb_win * sqrt(1/energy_Cwin);

% Obtain windowed version of input signal
x_Cwin = x .* Cheb_win_scaled';
X2 = fft(x_Cwin, N);

subplot(2,1,2);
plot(freqs, abs(X2));      
title("512-point DFT of 6MHz signal using a Chebyshev window");  
xlabel("Frequency [Hz]");
ylabel("Magnitude [dB]");

% Create Superimposed graphs around the peak index:
location = LOCS(1);
%superimposed_range = location + 1 + N/2;
range_shift = (location - 10):(location + 10);
freq_range_shift = (range_shift) * bin_spacing;

figure;
plot(freq_range_shift, 20*log10(abs(X1(range_shift))), '-r', ...
     freq_range_shift, 20*log10(abs(X2(range_shift))), '-b');
title("Superimposed graphs of DFT of signal using different windows");
ylabel("|X(f)| [dB]");
xlabel("Frequency [Hz]");
legend("Rectangular Window", "Chebyshev Window");

% Question 2d:
SNR = 20;               % SNR of 20dB to corrupt the input signal
signalpower = mean(x.^2);
noisepower = signalpower / (10^(SNR / 10));
x_corrupted = x + sqrt(noisepower) * randn(1, n);

% Window the corrupted signal using a rectangular window and take 512-point
% DFT:
x_corrupted_Rwin = x_corrupted .* Rect_win_scaled';
X_corrupted_Rwin = fft(x_corrupted_Rwin, N);
X_corrupted_Rwin_shifted = fftshift(X_corrupted_Rwin);

% Create plot of the magnitude of the DFT of X in the range -Fs/2 to Fs/2:
shift_interval = (-N/2:(N-1)/2) * bin_spacing;
figure;
subplot(2,1,1);
plot(shift_interval, 20*log10(abs(X_corrupted_Rwin_shifted)));
xlabel("Frequency [Hz]");
ylabel("|X(f)| in dB");
title("DFT of 6MHz signal corrupted with 20dB SNR using a Rectangular window");

% Window the corrupted signal using the Chebyshev window and take 512-point
% DFT:
subplot(2,1,2);
x_corrupted_Cwin = x_corrupted .* Cheb_win_scaled';
X_corrupted_Cwin = fft(x_corrupted_Cwin, N);
X_corrupted_Cwin_shifted = fftshift(X_corrupted_Cwin);

% Create plot of the magnitude of the DFT of X in the range -Fs/2 to Fs/2:

plot(shift_interval, 20*log10(abs(X_corrupted_Cwin_shifted)));
xlabel("Frequency [Hz]");
ylabel("|X(f)| in dB");
title("DFT of 6MHz signal corrupted with 20dB SNR using a Chebyshev window");

% Create superimposed graphs of the two cases:
k0 = LOCS(1);
k0_range = k0 + 1 + N/2;
k0_range_shift = (k0_range - 10):(k0_range + 10);
f0_range_shift = (k0_range_shift - N/2 - 1) * bin_spacing;

figure;
plot(f0_range_shift, 20*log10(abs(X_corrupted_Rwin_shifted(k0_range_shift))), '-r', ...
     f0_range_shift, 20*log10(abs(X_corrupted_Cwin_shifted(k0_range_shift))), '-b');
title("Superimposed graph of DFT of noisy signal using different windows");
ylabel("|X(f)| [dB]");
xlabel("Frequency [Hz]");
legend("Rectangular Window", "Chebyshev Window");

%% Problem 3
Fs = 44.1e3;                        % Sample frequency [Hz] (Common for audio applications)
N0 = Fs/2;                          % Given 2Hz bin spacing, this implies f_s/2Hz data samples
N = 2.^15;                          % Minimum power of 2 that makes DFT efficient (32,768)
tempo = 100;                        % 100 beats per min tempo
BPS = 100/60;                       % Beats per second
duration = 1/BPS;                   % Duration of each beat during one second
frequencies = [392 440 587.33];     % List frequencies to randomly choose from
T = 1/Fs;                           % Sampling period

% Calculate samples per note considering the duration of a beat:
samples_per_note = round(Fs * duration);

% Calculate the number of total samples needed for 100 blocks 50% overlap criterion
total_samples = (101*N0)/2;

% Calculate total number of notes to generate based on the total samples criterion:
num_notes = round(total_samples / samples_per_note);

% Synthesize signal with requirement of 100 blocks:
outputSignal = zeros(1,total_samples);
for i = 1:num_notes
    % Randomly choose the 2 tones of the 3 total tones to combine
    rand_freq = randsample(frequencies, 2);
    % Create time vector for each note
    time = (0:T:(duration - T));  
    % Generate two sinusoids with the randomly chosen frequencies
    x1 = sin(2*pi*rand_freq(1)*time);
    x2 = sin(2*pi*rand_freq(2)*time);
    
    % Create quarter note from sum of two sinusoids
    note = x1 + x2;
    
    % List indices for concatenation
    start_idx = (i - 1) * samples_per_note + 1;
    end_idx = i * samples_per_note;

    % Assign the output signal the quarter note
    outputSignal(start_idx:end_idx) = note;   
end

% Add noise with SNR of 40dB to output signal:
power_signal = sum(abs(outputSignal).^2) / total_samples;
SNR_linear_scale = 10^(40 / 10);
awgn = sqrt(power_signal/SNR_linear_scale) * randn(size(outputSignal));
noisy_signal = outputSignal + awgn;

% Create hamming window:
HamWin = hamming(N0);

% Estimate PSD using pwelch:
% Window output signal with 2Hz window bin spacing
[pxx, W] = pwelch(noisy_signal, HamWin, N0/2, N, Fs);         

% Plot PSD from 0 to Fs/2
figure;
plot(W, 10*log10(pxx));
xlabel("Frequency [Hz]");
ylabel("Power spectrum (dB/Hz)");
title("Periodogram from 0 to f_s/2");
xlim([0 Fs/2]);

% Plot PSD from 300 Hz to 600 Hz
figure;
plot(W, 10*log10(pxx));
xlabel("Frequency [Hz]");
ylabel("Power spectrum (dB/Hz)");
title("Zoomed in periodogram from 300 to 600 Hz");
xlim([300 600]);

% Problem 3d:
% Create spectrogram plot to see frequency spectrum in "real time"
figure;
spectrogram(noisy_signal, HamWin, N0/2, N, Fs);
title("Spectrogram of output signal (frequency as time evolves)");
% Limit frequency band from 300 Hz to 700 Hz 
xlim([0.3 0.7]);