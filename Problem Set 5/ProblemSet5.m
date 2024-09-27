%% Alexander Faust
%
% ECE-310 Problem Set 5 Multidimensional
%
% November 27, 2023
clc; clear; close all;

% Load images:
load("LilyImg.mat");

% All so powerful Rodan! - sorry Rodan I will get back to work
load("Rodanimg.mat");

Lily_img = Lilyx;
Rodan_img = Rodanx;

%% Question 2 - Laplacian
% Part a) - Specify 2D Laplacian estimation operator and compute freqz2:
h_lap = (1/6) * [1, 4, 1, ; 4, -20, 4 ; 1, 4, 1];
[H_lap,fx,fy]=freqz2(h_lap);

fx = fx*pi;
fy = fy*pi;

% Specify the horizontal and vertical sobel operators:
const = 1/2;
hx = const * [-1, 0, 1 ; 
              -2, 0, 2 ; 
              -1, 0, 1];

hy = const * [-1, -2, -1 ;
               0,  0,  0 ;
               1,  2,  1];

% Specify the laplacian kernel in terms of the Sobel operators:
h_lap_sobel = conv2(hx, hx, 'full') + conv2(hy, hy, 'full');

% Part b) - Compute frequency response of the laplacian and laplacian
%           estimation when using the sobel operators:
[H_lap_sobel, ~, ~] = freqz2(h_lap_sobel);

% Part c) - Obtain surface plots and contour plots for each filter
figure;
surf(abs(H_lap));
title("Surface plot of H_{Lap}");
ylabel("Frequency k_{y}");
xlabel("Frequency k_{x}");
zlabel("Amplitude");

figure;
surf(abs(H_lap_sobel));
title("Surface plot of H_{LapSob}");
ylabel("Frequency k_{y}");
xlabel("Frequency k_{x}");
zlabel("Amplitude");

figure;
contour(abs(H_lap));
title("Contour plot of H_{Lap}");
ylabel("Frequency k_{y}");
xlabel("Frequency k_{x}");

figure;
contour(abs(H_lap_sobel));
title("Surface plot of H_{LapSob}");
ylabel("Frequency k_{y}");
xlabel("Frequency k_{x}");

% Part d) - Use the Rodan and Lily images to create gray-scale images
%           & the result of applying h_Lap and H_LapSob to each

% Show the un-altered gray-scale images:
figure;
image(Lily_img);
colormap("gray");

figure;
image(Rodan_img);
colormap("gray");

% Filter Rodan and Lily using Laplacian filter first:
% Convert to double before filtering:
Lily_gray = im2gray(Lily_img);
Lily_gray = double(Lily_gray);

Rodan_gray = im2gray(Rodan_img);
Rodan_gray = double(Rodan_gray);

% Apply Laplacian filter:
LapFilteredRodan = filter2(h_lap, Rodan_gray);
LapFilteredLily = filter2(h_lap, Lily_gray);

% Apply Laplacian filter with the Sobel filters:
LapSobFilteredRodan = filter2(h_lap_sobel, Rodan_gray);
LapSobFilteredLily = filter2(h_lap_sobel, Lily_gray);

% Create image plots of the filtered images
figure;
subplot(1, 2, 1);
image(LapFilteredRodan);
title("Laplacian Kernel");
colormap("gray");

subplot(1, 2, 2);
image(LapSobFilteredRodan);
title("Laplacian with Sobel Approximation");
colormap("gray");

figure;
subplot(1, 2, 1);
image(LapFilteredLily);
title("Laplacian Kernel");
colormap("gray");

subplot(1, 2, 2);
image(LapSobFilteredLily);
title("Laplacian with Sobel Approximation");
colormap("gray");

%% Question 3 - Filtering Continued
% Part a/b) - Show upsampled versions of Lily and Rodan images using the
%           upsampling function
[upsampledLily, originalLily] = upsampleImage(Lily_gray);

[upsampledRodan, originalRodan] = upsampleImage(Rodan_gray);

% Part c) - Compute 2D DFT of images before and after upsampling

% DFT of images before upsampling:
LilyDFT = fftshift(fft2(Lily_gray));
RodanDFT = fftshift(fft2(Rodan_gray));

% DFT of images after upsampling:
upsampledLilyDFT = fftshift(fft2(upsampledLily));
upsampledRodanDFT = fftshift(fft2(upsampledRodan));

% Create the Lily images:
figure;
subplot(1, 2, 1);
imshow(20*log10(1 + abs(LilyDFT)), []);
colormap("gray");
title("Lily DFT before upsampling");
axis on;

subplot(1, 2, 2);
imshow(20*log10(1 + abs(upsampledLilyDFT)), []);
colormap("gray");
title("Lily DFT after upsampling");
axis on;

% Create the Rodan images:
figure;
subplot(1, 2, 1);
imshow(20*log10(1 + abs(RodanDFT)), []);
colormap("gray");
title("Rodan DFT before upsampling");
axis on;

subplot(1, 2, 2);
imshow(20*log10(1 + abs(upsampledRodanDFT)), []);
colormap("gray");
title("Rodan DFT after upsampling");
axis on;

%% Comment on DFT images
%
% What is being observed in the upsampled spectrum is aliasing distortion.
% Much like the classical example of aliasing from filtering, we can
% observe an aliased version of the image due to taking a DFT.




%% Functions Created
% Function to upsample an image 2D matrix
function [upsampled_image, u] = upsampleImage(u)
    % Declare the upsample matrix:
    M = [2, 0 ;
         0, 2];

    upsampled_image = kron(u, M);
    
    % Create figure comparing the upsampled matrix and the original:
    figure;
    subplot(1, 2, 1);
    imshow(u, []);
    colormap("gray");
    title("Image Before Upsampling");

    subplot(1, 2, 2);
    imshow(upsampled_image, []);
    colormap("gray");
    title("Image After Upsampling");

end









