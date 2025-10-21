clearvars; close all; clc;

%% Parameters
fs = 250;                % Sampling frequency
N = 3000;                 % Segment length (adjust as needed)             % Discrete time axis
n = 1:3*N;

x = zeros(size(n)); 
x(n >= 1 & n < (3*N)/2) = sin((1/fs)*0.5 * pi * n(n >= 1 & n < (3*N)/2));
x(n >= (3*N)/2 & n < 3*N) = sin((1/fs)*1.5 * pi * n(n >= (3*N)/2 & n < 3*N));

figure;
plot(n/fs, x, 'w','LineWidth',0.2);
xlabel('Time (s)'), ylabel('x[n]');
title('Original Signal');
xlim([0 35]);
grid on;


 
s_list = 0.01:0.01:2;      % scale resolution


t = (-N:N)/fs;
%% Continuous wavelet transform (convolution)
coeffs = zeros(numel(s_list), length(x));

for k = 1:numel(s_list)
    s = s_list(k);
    
    % scaled daughter wavelet
    psi_values = (2 / (sqrt(3*s) * pi^(-1/4)))*(1 - (t/s).^2).* exp(-1*((t/s).^2 )/ 2);
    wavelet = psi_values;
    %psi_s = psi_s / sqrt(trapz(tau_s, psi_s.^2));  % normalize
    % convolution
    coeffs(k,:) = conv(x, wavelet, 'same');
end

figure;
h = pcolor(n, s_list, coeffs);
set(h,'EdgeColor','none');
colormap jet; colorbar;
xlabel('Time (s)'), ylabel('Scale (s)');
title('Mexican Hat Wavelet Coefficients (Scalogram)');

