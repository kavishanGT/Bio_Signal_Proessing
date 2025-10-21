clearvars; close all; clc;

fs = 512;                % Sampling frequency
N = 1024;                 % Segment length (adjust as needed)             % Discrete time axis
n = 0:N-1;         % Discrete time index (0 to 1023)
x_1 = zeros(1, N); % Preallocate signal 
x_1(n >= 0 & n < 511) = 2*sin((1/fs)*20 * pi * n(n >= 0 & n < 511)) + sin((1/fs)*80 * pi * n(n >= 0 & n < 511));
x_1(n >= 512 & n < 1024) = 0.5*sin((1/fs)*40 * pi * n(n >= 512 & n < 1024)) + sin((1/fs)*60 * pi * n(n >= 512 & n < 1024));



x2 = zeros(1, N);     % preallocate

% Apply piecewise definitions
x2(n >= 0   & n < 64)   = 1;
x2(n >= 192 & n < 256)  = 2;
x2(n >= 128 & n < 512)  = -1;
x2(n >= 512 & n < 704)  = 3;
x2(n >= 704 & n < 960)  = 1;


% Add AWGN
y1 = awgn(x_1, 10, 'measured');
y2 = awgn(x2,10,'measured');

%y1_haar
[C_y1_haar, L_y1_haar] = wavedec(y1, 10, 'haar');
coeffs_y1_haar = abs(C_y1_haar);
coeffs_sorted_y1_haar = sort(coeffs_y1_haar,'descend');

%y2_haar
[C_y2_haar, L_y2_haar] = wavedec(y2, 10, 'haar');
coeffs_y2_haar = abs(C_y2_haar);
coeffs_sorted_y2_haar = sort(coeffs_y2_haar,'descend');

%y1_db9
[C_y1_db9, L_y1_db9] = wavedec(y1, 10, 'db9');
coeffs_y1_db9 = abs(C_y1_db9);
coeffs_sorted_y1_db9 = sort(coeffs_y1_db9,'descend');

%y2_db9
[C_y2_db9, L_y2_db9] = wavedec(y2, 10, 'db9');
coeffs_y2_db9 = abs(C_y2_db9);
coeffs_sorted_y2_db9 = sort(coeffs_y2_db9,'descend');

figure;

% y1_haar
subplot(2,2,1);
stem(coeffs_sorted_y1_haar, 'Marker','none');
title('y1 - Haar');
xlabel('Coefficient Index');
ylabel('Magnitude');

% y2_haar
subplot(2,2,2);
stem(coeffs_sorted_y2_haar, 'Marker','none');
title('y2 - Haar');
xlabel('Coefficient Index');
ylabel('Magnitude');

% y1_db9
subplot(2,2,3);
stem(coeffs_sorted_y1_db9, 'Marker','none');
title('y1 - db9');
xlabel('Coefficient Index');
ylabel('Magnitude');

% y2_db9
subplot(2,2,4);
stem(coeffs_sorted_y2_db9, 'Marker','none');
title('y2 - db9');
xlabel('Coefficient Index');
ylabel('Magnitude');

sgtitle('Sorted Wavelet Coefficients (10-level decomposition)');

%y1_haar
theta_y1_haar = 0.2 * max(coeffs_y1_haar); % threshold chosen by observation
fprintf('threshold for y1_haar: %.2f\n', theta_y1_haar);

C_thresh_y1_haar = C_y1_haar;
C_thresh_y1_haar(abs(C_thresh_y1_haar)<theta_y1_haar) = 0;

y_denoised_haar = waverec(C_thresh_y1_haar,L_y1_haar,'haar');

rmse_y1_haar = sqrt(mean((x_1 - y_denoised_haar).^2));
fprintf('rmse between original signal x1 and denoised y1 signal using haar: %.14f\n', rmse_y1_haar);

%y1_db9
theta_y1_db9 = 0.2 * max(coeffs_y1_db9); % threshold chosen by observation
fprintf('threshold for y1_db9: %.2f\n', theta_y1_db9);

C_thresh_y1_db9 = C_y1_db9;
C_thresh_y1_db9(abs(C_thresh_y1_db9)<0.9) = 0;

y_denoised_db9 = waverec(C_thresh_y1_db9,L_y1_db9,'db9');

rmse_y1_db9 = sqrt(mean((x_1 - y_denoised_db9).^2));
fprintf('rmse between original signal x1 and denoised y1 signal using db9: %.14f\n', rmse_y1_db9);

figure;

subplot(2,1,1);

plot(x_1, 'b', 'DisplayName', 'Original Signal');hold on;
plot(y_denoised_haar , 'r', 'DisplayName', 'Denoised Signal'); 
xlabel('Time (s)'), ylabel('Amplitude');
title('Signal x1[n] and denoisy y1[n]');
legend;
grid on;

subplot(2,1,2);
 
plot(x_1, 'b', 'DisplayName', 'Original Signal');hold on;
plot(y_denoised_db9 , 'r', 'DisplayName', 'Denoised Signal'); 
xlabel('Time (s)'), ylabel('Amplitude');
title('Signal x1[n] and denoisy y1[n]');
legend;
grid on;

%y2_haar
theta_y2_haar = 0.2 * max(coeffs_y2_haar); % threshold chosen by observation
fprintf('threshold for y2_haar: %.2f\n', theta_y2_haar);

C_thresh_y2_haar = C_y2_haar;
C_thresh_y2_haar(abs(C_thresh_y2_haar)<3.0) = 0;

y2_denoised_haar = waverec(C_thresh_y2_haar,L_y2_haar,'haar');

rmse_y2_haar = sqrt(mean((x2 - y2_denoised_haar).^2));
fprintf('rmse between original signal x2 and denoised y2 signal using haar: %.14f\n', rmse_y2_haar);

%y1_db9
theta_y2_db9 = 0.2 * max(coeffs_y2_db9); % threshold chosen by observation
fprintf('threshold for y2_db9: %.2f\n', theta_y2_db9);

C_thresh_y2_db9 = C_y2_db9;
C_thresh_y2_db9(abs(C_thresh_y2_db9)<1.2) = 0;

y2_denoised_db9 = waverec(C_thresh_y2_db9,L_y2_db9,'db9');

rmse_y2_db9 = sqrt(mean((x2 - y2_denoised_db9).^2));
fprintf('rmse between original signal x2 and denoised y2 signal using db9: %.14f\n', rmse_y2_db9);

figure;

subplot(2,1,1);

plot(x2, 'b', 'DisplayName', 'Original Signal');hold on;
plot(y2_denoised_haar , 'r', 'DisplayName', 'Denoised Signal'); 
xlabel('Time (s)'), ylabel('Amplitude');
title('Signal x2[n] and denoisy y2[n]');
legend;
grid on;

subplot(2,1,2);
 
plot(x2, 'b', 'DisplayName', 'Original Signal');hold on;
plot(y2_denoised_db9 , 'r', 'DisplayName', 'Denoised Signal'); 
xlabel('Time (s)'), ylabel('Amplitude');
title('Signal x2[n] and denoisy y2[n]');
legend;
grid on;