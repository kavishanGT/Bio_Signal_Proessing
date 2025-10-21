% 2.4. Signal Compression with DWT

clearvars;
close all;
clc;

% Load ECG signal
load('ECGsig.mat'); % variable: ecg_signal

figure;
plot(aVR, 'y','LineWidth',0.2);title('Original Signal');
xlabel('Sample Index');
ylabel('Amplitude');
grid on;

% Decomposition
[C_db9,L_db9] = wavedec(aVR,10,'db9');
[C_haar,L_haar] = wavedec(aVR,10,'haar');

% Sort coefficients by magnitude
coeffs_db9 = sort(abs(C_db9),'descend');
coeffs_haar = sort(abs(C_haar),'descend');

E_total_haar = sum(coeffs_haar.^2);
E_total_db9 = sum(coeffs_db9.^2);

E_cum_haar = cumsum(coeffs_haar.^2);
E_cum_db9 = cumsum(coeffs_db9.^2);

figure;

% y1_haar
subplot(1,2,1);
stem(coeffs_haar, 'Marker','none');
title('ecg - Haar');
xlabel('Coefficient Index');
ylabel('Magnitude');

% y2_haar
subplot(1,2,2);
stem(coeffs_db9, 'Marker','none');
title('ecg - db9');
xlabel('Coefficient Index');
ylabel('Magnitude');

% Find number of coefficients for 99% energy
k_db9 = find(E_cum_db9 >= 0.99*E_total_db9, 1);
k_haar = find(E_cum_haar >= 0.99*E_total_haar, 1);

fprintf('db9: %d coefficients for 99%% energy\n', k_db9);
fprintf('Haar: %d coefficients for 99%% energy\n', k_haar);

% Compression using db9
C_db9_comp = C_db9;
C_db9_comp(abs(C_db9) < coeffs_db9(k_db9)) = 0;
y_compressed_db9 = waverec(C_db9_comp,L_db9,'db9');

% Compression using Haar
C_haar_comp = C_haar;
C_haar_comp(abs(C_haar) < coeffs_haar(k_haar)) = 0;
y_compressed_haar = waverec(C_haar_comp,L_haar,'haar');

% Compression ratio
CR_db9 = length(C_db9)/k_db9;
CR_haar = length(C_haar)/k_haar;
fprintf('Compression ratio db9: %.2f\n', CR_db9);
fprintf('Compression ratio Haar: %.2f\n', CR_haar);

% Plot the original and reconstructed signals
figure;
subplot(2,1,1);
plot(aVR, 'b', 'DisplayName', 'Original Signal');
hold on;
plot(y_compressed_db9, 'r', 'DisplayName', 'Reconstructed (db9)');
title('Original vs Reconstructed Signal (db9)');
xlabel('Sample Index');
ylabel('Amplitude');
legend;

subplot(2,1,2);
plot(aVR, 'b', 'DisplayName', 'Original Signal');
hold on;
plot(y_compressed_haar, 'r', 'DisplayName', 'Reconstructed (haar)');
title('Original vs Reconstructed Signal (haar)');
xlabel('Sample Index');
ylabel('Amplitude');
legend;