clearvars; close all; clc;

fs = 512;                % Sampling frequency
N = 1024;                 % Segment length (adjust as needed)             % Discrete time axis
n = 0:N-1;         % Discrete time index (0 to 1023)
x_1 = zeros(1, N); % Preallocate signal 
x_1(n >= 0 & n < 511) = 2*sin((1/fs)*20 * pi * n(n >= 0 & n < 511)) + sin((1/fs)*80 * pi * n(n >= 0 & n < 511));
x_1(n >= 512 & n < 1024) = 0.5*sin((1/fs)*40 * pi * n(n >= 512 & n < 1024)) + sin((1/fs)*60 * pi * n(n >= 512 & n < 1024));

figure;
plot(n/fs, x_1, 'y','LineWidth',0.2);
xlabel('Time (s)'), ylabel('x[n]');
title('Original Signal');
xlim([0 2]);
ylim([-3.5 3.5]);
grid on;


x2 = zeros(1, N);     % preallocate

% Apply piecewise definitions
x2(n >= 0   & n < 64)   = 1;
x2(n >= 192 & n < 256)  = 2;
x2(n >= 128 & n < 512)  = -1;
x2(n >= 512 & n < 704)  = 3;
x2(n >= 704 & n < 960)  = 1;

% Time axis
t = n/fs;

% Plot
figure;
plot(t, x2, 'LineWidth',0.2);
xlabel('Time (s)'), ylabel('x_2[n]');
title('Piecewise Signal x_2[n]');
grid on;
xlim([0 N/fs]);


% Add AWGN
y1 = awgn(x_1, 10, 'measured');
y2 = awgn(x2,10,'measured');

figure;

subplot(2,1,1);

plot(t, y1, 'r--');hold on;
plot(t, x_1, 'w', 'LineWidth', 0.2); 
xlabel('Time (s)'), ylabel('Amplitude');
title('Signal x1[n] and Noisy y1[n]');
legend('y1[n] (AWGN 10 dB)','x1[n] (clean)');
grid on;

subplot(2,1,2);
 
plot(t, y2, 'r','LineWidth', 1.2);hold on;
plot(t, x2, 'w', 'LineWidth', 1.2);
xlabel('Time (s)'), ylabel('Amplitude');
title('Signal x2[n] and Noisy y2[n]');
legend('y2[n] (AWGN 10 dB)','x2[n] (clean)');
grid on;

% ii. Observe the morphology of the wavelet and scaling functions of Haar Daubechies tap 9

[phi_haar,psi_haar,x_haar] = wavefun('haar',10);   % level = 10 for resolution
figure;
subplot(2,1,1);
plot(x_haar, phi_haar, 'y','LineWidth',1.2);
title('Haar Scaling Function \phi');
xlabel('t'); ylabel('\phi(t)'); grid on;

subplot(2,1,2);
plot(x_haar, psi_haar, 'y','LineWidth',1.2);
title('Haar Wavelet Function \psi');
xlabel('t'); ylabel('\psi(t)'); grid on;

% Daubechies 9 wavelet
[phi_db9,psi_db9,x_db9] = wavefun('db9',10);
figure;
subplot(2,1,1);
plot(x_db9, phi_db9, 'y','LineWidth',1.2);
title('Daubechies 9 Scaling Function \phi');
xlabel('t'); ylabel('\phi(t)'); grid on;

subplot(2,1,2);
plot(x_db9, psi_db9, 'y','LineWidth',1.2);
title('Daubechies 9 Wavelet Function \psi');
xlabel('t'); ylabel('\psi(t)'); grid on;


%% --- 10-level Wavelet Decomposition using Haar ---
% y1 with haar
[C_y1_haar,L_y1_haar] = wavedec(y1,10,'haar');

% % % For y1 with 'haar' wavelet
[a10_y1_haar, d1_y1_haar, d2_y1_haar, d3_y1_haar, d4_y1_haar, d5_y1_haar, ...
d6_y1_haar, d7_y1_haar, d8_y1_haar, d9_y1_haar, d10_y1_haar] = wavelet_coeff_calculation(y1, 'haar', "y1_haar_");

D_y1_haar = cell(1, 10);
A_y1_haar = cell(1, 1);

for i = 1:10
    D_y1_haar{i} = wrcoef('d', C_y1_haar, L_y1_haar, 'haar', i);
end
A_y1_haar{1} = wrcoef('a', C_y1_haar, L_y1_haar, 'haar', 10);

figure('Name', 'y1_haar_det_approx_fn', 'NumberTitle', 'off');
subplot(11, 1, 1);
plot(A_y1_haar{1});
title(['A^' num2str(10)]);

for i = 1:10
    subplot(11, 1, i+1);
    plot(D_y1_haar{i});
    title(['D^' num2str(i)]);
end

y1_haar = 0;
for i = 1:10
    y1_haar = y1_haar + D_y1_haar{i};
end
y1_haar = y1_haar + A_y1_haar{1};


energy_diff_y1_haar = sum(y1.^2)- sum(y1_haar.^2);
fprintf('energy differnece between y1 and reconstructed y1 using haar: %.14f\n', abs(energy_diff_y1_haar));

% y1 with db9
[C_y1_db9,L_y1_db9] = wavedec(y1,10,'db9');
D_y1_db9 = cell(1, 10);
A_y1_db9 = cell(1, 1);

for i = 1:10
    D_y1_db9{i} = wrcoef('d', C_y1_db9, L_y1_db9, 'db9', i);
end
A_y1_db9{1} = wrcoef('a', C_y1_db9, L_y1_db9, 'db9', 10);

figure('Name', 'y1_haar_det_approx_fn', 'NumberTitle', 'off');
subplot(11, 1, 1);
plot(A_y1_db9{1});
title(['A^' num2str(10)]);

for i = 1:10
    subplot(11, 1, i+1);
    plot(D_y1_db9{i});
    title(['D^' num2str(i)]);
end

y1_db9 = 0;
for i = 1:10
    y1_db9 = y1_db9 + D_y1_db9{i};
end
y1_db9 = y1_db9 + A_y1_db9{1};


energy_diff_y1_db9 = sum(y1.^2)- sum(y1_db9.^2);
fprintf('energy differnece between y1 and reconstructed y1 using db9: %.14f\n', abs(energy_diff_y1_db9));

figure;
subplot(2,1,1);
plot(y1_haar, 'y','LineWidth',1.2);
title('Reconstructed y1 with haar wavelet');
xlabel('t'); ylabel('y1_haar'); grid on;

subplot(2,1,2);
plot(y1_db9, 'y','LineWidth',1.2);
title('Reconstructed y1 with db9 wavelet');
xlabel('t'); ylabel('y1_db9'); grid on;


% y2 with haar
[C_y2_haar,L_y2_haar] = wavedec(y2,10,'haar');
D_y2_haar = cell(1, 10);
A_y2_haar = cell(1, 1);

for i = 1:10
    D_y2_haar{i} = wrcoef('d', C_y2_haar, L_y2_haar, 'haar', i);
end
A_y2_haar{1} = wrcoef('a', C_y2_haar, L_y2_haar, 'haar', 10);

figure('Name', 'y2_haar_det_approx_fn', 'NumberTitle', 'off');
subplot(11, 1, 1);
plot(A_y2_haar{1});
title(['A^' num2str(10)]);

for i = 1:10
    subplot(11, 1, i+1);
    plot(D_y2_haar{i});
    title(['D^' num2str(i)]);
end

y2_haar = 0;
for i = 1:10
    y2_haar = y2_haar + D_y2_haar{i};
end
y2_haar = y2_haar + A_y2_haar{1};


energy_diff_y2_haar = sum(y2.^2)- sum(y2_haar.^2);
fprintf('energy differnece between y2 and reconstructed y2 using haar: %.14f\n', abs(energy_diff_y2_haar));

% y2 with db9
[C_y2_db9,L_y2_db9] = wavedec(y2,10,'db9');
D_y2_db9 = cell(1, 10);
A_y2_db9 = cell(1, 1);

for i = 1:10
    D_y2_db9{i} = wrcoef('d', C_y2_db9, L_y2_db9, 'db9', i);
end
A_y2_db9{1} = wrcoef('a', C_y2_db9, L_y2_db9, 'db9', 10);

figure('Name', 'y2_db9_det_approx_fn', 'NumberTitle', 'off');
subplot(11, 1, 1);
plot(A_y2_db9{1});
title(['A^' num2str(10)]);

for i = 1:10
    subplot(11, 1, i+1);
    plot(D_y2_db9{i});
    title(['D^' num2str(i)]);
end

y2_db9 = 0;
for i = 1:10
    y2_db9 = y2_db9 + D_y2_db9{i};
end
y2_db9 = y2_db9 + A_y2_db9{1};


energy_diff_y2_db9 = sum(y2.^2)- sum(y2_db9.^2);
fprintf('energy differnece between y2 and reconstructed y2 using db9: %.14f\n', abs(energy_diff_y2_db9));

figure;
subplot(2,1,1);
plot(y2_haar, 'y','LineWidth',1.2);
title('Reconstructed y2 with haar wavelet');
xlabel('t'); ylabel('y2_haar'); grid on;

subplot(2,1,2);
plot(y2_db9, 'y','LineWidth',1.2);
title('Reconstructed y2 with db9 wavelet');
xlabel('t'); ylabel('y2_db9'); grid on;