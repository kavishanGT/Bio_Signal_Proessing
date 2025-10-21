clearvars; close all; clc;

fs = 250;           % Sample frequency
N = 3000;           % Half data length -> total length = 2N+1
t = (-N:N)/fs;      % Time scale
dt = 1/fs;
s_list = 0.01:0.1:2; % Values of scaling factor

% --- Step 1: Define base Mexican hat (unnormalized mother) ---
a = 1/sqrt(2*pi);                            % Gaussian prefactor
m_t = a .* (1 - t.^2) .* exp(-t.^2/2);       % Unnormalized Mexican hat
E_m = trapz(t, m_t.^2);                      % Energy of unnormalized
A = 1/sqrt(E_m);                             % Normalization constant


means = zeros(1,numel(s_list));
energies = zeros(1,numel(s_list));
% --- Step 2: Generate daughter wavelets for each scale ---
wavelt = zeros(numel(s_list), length(t));    % Store all daughter wavelets
for k = 1:numel(s_list)
    s = s_list(k);
    psi_s = (1/sqrt(s)) .* A .* a .* (1 - (t./s).^2) .* exp(-(t./s).^2/2);
    wavelt(k,:) = (1/sqrt(s)) .* A .* a .* (1 - (t./s).^2) .* exp(-(t./s).^2/2);
    means(k) = trapz(t, psi_s);
    energies(k) = trapz(t, psi_s.^2);
end

T = table(s_list.', means.', energies.', ...
          'VariableNames', {'Scale','Mean','Energy'});
disp(T);

% --- Step 3: Plot time-domain wavelets (subset for clarity) ---
figure;
plot(t, wavelt.')
xlim([-6 6])
ylim([-2 3])
xlabel('Time (s)'), ylabel('\psi_s(t)')
title('Mexican hat daughter wavelets (all scales)')
grid on
legend(arrayfun(@(x) sprintf('s=%.2f',x), s_list, 'UniformOutput', false), ...
       'Location','eastoutside')



% --- Step 4: Generate and plot spectrum of a chosen wavelet (e.g. s=1) ---
%idx = find(abs(s_list-1)<1e-6, 1);   % choose scale s=1
figure;
hold on
for k = 1:numel(s_list)
    Fwavelt = fft(wavelt(k,:)) / length(t);    % spectrum for scale k
    hz = linspace(0, fs/2, floor(length(t)/2)+1);
    plot(hz, 2*abs(Fwavelt(1:length(hz))))
end
hold off
xlim([0 5])
xlabel('Frequency (Hz)'), ylabel('Amplitude')
title('Spectra of Mexican hat daughter wavelets (all scales)')
grid on
legend(arrayfun(@(x) sprintf('s=%.2f',x), s_list, 'UniformOutput', false), ...
       'Location','northeastoutside')   % show subset of scales in legend
% --- Step 5: Spectrogram-like image (scale vs time) ---
figure;
imagesc(t, s_list, wavelt);
set(gca,'YDir','normal')
xlabel('Time (s)'), ylabel('Scale (s)')
title('Spectrogram-like representation of daughter wavelets')
colormap jet; colorbar;
