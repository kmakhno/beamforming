close all;
clear;
clc;

c = physconst('LightSpeed');
f0 = 24e9;
B = 100e6;
T = 100e-6;
alpha = B/T;
fs = 2*B;
Np = 128;
%coherent processign interval
Tcpi = Np*T;
%compression rate
comp = B*T;

%LFM signal
t = 0:1/fs:Np*T-1/fs;
N = floor(fs*T*Np);
s_tx = exp(1i*2*pi*(-B/2*mod(t,T) + .5*alpha*mod(t,T).^2));

%Distance to source
r = [500 3000];
%Speed of targets
v = [0 10];
mov = v'*t;
N_src = length(r);
%different delay for each channel
td = 2*(r'-mov)/c;
%Received signal without noize
s = {};
for i = 1:N_src
    s{i} = exp(1i*2*pi*(-B/2*mod(t-td(i,:),T) + .5*alpha*mod(t-td(i,:),T).^2)).*exp(1i*2*pi*f0*td(i,:));
end

%Find total signal at each channel from each source
s_rx = zeros(1, N);
for i = 1:N_src
    s_rx = s_rx + s{i};
end

SNR = 15;
%Power of transmit signal
P_tx = s_tx*s_tx'/N;
%Noise power
P_n = P_tx/10^(SNR/10);
variance = P_n;
noise = sqrt(variance)*(randn(1,N) + 1i*randn(1,N))/sqrt(2);
%Received signal with noize
s_rx_n = s_rx;% + noise;

%% Target parameters estimation
%Range processing
df = 4;
rmax = c/2/B*(round(fs*T/df));
%matched filtering 
s_rx_0 = conj(s_rx_n) .* s_tx;
s_rx_decim = decimate(s_rx_0, df);
S_RX = reshape(s_rx_decim, N/Np/df, Np);

%Double delay-line canceler(needed to remove stationary clutter)
%Richards 5.2.1 Pulse Cancellers
% H(z) = 1 - 2*z^-1 + z^-2
b = [1 -2 1];
a = 1;
S_RX = filter(b, a, S_RX, [], 2);

%Range processing
Srange = fft(S_RX, [], 1);
%I don't know that's true...
Srange_sum = sum(Srange,2);
[mx, indx] = max(abs(Srange_sum),[],1);
%computed distance to target
Nrfft = length(Srange_sum);
ff = (-Nrfft/2:Nrfft/2-1)*fs/Nrfft/df;
dist = ff(indx)*T*c/(2*B);

%Doppler processing
Ndfft = Np;
Sdop = fftshift(fft(Srange, Ndfft, 2),2);
Sdop_sum = sum(Sdop,1)/Ndfft;

%% Plots
dd = ff*T*c/(2*B);
figure
stem(dd, abs(fftshift(Srange_sum)))
title('Range processing');
xlabel('Distance, m');
grid on

fd = (-Ndfft/2:Ndfft/2-1)*1/Ndfft/T;
vv = fd*c/(2*f0);
figure
stem(vv, abs(Sdop_sum));
title('Doppler processing');
xlabel('Velocity, m/s');
grid on

figure
surf(vv, dd, 10*log10(abs(fftshift(Sdop,1))), 'EdgeColor', 'none');
title('Range-Doppler map');
ylabel('Range, m');
xlabel('Velocity, m/s');
xlim([vv(1) vv(end)]);
ylim([dd(1) dd(end)]);
view(2)