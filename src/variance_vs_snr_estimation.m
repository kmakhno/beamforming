close all;
clear;
clc;

c = physconst('LightSpeed');
f0 = 36e9;
B = 100e6;
T = 100e-6;
alpha = B/T;
fs = 2*B;
Np = 2;
%ULA parameters
N_tx = 1;
N_rx = 4;
lambda = c/(f0+B/2);
d = lambda/2;
%signal to noise ratio in dB
SNR = 2;

%LFM signal
t = 0:1/fs:Np*T-1/fs;
N = floor(fs*T*Np);
s_tx = exp(1i*2*pi*(-B/2*mod(t,T) + .5*alpha*mod(t,T).^2));

th_min = deg2rad(-90);
th_max = deg2rad(90);
theta = linspace(th_min, th_max, 18000);
phi = 2*pi*(0:N_rx-1)'*d.*sin(theta)/lambda;
a = exp(-1i*phi);

%Source angle
theta0 = [deg2rad(0)];
%Distance to source
r = [1000];
N_src = length(theta0);
td = 2*r/c;
a0 = exp(-1i*2*pi*(0:N_rx-1)'.*d*sin(theta0)/lambda);
%Received signal without noize
s = {1,N_src};
for i = 1:N_src
    s{i} = exp(1i*2*pi*(-B/2*mod(t-td(i),T) + .5*alpha*mod(t-td(i),T).^2))*exp(1i*2*pi*f0*td(i)).*a0(:,i);
end

%Find total signal at each channel from each source
s_rx = zeros(N_rx, N);
for i = 1:N_src
    s_rx = s_rx + s{i};
end

%Power of transmit signal
P_tx = s_tx*s_tx'/N;
%Noise power
P_n = P_tx/10^(SNR/10);
variance = P_n;
noise = sqrt(variance)*randn(1,N);
%Received signal with noize
s_rx_n = s_rx;% + noise;

%Calculate covariance matrix of received signal
mn = sum(s_rx_n,2)/N;
s_rx_n = s_rx_n - mn;
R = s_rx_n*s_rx_n'/N;

%% MUSIC algorithm
[V,D] = eig(R);

En = V(:,1:N_rx-N_src);

for i = 1:length(a)
    pmu(i) = (a(:,i)'*a(:,i))/(((a(:,i)'*En)*En')*a(:,i));
end
%subspace spectrum
pmu = 10*log10(pmu/max(pmu));

%simple sum of rx signal from each channel
s_rx_sum = sum(s_rx_n,1);

rs = xcorr(s_tx, s_tx);

%% Plots
figure
for i = 1:N_rx
    subplot(N_rx, 1, i);
    plot(real(s_rx_n(i,:)));
end

figure
plot(rad2deg(theta), pmu);
title('MUSIC algorithm');
grid on 

figure
plot(real(s_rx_sum));
title('Simple rx signals sum');
grid on

figure
plot(real(rs));
title('Correlation');
grid on