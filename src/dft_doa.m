close all;
clear;
clc

f0 = 9.6e9;
T = 1/f0;

tt = linspace(0,1*T, 1*1024);
fs = 1024*f0;
ts = 1/fs;

s = exp(1i*2*pi*f0*tt);

N_rx = 8;
c = 3e8;
lambda = c/f0;
d = lambda/2;

theta = [-40 30 -5];
N_src = length(theta);
theta_rad = deg2rad(theta);

P_tx = s*s'/length(s);
SNR = 100;

P_n = P_tx/(ts*10^(SNR/10));
variance = 1;

n = sqrt(variance)*(randn(1, length(s)) + 1i*randn(1, length(s)))/sqrt(2);
n = repmat(n, N_rx, 1);
x = zeros(N_rx, length(s));
xn = {};
for i = 1:N_src
    xn{i} = s.*exp(-1i*2*pi*(0:N_rx-1)'*d*sin(theta_rad(i))/lambda);
end

for i = 1:N_src
    x = x + xn{i};
end

% x = x + n;
% mn = sum(x,2)/length(x);
% x = x - mn;

Nfft = N_rx*8;
X = fftshift(fft(x, Nfft, 1));

Xs = sum(X,2);

thm = asin(2*(-Nfft/2:Nfft/2-1)/Nfft);

figure
subplot(2,1,1)
stem((-Nfft/2:Nfft/2-1), abs(Xs).^2)
grid on
subplot(2,1,2)
stem((-Nfft/2:Nfft/2-1), rad2deg(thm))
grid on

%MUSIC
th = linspace(-90, 90, 18000);
th_rad = deg2rad(th);
a = exp(-1i*2*pi*(0:N_rx-1)'*d.*sin(th_rad)/lambda);

R = x*x'/length(x);
[V,D] = eig(R);
[vi, di] = sort(diag(D),'descend');
VV = V(:,di);
En = VV(:, N_src+1:end);

for i = 1:length(a)
    P(i) = 1/(a(:,i)'*En*En'*a(:,i));
end
Pm = abs(P);
P_max = max(abs(P));
% figure
% plot(th, 10*log10(Pm/P_max))