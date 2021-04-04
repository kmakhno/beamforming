close all;
clear;
clc

f0 = 9.6e9;
T = 1/f0;

tt = linspace(0,1*T, 1*1024);
fs = 1024*f0;
ts = 1/fs;

s = exp(1i*2*pi*f0*tt);

N_rx = 4;
c = 3e8;
lambda = c/f0;
d = lambda/2;

theta = [-17 7];
N_src = length(theta);
theta_rad = deg2rad(theta);

P_tx = s*s'/length(s);
SNR = 100;

P_n = P_tx/(ts*10^(SNR/10));
variance = 0.01;

n = sqrt(variance)*(randn(1, length(s)) + 1i*randn(1, length(s)))/sqrt(2);
n = repmat(n, N_rx, 1);
x = zeros(N_rx, length(s));
xn = {};
for i = 1:N_src
    xn{i} = s.*exp(1i*2*pi*(0:N_rx-1)'*d*sin(theta_rad(i))/lambda);
end

for i = 1:N_src
    x = x + xn{i};
end
% x = x/N_src;


%ESPRIT
mn = sum(x,2)/length(x);
x = x - mn;
R = x*x'/length(x);
[V,D] = eig(R);
[vi, di] = sort(diag(D),'descend');
VV = V(:,di);
Es = VV(:, 1:N_src);
