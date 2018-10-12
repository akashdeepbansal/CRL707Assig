%Assignment 4 CRL707
clc;
clear all;
close all;

%Reading Input signal
[A,Fs] = audioread('akashdeep.wav');
A = A./max(abs(A));

%Sampling frequency reduced by half
j=1;
for i=1:2:length(A)
    a(j) = A(i);
    j=j+1;
end

%length(A)
%length(a)
%plot(a)

%Spectrum of the voiced segment
vs = a(39501:39660);%38401:38560);%47161:47320);
N = length(vs);
%w = hamming(N);
%vs = vs.*w;
vs_fft = 10*log10(abs(fft(vs)));
f = Fs/2*(0:N-1)/N;
figure(1)
plot(f,vs_fft)

%Linear Prediction Filter
for i=1:11
    s1 = vs(i:N-i);
    s2 = vs(1:N-2*i+1);
    R(i) = s1*s2';
end

%% RL and RR Matrix
RR = R(2:11);
for i=1:10
    for j=1:10
        RL(i,j) = R(abs(i-j)+1);
    end
end

a_c = inv(RL)*RR';
b(1) = 1;
b(2:11) = -a_c;

% Spectrum of the filter function
[h,w] = freqz(1,b);
figure(2)
plot(w/pi,20*log10(abs(h)));

%% Residual Signal
for i = 11:N
    e(i) = a(i) - a_c'*a(i-1:-1:i-10)';
end
figure(3)
plot(e)