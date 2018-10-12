%% Assignment 2 CRL707
clc;
clear all;
close all;

%% Reading audio file in wave format
[A, Fs] = audioread('deepika.wav');
plot(A)

%% Zero
zero = A(2.18e5:2.22e5);
%plot(zero)
xlabel('time(sec)')
ylabel('amplitude')
title('nine')
%spectrogram(zero,hamming(960),180,f,fs,'yaxis');
 %plot(A(2.85e4:3.5e4))

%%Rectangular window
% r = rectwin(320)
% plot(r)

N = length(zero);
nf = ceil(N/160)-1;
for i=1:nf 
    E(i) = sum((abs(zero(160*(i-1)+1:160*(i-1)+160))).^2);
end
%plot(E)
ylabel('Energy')
xlabel('Frame')
title('nine')