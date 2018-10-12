%% Assignment 2 CRL707
clc;
clear all;
close all;

%% Reading audio file in wave format
[A, Fs] = audioread('akashdeep.wav');
% plot(A)

%% Zero
zero = A(2.85e4:3.5e4);
% plot(A(2.85e4:3.5e4))

%%Rectangular window
% r = rectwin(320)
% plot(r)

N = 3.5e4-2.85e4;
nf = ceil(N/160)-1;
for i=1:nf 
    E(i) = sum((abs(zero(160*(i-1)+1:160*(i-1)+160))).^2);
end
plot(E)