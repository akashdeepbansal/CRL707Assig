%% CRL707 Assignment 1

clc;
clear all;
close all;
%% Reading audio file in wave format
[A, Fs] = audioread('akashdeep.wav');

%plot(A)

%% Zero
zero = A(2.85e4:3.5e4);
%plot(A(2.85e4:3.5e4))
Ez = sum((abs(zero)).^2);
Pz = sum((abs(zero)).^2)*Fs/length(zero)

noise = A(1:2.8499e4);
En = sum((abs(noise)).^2);
Pn = sum((abs(noise)).^2)*Fs/length(noise)

SNR = Pz/Pn

%% Spectrogram
%s = spectrogram(zero)
f = 0:4e3;
%spectrogram(zero,hamming(160),80,f,Fs)

%% One
% one = A(4.48e4:5.02e4);
% % plot(A(4.48e4:5.02e4))
% spectrogram(one,hamming(160),80,f,Fs)
% 
% %% Two
% two = A(6.14e4:6.57e4);
% % plot(A(6.14e4:6.57e4))
% spectrogram(two,hamming(160),80,f,Fs)
% 
% %% Three
% three = A(7.64e4:8.18e4);
% % plot(A(7.64e4:8.18e4))
% spectrogram(three,hamming(160),80,f,Fs)
% 
% %% Four
% four = A(9.17e4:9.82e4);
% % plot(A(9.17e4:9.82e4))
% spectrogram(four,hamming(160),80,f,Fs)
Max4 = A(94700:94860);
% plot(Max4)
Emax4 = sum((abs(Max4)).^2)
SPLmax4 = 20*log(Emax4/20e-8)

Min4 = A(92700:92860);
% plot(Min4)
Emin4 = sum((abs(Min4)).^2)
SPLmin4 = 20*log(Emin4/20e-8)

% 
% %% Five
% five = A(10.9e4:11.4e4);
% % plot(A(10.9e4:11.4e4))
% spectrogram(five,hamming(160),80,f,Fs)
Max5 = A(11e4:110160);
% plot(Max5)
Emax5 = sum((abs(Max5)).^2)
SPLmax5 = 20*log(Emax5/20e-8)
% 

Min5 = A(112800:112960);
% plot(Min5)
Emin5 = sum((abs(Min5)).^2)
SPLmin5 = 20*log(Emin5/20e-8)

% %% Six
% six = A(12.2e4:13.25e4);
% % plot(A(12.2e4:13.25e4))
% spectrogram(six,hamming(160),80,f,Fs)
% 
% %% Seven
% seven = A(13.8e4:14.7e4);
% % plot(A(13.8e4:14.7e4))
% spectrogram(seven,hamming(160),80,f,Fs)
% 
% %% Eight
% eight = A(15.6e4:16.05e4);
% % plot(A(15.6e4:16.05e4))
% spectrogram(eight,hamming(160),80,f,Fs)
% 
% %% Nine
% nine = A(17.05e4:17.8e4);
% % plot(A(17.05e4:17.8e4))
% spectrogram(nine,hamming(160),80,f,Fs)
