%% Assignment 3
clc;
clear all;
close all;

%% Reading audio file in wave format
[A, Fs] = audioread('akashdeep.wav');
% plot(A)

%% Zero
zero = A(2.85e4:3.5e4);
c = cceps(zero);
plot(c)