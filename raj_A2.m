clc;
clear all;
close all;
[a,fs]=audioread('akashdeep.wav');
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
digit2=a(6.14e4:6.57e4); %two = A(6.14e4:6.57e4)
figure(1)
plot(digit2);
t=1:length(digit2);
plot(t/16000,digit2);
xlabel('Time(msec)');
ylabel('Amplitude');
title('Digit 2');
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
window_length=0.010;
window_overlap=0.005;
 
frame_size=fs*window_length;
frame_overlap=fs*window_overlap;
energy=[];
number_of_frames=floor(length(digit2)/frame_overlap)-1;
for i=1:number_of_frames
    signal=digit2((i-1)*frame_overlap+1:(i-1)*frame_overlap+frame_size);
    ener=0;
    for j=1:length(signal)
        ener=ener+signal(j)*signal(j);
    end
    energy=[energy;ener];
end
figure(2)
stem(energy);
xlabel('Frame Number');
ylabel('Energy');
title('Short-Time Energy of Digit 2');
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3);
window_length=0.003;
window_overlap=0.0015;
 
frame_size=fs*window_length;
frame_overlap=fs*window_overlap;
t=1:length(digit2);
f=0:8000;
plot(t/16000,digit2);
figure(3);
%spectrogram(digit2,hamming(frame_size),frame_overlap,f,fs,'yaxis');
colormap(flipud(gray));
title('Digit 2 Spectrogram');
