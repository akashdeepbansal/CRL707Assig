clear all;
clc;
close all;
[y,fs]= audioread('deepika.wav');

l1 = 7.6e4-6.7e4;
nf = ceil(l1/480);
yn = y(1:480);
snr1=[];
c=6.7e4;
for i=1:nf 
    y1= y(c:c+479);
    snr1(i)=snr(y1,yn) ;
    
    c= c+480 ;
    
end     
   
%yn = y(1:480);
%snr=[];
%for i=1:nf 
   % snr(i)=snr(y1(i),yn) ;
%end





%l=length(y);
%x=(0:(1/fs):((l-1)/fs))
figure(1)
xlabel('time(sec)')
ylabel('amplitude')
plot(snr1)
xlabel('No. of Frames ')
ylabel('Signal to Noise Ratio(zero)')

%for narrow band spectrogram 
zero = y(2.03e5:2.06e5);
f = 0:400;
spectrogram(zero,hamming(960),180,f,fs);

Max5 = y(2.038e5:203960);
% plot(Max5)
Emax5 = sum((abs(Max5)).^2)
SPLmax5 = 20*log(Emax5/20e-8)
% 

Min5 = y(2.03e5:203160);
% plot(Min5)
Emin5 = sum((abs(Min5)).^2)
SPLmin5 = 20*log(Emin5/20e-8)

% %% Six
% six = A(12.2e4:13.25e4);
% % plot(A(12.2e4:13.25e4))
% spectrogram(six,hamming(160),80,f,Fs)
% 



