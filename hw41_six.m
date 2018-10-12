clc;
clear all;
close all;

%% Part-1:Plot the magnitude spectrum of voiced and unvoiced downsampled
% segments

%\\ Part-1.a :Plotting the time waveform for the utterance six 
%[six1,fs1,nbits1]=wavread('pratik_six'); %% \\CHANGE HERE %%
[six1,fs1]=audioread('6.wav');
N1=length(six1);
m1=0:N1-1;
t1=m1/fs1;

fs2=8000;% downsampled frequency
six2=resample(six1,fs2,fs1);
N2=length(six2);
m2=0:N2-1;
t2=m2/fs2;

figure();
subplot(2,1,1);plot(t1,six1,'k');
xlabel('time(in sec)');xlim([min(t1) max(t1)]);ylabel('Amplitude');title('Time waveform for the original utterance - six');
subplot(2,1,2);plot(t2,six2,'k');
xlabel('time(in sec)');xlim([min(t2) max(t2)]);ylabel('Amplitude');title('Time waveform for the downsampled utterance - six');

%\\Part-1.b :extracting 20 ms of voiced and unvoiced speech and plotting them 

%\\\ CHANGE HERE
n_unvoiced=1201:1360;% 20ms segment of 8 KHz speech means 160 samples
unvoiced_six=six2(n_unvoiced);
t_unvoiced=n_unvoiced/fs2;

%\\\ CHANGE HERE
n_voiced=1841:2000; % 20ms segment of 8 KHz speech means 160 samples
voiced_six=six2(n_voiced);
t_voiced=n_voiced/fs2;

figure();
subplot(2,1,1);
plot(t_voiced,voiced_six,'k');
xlabel('time(in sec)');xlim([min(t_voiced) max(t_voiced)]);ylabel('Amplitude');
title('Voiced Segment for the utterance - six');%ylim([-0.03 0.03]);
subplot(2,1,2);
plot(t_unvoiced,unvoiced_six,'k');
xlabel('time(in sec)');xlim([min(t_unvoiced) max(t_unvoiced)]);ylabel('Amplitude');
title('Unvoiced Segment for the utterance - six');%ylim([-0.005 0.005]);

%\\Part-1.c:Magnitude spectrum of voiced and unvoiced segments

nfft=1024;%2^nextpow2(length(n_voiced))
wd=length(voiced_six);
hvoiced_six=voiced_six.*hamming(wd);
hunvoiced_six=unvoiced_six.*hamming(wd);


f=fs2/2*linspace(0,1,nfft/2);
y1=(fft(hvoiced_six,nfft));
yv=20*log10(abs(y1(1:nfft/2)));
% [ueyv,leyv]=envelope(yv);
y2=(fft(hunvoiced_six,nfft));
yuv=20*log10(abs(y2(1:nfft/2)));
% [ueyuv,leyuv]=envelope(yuv);

figure()
subplot(2,1,1);plot(f,yv,'k'); xlabel('Frequency(Hz)');ylabel('Magnitude(dB)');ylim([-82 10]);
title('Magnitude spectrum of Voiced Segment of the utterance - six');
subplot(2,1,2);plot(f,yuv,'k');xlabel('Frequency(Hz)');ylabel('Magnitude(dB)');ylim([-82 0]);
title('Magnitude spectrum of Unvoiced Segment of the utterance - six');


%% Part-2:linear predictor model

%\\ Part-2.a:10th order linear predictor Ai(z) using the Autocorrelation method
% and the gain G

N=160;
p=10;

% LPC for voiced segment
[va,vg]   = lpc(voiced_six,p);
va1=va(2:end);
% LPC for unvoiced segment
[uva,uvg] = lpc(unvoiced_six,p);
uva1=uva(2:end);

% f=fs2/2*linspace(0,1,nfft/2);
ey1=(fft(va,nfft));
eyv=20*log10(abs(ey1(1:nfft/2)));
ey2=(fft(uva,nfft));
eyuv=20*log10(abs(ey2(1:nfft/2)));

figure()
subplot(2,1,1);plot(f,yv,'k',f,-eyv-20,'k--');% adjust 20 to superimpose on original signal
xlabel('Frequency(Hz)');ylabel('Magnitude(dB))');ylim([-85 10]);
title('Magnitude spectrum of Voiced Segment of the utterance - six');
legend('Original','LPC estimate');
subplot(2,1,2);plot(f,yuv,'k',f,-eyuv-40,'k--');% adjust 20 to supermipose on original signal
xlabel('Frequency(Hz)');ylabel('Magnitude(dB))');ylim([-85 10]);
title('Magnitude spectrum of Unvoiced Segment of the utterance - six');
legend('Original','LPC estimate');

%\\ Part-2.b: Display LPC and Gain for voiced and unvoiced segments

% est_voiced_six = filter([0 -va(2:end)],1,voiced_six);
s=voiced_six;
for n=0+1:N-1+p+1
sum=0;
for i=1:p
    m=n-i;
    if((0< m)&&(m<= N))
    sum=sum+va(i+1)*s(n-i);
    end
end   
est_s(n)=sum;
end
est_voiced_six=est_s';
voiced_six_pad=padarray(voiced_six,[p 0],'post');% make length equal
e_voiced_six = voiced_six_pad- (-est_voiced_six);% Residual signal for voice segment
he_voiced_six=e_voiced_six .*hamming(170);

R0_v=max(xcorr(voiced_six));
R0_ev=max(xcorr(est_voiced_six));
% p1=e_voiced_six.^2;
% energy_v=sum(p1,2);
energy_v=mean(e_voiced_six.^2);
Gv=sqrt(energy_v);% gain for voiced segment
Gv1=sqrt(R0_v-R0_ev);

% est_unvoiced_six = filter([0 -uva(2:end)],1,unvoiced_six);
N=160;
s=unvoiced_six;
for n=0+1:N-1+p+1
sum=0;
for i=1:p
    m=n-i;
    if((0< m)&&(m<= N))
    sum=sum+uva(i+1)*s(n-i);

    end

end   
est_s(n)=sum;
end
est_unvoiced_six=est_s';
unvoiced_six_pad=padarray(unvoiced_six,[p 0],'post');
e_unvoiced_six = unvoiced_six_pad- (-est_unvoiced_six);% Residual signal for unvoice segment
he_unvoiced_six=e_unvoiced_six .*hamming(170);

energy_uv=mean(e_unvoiced_six.^2);
Guv=sqrt(energy_uv);
R0_uv=max(xcorr(unvoiced_six));
R0_euv=max(xcorr(est_unvoiced_six));
Guv2=sqrt(R0_uv-R0_euv);

display(va);
display(Gv);
display(uva);
display(Guv);


%\\ Part-3.a: the linear prediction residual signal for each segment 
% the spectral magnitude

% Residul signal Analysis
nfft=1024;
w1=fftshift(fft(he_voiced_six,nfft));
wv=20*log10(abs(w1(1:nfft/2)));
w2=fftshift(fft(he_unvoiced_six,nfft));
wuv=20*log10(abs(w2(1:nfft/2)));

%%Estimated signal analysis
% wd=length(est_voiced_six);
% hest_voiced_six=est_voiced_six.*hamming(wd);
% hest_unvoiced_six=est_unvoiced_six.*hamming(wd);

% f=fs2/2*linspace(0,1,nfft/2);
% z1=fftshift(fft(hest_voiced_six,nfft));
% zv=20*log10(abs(z1(1:nfft/2)));
% z2=fftshift(fft(hest_unvoiced_six,nfft));
% zuv=20*log10(abs(z2(1:nfft/2)));


figure()
subplot(2,1,1);plot(1:length(e_voiced_six),e_voiced_six,'k');xlim([1 length(e_voiced_six)]);
xlabel('Sample Number');ylabel('Amplitude');title('Residual signal for Voiced Segment');
subplot(2,1,2);plot(f,wv,'k'); xlabel('Frequency(Hz)');ylabel('Magnitude(dB))');ylim([-85 0]);
title('Magnitude Spectrum of Residual signal for Voiced Segment'); 

figure()
subplot(2,1,1);plot(1:length(e_unvoiced_six),e_unvoiced_six,'k');xlim([1 length(e_unvoiced_six)]);
xlabel('Sample Number');ylabel('Amplitude');title('Residual signal for Unvoiced Segment');
subplot(2,1,2);plot(f,wuv,'k'); xlabel('Frequency(Hz)');ylabel('Magnitude(dB))');ylim([-85 0]);
title('Magnitude Spectrum of Residual signal for Unvoiced Segment'); 

%\\ Part-4.a: autocorrelation function for speech segment and residual signal

[vacs,vlags] = xcorr(voiced_six);
[veacs,velags] = xcorr(e_voiced_six);
% [uvacs,uvlags] = xcorr(unvoiced_six);
% [uveacs,uvelags] = xcorr(e_unvoiced_six);

figure();
subplot(3,1,1);plot(1:length(n_voiced),voiced_six,'k');xlabel('Sample Number');ylabel('Amplitude');xlim([1 160]); title('Original voiced segment');
subplot(3,1,2);plot(vlags,vacs,'k');xlabel('Lags');ylabel('Amplitude');xlim([min(vlags) max(vlags)]); title('Autocorrelation of voiced segment');
subplot(3,1,3);plot(velags,veacs,'k');xlabel('Lags');ylabel('Amplitude');xlim([min(velags) max(velags)]);title('Autocorrelation of voiced residual segment');

