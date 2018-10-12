%Assignment 3 CRL707
clc;
clear all;
close all;

[A, Fs] = audioread('deepika.wav');
plot(A);
Seven = A(1.16e5:1.24e5);
plot(Seven)
xlabel('sample value');
ylabel('amplitude');
title('Seven'),
A=A./(1.01*abs(max(A)));

% %Unvoiced Segment
% uvs = Seven(200:520);
% c1 = rceps(uvs);
% t = 1:length(uvs);
% figure(1)
% plot(t/16000,c1);
% xlabel('Quefrency (sec) ---->');
% ylabel('Amplitude ---->')
% title('Cepstrum for Unvoiced Segment of 20ms');
% % High Liftering
% c1l = c1(1:160);
% L=zeros(1,length(c1l));
% L(15:length(L))=1;
% c1lh = real(c1l.*L');
% [val id] = max(c1lh);
% pitch_period = id;
% pitch_freq = (1/pitch_period)*Fs
% %High liftering log spectrum
% c1lhf = fft(c1l(15:length(c1l)),16000);
% c1lhf = real(c1lhf(1:8000));
% figure(2)
% plot(c1lhf)
% xlabel('Frequency (Hz) ----->')
% ylabel('Magnitude (in dB) ----->')
% title('Log Spectrum of Unvoiced Segment with High pass liftering')
% 
% 
% % Low pass liftering
% c1ll = fft(c1l(1:15),16000);
% c1ll = real(c1ll(1:8000));
% figure(3)
% plot(c1ll)
% xlabel('Frequency (Hz) ----->')
% ylabel('Magnitude (in dB) ----->')
% title('Log Spectrum of Unvoiced Segment with Low pass liftering')
% 
% %Formant frequencies
% k =1;
% for i=2:length(c1ll)-1
%     if(c1ll(i-1)<c1ll(i)&&c1ll(i)>c1ll(i+1))
%         Formant(k) = i;
%         k=k+1;
%     end
% end

%Voiced Segment
vs = Seven(5000:5320);

c1 = fftshift(ifft(log(abs(fft(vs)))));
t = 1:length(vs);
figure(1)
plot(t/16000,c1);
xlabel('Quefrency (sec) ---->');
ylabel('Amplitude ---->')
title('Cepstrum for Voiced Segment of 20ms');
% High Liftering
c1l = c1(1:160);
L=zeros(1,length(c1l));
L(15:length(L))=1;
c1lh = real(c1l.*L');
[val id] = max(c1lh);
pitch_period = id;
pitch_freq = (1/pitch_period)*Fs

% % Low pass liftering
% c1ll = fft(c1l(1:15),16000);
% c1ll = real(c1ll(1:8000));
% figure(2)
% plot(c1ll)
% xlabel('Frequency (Hz) ----->')
% ylabel('Magnitude (in dB) ----->')
% title('Log Spectrum of Voiced Segment with Low pass liftering')
% 
% %Formant frequencies
% k =1;
% for i=2:length(c1ll)-1
%     if(c1ll(i-1)<c1ll(i)&&c1ll(i)>c1ll(i+1))
%         Formant(k) = i;
%         k=k+1;
%     end
% end Formant