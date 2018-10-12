[y,fs]= audioread('deepika.wav');
plot(y);
q = y(7.6e4:6.7e4);
f=[0:2000];
spectrogram(q,160,80,f,fs);