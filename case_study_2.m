%% Case Study 2: PAM Communication

Ts = .1; % symbol time
dt = .01; % sample period
t = -5*Ts : dt : 5*Ts; % time vector
L = 1024; % length of the signal

fs = 1/Ts; % sampling frequency
f = -fs/2:fs/L:fs/2-fs/L; % frequency range 

% 1st pulse shape: truncated sinc
p1 = sinc(t/Ts) .* double(abs(t) <= 3*Ts);
P1 = fft(p1,L);

figure;
subplot(2,1,1), plot(t,p1);
title('p1(t): truncated sinc');
xlabel('time(s)'); ylabel('p1(t)');
subplot(2,1,2), plot(f,fftshift(abs(P1/N)));
title('fourier transform of p1(t)');
xlabel('frequency(Hz)'); ylabel('|P1|');

% 2nd pulse shape: rectangular pulse with Gaussian apodization
p2 = exp(-pi*t.^2/Ts.^2) .* double(abs(t) <= Ts);
P2 = fft(p2,L);

figure;
subplot(2,1,1), plot(t,p2);
title('p2(t): rectangular pulse with Gaussian apodization');
xlabel('time(s)'); ylabel('p1(t)');
subplot(2,1,2), plot(f,fftshift(abs(P2/N)));
title('fourier transform of p2(t)');
xlabel('frequency(Hz)'); ylabel('|P2|');
