%% Case Study 2: PAM Communication

Ts = .1; % symbol time
dt = .01; % sample period
t = -5*Ts : dt : 5*Ts; % time vector
L = 1024; % length of the signal in frequency domain

fs = 1/dt; % length in frequency domain 
f = -fs/2:fs/L:fs/2-fs/L; % frequency range 

% 1st pulse shape: truncated sinc
p1 = sinc(t/Ts) .* double(abs(t) <= 3*Ts); % this truncating time can be varied
P1 = fft(p1,L);

figure;
subplot(2,1,1), plot(t,p1);
title('p1(t): truncated sinc');
xlabel('time(s)'); ylabel('p1(t)');
subplot(2,1,2), plot(f,fftshift(abs(P1/L)));
title('fourier transform of p1(t)');
xlabel('frequency(Hz)'); ylabel('|P1|');

% 2nd pulse shape: rectangular pulse with Gaussian apodization
p2 = exp(-pi*t.^2/Ts.^2) .* double(abs(t) <= Ts);
P2 = fft(p2,L);

figure;
subplot(2,1,1), plot(t,p2);
title('p2(t): rectangular pulse with Gaussian apodization');
xlabel('time(s)'); ylabel('p2(t)');
subplot(2,1,2), plot(f,fftshift(abs(P2/L)));
title('fourier transform of p2(t)');
xlabel('frequency(Hz)'); ylabel('|P2|');

% 3rd pulse shape: raised cosine filter
p3 = rcosdesign(0.25, 7 , 14);
t3 = - length(p3)/2 * dt : dt: length(p3)/2 * dt-dt; 
P3 = fft(p3,L);

figure;
subplot(2,1,1), plot(t3,p3);
title('p3(t): raised cosine filter');
xlabel('time(s)'); ylabel('p3(t)');
subplot(2,1,2), plot(f,fftshift(abs(P3/L)));
title('fourier transform of p3(t)');
xlabel('frequency(Hz)'); ylabel('|P3|');

%% Test  Nyquist filtering Criterion

%  Test first pulse 
num_of_sample = int32(length(p1)/(Ts/dt))+1;  
sum1 = zeros(1, num_of_sample); % initialize the container

sum1(1)= abs(P1(1));

for k = 2: num_of_sample
    sum1(k) = sum1(k)+abs(P1(k* (Ts/dt))); 
end
figure; 
indexing =  -5: 1 : 5; % indexing vector
plot(indexing, sum1); 
title('Sum of p1(t)  at integer multiple of Ts');
xlabel('Integer multiple of Ts '); ylabel('|P1|');


%  Test Second pulse 
num_of_sample = int32(length(p2)/(Ts/dt))+1;  
sum2 = zeros(1, num_of_sample); % initialize the container
sum2(1)= abs(P2(1));

for k = 2: num_of_sample
    sum2(k) = sum2(k)+abs(P2(k* (Ts/dt))); 
end

figure; 
plot(indexing, sum2); 
title('Sum of p2(t)  at integer multiple of Ts');
xlabel('Integer multiple of Ts '); ylabel('|P2|');


%  Test Third pulse 
midpoint = int32(length(p3)/2); 
interval = Ts/dt; 
half_length = int32(length(p3(1: midpoint))-1); 
half_num_sample = idivide(half_length, interval, 'floor'); 
num_of_sample = 1+2* half_num_sample; 
first_index = midpoint - interval * half_num_sample; 

sum3 = zeros(1, num_of_sample); % initialize the container
sum3(1)= abs(P3(first_index));

for k = 1: num_of_sample-1
    sum3(k) = sum3(k+1)+abs(P3(first_index + k* (Ts/dt))); 
end

indexing = -half_num_sample : half_num_sample; 
figure; 
plot(indexing, sum3); 
title('Sum of p3(t)  at integer multiple of Ts');
xlabel('Integer multiple of Ts '); ylabel('|P3|');
