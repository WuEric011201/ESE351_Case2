%% Case Study 2: PAM Communication

% fixed variables
Ts = .1; % symbol time
dt = Ts/50; % sample period
L = 1000; % length of the signal in frequency domain
interval = Ts/dt; 

fs = 1/dt; % length in frequency domain 
f = -fs/2 : fs/L : fs/2-fs/L; % frequency range 

xn1 = 2*((rand(1,11)>0.5)-0.5); % a 1 by 11 message for pulse analysis
xn2 = 2*((rand(1,10001)>0.5)-0.5); % a 1 by 10001 message for error analysis

% 1st pulse shape: truncated sinc
Tp1 = 5*Ts;
tp1 = -Tp1 : dt : Tp1; % time vector;
p1 = sinc(t/Ts) .* double(abs(t) <= 3*Ts); % this truncating time can be varied
P1 = fft(p1,L);

figure;
subplot(2,1,1), plot(t,p1);
title('p1(t): truncated sinc');
xlabel('time(s)'); ylabel('p1(t)');
subplot(2,1,2), plot(f,fftshift(abs(P1/L)));
title('Fourier transform of p1(t)');
xlabel('frequency(Hz)'); ylabel('|P1|');

pulse_analysis(p1, dt, Tp1, Ts, xn1, 20, 0.5);
error_rate(p1, dt, Tp1, Ts, xn2, 20, 1);

% % 2nd pulse shape: raised cosine filter
% tp2 = - length(p2)/2 * dt : dt: length(p2)/2 * dt-dt; 
% p2 = rcosdesign(0, 10 , Ts/dt);
% P2 = fft(p2,L);
% 
% figure;
% subplot(2,1,1), plot(tp2,p2);
% title('p2(t): raised cosine filter');
% xlabel('time(s)'); ylabel('p2(t)');
% subplot(2,1,2), plot(f,fftshift(abs(P2/L)));
% title('Fourier transform of p2(t)');
% xlabel('frequency(Hz)'); ylabel('|P2|');
% 
% pulse_analysis(p2, dt, Tp1, Ts, xn, 20, sigma)

% 2nd pulse shape: longer triangular pulse with 0 at symbol periods
Tp2 = 2*Ts;
tp2 = -2*Ts : dt : 2*Ts;
p2 = 1-abs(tp2)./(2*Ts);
p2 = p2 - p2 .* double(abs(abs(tp2)-Ts) <= 1E-5);
P2 = fft(p2,L);

figure;
subplot(2,1,1), plot(tp2,p2);
title('p3(t): longer triangular pulse with 0 at symbol periods');
xlabel('time(s)'); ylabel('p2(t)');
subplot(2,1,2), plot(f,fftshift(abs(P2/L)));
title('Fourier transform of p2(t)');
xlabel('frequency(Hz)'); ylabel('|P2|');

pulse_analysis(p2, dt, Tp2, Ts, xn1, 20, 0.5);
error_rate(p2, dt, Tp2, Ts, xn2, 20, 1);

%% Test  Nyquist filtering Criterion

%  Two graphs are plotted, one with each value at integer multiple of Ts
%  and the other is the the sum of values at integer multiple of Ts of the
%  pulse in the time domain. 

%  Test first pulse 
num_of_sample = int32(length(p1)/interval )+1;  
sum1 = zeros(1, num_of_sample); % initialize the container --11
mid_index = 1+ idivide(int32(num_of_sample), 2, 'floor'); 
each_P1 = zeros(1, num_of_sample); 

each_P1(1) = abs(p1(1)); 
sum1(1)= abs(p1(1));

for k = 2: num_of_sample
     if k == mid_index
         sum1(k) = sum1(k-1); 
        continue; 
     end
     each_P1(k) = abs(p1(1+(k-1)* interval)); 
    sum1(k) = sum1(k-1)+abs(p1(1+(k-1)* interval)); 
end
figure; 
hold on
indexing =  -5: 1 : 5; % indexing vector
plot(indexing, each_P1); 
plot(indexing, sum1); 
legend('Value at integer multiple of Ts', 'Sum of the values'); 
title('Sum of p1(t)  at integer multiple of Ts');
xlabel('Integer multiple of Ts '); ylabel('|P1|');

%  Test Second pulse 
midpoint = idivide(int32(length(p2)), 2, 'ceil'); 
half_length = int32(length(p2(1: midpoint))-1); 
half_num_sample = idivide(half_length, interval, 'floor'); 
num_of_sample = 1+2* half_num_sample; 
first_index = midpoint - interval * half_num_sample; 
mid_index =1+ idivide(int32(num_of_sample), 2, 'floor'); 

sum2 = zeros(1, num_of_sample); % initialize the container
each_P = zeros(1, num_of_sample); % initialize the container
sum2(1)= abs(p2(first_index));
each_P(1) = abs(p2(first_index)); 

for k = 2: num_of_sample
    if k == mid_index
        sum2(k) = sum2(k-1); 
        continue; 
    end
    each_P(k) = abs(p2(first_index + (k-1)* interval)); 
    sum2(k) = sum2(k-1)+abs(p2(first_index + (k-1)* interval)); 
end

indexing = -half_num_sample : half_num_sample; 
figure; 
hold on
plot(indexing, each_P); 
plot(indexing, sum2); 
title('p2(t)  at integer multiple of Ts');
legend('Value at integer multiple of Ts', 'Sum of the values'); 
xlabel('Integer multiple of Ts '); ylabel('|p2|');

%% Building and testing the PAM system using three signal inputs

N = 101; % choose a large N for accuracy
xn1 = zeros(3, N); 

for i = 1: 3 % Populate xn
    xn1(i, : ) = 2*((rand(1,N)>0.5)-0.5); % binary message
end

% Using functions to go through the PAM process and test the error
[tImp, r, y, y_total, y_up, y_rec, xn_est]= pam(p2, xn1, dt, 2*Ts, Ts, [20, 30, 40], 0.1, 3); 
graphing(tImp, fs, y, y_up, y_total, r, y_rec, xn1, xn_est, 3); 
error_rate(p2, dt, 2*Ts, Ts, xn1, [20, 30, 40], 3)

%% Testing Other Forms of Messages

% generate the text message and preprocess
message = 'case study 2 finished';  
binary = str2num(reshape(dec2bin(message)',1,[])'); 
binary = -1 * double(binary == 0) + binary;

% copy the same message 3 times
N = length(binary);
xn1 = zeros(3, N);
for i = 1:3
    xn1(i, : ) = binary'; 
end

% with sigma = 1
[~, ~, ~, ~, ~, ~, xn_est]= pam(p1, xn1, dt, 5*Ts, Ts, [20, 30, 40], 1, 3); 

% process the output and display
for i = 1:3
    answer = xn_est(i, : )'; 
    answer(answer <0) = 0; 
    
    a = num2str(reshape(answer,7,[]))'; 
    messageOut = char(bin2dec(a))'; 
    disp(messageOut);
end