%% Case Study 2: PAM Communication

Ts = .1; % symbol time
dt = Ts/50; % sample period
t = -5*Ts : dt : 5*Ts; % time vector
L = 1024; % length of the signal in frequency domain
interval = Ts/dt; 

fs = 1/dt; % length in frequency domain 
f = -fs/2 : fs/L : fs/2-fs/L; % frequency range 

% 1st pulse shape: truncated sinc
p1 = sinc(t/Ts) .* double(abs(t) <= 3*Ts); % this truncating time can be varied
P1 = fft(p1,L);

figure;
subplot(2,1,1), plot(t,p1);
title('p1(t): truncated sinc');
xlabel('time(s)'); ylabel('p1(t)');
subplot(2,1,2), plot(f,fftshift(abs(P1/L)));
title('Fourier transform of p1(t)');
xlabel('frequency(Hz)'); ylabel('|P1|');

% 2rd pulse shape: raised cosine filter
p2 = rcosdesign(0, 10 , Ts/dt);
t3 = - length(p2)/2 * dt : dt: length(p2)/2 * dt-dt; 
P3 = fft(p2,L);

figure;
subplot(2,1,1), plot(t3,p2);
title('p2(t): raised cosine filter');
xlabel('time(s)'); ylabel('p2(t)');
subplot(2,1,2), plot(f,fftshift(abs(P3/L)));
title('Fourier transform of p2(t)');
xlabel('frequency(Hz)'); ylabel('|P2|');

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
xn = zeros(3, N); 

for i = 1: 3 % Populate xn
    xn(i, : ) = 2*((rand(1,N)>0.5)-0.5); % binary message
end

% Using functions to go through the PAM process and test the error
[tImp, r, y, y_total, y_up, y_rec, xn_est]= pam(p1, xn, dt, 5*Ts, Ts, [20, 30, 40], 0.1, 3); 
graphing(tImp, fs, y, y_up, y_total, r, y_rec, xn, xn_est, 3); 
error_rate(p1, dt, 5*Ts, Ts, xn, [20, 30, 40], 3)

%% Testing Other Forms of Messages

message = 'I love ESE 232';  
binary = str2num(reshape(dec2bin(message)',1,[])'); 

N = length(binary);
xn = zeros(3, N);
xn(1, : ) = binary'; 

[tImp, r, y, y_total, y_up, y_rec, xn_est]= pam(p1, xn, dt, 5*Ts, Ts, [20, 30, 40], 0.1, 1); 
% graphing(tImp, fs, y, y_up, y_total, r, y_rec, xn, xn_est, 3); 
% error_rate(p1, dt, 5*Ts, Ts, xn, [20, 30, 40], 1)

answer = xn_est(1, : )'; 
answer(answer <0) = 0; 

a = num2str(reshape(answer,7,[]))'; 
messageOut = char(bin2dec(a))'; 
disp(messageOut);