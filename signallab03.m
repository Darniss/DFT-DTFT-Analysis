
%                           DSP 
% REG : AMENP2WNA21003
% Name : Darniss 





%% Generate and plot a complex exponential discrete-time signal 𝑥(𝑛) = 𝑒(−𝛼+𝑗𝜔)𝑛[𝑢(𝑛) − 𝑢(𝑛 − 25)]
% where α = 0.1, ω = π/6. Plot the absolute and angle parts and the real and imaginary parts 
% of x(n) against n. plot all the four plots as 2x2 plots with the title exactly as given 
% above and appropriate x and y labels.


close all;
clc;
clear all;
n=-20:20;
uu =  heaviside(n)-heaviside(n-25);
angl= pi/6;
alfa=-0.1+angl*j;
x_n=exp(alfa*n).*uu;

sgtitle('Complex Exponential discrete-time signal') 
subplot(2,2,1);
stem(n,real(x_n));
xlabel('Time Samples');
ylabel('Real Part');
title('Real part');

subplot(2,2,2);
stem(n,imag(x_n));
xlabel('Time Samples');
ylabel('Imaginary Part');
title('Imaginary part');

subplot(2,2,3);
stem(n,abs(x_n));
xlabel('Time Samples');
ylabel('Magnitude');
title('Absolute part');

subplot(2,2,4);
stem(n,angle(x_n));
xlabel('Time Samples');
ylabel('Phase');
title('Angle part');


%%  Generate and plot a continuous time sinusoidal signal 𝑦(𝑡) = 𝐴 sin (2𝜋𝐹𝑡 + 𝜃)
% where F is the frequency in Hz, θ is the phase angle in rad/s and t is the time in s and A
% is the maximum amplitude. If F = 10 Hz, θ = 30° and A = 3, generate y(t) for 5 cycles,
% starting with t = 0. Choose suitable increment in time to get a smooth plot. Give a 
% title exactly as given in the question with the variable values given in the problem. 
% Give appropriate x and y labels.

close all;
clc;
clear all;

F=10; %frequency [Hz]
t=(0:1/(F*100):1);
A=3;    %amplitude 
phi=30;  %phase
y=A*sin(2*pi*F*t+phi);
plot(t,y)
title('Sinusoidal wave');
xlabel('time(s)')
ylabel('amplitude(V)')




%% Sample the above signal with a sampling frequency 10 times of the maximum frequency 
% content of the above signal. Let the number of samples be the number of samples of 5 
% cycles of the above signal y(t). Generate and stem plot y(n).Let the title be 
% 𝑦(𝑛) = 𝐴 sin (𝜔𝑛 + 𝜃) with the actual values of ω and θ and appropriate x and y labels.

close all;
clc;
clear all;

F=10; %frequency [Hz]
t=(0:1/(F*100):1);
A=3;    %amplitude 
phi=30;  %phase
y=A*sin(2*pi*F*t+phi);
subplot(2,1,1);
plot(t,y)
title('Sinusoidal wave');
xlabel('time(s)')
ylabel('amplitude(V)')

sampleRate = 100;
samplePeriod = 1/sampleRate;
nT = 0:samplePeriod:1;
ys= A*sin(2*pi*F*nT+phi);

subplot(2,1,2);
stem(nT,ys)
title('Sampled signal of Sinusoidal wave');
xlabel('time(s)')
ylabel('amplitude(V)')

%% A random noise corrupted signal is generated by the Matlab code shown below:
% clear
% clc
% N = 50;
% n = 0:N-1;
% d = (rand(1,N)-0.5)*1.2; s = 2*n.*(0.9).^n;
% x = s+d;
% The noise can be reduced considerably by a three-point moving average filter which can be expressed as
% 𝑦(𝑛) = 1 [𝑥(𝑛 − 1) + 𝑥(𝑛) + 𝑥(𝑛 − 1)]. 3
% Write code for the above filter and remove the noise. Plot (using stem plot) the noisy signal and 
% the smoothed signal as 2x1 subplots. Let the range of the output plot be the same as that of input.


close all;
clc;
clear all;

N = 50;
n = 0:N-1;
d = (rand(1,N)-0.5)*1.2; s = 2*n.*(0.9).^n;
x = s+d;
% random noise
% uncorrupted signal
% random noise corrupted signal

% y= 1/3.*(x(n-1)+x(n)+x(n-1));

num= (1/3)*[1 1 1];
den= [1];
y= filter(num,den,x);

sgtitle('3 Point Moving Filter') 
subplot(2,1,1);
plot(x);
xlabel('N');
ylabel('x(n)');
title('Noisy Signal');

subplot(2,1,2);
stem(y);
xlabel('N');
ylabel('y(n)');
title('Smoothed Signal');




%% Let the impulse response h(n) of a causal FIR system be h(𝑛) = {3,2,1, −2,1,0, −4,0,3}. 
% To this system a causal sequence 𝑥(𝑛) = {1, −2,3, −4,3,2,1} is given. Compute the output
% using the convolution function as well as filter function. Plot the output by both methods
% and see if they produce the same result. Give appropriate titles and labels.

close all;
clc;
clear all;

h=[3,2,1, -2,1,0, -4,0,3];
x=[1, -2,3, -4,3,2,1];

filter_op= filter(h,1,x);
convolution_op= conv(h,x);
yf2 = filter(h,1,[x,0,0,0,0,0,0,0,0]);
disp(convolution_op);
disp(filter_op)
disp(yf2)
figure
sgtitle('Convolution Output') 
subplot(2,1,1);
stem(yf2);
title('Using Filter Function');

subplot(2,1,2);
stem(convolution_op);
title('Using Convolution Function');