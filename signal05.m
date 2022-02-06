%                           DSP
% REG : AMENP2WNA21003
% Name : Darniss 

%% 1. Obtain the DTFT and DFT for the causal signal 𝑥(𝑛) = {3,2,1,0,1,2}
% and plot the absolute and angle values in the frequency range 0 ≤ 𝜔 ≤ 𝜋 as
% 2x2 subplots with appropriate titles and labels.

close all;
clc;
clear all;


x = [3 2 1 0 1 2];  % input sequence
N = length(x);
X = zeros(6,1);
w=0:0.01:pi;
% DFT
for k = 0:N-1
    for n = 0:N-1
        X(k+1) = X(k+1) + x(n+1)*exp(-j*pi/2*n*k);  % dft algorithm
    end
end
t = 0:N-1;
% plot
subplot(221)
stem(t,X)
xlabel('Samples');
ylabel('|X(k)|');
title('DFT - Magnitude response')

subplot(222)
stem(t,angle(X))
xlabel('Samples');
ylabel('Phase');
title('DFT - Phase response')

clear
% DTFT
disp("yes")
x = [3 2 1 0 1 2];  % input sequence
n = 0:5;
w=0:0.01:pi;

for i=1:length(w);
    X(i)=0;
    for k=1:length(x);
          X(i)=X(i)+x(k).*exp(-j.*w(i).*n(k));
      end
end

% plot
subplot(223)
plot(w,X)
xlabel('Frequency');
ylabel('|X(k)|');
title('DTFT - Magnitude response')

subplot(224)
plot(w,angle(X))
xlabel('Frequency');
ylabel('Phase');
title('DTFT - Phase response')





%% 2. Repeat the above by plotting the absolute and angle values in the frequency
% range 0 ≤ 𝜔 ≤ 𝜋 of both the DTFT and DFT in the same plot as 2x1 subplots
% with appropriate titles, labels and legends.


x = [3 2 1 0 1 2];

n1 = 0:5;
y1 = x(1:6);
Y1 = dft(y1,length(x));
magY1 = abs(Y1);
angY1= angle(Y1);

k1 = 0:1:length(x)-1;
N = length(x);
w1 = (2*pi/N)*k1;%linspace(-pi,pi,N); % (2*pi/N)*k1;
w0= (2*pi/N)*k1;
disp(w0)
disp(w1)

%Discrete-time Fourier Transform
K = 500;
k = 0:1:K;
w = linspace(-pi,pi,500);%2*pi*k/K; %plot DTFT in [0,2pi];
X = y1*exp(-j*n1'*w);

magX = abs(X);
angX= angle(X);


% Plot
subplot(2,1,2)
stem(n1,y1);
stem(w1/pi,angY1);
title('Phase Spectra');
xlabel('Frequency in pi units');
ylabel('Angle');
hold on 
plot(w/pi,angX);
hold off
xlim([0 1])

disp(length(w1))
subplot(2,1,1);
stem(w1/pi,magY1);
title('Magnitude Spectra');
xlabel('Frequency in pi units');
ylabel('Magnitude');
hold on 
plot(w/pi,magX);
hold off
xlim([0 1])




%% 3. Artificially increase the length of x(n) to 64 and repeat the same of
% question 1. What do you infer from the plots?

close all;
clc;
clear all;

inp= randi([1 10],64);
disp(length(inp))

n = 0:64;
x =[inp];

n1 = 0:63;
y1 = x(1:64);
Y1 = dft(y1,length(x));
magY1 = abs(Y1);
angY1= angle(Y1);

k1 = 0:1:length(x)-1;%0:1:length(x)-1;
N = length(x);
w1 = (2*pi/N)*k1;%linspace(-pi,pi,N); % (2*pi/N)*k1;
w0= (2*pi/N)*k1;
disp(w0)
disp(w1)

%Discrete-time Fourier Transform
K = 500;
k = 0:1:K;
w = linspace(-pi,pi,500);%2*pi*k/K; %plot DTFT in [0,2pi];
X = y1*exp(-j*n1'*w);

magX = abs(X);
angX= angle(X);


% Plot
subplot(2,1,2)
stem(n1,y1);
stem(w1/pi,angY1);
title('Phase Spectra');
xlabel('Frequency in pi units');
ylabel('Angle');
hold on 
plot(w/pi,angX);
hold off
xlim([0 1])

disp(length(w1))
subplot(2,1,1);
stem(w1/pi,magY1);
title('Magnitude Spectra');
xlabel('Frequency in pi units');
ylabel('Magnitude');
hold on 
plot(w/pi,magX);
hold off
xlim([0 1])

disp("Comments:")
disp("Increasing the sampling rate increases the frequency of the waveform")


%% 4 Obtain the DTFT and DFT for the causal signal
% 𝑥(𝑛)=cos(0.48𝜋𝑛)+sin(0.52𝜋𝑛+𝜋), 0≤𝑛≤15 6
% and plot the absolute and angle values in the frequency range 0 ≤ 𝜔 ≤ 𝜋 in
% the same plot as 2x1 subplots with appropriate titles, labels and legends.

close all;
clc;
clear all;


n = 0:15;
x = cos(0.48*pi*n) + sin(0.52*pi*n+ pi/6);

n1 = 0:15;
y1 = x(1:16);
Y1 = dft(y1,length(x));
magY1 = abs(Y1);
angY1= angle(Y1);

k1 = 0:1:length(x)-1;%0:1:length(x)-1;
N = length(x);
w1 = (2*pi/N)*k1;%linspace(-pi,pi,N); % (2*pi/N)*k1;
w0= (2*pi/N)*k1;
disp(w0)
disp(w1)

%Discrete-time Fourier Transform
K = 500;
k = 0:1:K;
w = linspace(-pi,pi,500);%2*pi*k/K; %plot DTFT in [0,2pi];
X = y1*exp(-j*n1'*w);

magX = abs(X);
angX= angle(X);


% Plot
subplot(2,1,2)
stem(n1,y1);
stem(w1/pi,angY1);
title('Phase Spectra');
xlabel('Frequency in pi units');
ylabel('Angle');
hold on 
plot(w/pi,angX);
hold off
xlim([0 1])

disp(length(w1))
subplot(2,1,1);
stem(w1/pi,magY1);
title('Magnitude Spectra');
xlabel('Frequency in pi units');
ylabel('Magnitude');
hold on 
plot(w/pi,magX);
hold off
xlim([0 1])







%% 5. Repeat the above for samples 0 ≤ 𝑛 ≤ 127 and repeat the same of question 4.
% What do you infer from the plots?
close all;
clc;
clear all;



n = 0:127;
x = cos(0.48*pi*n) + sin(0.52*pi*n+ pi/6);

n1 = 0:127;
y1 = x(1:128);
Y1 = dft(y1,length(x));
magY1 = abs(Y1);
angY1= angle(Y1);

k1 = 0:1:length(x)-1;%0:1:length(x)-1;
N = length(x);
w1 = (2*pi/N)*k1;%linspace(-pi,pi,N); % (2*pi/N)*k1;
w0= (2*pi/N)*k1;
disp(w0)
disp(w1)

%Discrete-time Fourier Transform
K = 500;
k = 0:1:K;
w = linspace(-pi,pi,500);%2*pi*k/K; %plot DTFT in [0,2pi];
X = y1*exp(-j*n1'*w);

magX = abs(X);
angX= angle(X);


% Plot
subplot(2,1,2)
stem(n1,y1);
stem(w1/pi,angY1);
title('Phase Spectra');
xlabel('Frequency in pi units');
ylabel('Angle');
hold on 
plot(w/pi,angX);
hold off
xlim([0 1])

disp(length(w1))
subplot(2,1,1);
stem(w1/pi,magY1);
title('Magnitude Spectra');
xlabel('Frequency in pi units');
ylabel('Magnitude');
hold on 
plot(w/pi,magX);
hold off
xlim([0 1])


% comments
disp("Inference")
disp("The more samples that are taken, the more detail about where the waves rise and fall is recorded")


%% 6.  The first six values of the 10-point DFT of a real-valued sequence x(n) 
% are given as 𝑋(𝑘) = {10, −2 + 𝑗3,3 + 𝑗4,2 − 𝑗3,4 + 𝑗5,12}
% Compute x(n).

close all;
clc;
clear all;

N=10;
xk= [10 -2+3j 3+4j 2-3j 4+5j 12];
n= length(xk);
xk1=[];
for k=1:N-n
    xk1(k)=conj(xk(k));
    
end

XK=[xk flip(xk1)];
% disp(XK)
Y = idft(XK,N);

subplot(211)
stem(abs(XK))
title('DFT');
disp('DFT of input sequence is ');
disp(abs(XK))

subplot(212)
stem(abs(Y))
%stem(fftshift(Y))
title('IDFT');
disp('IDFT of input sequence is ');
disp(abs(Y))







%% 7. From the DFT of x(n) in question 6, and using the DFT properties, obtain the
% DFTs of the following sequences and verify the results:

close all;
clc;
clear all;

x = [4.2297 1.2700 1.7064 1.6798 3.0347 1.5264 0.8064 0.6967 3.1605 0.8943];
disp('Input sequence')
disp(x)
N = length(x); %data length
n=1:10;
W = exp(-j*pi*2/N);
X = fft(x);

disp("Welcome to Properties of DFT")
disp("1) x1(n)=x((n−2))")
disp("2)  x2(n) = x((n + 5))")
disp("3) x3(n)=x((2−n))")
disp("4) x4(n)=x((−n−3))");
disp("5) x5(n)=x(n)∙x((−n))")
disp("6) x6(n)=x(n)e4𝜋𝑛/5")

option=input('Enter Your Option : ');

% disp(option);

switch option
    case 1
        disp('a) x1(n)=x((n−2)) Time Shifting property')
        k = 0:N-1; % frequency index
        m = -2; % number of shift
        Y= W.^(-m*k) .* X;
        y= ifft(Y);
        subplot(211)
        stem(x)
        title('Original Signal');
        subplot(212)
        stem(y)
        title('Shifted Signal');
    case 2
        disp('b) x2(n) = x((n + 5)) Time Shifting property')
        k = 0:N-1; % frequency index
        m = 5; % number of shift
        Y= W.^(-m*k) .* X;
        y= ifft(Y);
        subplot(211)
        stem(x)
        title('Original Signal');
        subplot(212)
        stem(y)
        title('Shifted Signal');
    case 3
        disp('c) x3(n)=x((2−n)) Time Reversal property wrong')
        y2= [X(1) flip(X(2:end))];
        k = 0:N-1; % frequency index
        m = 2; % number of shift
        X2=fft(y2);
        Y= W.^(-m*k) .* y2;
        y= ifft(Y);        
        subplot(211)
        stem(x)
        title('Original Signal');
        subplot(212)
        stem(n,y)
        title('Time reversed and shifted Signal');        
    case 4
        disp('d) x4(n)=x((−n−3)) Time Reversal property wrong')
        y2=[X(1) flip(X(2:end))];
        k = 0:N-1; % frequency index
        m = -3; % number of shift
        X2=fft(y2);
        Y= W.^(-m*k) .* y2;
        y= ifft(Y);        
        subplot(211)
        stem(x)
        title('Original Signal');
        subplot(212)
        stem(n,y)
        title('Time reversed and shifted Signal');  
    case 5
        disp('e) x5(n)=x(n)∙x((−n)) Multiplication property')
        rev= [X(1) flip(X(2:end))];
        rev0=ifft(rev);
        con= cconv(rev,X,N)/N;
        an=ifft(con);
        subplot(211)
        stem(x)
        title('Original Signal');
        subplot(212)
        stem(n,an)
        title('Multiplied Signal');         

    case 6
        disp('f) x6(n)=x(n)e4𝜋𝑛/5 Frequency Shifting property')
        k = 0:N-1; % frequency index
        m=4;
        y=W.^(m*k) .* X;
        y2=ifft(y);
        disp(y2)
        subplot(211)
        stem(x)
        title('Original Signal');
        subplot(212)
        stem(n,y2)
        title('Frequency shifted Signal');       

    otherwise
        warning('Please Enter the correct number')
end





%% 8. From the DFT of x(n) in question 6, and using the DFT properties, obtain the following sequences:
% 𝑎) 𝑦1(𝑛)=10−𝑝𝑜𝑖𝑛𝑡𝑐𝑜𝑛𝑣𝑜𝑙𝑢𝑡𝑖𝑜𝑛𝑜𝑓𝑜𝑓𝑥(𝑛)𝑤𝑖𝑡h𝑥((−𝑛))10.
% 𝑏) 𝑦2(𝑛)=𝑙𝑖𝑛𝑒𝑎𝑟𝑐𝑜𝑛𝑣𝑜𝑙𝑢𝑡𝑖𝑜𝑛𝑜𝑓𝑥(𝑛)𝑤𝑖𝑡h𝑥((−𝑛))10.
% Plot the linear convolution obtained by the above method and obtained using the conv function.

close all;
clc;
clear all;

x = [4.2297 1.2700 1.7064 1.6798 3.0347 1.5264 0.8064 0.6967 3.1605 0.8943];
disp('Input sequence')
disp(x)
N = length(x); %data length
n=1:10;
W = exp(-j*pi*2/N);
X = fft(x);
rev= [X(1) flip(X(2:end))];
rev0=ifft(rev);

disp("Question 8")
option=input('Enter Your Option (1 or 2) : ');
%  disp(option);

switch option
    case 1
        circular= cconv(x,rev0,N);        
        mul=rev.*X;
        disp(ifft(mul))
        disp(circular)
        subplot(211)
        stem(circular)
        title('Circularly convoluted Signal');
        subplot(212)
        stem(ifft(mul))
        title('Multiplication of DFT of x(n) and x(-n)');          
        

    case 2
        
        linear= conv(x,rev0,"same");
        mul=rev.*X;

        disp(ifft(mul))
        disp(linear)   
        subplot(211)
        stem(linear)
        title('Linearly convoluted Signal');
        subplot(212)
        stem(ifft(mul))
        title('Multiplication of DFT of x(n) and x(-n)');    

    otherwise
        warning('Please Enter the correct number')
end





%% 9.  An analog signal 𝑥(𝑡) = 4 + 2 cos (150𝜋𝑡 + 𝜋) + 4𝑠𝑖𝑛(350𝜋𝑡) is sampled by a 3
% sampling frequency 𝑓 = 400 𝑠𝑎𝑚𝑝𝑙𝑒𝑠/𝑠 to obtain the discrete-time signal x(n). 𝑠
% Compute the signal x(n) and its corresponding DFT 𝑋(𝑘). Plot the magnitude
% of DFT in terms of π units as well as actual frequency (Hz) using stem plot.

close all;
clc;
clear all;


Fs = 400;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
Fn=Fs/2;
t=0:1/Fs:5;%(0:L-1)*T; % Time vector


x = 4+2*cos(150*pi*t + pi/3)+ 4*sin(350*pi*t);
L = length(x); % Length of signal
disp(L)
tim = linspace(0, L, L)*T; 
Y = fft(x);

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

% f = Fs*(0:(L/2))/L;
% disp(length(f))
% w= linspace(-pi,pi,L/2+.5);
Fvrs = linspace(0, 1, fix(L/2)+1)*pi; 
FvHz = linspace(0, 1, fix(L/2)+1)*Fn; 

figure(1)
plot(tim,x(1: length(tim)))
xlim([0 5])

figure(2)
subplot(211)
plot(FvHz,abs(P1(1:length(FvHz)))*2) 
title('DFT')
xlabel('f (Hz)')
ylabel('Magnitude')


subplot(212)
plot(Fvrs,abs(P1(1:length(Fvrs)))*2) 
title('DFT')
xlabel('f in pi units')
ylabel('Magnitude')








%% 10. Repeat the above for a sampling frequency 𝑓 = 200 𝑠𝑎𝑚𝑝𝑙𝑒𝑠/𝑠. What do you infer 𝑠
% from the plots of this and the previous question’s (question 9) plot? Comment on this.
 
close all;
clc;
clear all;


Fs = 200;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
Fn=Fs/2;
t=0:1/Fs:5;%(0:L-1)*T; % Time vector


x = 4+2*cos(150*pi*t + pi/3)+ 4*sin(350*pi*t);
L = length(x); % Length of signal
disp(L)
tim = linspace(0, L, L)*T; 
Y = fft(x);

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

% f = Fs*(0:(L/2))/L;
% disp(length(f))
% w= linspace(-pi,pi,L/2+.5);
Fvrs = linspace(0, 1, fix(L/2)+1)*pi; 
FvHz = linspace(0, 1, fix(L/2)+1)*Fn; 


figure(1)
plot(tim,x(1: length(tim)))
xlim([0 5])

figure(2)
subplot(211)
plot(FvHz,abs(P1(1:length(FvHz)))*2) 
title('DFT')
xlabel('f (Hz)')
ylabel('Magnitude')


subplot(212)
plot(Fvrs,abs(P1(1:length(Fvrs)))*2) 
title('DFT')
xlabel('f in pi units')
ylabel('Magnitude')


% comments
disp("Comments:")
disp("As the sampling frequency decreases, the signal separation also decreases.")








