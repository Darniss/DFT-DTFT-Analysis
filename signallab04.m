%                           DSP
% REG : AMENP2WNA21003
% Name : Darniss 




%% 1.Obtain the DTFT of the sequence h(𝑛) = {1, 2, 3, 4, 5}, −1 ≤ 𝑛 ≤ 3. 
% Try various methods to compute the DTFT in MATLAB.
close all;
clc;
clear all;

w=-pi:0.01:pi;
% n=-1:3;
% x=[1 2 3 4 5];
x = [3 2 1 0 1 2];  % input sequence
n = 0:5;
for i=1:length(w);
    X(i)=0;
    for k=1:length(x);
          X(i)=X(i)+x(k).*exp(-j.*w(i).*n(k));
      end
end
[Y wu]=dtft_m(x,n);
sgtitle('DTFT');
subplot(2,1,1);
plot(w,X);
xlabel('radian/pi')
ylabel("magnitude of DTFT");
title('DTFT without mydsp tool box');
subplot(2,1,2);
plot(wu,Y);
xlabel('radian/pi')
ylabel("magnitude of DTFT");
title('DTFT with mydsp tool box');


%% Use the dtft_m function to compute the DTFT of the following finite-duration 
% sequences over the range −𝜋 ≤ 𝜔 ≤ 𝜋. Plot DTFT magnitude and angle as a function 
% of ω and plot the real and imaginary parts as a function of ω in one figure window 
% as 2x2 subplots. Comment on the plots.

% a) 𝑥(𝑛) = 𝑛(0.9)𝑛[𝑢(𝑛 + 10) − 𝑢(𝑛 − 21)].
close all;
clc;
clear all;

n=-3:3;
x=((0.9).^n).*(heaviside(n+10)-heaviside(n-21));
disp(x)
 [Y wu]=dtft_m(x,n);
 sgtitle('DTFT');
 subplot(2,1,1);
 plot(wu,abs(Y));
 xlabel('radian/pi')
 ylabel("magnitude of DTFT");
 title('DTFT');
 subplot(2,1,2);
 plot(wu,angle(Y));
xlabel('radian/pi')
ylabel("angle of DTFT");
title('DTFT');

%% b) 𝑥(𝑛) = (0.6)𝑛[𝑢(𝑛 + 10) − 𝑢(𝑛 − 11)].
close all;
clc;
clear all;

n=-10:10;
x=((0.6).^n).*(heaviside(n+10)-heaviside(n-11));
disp(x)
 [Y wu]=dtft_m(x,n);
 sgtitle('DTFT');
 subplot(2,1,1);
 plot(wu,abs(Y));
 xlabel('radian/pi')
 ylabel("magnitude of DTFT");
 title('DTFT');
 subplot(2,1,2);
 plot(wu,angle(Y));
xlabel('radian/pi')
ylabel("angle of DTFT");
title('DTFT');


%% c
close all;
clc;
clear all;

n=-10:10;
for i=1:length(n)
if (n(i)==0)
    x=0;
elseif(n(i)<=-10 || n(i)<=-1)
    x=-(2).^n;
elseif(n(i)<=1 || n(i)<= 10)
   x= (1/2).^n;
end
end

 [Y wu]=dtft_m(x,n);
 sgtitle('DTFT');
 subplot(2,1,1);
 plot(wu,abs(Y));
 xlabel('radian/pi')
 ylabel("magnitude of DTFT");
 title('DTFT');
 subplot(2,1,2);
 plot(wu,angle(Y));
xlabel('radian/pi')
ylabel("angle of DTFT");
title('DTFT');

%% d
close all;
clc;
clear all;

n=-3:3;
x=(cos(0.5*pi*n)+ j*sin(0.5*pi*n)).*(heaviside(n+10)-heaviside(n-51));
 [Y wu]=dtft_m(x,n);
 sgtitle('DTFT');
 subplot(2,1,1);
 plot(wu,abs(Y));
 xlabel('radian/pi')
 ylabel("magnitude of DTFT");
 title('DTFT');
 subplot(2,1,2);
 plot(wu,angle(Y));
xlabel('radian/pi')
ylabel("angle of DTFT");
title('DTFT');


%% e
close all;
clc;
clear all;

n=0:7;
x=[4 3 2 1 1 2 3 4];
disp(x)
 [Y wu]=dtft_m(x,n);
 sgtitle('DTFT');
 subplot(2,1,1);
 plot(wu,abs(Y));
 xlabel('radian/pi')
 ylabel("magnitude of DTFT");
 title('DTFT');
 subplot(2,1,2);
 plot(wu,angle(Y));
xlabel('radian/pi')
ylabel("angle of DTFT");
title('DTFT');


%% f
close all;
clc;
clear all;

n=0:7;
x=[4 3 2 1 -1 -2 -3 -4];
disp(x)
 [Y wu]=dtft_m(x,n);
 sgtitle('DTFT');
 subplot(2,1,1);
 plot(wu,abs(Y));
 xlabel('radian/pi')
 ylabel("magnitude of DTFT");
 title('DTFT');
 subplot(2,1,2);
 plot(wu,angle(Y));
xlabel('radian/pi')
ylabel("angle of DTFT");
title('DTFT');


%% 3
close all;
clc;
clear all;

n=0:7;
x1=[1 2 2 1];
disp(x1(2))
for i=1:length(n)
if (n(i)<=0 || n(i)<=3)
   x2=x1(i);
elseif(n(i)<=4 || n(i)<=7)
    x2=x1(n(i)-4);
else
   x2=0;
end
end
 [Y wu]=dtft_m(x2,n);
 sgtitle('DTFT');
 subplot(2,1,1);
 plot(wu,abs(Y));
 xlabel('radian/pi')
 ylabel("magnitude of DTFT");
 title('DTFT');
 subplot(2,1,2);
 plot(wu,angle(Y));
xlabel('radian/pi')
ylabel("angle of DTFT");
title('DTFT');


%% 4

close all;
clc;
clear all;

w=-4*pi:0.01:4*pi;
n=-1:3;
x=[-1 4 -3 2 5];
for i=1:length(w);
    X(i)=0;
    for k=1:length(x);
          X(i)=X(i)+x(k).*exp(-j.*w(i).*n(k));
      end
end
[Y wu]=dtft_m(x,n);
sgtitle('DTFT');
subplot(2,2,1);
plot(wu,abs(Y));
xlabel('radian/pi')
ylabel("magnitude");
title('Magnitude part');


subplot(2,2,2);
plot(wu,angle(Y));
xlabel('radian/pi')
ylabel("phase");
title('Phase part');

subplot(2,2,3);
plot(wu,real(Y));
xlabel('radian/pi')
ylabel("real part");
title('Real part');


subplot(2,2,4);
plot(wu,imag(Y));
xlabel('radian/pi')
ylabel("imaginary part");
title('Imaginary Part');


%% 5) Also, obtain the DTFTs of the following discrete-time signals and plot 
% the real part and imaginary parts against ω in the range -π to π. Observe the 
% plots and make your comments:


% i)
close all;
clc;
clear all;

n=-5:5;
x=[6 -4 5 -3 -2 7 -2 -3 5 -4 6];
disp(x)
 [Y wu]=dtft_m(x,n);
 sgtitle('DTFT');
 subplot(2,1,1);
 plot(wu,real(Y));
 xlabel('radian/pi')
 ylabel("real >");
 title('Real part of DTFT');
 subplot(2,1,2);
 plot(wu,imag(Y));
xlabel('radian/pi')
ylabel("imagninary >");
title('Imaginary part of DTFT');

%% ii) 

close all;
clc;
clear all;

n=-5:5;
x=[6 -4 5 -3 -2 0 2 3 -5 4 -6];
disp(x)
 [Y wu]=dtft_m(x,n);
 sgtitle('DTFT');
 subplot(2,1,1);
 plot(wu,real(Y));
 xlabel('radian/pi')
 ylabel("real >");
 title('Real part of DTFT');
 subplot(2,1,2);
 plot(wu,imag(Y));
xlabel('radian/pi')
ylabel("imagninary >");
title('Imaginary part of DTFT');


%% 6 
close all;
clc;
clear all;

n=-5:5;
x1= (1/2).^n.*heaviside(n);
x2= (1/3).^n.*heaviside(n);
y1=conv(x1,x2,"same");
 [Y wu]=dtft_m(y1,n);
  [X1 wu1]=dtft_m(x1,n);  % dft of x1
  [X2 wu2]=dtft_m(x2,n);  % dft of x2
  Y0=X1.*X2;

  figure(1)
 sgtitle('DTFT');
 subplot(2,2,1);
 plot(wu,abs(Y));
 xlabel('radian/pi')
 ylabel("magnitude of DTFT");
 title('DTFT');
 subplot(2,2,2);
 plot(wu,angle(Y));
xlabel('radian/pi')
ylabel("angle of DTFT");
title('DTFT');

 subplot(2,2,3);
 plot(wu,abs(Y0));
 xlabel('radian/pi')
 ylabel("magnitude of DTFT");
 title('DTFT');
 subplot(2,2,4);
 plot(wu,angle(Y0));
xlabel('radian/pi')
ylabel("angle of DTFT");
title('DTFT');

%% 7
close all;
clc;
clear all;

M = [10 101]; 
R1 = ones(1,M(1)); 
R2 = ones(1,M(2));
N1 = 0:M(1)-1;
N2 = 0:M(2)-1;
% disp(N1)
% Rf1 = fft(R1,N1);
% Rf2 = fft(R2,N2);
% disp(Rf1)
Rec_window1= rectwin(length(N1)); % For M1
hann_window1 = 0.5 + 0.5* cos (2 * pi* N1./ (M(1)-1)).*Rec_window1; 
hamm_window1 = 0.54 + 0.46* cos (2 * pi* N1./ (M(1)-1)).*Rec_window1; 

Rec_window2= rectwin(length(N2)); % For M2
hann_window2 = 0.5 + 0.5* cos (2 * pi* N2./ (M(2)-1)).*Rec_window2; 
hamm_window2 = 0.54 + 0.46* cos (2 * pi* N2./ (M(2)-1)).*Rec_window2; 
 [Y1 wu1]=dtft_m(Rec_window1,N1);
 [Y2 wu2]=dtft_m(hann_window1,N1);
 [Y3 wu3]=dtft_m(hamm_window1,N1);
 [Y4 wu4]=dtft_m(Rec_window2,N2);
 [Y5 wu5]=dtft_m(hann_window2,N2);
 [Y6 wu6]=dtft_m(hamm_window2,N2);

 % plot
sgtitle("DFT of the Window");
subplot(3,2,1);
plot(wu1,abs(Y1));
xlabel('radian/pi')
ylabel("magnitude of DTFT");
ylim([0 1])
title("Rectangular Window DFT (M=10)")
subplot(3,2,5);
plot(wu2,abs(Y2));
xlabel('radian/pi')
ylabel("magnitude of DTFT");
ylim([0 1])
title("Hanning Window DFT (M=10)")
subplot(3,2,3);
plot(wu3,abs(Y3));
xlabel('radian/pi')
ylabel("magnitude of DTFT");
ylim([0 1])
title("Hamming Window DFT (M=10)")
subplot(3,2,2);
plot(wu4,abs(Y4));
xlabel('radian/pi')
ylabel("magnitude of DTFT");
ylim([0 1])
title("Rectangular Window DFT (M=101)")
subplot(3,2,6);
plot(wu5,abs(Y5));
xlabel('radian/pi')
ylabel("magnitude of DTFT");
ylim([0 1])
title("Hanning Window DFT (M=101)")
subplot(3,2,4);
plot(wu6,abs(Y6));
xlabel('radian/pi')
ylabel("magnitude of DTFT");
ylim([0 1])
title("Hamming Window DFT (M=101)")

%% 8
close all;
clc;
clear all;

n=-1:3;
xa=((1/2).^(n-1)).*heaviside(n-1);
xb= ((1/2).^n).*heaviside(n+3);
xc= dd(n+2)+dd(n-2);

 [Y1 wu1]=dtft_m(xa,n);
 [Y2 wu2]=dtft_m(xb,n);
 [Y3 wu3]=dtft_m(xc,n);

sgtitle('DTFT');
 subplot(3,2,1);
 plot(wu1,abs((Y1)));
 xlabel('radian/pi')
 ylabel("magnitude of DTFT");
 xlim([0 pi])
 title('xa (magnitude)');
 subplot(3,2,2);
 plot(wu1,angle(Y1));
xlabel('radian/pi')
ylabel("angle of DTFT");
 xlim([0 pi])
title('xa (angle)');

 subplot(3,2,3);
 plot(wu2,abs(Y2));
 xlabel('radian/pi')
 ylabel("magnitude of DTFT");
 xlim([0 pi])
 title('xb (magnitude)');
 subplot(3,2,4);
 plot(wu2,angle(Y2));
xlabel('radian/pi')
ylabel("angle of DTFT");
 xlim([0 pi])
title('xb (angle)');

 subplot(3,2,5);
 plot(wu3,abs(Y3));
 xlabel('radian/pi')
 ylabel("magnitude of DTFT");
 xlim([0 pi])
 title('xc (magnitude)');
 subplot(3,2,6);
 plot(wu3,angle(Y3));
xlabel('radian/pi')
ylabel("angle of DTFT");
 xlim([0 pi])
title('xc (angle)');


%% 9


close all;
clc;
clear all;

% L=1000;
% dw=2*pi/L;
% w = 0:dw:pi-dw;
% 𝑦(𝑛) − 1/6 𝑦(𝑛 − 1) − 1/6 𝑦(𝑛 − 2) = 𝑥(𝑛)
b= [1 -1/6 -1/6];
a= [1];
[h,t]=freqz(a,b);  % Frequency Response
[hh,tt]=impz(a,b,10);  % Frequency Response
% filter(b,a,[1 zeros(1,n-1)])
figure(1)
sgtitle("Frequency Response using freqz()");
subplot(2,1,1), plot(t/pi,mag2db(abs(h)));
xlabel('Frequency [Hz]');
ylabel('Magnitude (bB)');
subplot(2,1,2), plot(t/pi,angle(h));
xlabel('Frequency [Hz]');
ylabel('Phase');

figure(2)
sgtitle("Impulse Response using impulz()");
subplot(2,1,1), stem(tt,(abs(hh)));
xlabel('Samples');
ylabel('Magnitude');
subplot(2,1,2), plot(tt,(hh));
xlabel('Time');
ylabel('Magnitude');

%% 10

close all;
clc;
clear all;

n=1:10-1;
xn= heaviside(n)-heaviside(n-10);
b= [1 -1/6 -1/6];
a= [1 0 0 0 0 0 0 0 0 0 1];

yn=filter(a,b,[1 zeros(1,length(n)-1)]);
figure(1)
sgtitle("Impulse Response from filter()");
subplot(2,1,1), stem(xn);
xlabel('Samples');
ylabel('Magnitude');
title("Input");
subplot(2,1,2), stem(yn);
xlabel('Sample');
ylabel('Magnitude');
title("Output");
