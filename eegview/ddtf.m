% Pr7_2.m
% Demo DTF based on signals generated with function Simulated_Signal,
% These signals are saved in File test.mat
% !! the function regres must be in the directory to detrend the data !!
% This program steps with 20 s windows (duration) through
% each of three 40 s signals (el1, el2, el3).
% These 20 second windows move ahead with steps of 5 s (increment).
% So there is 15 s overlap between the 20 s analysis windows)
% Within each window of 20 seconds, the average (cross)spectra
% are computed from fft-analysis epochs of 128 points (step).
% NOTE: THIS PROGRAM IS NOT OPTIMIZED.
clear;
load test
% Parameters
cmd0=[‘N=length(el1)’];
eval([cmd0 ‘;’])
sample_rate=400;
duration=20;
step=128;
increment=5;
dt=1/sample_rate;
fNyq=sample_rate/2;
df=1/(step*dt);
f=0:df:fNyq;
% load the data with 40 s input traces el1 - el3
% Determine the length of the signal
% 400 Hz sample rate
% duration of the total analysis window in
% seconds
% # of points in the FFT analysis window
% steps of the analysis window in seconds
% sample interval
% Nyquist frequency
% Frequency step for the FFT
% Frequency axis for the FFT
% Plot the three signals el1 - el3 in the top panels
figure
subplot(4,3,1);
plot(el1);hold;axis([0 N min(el1) max(el1)]);
t=[‘el1’];title(t);
axis(‘off’);
subplot(4,3,2);
plot(el2);hold;axis([0 N min(el2) max(el2)]);
t=[‘el2’];title(t);
axis(‘off’);
subplot(4,3,3);
plot(el3);hold;axis([0 N min(el3) max(el3)]);
t=[‘el3’];title(t);
axis(‘off’);
Causality
% MAIN LOOP: STEPPING THROUGH THE DATA
% AND COMPUTING THE (CROSS)SPECTRA
count=0;
for w=1:increment*sample_rate:N-duration*sample_rate
% Move data window into x, y, z
x=el1(w:w+duration*sample_rate-1);
y=el2(w:w+duration*sample_rate-1);
z=el3(w:w+duration*sample_rate-1);
% Initialize the Cross-Spectral arrays for averaging
Sxx=zeros(1,step);
Syy=Sxx;
Szz=Sxx;
Sxy=Sxx;
Sxz=Sxx;
Syz=Sxx;
% SECOND LOOP TO COMPUTE AVERAGE (CROSS)SPECTRA
xtemp=0:step-1;
for i=1:step:sample_rate*duration-step;
% pre-processing x
[m,b,r]=regres(xtemp,x(i:i+step-1)); % Use regression to compute trend
trend=m*xtemp+b;
x(i:i+step-1)=x(i:i+step-1)-trend;
% detrend
x(i:i+step-1)=x(i:i+step-1)-mean(x(i:i+step-1));
% demean
fx=fft(x(i:i+step-1).*hann(step)’);
% windowed fft
% pre-processing y
[m,b,r]=regres(xtemp,y(i:i+step-1));
trend=m*xtemp+b;
y(i:i+step-1)=y(i:i+step-1)-trend;
y(i:i+step-1)=y(i:i+step-1)-mean(y(i:i+step-1));
fy=fft(y(i:i+step-1).*hann(step)’);
% pre-processing z
[m,b,r]=regres(xtemp,z(i:i+step-1));
trend=m*xtemp+b;
z(i:i+step-1)=z(i:i+step-1)-trend;
z(i:i+step-1)=z(i:i+step-1)-mean(z(i:i+step-1));
fz=fft(z(i:i+step-1).*hann(step)’);
% compute all 9 spectra which are proportinal with jHj^2, Eq (7.6c)
Sxx=Sxx+fx.*conj(fx);
Syy=Syy+fy.*conj(fy);
Szz=Szz+fz.*conj(fz);
Sxy=Sxy+fx.*conj(fy);
Sxz=Sxz+fx.*conj(fz);
Syz=Syz+fy.*conj(fz);
Syx=conj(Sxy);
Szx=conj(Sxz);
Szy=conj(Syz);
end;
% Compute the power
S11=abs(Sxx).^2;
S12=abs(Sxy).^2;
S13=abs(Sxz).^2;
S21=abs(Syx).^2;
S22=abs(Syy).^2;
S23=abs(Syz).^2;
S31=abs(Szx).^2;
S32=abs(Szy).^2;
S33=abs(Szz).^2;
% Normalize
NS11=S11./max(S11);
% on diagonal the normalized power spectrum
NS12=S12./(S11+S12+S13); % Eq (7.7)
NS13=S13./(S11+S12+S13); % Eq (7.7)
NS21=S21./(S21+S22+S23); % Eq (7.7)
NS22=S22./max(S22);
% on diagonal the normalized power spectrum
NS23=S23./(S21+S22+S23); % Eq (7.7)
NS31=S31./(S31+S32+S33); % Eq (7.7)
NS32=S32./(S31+S32+S33); % Eq (7.7)
NS33=S33./max(S33);
% on diagonal the normalized power spectrum
count=count+1;
% Plot the results in the corresponding panels and
% superimpose the results for different epochs
% Titles for the panels
ttle1=[‘ ’ num2str(count) ‘ ’ ‘Spectrum el1’];
ttle2=‘ el2 - . el1’;
ttle3=‘ el3 - . el1’;
ttle4=‘ el1 - . el2’;
ttle5=‘ Spectrum el2’;
ttle6=‘ el3 - . el2’;
ttle7=‘ el1 - . el3’;
ttle8=‘ el2 - . el3’;
ttle9=‘ Spectrum el3’;
Causality
173
% Draw a red horizontal line for each 20 s analysis window
Y=[0 0];
X=[w w+duration*sample_rate];
XP=[w+1-increment*sample_rate w];
subplot(4,3,1);plot(X,Y,‘r’);if (count . 1);plot(XP,Y);end;
subplot(4,3,2);plot(X,Y,’r’);if (count . 1);plot(XP,Y);end;
subplot(4,3,3);plot(X,Y,’r’);if (count . 1);plot(XP,Y);end;
% Plot the (cross)spectral information in the lower 3 3 3 panels
subplot(4,3,4);hold on; plot(f(1:step/4),NS11(1:step/4),’k’);axis([0 60 0 1]);
title(ttle1);
subplot(4,3,5);hold on; plot(f(1:step/4),NS12(1:step/4),’k’);axis([0 60 0 1]);
title(ttle2);
subplot(4,3,6);hold on; plot(f(1:step/4),NS13(1:step/4),’k’);axis([0 60 0 1]);
title(ttle3);
subplot(4,3,7);hold on; plot(f(1:step/4),NS21(1:step/4),’k’);axis([0 60 0 1]);
title(ttle4);
subplot(4,3,8);hold on; plot(f(1:step/4),NS22(1:step/4),’k’);axis([0 60 0 1]);
title(ttle5);
subplot(4,3,9);hold on; plot(f(1:step/4),NS23(1:step/4),’k’);axis([0 60 0 1]);
title(ttle6);
subplot(4,3,10);hold on; plot(f(1:step/4),NS31(1:step/4),’k’);axis([0 60 0 1]);
title(ttle7);
subplot(4,3,11);hold on; plot(f(1:step/4),NS32(1:step/4),’k’);axis([0 60 0 1]);
title(ttle8);
subplot(4,3,12);hold on; plot(f(1:step/4),NS33(1:step/4),’k’);axis([0 60 0 1]);
title(ttle9);
% Force the script to draw the plots and pause for 1 second
drawnow;
pause(1);
% END MAIN LOOP
end;
