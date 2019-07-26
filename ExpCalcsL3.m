%Experimental Calculations

close all
clear all
clc

%% EXP Output & Input 

% Impulse Output Data 
% Amplitude in [v*10^2]
f19 = load('10_21_2014test1_tuesday1230_2014dat.txt');
f21 = load('10_21_2014test2_tuesday1230_2014dat.txt');
f23 = load('10_21_2014test3_tuesday1230_2014dat.txt');
t19 = load('10_21_2014test1_tuesday1230_2014fft.txt');
t21 = load('10_21_2014test2_tuesday1230_2014fft.txt');
t23 = load('10_21_2014test3_tuesday1230_2014fft.txt');

tout19 = t19(:,3); tt19 = t19(:,1);
tout21 = t21(:,3); tt21 = t21(:,1);
tout23 = t23(:,3); tt23 = t23(:,1);
tin19 = t19(:,2); 
tin21 = t21(:,2);
tin23 = t23(:,2);

aout19 = f19(:,2); ft19 = f19(:,1);
aout21 = f21(:,2); ft21 = f21(:,1);
aout23 = f23(:,2); ft23 = f23(:,1);
ain19 = f19(:,3);
ain21 = f21(:,3);
ain23 = f23(:,3);

%% Impulse plots 
figure(1)
subplot(3,1,1)
plot(tt19,tout19); title('Output Impulse Test 19in'); 
xlabel('Time (s)'); 
ylabel('Amplitude'); 
xlim([0.55 .75]); grid on 

subplot(3,1,2)
plot(tt21,tout21); title('Impulse Test 21in'); 
xlabel('Time (s)'); 
ylabel('Amplitude'); 
xlim([0.7 1]); grid on 

subplot(3,1,3)
plot(tt23,tout23); title('Impulse Test 23in'); 
xlabel('Time (s)'); 
ylabel('Amplitude'); 
xlim([0.9 1.3]); grid on 

figure(2)
subplot(3,1,1)
plot(tt19,tin19); title('Input Impulse Test 19in'); 
xlabel('Time (s)'); 
ylabel('Amplitude'); 
xlim([0.55 .65]); grid on 

subplot(3,1,2)
plot(tt21,tin21); title('Impulse Test 21in'); 
xlabel('Time (s)'); 
ylabel('Amplitude'); 
xlim([0.7 .75]); grid on 

subplot(3,1,3)
plot(tt23,tin23); title('Impulse Test 23in'); 
xlabel('Time (s)'); 
ylabel('Amplitude'); 
xlim([0.9 1]); grid on 

figure(3)
subplot(3,1,1)
plot(ft19,ain19); title('Input Impulse Test 19in'); 
xlabel('Frequency (Hz)'); 
ylabel('Amplitude'); 
xlim([0 50]); grid on 

subplot(3,1,2)
plot(ft21,ain21); title('Impulse Test 21in'); 
xlabel('Frequency (Hz)'); 
ylabel('Amplitude'); 
xlim([0 50]); grid on 

subplot(3,1,3)
plot(ft23,ain23); title('Impulse Test 23in'); 
xlabel('Frequency (Hz)'); 
ylabel('Amplitude'); 
xlim([0 50]); grid on 

figure(4)
subplot(3,1,1)
plot(ft19,aout19); title('Output Impulse Test 19in'); 
xlabel('Frequency (Hz)'); 
ylabel('Amplitude'); 
xlim([0 1500]); grid on 

subplot(3,1,2)
plot(ft21,aout21); title('Impulse Test 21in'); 
xlabel('Frequency (Hz)'); 
ylabel('Amplitude'); 
xlim([0 1500]); grid on 

subplot(3,1,3)
plot(ft23,aout23); title('Impulse Test 23in'); 
xlabel('Frequency (Hz)'); 
ylabel('Amplitude'); 
xlim([0 1500]); grid on 

%% Calculations

t = 6.35; %mm
w = 25.4; %mm

Length  = [ 482.4, 533.4, 584.2 ]; % Beam lengths (mm)
maccel = 3.3*10^-3;             % mass of accelerometer (kg)
rho  = 7.87*10^-6;            % density of Carbon Steel 1018 (kg/mm^3)
E    = 200*10^9;              % Young's Modulus (Pa)
I    = (1/12)*w*t^3*10^-12;   % moment of inertia of cross section (m^4)
z = 0;

massb  = (rho)*(t*w)*(Length);     % mass of beam (kg)
masseq = (33/140)*massb + maccel;   % mass equivalent (kg)
springeq = 3*E*I./(Length*10^-3).^3; % Spring constant equivalnet (N/m)

natfreq   = sqrt(springeq./masseq);   % Natural Frequency in (rad/s)
natf   = natfreq/(2*pi);       % Natural Frequency in (Hz)

resf = natf.*sqrt(1-(2*z^2)); %resonant freq Hz
dmpf = natf.*sqrt(1-(z^2));    %damped freq Hz

% Zeta 19in
y19out1  = 194.2;
y19out2 = 111.5; 

% Zeta 21in 
y21out1  = 176.4;
y21out2 = 130; 

% Zeta 23in 
y23out1  = 155.6;
y23out2 = 96.2; 

dampratio19 = (1/2*pi)/log(y19out1/y19out2);
dampratio21 = (1/2*pi)/log(y21out1/y21out2);
dampratio23 = (1/2*pi)/log(y23out1/y23out2);




%% Magnitude of Response

TransformFunc19in = fft(tin19);
time19in          = tt19(2) - tt19(1);  %Making vectors like sizes
lt19in            = length(TransformFunc19in);
freq19in          = 1/time19in;
freqint19in       = (freq19in/2)*linspace(0,1,(lt19in/2 + 1));
complexFFT19in    = TransformFunc19in(1:(lt19in/2 + 1));
Magres19in = sqrt(real(complexFFT19in).^2+imag(complexFFT19in).^2);

TransformFunc21in = fft(tin21);
time21in          = tt21(2) - tt21(1);  %Making vectors like sizes
lt21in            = length(TransformFunc21in);
freq21in          = 1/time21in;
freqint21in       = (freq21in/2)*linspace(0,1,(lt21in/2 + 1));
complexFFT21in    = TransformFunc21in(1:(lt21in/2 + 1));
Magres21in = sqrt(real(complexFFT21in).^2+imag(complexFFT21in).^2);

TransformFunc23in = fft(tin23);
time23in          = tt23(2) - tt23(1);  %Making vectors like sizes
lt23in            = length(TransformFunc23in);
freq23in          = 1/time23in;
freqint23in       = (freq23in/2)*linspace(0,1,(lt23in/2 + 1));
complexFFT23in    = TransformFunc23in(1:(lt23in/2 + 1));
Magres23in = sqrt(real(complexFFT23in).^2+imag(complexFFT23in).^2);

figure(5)
subplot(3,1,1)
plot(freqint19in,Magres19in); title('Output Magnitude of Response 19in'); 
xlabel('Frequency (Hz)'); 
ylabel('Mag of Res');  
xlim([0 30])

subplot(3,1,2)
plot(freqint21in,Magres21in); title('Magnitude of Response 21in'); 
xlabel('Frequency (Hz)'); 
ylabel('Magnitude of Response'); 
xlim([0 30])

subplot(3,1,3)
plot(freqint23in,Magres23in); title('Magnitude of Response 23in'); 
xlabel('Frequency (Hz)'); 
ylabel('Mag of Res'); 
xlim([0 30])

TransformFunc19out = fft(tout19);
time19out          = tt19(2) - tt19(1);  %Making vectors like sizes
lt19out            = length(TransformFunc19out);
freq19out          = 1/time19out;
freqint19out       = (freq19out/2)*linspace(0,1,(lt19out/2 + 1));
complexFFT19out    = TransformFunc19out(1:(lt19out/2 + 1));
Magres19out = sqrt(real(complexFFT19out).^2+imag(complexFFT19out).^2);

TransformFunc21out = fft(tout21);
time21out          = tt21(2) - tt21(1);  %Making vectors like sizes
lt21out            = length(TransformFunc21out);
freq21out          = 1/time21out;
freqint21out       = (freq21out/2)*linspace(0,1,(lt21out/2 + 1));
complexFFT21out    = TransformFunc21out(1:(lt21out/2 + 1));
Magres21out = sqrt(real(complexFFT21out).^2+imag(complexFFT21out).^2);

TransformFunc23out = fft(tout23);
time23out          = tt23(2) - tt23(1);  %Making vectors like sizes
lt23out            = length(TransformFunc23out);
freq23out          = 1/time23out;
freqint23out       = (freq23out/2)*linspace(0,1,(lt23out/2 + 1));
complexFFT23out    = TransformFunc23out(1:(lt23out/2 + 1));
Magres23out = sqrt(real(complexFFT23out).^2+imag(complexFFT23out).^2);

figure(6)
subplot(3,1,1)
plot(freqint19out,Magres19out); title('Input Magnitude of Response 19in'); 
xlabel('Frequency (Hz)'); 
ylabel('Mag of Res');  
xlim([0 30])

subplot(3,1,2)
plot(freqint21out,Magres21out); title('Magnitude of Response 21in'); 
xlabel('Frequency (Hz)'); 
ylabel('Magnitude of Response'); 
xlim([0 30])

subplot(3,1,3)
plot(freqint23out,Magres23out); title('Magnitude of Response 23in'); 
xlabel('Frequency (Hz)'); 
ylabel('Mag of Res');
xlim([0 30])

%% Transfer Phase

Transph19 = complexFFT19out/complexFFT19in;
T19r = real(Transph19(:,2));
T19i = imag(Transph19(:,2));
phaselag19 = atan(T19i./T19r);

figure(7)
subplot(3,1,1)
plot(freqint19out,phaselag19); title('Phase Lag 19in'); 
xlabel('Frequency (Hz)'); 
ylabel('Phase Lag');
xlim([0 70])

Transph21 = complexFFT21out/complexFFT21in;
T21r = real(Transph21(:,2));
T21i = imag(Transph21(:,2));
phaselag21 = atan(T21i./T21r);

subplot(3,1,2)
plot(freqint21out,phaselag21); title('Phase Lag 21in'); 
xlabel('Frequency (Hz)'); 
ylabel('Phase Lag'); 
xlim([0 70])

Transph23 = complexFFT23out/complexFFT23in;
T23r = real(Transph23(:,3));
T23i = imag(Transph23(:,3));
phaselag23 = atan(T23i./T23r);
 
subplot(3,1,3)
plot(freqint23out,phaselag23); title('Phase Lag 23in'); 
xlabel('Frequency (Hz)'); 
ylabel('Phase Lag');
xlim([0 70])

















%% Shaker Plots

shaker19 = load('shaker19in.txt');

shaker21 = load('shaker21in.txt');

shaker23 = load('shaker23in.txt');

figure(10)
errorbar(shaker19(:,2),shaker19(:,3),shaker19(:,3)*.3,'k'); hold on
errorbar(shaker21(:,2),shaker21(:,3),shaker21(:,3)*.3,'b'); hold on
errorbar(shaker23(:,2),shaker23(:,3),shaker23(:,3)*.3,'r'); hold on

title('Shaker Plots')
xlabel('Frequency (Hz)')
ylabel('Amplitude (V)')
legend('19in','21in','23in')






