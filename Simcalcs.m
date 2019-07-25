%% Simulation Calculations
%Ben Ausburn

clear all
close all
clc

%% Sine Wave
sint = load('sint.txt');
sinf = load('sinf.txt');

%Fourier Integration
syms t A f n1 

x0   = A*sin(2*pi*f*t);
a0t  = int(x0*f,t,-1/(2*f),1/(2*f));
bt = int(2*f*x0*sin(2*pi*n1*f*t),t,-1/(2*f),1/(2*f));
at = int(2*f*x0*cos(2*pi*n1*f*t),t,-1/(2*f),1/(2*f));

f = 500;            %frequency (Hz)
A = max(sint(:,2)); %Amplitude
v = sint(:,2);      %voltage (volts)
t = sint(:,1);      %time
dif = 83;
n = sinf(:,1)/f;
n1 = 1;
a0t = subs(a0t); 
at = subs(at); 
bt = subs(bt);
T = 1/f; %Period

%Fourier Application
x1 = a0t;
for i = length(dif)
    k = dif(i);
    nk = n(k);
    x1 = x1 + at*cos(2*pi*n1*f*t)+ bt*sin(2*pi*n1*f*t);
    
end

dif2 = find(t <= 1/f,1,'last');
a0 = 0;
for i = 1:(dif2-1)
    k = i+1;
    d1    = t(i);
    c    = t(k);
    fa  = v(i);
    fb  = v(k);
    a0  = a0 + (c-d1)*(fa+fb)/2;
end

syms n2
an = 0;
for i = 1:(dif2-1)
    k = i+1;
    d1    = t(i);
    c    = t(k);
    fa  = (2*f)*v(i)*cos(2*pi*n2*f*t(i));
    fb  = (2*f)*v(k)*cos(2*pi*n2*f*t(k));
    an  = an + (c-d1)*(fa+fb)/2;
end

bn = 0;
for i = 1:(dif2-1)
    k = i+1;
    d1    = t(i);
    c    = t(k);
    fa  = (2*f)*v(i)*sin(2*pi*n2*f*t(i));
    fb  = (2*f)*v(k)*sin(2*pi*n2*f*t(k));
    bn  = bn + (c-d1)*(fa+fb)/2;
end

n2 = n1;
an = subs(an); bn = subs(bn);

x2 = a0;
for i = length(dif)
    k = dif(i);
    nk = n(k);
    x2 = x2 + an*cos(2*pi*nk*f*t)+ bn*sin(2*pi*nk*f*t);
end 

%% Sine Plots
figure(1)
subplot(2,2,1)
plot(sinf(:,1),sinf(:,2))
xlabel('Frequency (Hz)')
ylabel('V(f)')
grid on
title('Simulation Sine Freq')

subplot(2,2,2)
plot(sint(:,1),sint(:,2))
title('Simulation Sine Time')
axis([0 T -10 10])   %-10 to 10 b/c bipolar range 
grid on
xlabel('Time (s)')
ylabel('V(t)')

subplot(2,1,2)
plot(t,x1)
title('Fourier Sine Wave')
axis([0 T -10 10])
grid on 
xlabel('Time (s)')
ylabel('V(t)')

figure(2)
subplot(2,1,1)
plot(t,x2)
grid on 
title('Reconstruction Sine Sim')
xlabel('Time (s)')
ylabel('V(t)')
axis([ 0 T -10 10])

subplot(2,1,2)
plot(sint(:,1),sint(:,2))
title('Simulation Sine Time')
axis([0 T -10 10])
grid on
xlabel('Time (s)')
ylabel('V(t)')


%% Square Wave
close all
squaret = load('squarefilteroffdat.txt');
squaref = load('squarefilterofffft.txt');

syms t ASqu f ni
xsq = ASqu;
asq = 0;
ansq = 0;
bnsq = 4*ASqu/(pi*ni)



freqsqu = 500;
ASqu = max(squaref(:,2));
vsq = squaret(:,2);
tisq = squaret(:,1);
nSqu = squaref(:,1)/freqsqu;
ni = [1 3 5 7 9 11 13];
asq = subs(asq); ansq = subs(ansq); bnsq = subs(bnsq);

rsqu = asq;
for i = 1:length(ni)
    ng = ni(i);
    rsqu = rsqu + ansq*cos(2*pi*ng*freqsqu*tisq)+ bnsq(i)...
        *sin(2*pi*ng*freqsqu*tisq);
end

tsid2 = find(tisq <= 1/freqsqu,1,'last');
a0squ = 0;
for i = 1:(tsid2-1)
    k = i+1;
    r1    = tisq(i);
    b1    = tisq(k);
    fa  = vsq(i);
    fb  = vsq(k);
    a0squ  = a0squ + (b1-r1)*(fa+fb)/2;
end

syms nt
asqu = 0;
for i = 1:(tsid2-1)
    k = i+1;
    r1    = tisq(i);
    b1    = tisq(k);
    fa  = (2*freqsqu)*vsq(i)*cos(2*pi*nt*freqsqu*tisq(i));
    fb  = (2*freqsqu)*vsq(k)*cos(2*pi*nt*freqsqu*tisq(k));
    asqu  = asqu + (b1-r1)*(fa+fb)/2;
end

bnsqu = 0;
for i = 1:(tsid2-1)
    k = i+1;
    r1    = tisq(i);
    b1    = tisq(k);
    fa  = (2*freqsqu)*vsq(i)*sin(2*pi*nt*freqsqu*tisq(i));
    fb  = (2*freqsqu)*vsq(k)*sin(2*pi*nt*freqsqu*tisq(k));
    bnsqu  = bnsqu + (b1-r1)*(fa+fb)/2;
end

nt = ni;
asqu = subs(asqu); bnsqu = subs(bnsqu);

zrsqu = a0squ;
for i = 1:length(ni)
    ng = ni(i);
    zrsqu = zrsqu + asqu(i)*cos(2*pi*ng*freqsqu*tisq)+ bnsqu(i)...
        *sin(2*pi*ng*freqsqu*tisq);
end


%% Plots
figure(3)
subplot(2,2,1)
plot(squaref(:,1),squaref(:,2))
xlabel('Frequency (Hz)')
ylabel('V(f)')
grid on
title('Sim Square Wave Freq')

subplot(2,2,2)
plot(squaret(:,1),squaret(:,2))
title('Sim Square Wave Time')
axis([0 1/freqsqu -10 10])
grid on
xlabel('Time (s)')
ylabel('V(t)')

subplot(2,1,2)
plot(tisq,rsqu)
title('Reconsturction SQ')
axis([0 1/freqsqu -10 10])
grid on 
xlabel('Time (s)')
ylabel('V(t)')

figure (4)
subplot(2,1,1)
plot(tisq,zrsqu)
grid on 
title('Reconstruction Sim')
xlabel('Time (s)')
ylabel('V(t)')
axis([ 0 1/freqsqu -10 10])

subplot(2,1,2)
plot(squaret(:,1),squaret(:,2))
title('Sim Square Wave Time')
axis([0 1/freqsqu -10 10])
grid on
xlabel('Time (s)')
ylabel('V(t)')


%% Harmonics SQ  
n1  = bnsqu(1)*sin(2*pi*ni(1)*freqsqu*tisq);
n3  = bnsqu(2)*sin(2*pi*ni(2)*freqsqu*tisq);
n5  = bnsqu(3)*sin(2*pi*ni(3)*freqsqu*tisq);
n7  = bnsqu(4)*sin(2*pi*ni(4)*freqsqu*tisq);
n9  = bnsqu(5)*sin(2*pi*ni(5)*freqsqu*tisq);
n11 = bnsqu(6)*sin(2*pi*ni(6)*freqsqu*tisq);
n13 = bnsqu(7)*sin(2*pi*ni(7)*freqsqu*tisq);

figure (7)
plot(tisq,n1,tisq,n3,tisq,n5,tisq,n7,tisq,n9,tisq,n11,tisq,n13)
axis([0 1/freqsqu -10 10])
grid on 
title('Harmonics Square')
xlabel('Time (s)')
ylabel('V(t)')




%% Sawtooth 

close all
sawt = load('sawtoothoffdat.txt');
sawf = load('sawtoothofffft.txt');


% Calculate Fourier 
syms saw t adsaw ni 
zsaw = 2*saw*adsaw*t;
a0tsaw = int(zsaw,t,-1/(2*saw),1/(2*saw))
adntsaw = int((2*saw)*zsaw*cos(2*pi*ni*saw*t),t,-1/(2*saw),1/(2*saw))
bdntsaw = int((2*saw)*zsaw*sin(2*pi*ni*saw*t),t,-1/(2*saw),1/(2*saw))


saw = 500;
T = 1/saw;
adsaw = max(sawf(:,2));
vSaw = sawt(:,2);
tSaw = sawt(:,1);
nSaw = sawf(:,1)/saw;
ni = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 ];
a0tsaw = subs(a0tsaw); adntsaw = subs(adntsaw); bdntsaw = subs(bdntsaw);

z2saw = a0tsaw;
for i = 1:length(ni)
    ng = ni(i);
    z2saw = z2saw + adntsaw*cos(2*pi*ng*saw*tSaw)+ bdntsaw(i)...
        *sin(2*pi*ng*saw*tSaw);
end

fit2 = find(tSaw <= 1/saw,1,'last');
a0saw = 0;
for i = 1:(fit2-1)
    k = i+1;
    a1    = tSaw(i);
    b1    = tSaw(k);
    f_a1  = vSaw(i);
    f_b1  = vSaw(k);
    a0saw  = a0saw + (b1-a1)*(f_a1+f_b1)/2;
end

syms nt
ansaw = 0;
for i = 1:(fit2-1)
    k = i+1;
    a2    = tSaw(i);
    b2    = tSaw(k);
    f_a2  = (2*saw)*vSaw(i)*cos(2*pi*nt*saw*tSaw(i));
    f_b2  = (2*saw)*vSaw(k)*cos(2*pi*nt*saw*tSaw(k));
    ansaw  = ansaw + (b2-a2)*(f_a2+f_b2)/2;
end

bnsaw = 0;
for i = 1:(fit2-1)
    k = i+1;
    a    = tSaw(i);
    b    = tSaw(k);
    f_a  = (2*saw)*vSaw(i)*sin(2*pi*nt*saw*tSaw(i));
    f_b  = (2*saw)*vSaw(k)*sin(2*pi*nt*saw*tSaw(k));
    bnsaw  = bnsaw + (b-a)*(f_a+f_b)/2;
end

nt = ni;
ansaw = subs(ansaw); bnsaw = subs(bnsaw);

x_rSqu = a0saw;
for i = 1:length(ni)
    ng = ni(i);
    x_rSqu = x_rSqu + ansaw(i)*cos(2*pi*ng*saw*tSaw)+ bnsaw(i)...
        *sin(2*pi*ng*saw*tSaw);
end


%% Sawtooth Plots
figure(5)
subplot(2,2,1)
plot(sawf(:,1),sawf(:,2))
xlabel('Frequency (Hz)')
ylabel('V(f)')
grid on
title('Sim Sawtooth Wave Freq')

subplot(2,2,2)
plot(sawt(:,1),sawt(:,2))
title('Sim Sawtooth Wave Time')
axis([0 1.4*10e-4 -10 10])
grid on
xlabel('Time (s)')
ylabel('V(t)')

subplot(2,1,2)
plot(tSaw,z2saw)
title('Reconsturction Saw')
axis([0 1/saw -10 10])
grid on 
xlabel('Time (s)')
ylabel('V(t)')

figure (6)
subplot(2,1,2)
plot(tSaw,x_rSqu)
grid on 
title('Reconstruction')
xlabel('Time (s)')
ylabel('V(t)')
axis([ 0 1.4*10e-4 -10 10])

subplot(2,1,1)
plot(sawt(:,1),sawt(:,2))
title('Sim Sawtooth Wave Time')
axis([0 1.4*10e-4 -10 10])
grid on
xlabel('Time (s)')
ylabel('V(t)')

% Harmonics of Sawtooth Wave 
n1  = bnsaw(1)*sin(2*pi*ni(1)*saw*tSaw);
n2  = bnsaw(2)*sin(2*pi*ni(2)*saw*tSaw);
n3  = bnsaw(3)*sin(2*pi*ni(3)*saw*tSaw);
n4  = bnsaw(4)*sin(2*pi*ni(4)*saw*tSaw);
n5  = bnsaw(5)*sin(2*pi*ni(5)*saw*tSaw);
n6  = bnsaw(6)*sin(2*pi*ni(6)*saw*tSaw);
n7  = bnsaw(7)*sin(2*pi*ni(7)*saw*tSaw);
n8  = bnsaw(8)*sin(2*pi*ni(8)*saw*tSaw);
n9  = bnsaw(9)*sin(2*pi*ni(9)*saw*tSaw);
n10 = bnsaw(10)*sin(2*pi*ni(10)*saw*tSaw);
n11 = bnsaw(11)*sin(2*pi*ni(11)*saw*tSaw);
n12 = bnsaw(12)*sin(2*pi*ni(12)*saw*tSaw);
n13 = bnsaw(13)*sin(2*pi*ni(13)*saw*tSaw);
n14 = bnsaw(14)*sin(2*pi*ni(14)*saw*tSaw);
figure (10)
plot(tSaw,n1,tSaw,n2,tSaw,n3,tSaw,n4,tSaw,n5,tSaw,n6,tSaw,n7,tSaw,n8...
    ,tSaw,n9,tSaw,n10,tSaw,n11,tSaw,n12,tSaw,n13,tSaw,n14)
axis([0 1/saw -10 10])
grid on 
title('Harmonics Sawtooth')
xlabel('Time (s)')
ylabel('V(t)')

%% Low 4bit sampling graphs
%% Sine low sample

sintl = load('2.3.1lowHresdat.txt');
sinfl = load('2.3.1lowFresfft.txt');

f = 500;
T = 1/f;

figure(1)
subplot(2,1,1)
plot(sinfl(:,1),sinfl(:,2))
xlabel('Frequency (Hz)')
ylabel('V(f)')
grid on
title('Simulation Sine Freq Low')

subplot(2,1,2)
plot(sintl(:,1),sintl(:,2))
title('Simulation Sine Time Low')
axis([0 T -10 10])   %-10 to 10 b/c bipolar range 
grid on
xlabel('Time (s)')
ylabel('V(t)')

%% Saw low sample

sawtl = load('sawtoothondat.txt');
sawfl = load('sawtoothonfft.txt');

f = 500;
T = 1/f;

figure(1)
subplot(2,1,1)
plot(sawfl(:,1),sawfl(:,2))
xlabel('Frequency (Hz)')
ylabel('V(f)')
grid on
title('Simulation Saw Freq Low')

subplot(2,1,2)
plot(sawtl(:,1),sawtl(:,2))
title('Simulation Saw Time Low')
axis([0 T -10 10])   %-10 to 10 b/c bipolar range 
grid on
xlabel('Time (s)')
ylabel('V(t)')



squaretl = load('squarefilterondat.txt');
squarefl = load('squarefilteronfft.txt');

f = 500;
T = 1/f;

figure(1)
subplot(2,1,1)
plot(squarefl(:,1),squarefl(:,2))
xlabel('Frequency (Hz)')
ylabel('V(f)')
grid on
title('Simulation Saw Freq Low')

subplot(2,1,2)
plot(squaretl(:,1),squaretl(:,2))
title('Simulation Saw Time Low')
axis([0 T -10 10])   %-10 to 10 b/c bipolar range 
grid on
xlabel('Time (s)')
ylabel('V(t)')




















