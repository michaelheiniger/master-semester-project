clear;
clc;

f = 10; % Hz
fs = 400; % samples/s
ts = 1/fs;
d = 1; % s
t = (0:(fs*d-1))*ts;
x = cos(2*pi*f*t);

deltaT = 1/40;
% x_s = cos(2*pi*f*(t+deltaT));
% phase_shift = exp(1j*2*pi*f*deltaT);
x_s = phase_shift*cos(2*pi*f*t);

plot(t,x,'r');
hold on;
plot(t,x_s,'b');
