%
% Design of variable fractional-delay (VFD) FIR digital filters
%
clear all;
clc;
N = 50;
M = 7;
wp1 = 0.3 * pi;
wp2 = 0.6 * pi;
wt = 0.1 * pi;
sampling_w = 200;
sampling_p = 60;
%
%
NH = N / 2;
nma = (M+1) * (NH+1);
deltaw = pi / sampling_w;
deltap = 1 / sampling_p;
%
%
