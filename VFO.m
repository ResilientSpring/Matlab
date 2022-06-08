%
% Design of variable fractional-order FIR differintegration
%

clear all;
clc;
N = 40;
M = 5;
w1 = 0.05 * pi;
w2 = 0.95 * pi;
p1 = -1;
p2 = 1;
pointw = 200;
pointp = 60;
%
%
NH = N / 2;
deltaw = (w2 - w1) / pointw;
deltap = (p2 - p1) / pointp;
NMA = (NH+1) * (M+1);
NMB = NH * (M+1);
point = (pointw + 1) * (pointp+1);
%
%
