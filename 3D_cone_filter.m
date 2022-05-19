%
% Design of 3D cone filters
%
%
% Step 1: Fan-type design
%

clc;
clear all;
N = 20;
theta = 70;
wt = 0.2 * pi;
Nsamp = 200;
%
%
theta_w = theta * pi / 180;

if theta >= 45
    wup = pi;
else
    wup = pi * tan(theta_w);
end

deltaw_wup = wup / Nsamp;
deltaw_pi = pi / Nsamp;
deltaw_pih = 0.5 * pi / Nsamp;
%
%
r = zeros(3, 1);
Q = zeros(3, 3);
for iw = 0:Nsamp
    w2 = iw * deltaw_wup;
    w1 = w2 / tan(theta_w);
    c = [1; -1-cos(w1)*cos(w2); -cos(w1)-cos(w2)];
    r = r - 2 * cos(w1) * c;
    Q = Q + c * c';
end