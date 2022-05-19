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
wc = acos(T(1));
t11 = T(2);
t01 = T(3);
t00 = t11;
t10 = 1 + t01;
%
% Design of prototype FIR lowpass filter
%
N_h = 2 * N + 1;
wp = wc;

ws = wp + wt;
%
%
NH = (N_h - 1) / 2;
Nsam_p = round(wp / deltaw_pi);
Nsam_s = round((pi - ws) / deltaw_pi);
NV = (0 :NH)';
%
%
P = zeros(NH+1, 1);
Qp = zeros(NH+1, NH+1);
for iw = 0:Nsam_p
    w = iw * deltaw_pi;
    P = P - 2 * cos(w * NV);
    Qp = Qp + cos(w * NV) * (cos(w * NV))';
end

P = wp * P / (Nsam_p + 1);
Qp = wp * Qp / (Nsam_p + 1);
%
Qs = zeros(NH+1, NH+1);
for iw = 0:Nsam_s
    w = ws + iw * deltaw_pi;
    Qs = Qs + cos(w * NV) * (cos(w * NV))';
    
end
Qs = (pi -ws) * Qs / (Nsam_p+1);
%
%
A = -0.5 * inv(Qp + Qs) * P;
h = zeros(N_h, 1);
h(NH+1) = A(1);
h(1:NH) = 0.5 * flipud(A(2:NH+1));
h(NH+2 : N_h) = 0.5 * A(2:NH+1);
%
%
subplot(2, 2, 1);

