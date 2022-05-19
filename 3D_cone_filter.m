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
r = r / (Nsamp + 1);
Q = Q / (Nsamp + 1);
T = -0.5 * inv(Q) * r;
%
%
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
FR = freqz(h, 1, 0:deltaw_pi : pi);
plot(0:deltaw_pi/pi : 1, abs(FR));
axis([0, 1, 0, 1.1]);
xlabel('Normalized frequency (\omega/\pi)');
ylabel('Magnitude response');
title('Prototype filter');
%
%
Ta_b = zeros(N+1, N+1);
Ta_b(1, 1) = 1; 
Ta_b(2, 2) = 1;
B = A(1) * Ta_b(1, :);
B = B + A(2) * Ta_b(2, :);

for in = 2:N
    Ta_b(in+1, 2:in+1) = 2 * Ta_b(in, 1:in);
    Ta_b(in+1, 1:in+1) = Ta_b(in+1, 1:in+1) - Ta_b(in-1, 1:in+1);
    B = B+A(in+1) * Ta_b(in+1, :);
end
%
%
F = zeros(3, 3);
F(1, 1) = 0.25 * t11;
F(1, 2) = 0.5 * t10;
F(1, 3) = 0.25 * t11;
F(2, 1) = 0.5 * t01;
F(2, 2) = t00;
F(2, 3) = 0.5 * t01;
F(3, 1) = 0.25 * t11;
F(3, 2) = 0.5 * t10;
F(3, 3) = 0.25 * t11;
%
%
h = zeros(2 * N + 1, 2*N+1);
h(N+1, N+1) = B(1);
h(N:N+2, N:N+2) = h(N:N+2, N:N+2) + B(2) *F;
FN = F;
for in = 2:N
    FN = conv2(FN, F);
    h(N+1-in:N+1+in, N+1-in : N+1+in) = h(N+1-in:N+1+in, N+1-in:N+1+in) + B(in+1)*FN;
end
%
%
FR = freqz2(h, -1:2/100:1, -1:2/100:1);
XX = zeros(101, 101);
for i = 0:100
    XX(:, i+1) = (-1:2/100:1)';
end

YY = XX';
subplot(2, 2, 2);
plot3(XX, YY, abs(FR));
axis([-1, 1, -1, 1, 0, 1.1]);
xlabel('\omega_{12}/\pi');
ylabel('\omega_3/\pi');
zlabel('Magnitude response');

%
% Step 2: Circular design
%
wr = wup;
r = zeros(2, 1);
Q = zeros(2, 2);

for iw = 0:Nsamp
    w1 = wr * cos(iw * deltaw_pih);
    w2 = wr * sin(iw * deltaw_pih);
    c = [1; 1-cos(w1)*cos(w2)];
    r = r - 2 * (0.5 * cos(w1) + 0.5 * cos(w2))*c;
    Q = Q + c * c';
end

r = r / (Nsamp + 1);
Q = Q / (Nsamp + 1);
T = -0.5 * inv(Q) * r;
%
%
r11 = T(2);
r00 = -r11;
r10 = 0.5;
r01 = 0.5;
%
% Design of prototype FIR lowpass filter
%
N_h = 2 * N + 1;
wp = wc;
ws = wp + wt;
%
%
NH = (N_h - 1) / 2;
Nsam_s = round((pi - ws) / deltaw_pi);
NV = (0:NH)';
%
%
P = zeros(NH+1, 1);
Qp = zeros(NH+1, NH+1);
for iw = 0:Nsam_p
    w = iw * deltaw_pi;
    P = P - 2*cos(w*NV);
    Qp = Qp + cos(w*NV)*(cos(w*NV))';
end

P = wp*P / (Nsam_p + 1);
Qp = wp * Qp / (Nsam_p + 1);
%
Qs = zeros(NH+1, NH+1);

for iw = 0:Nsams
    w = ws + iw * deltaw_pi;
    Qs = Qs + cos(w * NV) * (cos(w*NV))';
end
Qs = (pi - ws) * Qs / (Nsam_p+1);
%
%
A = -0.5 * inv(Qp + Qs) * P;
%
%
Ta_b = zeros(N+1, N+1);
Ta_b(1, 1) = 1;
Ta_b(2, 2) = 1;
B_c = A(1) * Ta_b(1, :);
B_c = B_c + A(2) * Ta_b(2, :);
for in = 2:N
    Ta_b(in+1, 2:in+1) = 2 * Ta_b(in, 1 : in);
    Ta_b(in+1, 1:in+1) = Ta_b(in+1, 1 : in+1) - Ta_b(in-1, 1:in+1);
    B_c = B_c + A(in + 1) * Ta_b(in+1, :);
end
%
%
F = zeros(3, 3);
F(1, 1) = 0.25 * r11;
F(1, 2) = 0.5 * r10;
F(1, 3) = 0.25 * r11;
F(2, 1) = 0.5 * r01;
F(2, 2) = r00;
F(2, 3) = 0.5 * r01;
F(3, 1) = 0.25 * r11;
F(3, 2) = 0.5 * r10;
F(3, 3) = 0.25 * r11;
%
%
h = zeros(2*N+1, 2*N+1);
h(N+1, N+1) = B_c(1);
h(N:N+2, N:N+2) = h(N:N+2, N:N+2) + B_c(2) * F;
FN = F;
for in = 2:N
    FN = conv2(FN, F);
    h(N+1-in:N+1+in, N+1-in:N+1+in) = h(N+1 - in:N+1+in, N+1-in:N+1+in) + B_c(in+1) * FN;
end
%
%
FR = freqz2(h, -1:2/100:1, -1:2/100:1);
XX = zeros(101, 101);
for i = 0:100
    XX(:, i+1) = (-1:2/100:1)';
end
YY = XX';
subplot(2, 2, 3);
plot3(XX, YY, abs(FR));
axis([-1, 1, -1, 1, 0, 1.1]);
xlabel('\omega_1/\pi');
ylabel('\omega_2/\pi');
zlabel('Magnitude response');
