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

for iw = 0:Nsam_p
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
%
% Step 2:
%
t000 = t00 + t10 * r00;
t100 = t10 * r10;
t010 = t10 * r01;
t001 = t01+t11*r00;
t110 = t10 * r11;
t101 = t11 * r10;
t011 = t11 * r01;
t111 = t11 * r11;
%
%
F = zeros(3, 3, 3);
F(1, 1, 1) = 0.125 * t111;
F(1, 1, 2) = 0.25 * t110;
F(1, 1, 3) = 0.125 * t111;
F(1, 2, 1) = 0.25 * t101;
F(1, 2, 2) = 0.5 * t100;
F(1, 2, 3) = 0.25 * t101;
F(1, 3, 1) = 0.125 * t111;
F(1, 3, 2) = 0.25 * t100;
F(1, 3, 3) = 0.125 * t111;
F(2, 1, 1) = 0.25 * t011;
F(2, 1, 2) = 0.5 * t010;
F(2, 1, 3) = 0.25 * t011;
F(2, 2, 1) = 0.5 * t001;
F(2, 2, 2) = t000;
F(2, 2, 3) = 0.5 * t001;
F(2, 3, 1) = 0.25 * t011;
F(2, 3, 2) = 0.5 * t010;
F(2, 3, 3) = 0.25 * t011;
F(3, 1, 1) = 0.125 * t111;
F(3, 1, 2) = 0.25 * t110;
F(3, 1, 3) = 0.125 * t111;
F(3, 2, 1) = 0.25 * t101;
F(3, 2, 2) = 0.5 * t100;
F(3, 2, 3) = 0.25 * t101;
F(3, 3, 1) = 0.125 * t111;
F(3, 3, 2) = 0.25 * t110;
F(3, 3, 3) = 0.125 * t111;
%
%
h3 = zeros(2*N + 1, 2*N + 1, 2*N + 1);
h3(N+1, N+1, N+1) = B(1);
h3(N:N+2, N:N+2, N:N+2) = h3(N:N+2, N:N+2, N:N+2) + B(2) * F;
FN = F;

for in = 2:N
    FN = convn(FN, F);
    h3(N+1-in:N+1+in, N+1-in:N+1+in, N+1-in:N+1+in) = h3(N+1-in:N+1+in, N+1-in:N+1+in, N+1-in:N+1+in) + B(in+1) * FN;
end
%
%
A3 = h3(N+1:2*N+1, N+1:2*N+1, N+1:2*N+1);
A3(2:N+1, 1, 1) = 2 * A3(2:N+1, 1, 1);
A3(1, 2:N+1, 1) = 2 * A3(1, 2:N+1, 1);
A3(1, 1, 2:N+1) = 2 * A3(1, 1, 2:N+1);
A3(2:N+1, 2:N+1, 1) = 4 * A3(2:N+1, 2:N+1, 1);
A3(2:N+1, 1, 2:N+1) = 4 * A3(2:N+1, 1, 2:N+1);
A3(1, 2:N+1, 2:N+1) = 4 * A3(1, 2:N+1, 2:N+1);
A3(2:N+1, 2:N+1, 2:N+1) = 8 * A3(2:N+1, 2:N+1, 2:N+1);
%
%
deltaw_3d = pi / 32;
FR3 = zeros(65, 65, 65);
XXX = zeros(65, 65, 65);
YYY = zeros(65, 65, 65);
ZZZ = zeros(65, 65, 65);

for i1 = 0:64
    i1
    w1 = -pi + i1 * deltaw_3d;

    for i2 = 0:64
        w2 = -pi + i2 * deltaw_3d;
        for i3 = 0:64
            w3 = -pi + i3 * deltaw_3d;
            XXX(i1+1, i2+2, i3+1) = w1/pi;
            YYY(i1+1, i2+2, i3+1) = w2/pi;
            ZZZ(i1+1, i2+2, i3+1) = w3/pi;

            for n1 = 0:N
                for n2 = 0:N
                    for n3 = 0:N
                        FR3(i1+1, i2+1, i3+1) = FR3(i1+1, i2+1, i3+1) + A3(n1+1, n2+1, n3+1) * cos(n1 * w1) * cos(n2 * w2) * cos(n3 * w3);
                    end
                end
            end
        end        
    end
end
%
%
[F, V] = isosurface(XXX, YYY, ZZZ, abs(FR3), 0.98);
%
[rowV, colV] = size(V);
%
[Z1 ind] = sort(V(:, 3));
X1 = V(ind, 1);
Y1 = V(ind, 2);
%
Zu = Z1(rowV/2+1:rowV, 1);
Xu = X1(rowV/2+1:rowV, 1);
Yu = Y1(rowV/2+1:rowV, 1);
Zd = Z1(1:rowV/2, 1);
Xd = X1(1:rowV/2, 1);
Yd = Y1(1:rowV/2, 1);
%
ti = -1:.0025:1;
[XI, YI] = meshgrid(ti, ti);
ZU = griddata(Xu, Yu, Zu, XI, YI);
ZD = griddata(Xd, Yd, Zd, XI, YI);
%
subplot(2, 2, 4);
contour3(XI, YI, ZU, 35, 'k');
hold;
contour3(XI, YI, ZD, 35, 'k');
%
view(3);
xlabel('\omega_1 / \pi', 'FontSize', 16);
set(gca, 'xtick', linspace(-1, 1, 5));
ylabel('\omega_2 / \pi', 'FontSize', 16);
set(gca, 'ytick', linspace(-1, 1, 5));
zlabel('\omega_3 / \pi', 'FontSize', 16);
set(gca, 'ztick', linspace(-1, 1, 5));
