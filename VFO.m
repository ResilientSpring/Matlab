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
ra = zeros(NMA, 1);
Qa = zeros(NMA, NMA);
rb = zeros(NMB, 1);
Qb = zeros(NMB, NMB);
for ip = 0:pointp
    p = p1 + ip * deltap;
    for iw = 0:pointw
        w = w1 + iw * deltaw;
        cwp = zeros(NMA, 1);
        for inm = 0:NMA - 1
            n = mod(inm, NH+1);
            m = floor(inm/(NH+1));
            cwp(inm+1) = p ^ m * cos(n * w);
        end
        ra = ra - 2 * w ^ p * cos(p * pi / 2) * cwp;
        Qa = Qa + cwp * cwp';

        swp = zeros(NMB, 1);
        for inm = 0:NMB - 1
            n = mod(inm, NH) + 1;
            m = floor(inm / NH);
            swp(inm + 1) = p ^ m * sin(n * w);
        end
        rb = rb - 2 * w ^ p * sin(p*pi/2) * swp;
        Qb = Qb + swp * swp';
    end
end
ra = (p2 - p1) * (w2 - w1) * ra / point;
Qa = (p2 - p1) * (w2 - w1) * Qa / point;
rb = (p2 - p1) * (w2 - w1) * rb / point;
Qb = (p2 - p1) * (w2 - w1) * Qb / point;
%
%
a = - 0.5 * inv(Qa) * ra;
a2 = reshape(a, NH+1, M+1);
b = - 0.5 * inv(Qb) * rb;
b2 = reshape(b, NH, M+1);
%
he = zeros(N+1, M+1);
he(NH+1, :) = a2(1, :);
he(NH+2:N+1, :) = 0.5 * a2(2:NH+1, :);
he(NH:-1:1, :) = he(NH+2:N+1, :);
%
ho = zeros(N+1, M+1);
ho(NH+2:N+1, :) = -0.5 * b2(1:NH, :);
ho(NH:-1:1, :) = -ho(NH+2:N+1, :);
%
h = he + ho;
%
MAG = zeros(pointw+1, pointp+1);
XX = zeros(pointw+1, pointp+1);
YY = zeros(pointw+1, pointp+1);
for ip = 0:pointp
    p = p1 + ip * deltap;
    hnp = h(:, 1);
    for im = 1:M
        hnp = hnp + p^im*h(:, im+1);
    end
    MAG(:, ip+1) = abs(freqz(hnp, 1, w1:deltaw:w2));
    for iw = 0:pointw
        w = w1 + iw * deltaw;
        XX(iw+1, ip+1) = w / pi;
        YY(iw+1, ip+1) = p;
    end
end

plot3(XX, YY, MAG);
xlabel('Normalized frequency');
ylabel('Variable p');
zlabel('Variable magnitude response');
