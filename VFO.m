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

