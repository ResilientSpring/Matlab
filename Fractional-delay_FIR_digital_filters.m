%
% Design of variable fractional-delay FIR digital filters.
%

clear all;
clc;
N = 50;
M = 7;
wp = 0.9 * pi;
pointw = 200;
pointp = 60;
%
%
j = 0 + i;
NH = N / 2;

if mod(M, 2) == 0
    Mc = M / 2;
    Ms = M / 2;
else
    Mc = (M - 1) / 2;
    Ms = (M + 1) / 2;    
end

nma = (NH+1) * Mc;
nmb = NH * Ms;
deltaw = wp / pointw;
deltap = 1/ pointp;
point = (pointw + 1) * (pointp + 1);
%
%
ra = zeros(nma, 1);
Qa = zeros(nma, nma);

for iw = 0:pointw
    w = iw * deltaw;
    for ip = 0:pointp
        p = - 0.5 + ip * deltap;
        c = zeros(nma, 1);
        for ic = 0:nma - 1
            n = mod(ic, NH + 1);
            m = floor(ic/ (NH+1))+1;
            c(ic + 1) = p^(2*m)*cos(n*w);
        end
        ra = ra - 2 * (cos(p * w) - 1) * c;
        Qa = Qa + c * c';
    end
end
ra = wp * ra / point;
Qa = wp * Qa / point;
a = -0.5 * inv(Qa) * ra;
%
%
rb = zeros(nmb, 1);
Qb = zeros(nmb, nmb);

for iw = 0:pointw
    w = iw * deltaw;

    for ip = 0:pointp
        s = zeros(nmb, 1);
        for is = 0:nmb - 1
            n = mod(is, NH) + 1;
            m = floor(is/NH) + 1;
            s(is+1) = p^(2*m-1)*sin(n*w);
        end
        rb = rb + 2*sin(p * w) * s;
        Qb = Qb + s * s';
    end
end
rb = wp * rb / point;
Qb = wp * Qb / point;
b = - 0.5 * inv(Qb)*rb;
%
%
a2 = reshape(a, NH+1, Mc);
b2 = reshape(b, NH, Ms);
h = zeros(N+1, M+1);
h(NH+1, 1) = 1;
for im = 1:Mc
    h(NH+1, 2 * im + 1) = a2(1, im);
    h(1 : NH, 2 * im + 1) = 0.5 * a2(NH+1:-1:2, im);
    h(NH + 2 : NH + 1, 2 * im + 1) = 0.5 * a2(2 : NH + 1, im);
end
for im = 1:Ms
    h(1:NH, 2 * im) = 0.5 * b2(NH:-1:1, im);
    h(NH+2:N+1, 2 * im) = -0.5 * b2(1:NH, im);
end
%
%
MR = zeros(pointw+1, pointp+1);
GD = zeros(pointw+1, pointp+1);

for ip = 0:pointp
    p = -0.5 + ip * deltap;
    hnp = h(:, 1);

    for im = 1:M
        hnp = hnp + h(:, im+1) * p ^ im;
    end
    
end
