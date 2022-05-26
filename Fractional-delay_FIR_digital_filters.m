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

