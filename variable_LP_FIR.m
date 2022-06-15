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
sampling_pass = 0;
sampling_stop = 0;
ra = zeros(nma, 1);
Qp = zeros(nma, nma);
Qs = zeros(nma, nma);
for iw = 0:sampling_w
    w = iw * deltaw;
    for ip = 0:sampling_p
        p = -0.5 + ip * deltap;
        cwp = zeros(nma, 1);
        for inm = 0:nma-1
            n = mod(inm, NH+1);
            m = floor(inm/(NH+1));
            cwp(inm+1) = p^(m)*cos(n*w);
        end
    end
end
