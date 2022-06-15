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
        if p >= (w-wp1)/(wp2-wp1)-0.5
            sampling_pass = sampling_pass + 1;
            ra = ra - 2 * cwp;
            Qp = Qp + cwp*cwp';
        elseif p <= (w - wp1 - wt) / (wp2 - wp1) -0.5
            sampling_stop = sampling_stop + 1;
            Qs = Qs +cwp*cwp';
        else
             
        end
    end
end
ra = 0.5 * (wp1 + wp2) * ra / sampling_pass;
Qp = 0.5 * (wp1 + wp2) * Qp / sampling_pass;
Qs = 0.5 * (pi - wp1 - wt + pi - wp2 - wt) * Qs / sampling_stop;
a = -0.5 * inv(Qp + Qs) * ra;
a2 = reshape(a, NH+1, M+1);
%
%
h = zeros(N+1, M+1);
for im = 0:M
    h(NH+1, im+1) = a2(1, im+1);
    h(1:NH, im+1) = 0.5 * a2(NH+1:-1:2, im+1);
    h(NH+2:N+1, im+1) = 0.5 * a2(2:NH+1, im+1);
end
%
%
MR = zeros(sampling_w+1, sampling_p+1);
for ip = 0:sampling_p
    p = -0.5 +  ip * deltap;
    h1 = h(:, 1);
    for im = 1:M
        h1 = h1 + p^im*h(:, im+1);
    end
    MR(:, ip+1) = abs(freqz(h1, 1, 0:deltaw:pi));
end
XX = zeros(sampling_w+1, sampling_p+1);
YY = zeros(sampling_w+1, sampling_p+1);
for nw = 0:sampling_w
    w = nw * deltaw;
    XX(nw+1, :) = (w/pi)*ones(1, sampling_p+1);
end
for np = 0:sampling_p
    p = -0.5 + np * deltap;
    YY(:, np+1) = p*ones(sampling_w+1, 1);    
end

close all;
plot3(XX, YY, MR);
axis([0, 1, -0.5, 0.5, 0, 1.1]);
xlabel('Normalized frequency (\omega/\pi');
ylabel('Variable p');
zlabel('Magnitude response');
