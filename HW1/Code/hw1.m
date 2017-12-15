% -----DSP HW1-----%
% R06942106 Chen Kuan-Chun

% -----part 1----- %

% ---parameter---
a1 = 1.3711242;
a2 = 0.93906244;
b = [1+a2, -2*a1, 1+a2];
a = [1, -a1, a2];
a = 2*a;

% ---transform---
[H, w] = freqz(b, a);
w = w/pi;
H_mag = abs(H);
H_pha = unwrap(angle(H));

% ---plot---
figure;
plot(w, H_mag);
xlabel('Normalized Freq');
title('Magnitude');

figure;
plot(w, H_pha);
xlabel('Normalized Freq');
title('Phase');

index = find(H_mag < 1e-6);
f0 = w(index); % normalized freq

% -----part 2----- %

% ---parameter---
F0 = 1250;
T = 0.0001; % sampling time
I = 1200; % 1200 points
n_axis = (1:I);
plus = max(length(a), length(b)) - 1; % fill zeros in front of the signals to avoid the negative index error
n_axis_plus = (1:I+plus);

% ---signals---
x1 = [zeros(1, plus), sin(2*pi*F0*n_axis*T)];
y1 = zeros(1, I+plus);
x2 = [zeros(1, plus), sin(2*pi*2*F0*n_axis*T)];
y2 = zeros(1, I+plus);

% ---transform---
for n = n_axis_plus(3:end) % use differential equations to find y[n]
    y1(n) = ((1+a2)*x1(n)-2*a1*x1(n-1)+(1+a2)*x1(n-2))/2 - (-a1*y1(n-1)+a2*y1(n-2));
end

for n = n_axis_plus(3:end) % use differential equations to find y[n]
    y2(n) = ((1+a2)*x2(n)-2*a1*x2(n-1)+(1+a2)*x2(n-2))/2 - (-a1*y2(n-1)+a2*y2(n-2));
end

% ---plot---
figure;
plot(y1);
title('y1');
xlabel('n');
legend('F0 = 1250Hz');

figure;
plot(y2);
title('y2');
xlabel('n');
legend('2F0 = 2500Hz');

% -----part 3----- %

% ---signals---
x3 = x1 + x2; % sin(2pi*F0*nT) + sin(2pi*2F0*nT)
y3 = zeros(1, I+2);

% ---transform---
for n = n_axis_plus(3:end) % use differential equations to find y[n]
    y3(n) = ((1+a2)*x3(n)-2*a1*x3(n-1)+(1+a2)*x3(n-2))/2 - (-a1*y3(n-1)+a2*y3(n-2));
end

% ---plot---
figure;
plot(y3);
title('y3');
xlabel('n');
