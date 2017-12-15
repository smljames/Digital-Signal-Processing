% ----- DSP HW2 ----- %
% R06942106 Chen Kuan-Chun

% Real Subroutine 
%     FFT(X, Y, N, M, KIND) 
% is in the FFT.m

% ----- part 1 ----- %
clc; clear;
x.r = [1 1 2 2 1 0 0 0];   % real part of data
x.i = [1 -1 1 -1 2 4 5 8]; % imaginary part of data
N = length(x.r);           % total length
M = log2(N);               % the power of 2

% ----- FFT -----
y = fft(x.r + 1i*x.i, N); % uncomment to see fft result of matlab
[X.r, X.i] = FFTDIT(x.r, x.i, N, M, 1);

% ----- display -----
x.r + 1i*x.i;

% -----IFFT -----
[x.r, x.i] = FFTDIT(X.r, X.i, N, M, -1); % IFFT by using FFT

% ----- display -----
x.r + 1i*x.i;


% ----- part 2-1 ----- %
clc; clear;
x.r = [1 2 3 4 5 6 7 8];
x.i = [0 0 0 0 0 0 0 0];
y.r = [2 3 5 6 1 9 8 7];
y.i = [0 0 0 0 0 0 0 0];
N = length(x.r);
M = log2(N);

matlab_circular_conv = cconv(x.r+1i*x.i, y.r+1i*y.i, 8);
                              % uncomment to see circular convolution
                              % result of matlab

% ----- FFT -----
[X.r, X.i] = FFTDIT(x.r, x.i, N, M, 1);
[Y.r, Y.i] = FFTDIT(y.r, y.i, N, M, 1);

% ----- multiply in frequency domain -----
Z.r = X.r.*Y.r - X.i.*Y.i;
Z.i = X.r.*Y.i + Y.r.*X.i;

% ----- IFFT -----
[z.r, z.i] = FFTDIT(Z.r, Z.i, N, M, -1);

% ----- display -----
z.r + 1i*z.i;


% ----- part 2-2 ----- %
x.r = [1 2 3 3 2 2 1 2 0 0 0 0 0 0 0 0]; % zero-padding in x
x.i = zeros(1, 16);
y.r = [2 3 5 6 1 9 8 7 0 0 0 0 0 0 0 0]; % zero-padding in y
y.i = zeros(1, 16);
N = length(x.r);
M = log2(N);

matlab_conv = conv(x.r+1i*x.i, y.r+1i*y.i); % normal convolution
matlab_conv(1:16); % uncomment to see the result of matlab

% ----- FFT -----
[x.r, x.i] = FFTDIT(x.r, x.i, N, M, 1);
[y.r, y.i] = FFTDIT(y.r, y.i, N, M, 1);

% ----- multiply in frequency domain -----
z.r = x.r.*y.r - x.i.*y.i;
z.i = x.r.*y.i + x.i.*y.r;

% ----- IFFT -----
[z.r, z.i] = FFTDIT(z.r, z.i, N, M, -1);

% ----- display -----
z.r;


% ----- part 3 ----- %
x.r = [1 2 3 4 5 6 7 8];
x.i = [0 0 0 0 0 0 0 0];
N = length(x.r);
M = log2(N);

[X.r, X.i] = FFTDIFnbr(x.r, x.i, N, M, 1);
[x.r, x.i] = FFTDITnbr(X.r, X.i, N, M, -1);

% ----- display -----
x.r+1i*x.i;


% ----- part 4-1 ----- %
clc; clear;
T = 8;
N = 16;
M = log2(N);
dt = T/N;
df = 1/T;
n1 = 0:N-1;

% ----- sample points in time domain -----
x1.r = exp(-n1*dt);
x1.i = zeros(1, N);

% ----- correct results in frequency domain -----
Fx1 = 1./(1+1i*2*pi*n1*df);

% ----- FFT -----
[X1.r, X1.i] = FFTDIT(x1.r, x1.i, N, M, 1);
tmp = X1.r(1);
X1.r = X1.r./tmp;
X1.i = X1.i./tmp;

% ----- plot -----
figure;
subplot(2,2,1);
hold on;
plot(n1.*df, real(Fx1), n1.*df, X1.r);
ylabel('Real F(f)');

subplot(2,2,2);
hold on;
plot(n1.*df, imag(Fx1), n1.*df, X1.i);
ylabel('Image F(f)');

subplot(2,2,3);
hold on;
plot(n1.*df, abs(Fx1), n1.*df, abs(X1.r+1i*X1.i));
ylabel('Magnitude response');

subplot(2,2,4);
hold on;
plot(n1.*df, angle(Fx1), n1.*df, angle(X1.r+1i*X1.i));
ylabel('Phase response');

hL = legend('correct result', 'FFT N=16');
newPosition = [0.45 0.46 0.1 0.1];
newUnits = 'normalized';
set(hL, 'Position', newPosition, 'Units', newUnits);


% ----- part 4-2 ----- %
T = 8;
N = 32;
M = log2(N);
dt = T/N;
df = 1/T;
n2 = 0:N-1;

% ----- sample points in time domain -----
x2.r = exp(-n2*dt);
x2.i = zeros(1, N);

% ----- correct results in frequency domain -----
Fx2 = 1./(1+1i*2*pi*n2*df);

% ----- FFT -----
[X2.r, X2.i] = FFTDIT(x2.r, x2.i, N, M, 1);
tmp = X2.r(1);
X2.r = X2.r./tmp;
X2.i = X2.i./tmp;

% ----- plot -----
figure;
subplot(2,2,1);
hold on;
plot(n2.*df, real(Fx2), n2.*df, X2.r);
ylabel('Real F(f)');

subplot(2,2,2);
hold on;
plot(n2.*df, imag(Fx2), n2.*df, X2.i);
ylabel('Image F(f)');

subplot(2,2,3);
hold on;
plot(n2.*df, abs(Fx2), n2.*df, abs(X2.r+1i*X2.i));
ylabel('Magnitude response');

subplot(2,2,4);
hold on;
plot(n2.*df, angle(Fx2), n2.*df, angle(X2.r+1i*X2.i));
ylabel('Phase response');

hL = legend('correct result', 'FFT N=32');
newPosition = [0.45 0.46 0.1 0.1];
newUnits = 'normalized';
set(hL, 'Position', newPosition, 'Units', newUnits);


% ----- part 4-3 ----- %
T = 8;
N = 64;
M = log2(N);
dt = T/N;
df = 1/T;
n3 = 0:N-1;

% ----- sample points in time domain -----
x3.r = exp(-n3*dt);
x3.i = zeros(1, N);

% ----- correct results in frequency domain -----
Fx3 = 1./(1+1i*2*pi*n3*df);

% ----- FFT -----
[X3.r, X3.i] = FFTDIT(x3.r, x3.i, N, M, 1);
tmp = X3.r(1);
X3.r = X3.r./tmp;
X3.i = X3.i./tmp;

% ----- plot -----
figure;
subplot(2,2,1);
hold on;
plot(n3.*df, real(Fx3), n3.*df, X3.r);
ylabel('Real F(f)');

subplot(2,2,2);
hold on;
plot(n3.*df, imag(Fx3), n3.*df, X3.i);
ylabel('Image F(f)');

subplot(2,2,3);
hold on;
plot(n3.*df, abs(Fx3), n3.*df, abs(X3.r+1i*X3.i));
ylabel('Magnitude response');

subplot(2,2,4);
hold on;
plot(n3.*df, angle(Fx3), n3.*df, angle(X3.r+1i*X3.i));
ylabel('Phase response');

hL = legend('correct result', 'FFT N=64');
newPosition = [0.45 0.46 0.1 0.1];
newUnits = 'normalized';
set(hL, 'Position', newPosition, 'Units', newUnits);


% ----- plot all -----
figure;
subplot(2,2,1);
hold all
plot(n3.*df, real(Fx3));
plot(n1.*df, X1.r);
plot(n2.*df, X2.r);
plot(n3.*df, X3.r);
ylabel('Real F(f)');

subplot(2,2,2);
hold all
plot(n3.*df, imag(Fx3));
plot(n1.*df, X1.i);
plot(n2.*df, X2.i);
plot(n3.*df, X3.i);
ylabel('Image F(f)');

subplot(2,2,3);
hold all
plot(n3.*df, abs(Fx3));
plot(n1.*df, abs(X1.r+1i*X1.i));
plot(n2.*df, abs(X2.r+1i*X2.i));
plot(n3.*df, abs(X3.r+1i*X3.i));
ylabel('Magnitude response');

subplot(2,2,4);
hold all
plot(n3.*df, angle(Fx3));
plot(n1.*df, angle(X1.r+1i*X1.i));
plot(n2.*df, angle(X2.r+1i*X2.i));
plot(n3.*df, angle(X3.r+1i*X3.i));
ylabel('Phase response');

hL = legend('correct result', 'FFT N=16', 'FFT N=32', 'FFT N=64');
newPosition = [0.45 0.46 0.1, 0.1];
newUnits = 'normalized';
set(hL, 'Position', newPosition, 'Units', newUnits);


