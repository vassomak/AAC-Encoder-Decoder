function frameT = iFilterbank(frameF, frameType, winType)    
N = 2 * length(frameF);
assert(N==2048);

[long, short] = windows(N, winType);
frameT = zeros(N,2);

%IMDCT and Windows
if strcmp(frameType, 'OLS')
    s = imdct(frameF);
    Win = long;
    frameT = s .* Win;

elseif strcmp(frameType, 'LSS')
    s = imdct(frameF);
    Win = [long(1:1024,:); ones(448, 2); short(129:256,:); zeros(448,2)];
    frameT = s .* Win;
    
elseif strcmp(frameType, 'LPS')
    s = imdct(frameF);
    Win = [zeros(448,2); short(1:128,:); ones(448,2); long(1025:2048,:)];
    frameT = s .* Win;
    
elseif strcmp(frameType, 'ESH')
    width = N/8;
    Win = short;
    
    for j =1:8
        range1 = (j-1) * width/2  + 1 : j * width/2 ;  
        range2 = (j-1) * width/2 + 1 + 448 : (j-1) * width /2 + width + 448;
        subFrame = frameF(range1,:);
        s = imdct(subFrame);
        WinSubFrameT = s .* Win;
        frameT(range2,:) = frameT(range2,:) + WinSubFrameT;
    end 
end

end

function [long, short] = windows(N, winType)

if strcmp(winType, 'KBD')
    long = kbd(N, 6);
    short = kbd(N/8, 4);
elseif strcmp(winType, 'SIN')
    long = sinW(N);
    short = sinW(N/8);
end
long = [long long];
short = [short short];
end

function kbdWin = kbd(N,a)
w = kaiser(N/2+1, a*pi);
w2 = cumsum(w(1:N/2));
kbdWin = sqrt([w2/sum(w); w2(N/2:-1:1)/sum(w)]);
end

function sinWin = sinW(N)
         sinWin = zeros(N,1);
        for i = 0:N-1
            sinWin(i+1) = sin((pi/N)*(i+1/2));
        end
end

function y = imdct(x)
% Marios Athineos, marios@ee.columbia.edu
% http://www.ee.columbia.edu/~marios/
% Copyright (c) 2002 by Columbia University.
% All rights reserved.

[flen, fnum] = size(x);
% Make column if it's a single row
if (flen == 1)
    x = x(:);
    flen = fnum;
    fnum = 1;
end

% We need these for furmulas below
N = flen;
M = N / 2;
twoN = 2 * N;
sqrtN = sqrt(twoN);

% We need this twice so keep it around
t = (0:(M - 1)).';
w = diag(sparse(exp(-1i*2*pi*(t + 1 / 8)/twoN)));

% Pre-twiddle
t = (0:(M - 1)).';
c = x(2*t+1,:) + 1i * x(N-1-2*t+1,:);
c = (0.5 * w) * c;

% FFT for N/2 points only !!!
c = fft(c, M);

% Post-twiddle
c = ((8 / sqrtN) * w) * c;

% Preallocate rotation matrix
rot = zeros(twoN, fnum);

% Sort
t = (0:(M - 1)).';
rot(2*t+1,:) = real(c(t+1,:));
rot(N+2*t+1,:) = imag(c(t+1,:));
t = (1:2:(twoN - 1)).';
rot(t+1,:) = -rot(twoN-1-t+1,:);

% Shift
t = (0:(3 * M - 1)).';
y(t+1,:) = rot(t+M+1,:);
t = (3 * M:(twoN - 1)).';
y(t+1,:) = - rot(t-3*M+1,:);
end