function frameF = filterbank(frameT, frameType, winType)


%% Apply Windows and MDCT
N = length(frameT);
[long, short] = windows(N, winType);

frameF = zeros(N/2,2);

if strcmp(frameType, 'OLS')
    Win = long;
    frameWin = frameT .* Win;
    frameF = mdct(frameWin);
    
elseif strcmp(frameType, 'LSS')
    Win = [long(1:1024,:); ones(448, 2); short(129:256,:); zeros(448,2)];
    frameWin = frameT .* Win; 
    frameF = mdct(frameWin);
    
elseif strcmp(frameType, 'LPS')
    Win = [zeros(448,2); short(1:128,:); ones(448,2); long(1025:2048,:)];
    frameWin = frameT .* Win;
    frameF = mdct(frameWin);
    
elseif strcmp(frameType, 'ESH') 
    Win = short;
    width = N/8;   
   
   for j=1:8
       range1 = (j-1) * width/2 + 1 + 448 : (j-1) * width /2 + width + 448;
       range2 = (j-1) * width/2  + 1 : j * width/2 ;  
       subFrame = frameT(range1,:);   
       frameWin = subFrame .* Win;
       frameF(range2,:) = mdct(frameWin);
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
            sinWin(i+1) = sin(pi / N *(i+1/2));
        end
end

function y = mdct(x)
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
% Make sure length is multiple of 4
if (rem(flen, 4) ~= 0)
    error('MDCT4 defined for lengths multiple of four.');
end

% We need these for furmulas below
N = flen; % Length of window
M = N / 2; % Number of coefficients
N4 = N / 4; % Simplify the way eqs look
sqrtN = sqrt(N);

% Preallocate rotation matrix
% It would be nice to be able to do it in-place but we cannot
% cause of the prerotation.
rot = zeros(flen, fnum);

% Shift
t = (0:(N4 - 1)).';
rot(t+1,:) = -x(t+3*N4+1,:);
t = (N4:(N - 1)).';
rot(t+1,:) = x(t-N4+1,:);

% We need this twice so keep it around
t = (0:(N4 - 1)).';
w = diag(sparse(exp(-1i*2*pi*(t + 1 / 8)/N)));

% Pre-twiddle
t = (0:(N4 - 1)).';
c = (rot(2*t+1,:) - rot(N-1-2*t+1,:)) - 1i * (rot(M+2*t+1,:) - rot(M-1-2*t+1,:));
% This is a really cool Matlab trick ;)
c = 0.5 * w * c;

% FFT for N/4 points only !!!
c = fft(c, N4);

% Post-twiddle
c = (2 / sqrtN) * w * c;

% Sort
t = (0:(N4 - 1)).';
y(2*t+1,:) = real(c(t+1,:));
y(M-1-2*t+1,:) = - imag(c(t+1,:));
end
