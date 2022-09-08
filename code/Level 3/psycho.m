function SMR = psycho(frameT, frameType, frameTprev1, frameTprev2)

load('TableB219.mat')

if ~strcmp(frameType, 'ESH')
 bval = B219a(:,5);
 wLow = B219a(:,2);
 wHigh = B219a(:,3);
 qsthr = B219a(:,6);
 N = 2048;
else 
    bval = B219b(:,5);
    wLow = B219b(:,2);
    wHigh = B219b(:,3);
    qsthr = B219b(:,6);
    N = 256;
end
Nb = length(bval);
%% Step 1: Spreading Function
x = zeros(Nb, Nb);
for i = 0 : Nb-1
    for j = i+1 : Nb
x(i+1,j) = spreadingfun(i,j,bval);
    end
end

%% Step 2: Hanning Window - FFT
win = zeros(N,1);
for n = 0:N-1
win(n+1) = 0.5 - 0.5 * cos(pi*(n+0.5)/N);
end

if ~strcmp(frameType, 'ESH')
    
    sw = frameT .* win;
    [r, f] = FFT(sw);
    
    sw_1 = frameTprev1 .* win;
    [r_1, f_1] = FFT(sw_1);
    
    sw_2 = frameTprev2 .* win;
    [r_2, f_2] = FFT(sw_2);
   
else
    
frameT = reshape(frameT, [N,8]);   
frameTprev1 = reshape(frameTprev1, [N,8]);   
sw = zeros(N,8);
sw_1 = zeros(N,8);
sw_2 = zeros(N,8);

%Current Frame
subFrame = frameT(:,1);   
sw(:,1) = subFrame .* win;

%Previous Frame
subFrPrev1 = frameTprev1(:, 8);
sw_1(:,1) = subFrPrev1 .* win;

%Pre-previous Frame
subFrPrev2 = frameTprev1(:,7);
sw_2(:,1) = subFrPrev2 .* win;

for j=2:8
    
       subFrPrev2 = subFrPrev1;    %Pre-previous Frame
       sw_2(:,j) = subFrPrev2 .* win;
    
       subFrPrev1 = subFrame;       %Previous Frame
       sw_1(:,j) = subFrPrev1 .* win;
       
       subFrame = frameT(:,j);  %Current Frame
       sw(:,j) = subFrame .* win;             
end

[r, f] = FFT(sw);
[r_1, f_1] = FFT(sw_1);
[r_2, f_2] = FFT(sw_2);
end




%% Step 3: Predictions
r_pred  = 2*r_1 - r_2;
f_pred = 2*f_1 - f_2;

%% Step 4: Predictability
c = sqrt((r .* cos(f) - r_pred .* cos(f_pred)).^2 + (r .* sin(f) - r_pred .* sin(f_pred)).^2) ./ (r + abs(r_pred));

%% Step 5: Energy and (burdened) predictability
numOfCol = size(frameT,2);
e = zeros(Nb,numOfCol);
c2 = zeros(Nb,numOfCol);

for i = 1 : numOfCol
for b = 0 : Nb-1
    
    w1 = wLow(b+1);
    w2 = wHigh(b+1);
    r1 = r(w1+1:w2+1,i);
    c1 = c(w1+1:w2+1,i);
    e(b+1,i) = sum(r1.^2);
    c2(b+1,i) = sum(c1.*(r1.^2));
    
end
end
c = c2;
%% Step 6: Combine energy and predictability with the spreading function
ecb = zeros(Nb,numOfCol);
ct = zeros(Nb,numOfCol);
cb = zeros(Nb,numOfCol);
en = zeros(Nb,numOfCol);

for i = 1 : numOfCol
for b = 0 : Nb-1
    
    ecb(b+1,i) = sum(e(:,i).*x(:, b+1));
    ct(b+1,i) = sum(c(:,i).*x(:, b+1));
    
    %Normalise
    cb(b+1,i) = ct(b+1,i) ./ ecb(b+1,i);
    en(b+1,i) = ecb(b+1,i) ./ sum(x(:,b+1));
end
end
%% Step 7: Tonality Index
tb = -0.299 - 0.43 * log(cb);
margin = 0.0001;
tb(tb<0) = margin;
tb(tb>1) = 1-margin;

%% Step 8: SNR of each band
NMT = 6;
TMN = 18;
SNR = tb * TMN + (1-tb) * NMT;

%% Step 9: dB to Energy
bc = 10 .^ (-SNR/10);

%% Step 10: Threashold
nb = en .* bc;

%% Step 11: Noise Level
q_thr = eps * (N/2) * 10 .^ (qsthr./10);
npart = zeros(Nb, numOfCol);
for i = 1:numOfCol
npart(:,i) = max(nb(:,i), q_thr);
end
%% Step 12: SMR
SMR = e ./ npart;
end

function x = spreadingfun(i, j, bval)

if i>=j
    tmpx = 3*(bval(j) - bval(i+1));
else
    tmpx = 1.5*(bval(j) - bval(i+1));
end

tmpz = 8 * min((tmpx - 0.5)^2 - 2*(tmpx - 0.5),0);
tmpy = 15.811389 + 7.5*(tmpx + 0.474) - 17.5*sqrt(1 + (tmpx + 0.474)^2);

if tmpy<-100
    x = 0;
else
    x = 10^((tmpz + tmpy)/10);
end
end


function [r, f] = FFT(y)
r = zeros(size(y,1)/2, size(y,2));
f = zeros(size(y,1)/2, size(y,2));

Y = fft(y,[],1);
r = abs(Y(1:end/2,:));
f = angle(Y(1:end/2,:));
end
