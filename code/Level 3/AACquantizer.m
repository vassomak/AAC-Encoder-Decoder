function [S, sfc, G] = AACquantizer(frameF, frameType, SMR)

n = size(frameF,2);

load('TableB219.mat')

if ~strcmp(frameType, 'ESH')
 wLow = B219a(:,2);
 wHigh = B219a(:,3);
else 
    wLow = B219b(:,2);
    wHigh = B219b(:,3);
end

Nb = length(wLow);
P = zeros(Nb,1);
S = zeros(size(frameF));
Xq = zeros(size(frameF));
Pe = zeros(Nb,1);
G = zeros(1,n);
sfc = zeros(Nb,n);

%Acoustics Threshold
for i = 1:n
for b = 0 : Nb - 1 
    w1 = wLow(b+1);
    w2 = wHigh(b+1);
    X = frameF(w1+1 : w2+1,i);
    P(b+1) = sum(X.^2);
end
T = P ./ SMR(:,i);


%% Step 1: Scale Factor Gain approximation
MQ = 8191;
a = floor((16/3) * log2(max(frameF(:,i)).^(3/4)/MQ));
a = ones(Nb,1) * a;
%% Step 2: Final Scale Factor Gain
MagicNumber = 0.4054;

for b = 0 : Nb - 1
    
    w1 = wLow(b+1);
    w2 = wHigh(b+1);
    
    while 1  
    S(:,i) = sign(frameF(:,i)) .* fix((abs(frameF(:,i)) .* (2 .^ (-a(b+1)/4))) .^ (3/4) + MagicNumber);
    Xq(:,i) = sign(S(:,i)) .* (abs(S(:,i)) .^ (4/3)) .* (2 .^ (a(b+1)/4));
    Pe(b+1) = sum((frameF(w1+1 : w2+1,i) - Xq(w1+1 : w2+1,i)).^2);
    
    if (Pe(b+1) >= T(b+1) || max(abs(diff(a)))>60)
            break;
    else
        a(b+1) = a(b+1) + 1;
    end
    
    end
end

% Global Gain 
G(i) = a(1);

%Scale Factor of each band
sfc(:,i) = [a(1); diff(a)];
end
S = reshape(S, [1024,1]);

    end

