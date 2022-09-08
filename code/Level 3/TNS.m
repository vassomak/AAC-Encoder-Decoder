function [frameFout, TNScoeffs] = TNS(frameFin, frameType)

[m,n] = size(frameFin);
TNScoeffs = zeros(4,n);
frameFout = zeros(m,n);

%% Step 1: normalise MDCT
Xw = normaliseMDCT(frameFin, frameType);

p=4;
step = 0.1;
numberOfBits = 4;

for i = 1:n
%% Step 2: Linear prediction filter coefficients Ra = r
a = -lpc(Xw(:,i), p);
a = a(2:end)';

%% Step 3: Uniform Quantise
a = uniformQuantize(a, step, numberOfBits);

% Assert that the inverse filter is stable 
a = [1; -a];
a = assertInvertible(a);
TNScoeffs(:,i) = -a(2:end);

%% Step 4: Apply fir filter
frameFout(:,i) = filter(a, 1, frameFin(:,i));
end
end

function Xw = normaliseMDCT(frameFin, frameType)
load('TableB219.mat')
N = 1024;

if ~strcmp(frameType, 'ESH')
    Nb = length(B219a);
    b = B219a(:,2);
else
    Nb = length(B219b);
    b = B219b(:,2);
    frameFin = reshape(frameFin, [1024,1]);
end

P = zeros(Nb,1);
for j = 1 : Nb - 1 
X = frameFin(b(j)+1 : b(j+1));
P(j) = sum(X.^2);
end
P(j+1) = sum(frameFin(Nb:end).^2);

Sw = zeros(N,1);
for j = 1 : Nb - 1
    k = b(j): b(j+1)-1;
    Sw(k+1) = sqrt(P(j));
end
Sw(k(end)+1 : N) = sqrt(P(Nb));

for k = 1023 : -1 : 1
    Sw(k) = (Sw(k) + Sw(k+1))/2;
end

for k = 2 : 1024
    Sw(k) = (Sw(k) + Sw(k-1))/2;
end

Xw = frameFin ./ Sw;
if strcmp(frameType, 'ESH')
    Xw = reshape(Xw, [128, 8]);
end
end

function a = uniformQuantize(signal, step, numberOfBits)

L = 2^numberOfBits;
n = -L/2 + 1 : L/2;
n = n.*step;
a = round(10.*signal)/10;
a(a<n(1)) = n(1);
a(a>n(L)) = n(L);
end

function a = assertInvertible(a)
r = roots(a);
e = 0.0001;

r(r<-1) = -1 + e;
r(r>1) = 1 - e;
a = poly(r);
end

