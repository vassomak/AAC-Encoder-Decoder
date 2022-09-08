function frameFout = iTNS(frameFin, frameType, TNScoeffs)

[m,n] = size(frameFin);
frameFout = zeros(m,n);

a = [ones(1,n); -TNScoeffs];

for i = 1:n
    frameFout(:,i) = filter(1, a(:,i),frameFin(:,i));    
end
end