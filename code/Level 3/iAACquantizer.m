function frameF = iAACquantizer(S, sfc, G, frameType)

if(strcmp(frameType, 'ESH'))
    S = reshape(S, [128, 8]);
    Nb = 42;
else 
    Nb=69;
end

numOfCol = size(S,2);
frameF = zeros(size(S));
for i = 1:numOfCol
    for b = 0 : Nb-1
        a = sum(sfc(b+1:-1:1,i));
    frameF(:,i) = sign(S(:,i)) .* (abs(S(:,i)) .^ (4/3)) .* 2 .^ (a/4);
    end
end