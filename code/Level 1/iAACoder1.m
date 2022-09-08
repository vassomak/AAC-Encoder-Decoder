function x = iAACoder1(AACSeq1, fNameOut)

    frameWidth = 2048;
    [K, ~] = size(AACSeq1);
    M = (K + 1) * frameWidth /2;
    decoded = zeros(M,2);
    
for j = 1:K
    frameF = [AACSeq1(j).chl.frameF AACSeq1(j).chr.frameF];
    frameT = iFilterbank(frameF, AACSeq1(j).frameType, AACSeq1(j).winType);
    range = (j - 1) * frameWidth/2 + 1 : (j + 1) * frameWidth/2;
    decoded(range,:) = decoded(range,:) + frameT;
end

% Remove padding
decoded(1:frameWidth/2,:)=[];
decoded(end-frameWidth/2+1:end,:)=[];
fs=48000;
audiowrite(fNameOut, decoded, fs);

if nargout == 1
    x = decoded;
end

end

