function x = iAACoder2(AACSeq2, fNameOut)

 frameWidth = 2048;
 [K, ~] = size(AACSeq2);
 M = (K + 1) * frameWidth/2;
 decoded = zeros(M,2);
 
 for j = 1:K
    frameFout1 = iTNS(AACSeq2(j).chl.frameF, AACSeq2(j).frameType, AACSeq2(j).chl.TNScoeffs);
    frameFout2 = iTNS(AACSeq2(j).chr.frameF, AACSeq2(j).frameType, AACSeq2(j).chr.TNScoeffs);
    if strcmp(AACSeq2(j).frameType, 'ESH')
        frameFout1 = reshape(frameFout1, [frameWidth/2,1]);
        frameFout2 = reshape(frameFout2, [frameWidth/2,1]);
    end
    frameF = [frameFout1 frameFout2];
    frameT = iFilterbank(frameF, AACSeq2(j).frameType, AACSeq2(j).winType);
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