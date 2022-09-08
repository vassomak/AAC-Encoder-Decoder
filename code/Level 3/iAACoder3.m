function x = iAACoder3(AACSeq3, fNameOut)

frameWidth = 2048;
[K, ~] = size(AACSeq3);
M = (K + 1) * frameWidth/2;
decoded = zeros(M,2);
huffLUT = loadLUT();

for j = 1:K
    S1 = decodeHuff(AACSeq3(j).chl.stream, AACSeq3(j).chl.codebook, huffLUT);
    S2 = decodeHuff(AACSeq3(j).chr.stream, AACSeq3(j).chr.codebook, huffLUT);
    
    sfc1 = decodeHuff(AACSeq3(j).chl.sfc, 12, huffLUT);
    sfc2 = decodeHuff(AACSeq3(j).chr.sfc, 12, huffLUT);
    if strcmp(AACSeq3(j).frameType, 'ESH')
            sfc1 = reshape(sfc1, [42, 8]);
            sfc2 = reshape(sfc1, [42, 8]);
            S1 = reshape(S1, [128, 8]);
            S2 = reshape(S2, [128, 8]);
    else
            S1 = S1(:);
            S2 = S2(:);
            sfc1 = sfc1(:);
            sfc2 = sfc2(:);
    end
    
    frameF1 = iAACquantizer(S1, sfc1, AACSeq3(j).chl.G, AACSeq3(j).frameType);
    frameF2 = iAACquantizer(S2, sfc2, AACSeq3(j).chr.G, AACSeq3(j).frameType);
    
    frameFout1 = iTNS(frameF1, AACSeq3(j).frameType, AACSeq3(j).chl.TNScoeffs);
    frameFout2 = iTNS(frameF2, AACSeq3(j).frameType, AACSeq3(j).chr.TNScoeffs);
    
    if strcmp(AACSeq3(j).frameType, 'ESH')
        frameFout1 = reshape(frameFout1, [frameWidth/2,1]);
        frameFout2 = reshape(frameFout2, [frameWidth/2,1]);
    end
    frameF = [frameFout1 frameFout2];
    frameT = iFilterbank(frameF, AACSeq3(j).frameType, AACSeq3(j).winType);
    range = (j - 1) * frameWidth/2 + 1 : (j + 1) * frameWidth/2;
    decoded(range,:) = decoded(range,:) + frameT;
 end

 % Remove padding
decoded(1:frameWidth/2,:)=[];
decoded(end-frameWidth/2+1:end,:)=[];

%Normalise to avoid clipping
for j = 1:2
    decoded(:,j) = decoded(:,j)./(max(abs(decoded(:,j))));
end

fs=48000;
audiowrite(fNameOut, decoded, fs);

if nargout == 1
    x = decoded;
end

end