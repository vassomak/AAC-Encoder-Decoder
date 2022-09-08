function AACSeq3 = AACoder3(fNameIn, fnameAACoded)

y = audioread(fNameIn);
frameWidth = 2048;

% zero padding
signal = [zeros(frameWidth/2,2); y; zeros(frameWidth/2,2)];
b = (frameWidth/2) * ceil(length(signal) / (frameWidth/2)); % next multiple of 1024 >= length(signal)
signal = [signal; zeros(b-length(signal),2)];


overlap = 0.5;
[N,~] = size(signal);
K = (1/overlap) * ((N/frameWidth) - 1);

winTypes = {'KBD', 'SIN'};
r = randi([1, 2], 1); 
winType = winTypes(r);  % Use KBD or SIN randomly
fprintf('Chosen window: %s\n\n', cell2mat(winType));

AACSeq3 = struct('frameType', {}, 'winType', {}, ...
    'chl', struct('TNScoeffs', {}, 'T', {}, 'G', {}, 'sfc', {}, 'stream', {}, 'codebook', {}), ...
    'chr', struct('TNScoeffs', {}, 'T', {}, 'G', {}, 'sfc', {}, 'stream', {}, 'codebook', {}));

prevFrameType = 'OLS';
frameTprev1 = zeros(frameWidth,2);
frameTprev2 = zeros(frameWidth,2);

huffLUT = loadLUT();

for i = 0:K-1
    frameT = slice(i,signal,overlap,frameWidth);
    nextFrameT = slice(i+1,signal,overlap,frameWidth);
    
    prevFrameType = SSC(frameT, nextFrameT, prevFrameType);
    AACSeq3(i+1).frameType = prevFrameType;
    AACSeq3(i+1).winType = winType;
    
    frameF = filterbank(frameT, AACSeq3(i+1).frameType, AACSeq3(i+1).winType);
    
    if  strcmp(AACSeq3(i+1).frameType, 'ESH')
         frameFin1 = reshape(frameF(:,1), [128,8]);
         frameFin2 = reshape(frameF(:,2), [128,8]);
    
    else
        frameFin1 = frameF(:,1);
        frameFin2 = frameF(:,2);
    end
    
    [frameFout1, AACSeq3(i+1).chl.TNScoeffs] = TNS(frameFin1, AACSeq3(i+1).frameType);
    [frameFout2, AACSeq3(i+1).chr.TNScoeffs] = TNS(frameFin2, AACSeq3(i+1).frameType);
    
    SMR1 = psycho(frameT(:,1), AACSeq3(i+1).frameType, frameTprev1(:,1), frameTprev2(:,1));
    SMR2 = psycho(frameT(:,2), AACSeq3(i+1).frameType, frameTprev1(:,2), frameTprev2(:,2));
    
    AACSeq3(i+1).chl.T = threshold(SMR1, frameFout1);
    AACSeq3(i+1).chr.T = threshold(SMR2, frameFout2);
    
    [S1, sfc1, AACSeq3(i+1).chl.G] = AACquantizer(frameFout1, AACSeq3(i+1).frameType, SMR1);
    [S2, sfc2, AACSeq3(i+1).chr.G] = AACquantizer(frameFout2, AACSeq3(i+1).frameType, SMR2);
    
    AACSeq3(i+1).chl.sfc = encodeHuff(sfc1(:), huffLUT, 12);
    AACSeq3(i+1).chr.sfc = encodeHuff(sfc2(:), huffLUT, 12);

    [AACSeq3(i+1).chl.stream, AACSeq3(i+1).chl.codebook] = encodeHuff(S1, huffLUT); 
    [AACSeq3(i+1).chr.stream, AACSeq3(i+1).chr.codebook] = encodeHuff(S2, huffLUT);
   
    
    frameTprev2 = frameTprev1;
    frameTprev1 = frameT;
end

save(fnameAACoded, 'AACSeq3');

end

function frame = slice(i, signal, overlap, frameWidth)
startF = frameWidth * overlap * i + 1;
endF = startF + frameWidth - 1;
frame = signal(startF:endF,:);
end

function T = threshold(SMR, frameF)
    
    load('TableB219.mat')
    if size(SMR,2) == 1
        wLow = B219a(:,2);
        wHigh = B219a(:,3);
    else 
        wLow = B219b(:,2);
        wHigh = B219b(:,3);
    end
    
 [m,n] = size(SMR);
 P = zeros(length(wLow),1);
 T = zeros(m,n);
 
for i = 1 : n
for b = 0 : length(wLow) - 1 
    w1 = wLow(b+1);
    w2 = wHigh(b+1);
    X = frameF(w1+1 : w2+1,i);
    P(b+1) = sum(X.^2);
end
  T(:,i) = P ./ SMR(:,i);
end

T = reshape(T, [m*n,1]);
 
end