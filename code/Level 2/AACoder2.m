function AACSeq2 = AACoder2(fNameIn)
y = audioread(fNameIn);
frameWidth = 2048;

% zero padding
signal = [zeros(frameWidth/2,2); y; zeros(frameWidth/2,2)];
b = (frameWidth/2) * ceil(length(signal) / (frameWidth/2)); % next multiple of 1024 >= length(signal)
signal = [signal; zeros(b-length(signal),2)];

overlap = 0.5;
[N,~] = size(signal);
K = (1 / overlap) * (N / frameWidth - 1);

winTypes = {'KBD', 'SIN'};
r = randi([1, 2], 1); 
winType = winTypes(r);  % Use KBD or SIN randomly
fprintf('Chosen window: %s\n\n', cell2mat(winType));

AACSeq2 = struct('frameType', {}, 'winType', {}, 'chl', struct('TNScoeffs', {}, 'frameF', {}),...
    'chr', struct('TNScoeffs', {}, 'frameF', {}));
    
prevFrameType = 'OLS';

for i = 0:K-1
    frameT = slice(i,signal,overlap,frameWidth);
    nextFrameT = slice(i+1,signal,overlap,frameWidth);
    
    prevFrameType = SSC(frameT, nextFrameT, prevFrameType);
    AACSeq2(i+1).frameType = prevFrameType;
    AACSeq2(i+1).winType = winType;
    
    frameF = filterbank(frameT, AACSeq2(i+1).frameType, AACSeq2(i+1).winType);
    
    if  strcmp(AACSeq2(i+1).frameType, 'ESH')
         frameFin1 = reshape(frameF(:,1), [128,8]);
         frameFin2 = reshape(frameF(:,2), [128,8]);
    
    else
        frameFin1 = frameF(:,1);
        frameFin2 = frameF(:,2);
    end
    
    [AACSeq2(i+1).chl.frameF, AACSeq2(i+1).chl.TNScoeffs] = TNS(frameFin1, AACSeq2(i+1).frameType);
    [AACSeq2(i+1).chr.frameF, AACSeq2(i+1).chr.TNScoeffs] = TNS(frameFin2, AACSeq2(i+1).frameType);

end
end

function frame = slice(i, signal, overlap, frameWidth)
startF = frameWidth * overlap * i + 1;
endF = startF + frameWidth - 1;
frame = signal(startF:endF,:);
end