function AACSeq1 = AACoder1(fNameIn)

frameWidth = 2048;
y = audioread(fNameIn);

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

AACSeq1 = struct('frameType', {}, 'winType', {}, 'chl', struct('frameF', {}), 'chr', struct('frameF', {}));
    
prevFrameType = 'OLS';
for i = 0:K-1
    
    frameT = slice(i,signal,overlap,frameWidth);
    nextFrameT = slice(i+1,signal,overlap,frameWidth);
    
    prevFrameType = SSC(frameT, nextFrameT, prevFrameType);
    AACSeq1(i+1).frameType = prevFrameType;
    AACSeq1(i+1).winType = winType;
    frameF = filterbank(frameT, AACSeq1(i+1).frameType, AACSeq1(i+1).winType);
    AACSeq1(i+1).chl.frameF = frameF(:,1);
    AACSeq1(i+1).chr.frameF = frameF(:,2);
    
end
end

function frame = slice(i, signal, overlap, frameWidth)
startF = frameWidth * overlap * i + 1;
endF = startF + frameWidth - 1;
frame = signal(startF:endF,:);
end