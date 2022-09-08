function [SNR, bitrate, compression] = demoAAC3(fNameIn, fNameOut, frameAACoded)
tic
clc
[y, fs] = audioread(fNameIn);

fprintf('Encoding:\n');
AACSeq3 = AACoder3(fNameIn, frameAACoded);


fprintf('Decoding:\n');
x = iAACoder3(AACSeq3', fNameOut);

Length = min(length(x), length(y));
y = y(1:Length,:);
x = x(1:Length,:);
noise = y-x;
SNR = snr(x,noise);

NRMSE =NRMSE_fun(x,y);
fprintf('Normal Root Mean Square Error (NRMSE) for channel 1: %g\n', NRMSE(1));
fprintf('Normal Root Mean Square Error (NRMSE) for channel 2: %g\n', NRMSE(2));
fprintf('SNR: %g dB\n', SNR);

    uncompressedData = audioinfo(fNameIn);
    bitrate_1 = uncompressedData.BitsPerSample * fs;

    compressedData = audioinfo(fNameOut);
    bitrate = compressedData.BitsPerSample * fs;
    
    compression = bitrate/bitrate_1;

    fprintf('Bitrate: %g kbps\n', bitrate/1000);
    fprintf('Compression: %g\n', compression);
    toc
end

function NRMSE = NRMSE_fun(xest,x)
% x: the original signal
% xest: x estimated

a = zeros(2048,2);
for k=1:2048
     a(k,:) = (xest(k,:) - x(k,:)).^2;
end

% Root Mean Square Error
RMSE = sqrt(sum(a/2048)); 

% Normalise
NRMSE(1) = RMSE(1)/(max(x(:,1))-min(x(:,1)));
NRMSE(2) = RMSE(2)/(max(x(:,2))-min(x(:,2)));

end
