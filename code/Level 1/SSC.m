function frameType = SSC(frameT, nextFrameT, prevFrameType)  

% if the previous frame type is LSS or LPS, the current frame type is OLS
if (strcmp(prevFrameType, 'LSS'))
    frameType = 'ESH';
    return;
    
elseif(strcmp(prevFrameType, 'LPS'))
    frameType = 'OLS';
    return;
end

%% Step 1: Filtering
b = [0.7548 -0.7548];
a = [1 -0.5095];
nextFrameT = filter(b, a, nextFrameT);    

%% Step 2: Split to 8 parts and compute s^2 for each channel
seg = nextFrameT(577:1600,:);        

subFrames_ch1 = reshape(seg(:,1), [128,8]);   
subFrames_ch2 = reshape(seg(:,2), [128,8]);

s2_ch1(1:8) = sum(subFrames_ch1(:,1:8).^2);
s2_ch2(1:8) = sum(subFrames_ch2(:,1:8).^2);

%% Step 3: Compute the attack values
ds1 = zeros(1,8);
ds2 = zeros(1,8);

for i=2:8
ds1(i) = s2_ch1(i)/(sum(s2_ch1(1:i-1))/i);
ds2(i) = s2_ch2(i)/(sum(s2_ch2(1:i-1))/i);
end

%% Step 4: Determine if nextFrameType is 'ESH'
for i = 2:8
    if (ds1(i)>10 && s2_ch1(i)>10^(-3))
        nextFrameType{1} = 'ESH';
        break
    else
        nextFrameType{1} = [];
        
    end
end

for i = 2:8
    if (ds2(i)>10 && s2_ch2(i)>10^(-3))
        nextFrameType{2} = 'ESH';
        break
    else
        nextFrameType{2} = [];     
    end
end

for i=1:2
if(strcmp(prevFrameType,'OLS'))
    if strcmp(nextFrameType{i},'ESH')
         frameType{i} = 'LSS';
        else
         frameType{i} = 'OLS';
    end
 

elseif(strcmp(prevFrameType, 'ESH'))
         if strcmp(nextFrameType{i},'ESH')
        frameType{i} = 'ESH';
        else
        frameType{i} = 'LPS';
         end
end
end
  
% Decide the final type of the frame using the type of each channel
if strcmp(frameType{1}, frameType{2})
    frameType = frameType{1};

elseif(any(strcmp(frameType, 'OLS')) && any(strcmp(frameType, 'LSS')))
        frameType = 'LSS';
        
elseif (any(strcmp(frameType, 'OLS')) && any(strcmp(frameType, 'ESH')))
                frameType = 'ESH';
                
elseif (any(strcmp(frameType, 'OLS')) && any(strcmp(frameType, 'LPS')))
                frameType = 'LPS';
                
elseif (any(strcmp(frameType, 'LSS')) && any(strcmp(frameType, 'ESH')))          
                frameType = 'ESH';
                
elseif (any(strcmp(frameType, 'LSS')) && any(strcmp(frameType, 'LPS')))  
                frameType = 'ESH';
                
elseif (any(strcmp(frameType, 'LPS')) && any(strcmp(frameType, 'ESH')))
                frameType = 'ESH';
end
end


