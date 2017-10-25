clc; clear all;
addpath('dependencies/OpticalFlow/mex');
addpath('dependencies/gco-v3.0/matlab');
addpath('dependencies/IAT_v0.9.2');
iat_setup

%% Load images:
%    Assume that photo with largest magnification (closest focal plane) 
%    comes first in the stack.

N = 12;  % Specify total number of images
I = cell(N,1);
for i = 1:N

    % Read photos without alteration
%     I{i} = im2double(imread( ...
%                     sprintf('photos/external/largemotion/%i.jpg',i)));

    % Read photos and downsample
    I{i} = imresize(im2double(imread( ...
                sprintf('photos/personal/leaves/%i.jpg',i))), 0.2);

end
[rows,cols,chans] = size(I{1});
mkdir('results');

%% Apply affine realignment:
%    Apply affine alignemnt from image i+1 to image i to compensate for
%    magnification changes or rolling-shutter effects. Improves alignment
%    in low-texture regions.

% Compute affine transforms between subsequent images
affineTransforms = cell(N,1);
for i = 2:N
    
    % we want to compute transform from i --> i-1
    baseImage = rgb2gray(I{i-1});
    refImage = rgb2gray(I{i});
    
    % define parameters for inverse compositional alignment
    ICparams = struct('iterations', 50, ...
                      'levels', 1, ...
                      'transform', 'affine', ...
                      'initwarp', [eye(2) zeros(2,1)]);
    
    % compute affine transformation using Lucas Kanade IC algo
    ICwarp = iat_LucasKanadeIC(refImage, baseImage, ICparams);
    affineTransforms{i} = [ICwarp; 0 0 1];
    
end

% Apply the affine transformation(s) to each image
Iaff = cell(N,1);
Iaff{1} = I{1};
tform = eye(3);
for i = 2:N
    
    % propogate transform to current image
    tform = tform * affineTransforms{i};
    
    % apply transform
    [wImage,~] = iat_inverse_warping(I{i}, tform(1:2,:), ....
                                         'affine', 1:cols, 1:rows);
    
    Iaff{i} = wImage;
    
end

disp('Done applying affine realignment...');

mkdir('results','affine');
for i = 1:N
    imwrite(Iaff{i}, sprintf('results/affine/%02i.jpg',i));
end



%% Compute optical flow:
%    Align images 2...N to image 1

% Define F{i,j} to represent flow field from image i to image j
F = cell(N,N);
for i = 1:N
    F{i,i} = struct('x', zeros(rows,cols), ...
                    'y', zeros(rows,cols));
end

% Set parameters for optical flow computation
alpha = 0.03;               % 0.03
ratio = 0.75;               % 0.85
minWidth = 20;              % 20
nOuterFPIterations = 4;     % 4
nInnerFPIterations = 1;     % 1
nSORIterations = 40;        % 40

ofParams = [alpha, ...
            ratio, ...
            minWidth, ...
            nOuterFPIterations, ...
            nInnerFPIterations, ... 
            nSORIterations];

% Compute optical flow between consecutive frames
for i = N:-1:2
    
    % check time cost
    tic;
    [Vx,Vy,~] = Coarse2FineTwoFrames(Iaff{i-1}, Iaff{i}, ofParams);
    toc
    
    % save in struct format
    F{i,i-1} = struct('x', Vx, 'y', Vy);
    
end

% Recursively compute flow relative to image 1
for i = 2:N
    F{i,1} = concatFlow(F{i,i-1}, F{i-1,1}); 
end

disp('Done generating flow fields...');


%% Apply flow field alignment (aligned to image 1)
Ialigned = cell(N,1);
mkdir('results','aligned');
for i = 1:N
    Ialigned{i} = flowWarp(F{i,1}, Iaff{i});
    imwrite(Ialigned{i}, sprintf('results/aligned/%02i.jpg',i));
end


%% Compute all-in-focus image using MRF formulation:
%    Solved using alpha-expansion graph cuts algorithm. This formulation
%    allows us to apply a Total Variation (TV) prior which achieves
%    smooth results.

% Define model parameters
lambda = 0.13;      % regularization for smoothness term
mu = 3;             % standard deviation of Gaussian neighborhood

% Define scale factors (needed because GC solver requires integer values)
data_scale = 256;
smooth_scale = round(lambda * data_scale);

% Define DATA COST as sum of gradient magnitude over local neighborhood
psf = fspecial('gaussian', [2*mu-1, 2*mu-1], mu);
otf = psf2otf(psf, [rows,cols]);
Edata = int32(zeros([N rows*cols]));
for i = 1:N
    
    % Get gradient magnitude using Sobel operator
    [Gmag,~] = imgradient(rgb2gray(Ialigned{i}),'sobel');
    
    % Sum over local neighborhood by convolving with Gaussian
    Eimg = ifft2( fft2(exp(Gmag)) .* otf );
    
    % Say actual cost values are inverse of computed weight (range [0.1768,1])
    % (rough testing shows that this is better than using negative of weight)
    Edata(i,:) = int32(data_scale * 1./Eimg(:)');
    
end

% Define SMOOTHNESS COST as linear TV
[mX,mY] = meshgrid(1:N,1:N);
Esmooth = int32(smooth_scale * abs(mX-mY));

% Define 4-connected NEIGHBORS using giant sparse matrix
matSize = [rows,cols];
sparseI = zeros(rows*cols*4,1);
sparseJ = zeros(rows*cols*4,1);
counter = 1;
for i = 1:rows
    for j = 1:cols
        
        % Use sub2ind to get index of current pixel
        currentIdx = sub2ind(matSize,i,j);
        
        % Get indices of neighbors (handling all edge cases)
        neighborIdx = zeros(4,1);
        if i ~= 1
            neighborIdx(1) = sub2ind(matSize,i-1,j);
        end
        if i ~= rows
            neighborIdx(2) = sub2ind(matSize,i+1,j);
        end
        if j ~= 1
            neighborIdx(3) = sub2ind(matSize,i,j-1);
        end
        if j ~= cols
            neighborIdx(4) = sub2ind(matSize,i,j+1);
        end
            
        % Add neighbors into vectors used to construct sparse matrix
        neighborIdx = neighborIdx(neighborIdx>0);
        for k = 1:length(neighborIdx)
            sparseI(counter) = currentIdx;      % col idx (i)
            sparseJ(counter) = neighborIdx(k);  % row idx (j)
            counter = counter + 1;
        end
    
    end
end
    
% Shrink vectors appropriately
sparseI = sparseI(1:counter-1);
sparseJ = sparseJ(1:counter-1);
sparseV = ones(counter-1,1);      % all neighbor weights (v) will = 1

% Construct sparse neighbor matrix
neighbors = sparse(sparseI, sparseJ, sparseV);

disp('Done creating model for MRF formulation...');

% Solve model
MRFmodel = GCO_Create(rows*cols, N);
GCO_SetDataCost(MRFmodel, Edata);
GCO_SetSmoothCost(MRFmodel, Esmooth);
GCO_SetNeighbors(MRFmodel, neighbors);
% GCO_Swap(MRFmodel);
GCO_Expansion(MRFmodel);
gc_result = GCO_GetLabeling(MRFmodel);
gc_result = reshape(gc_result,[rows cols]);

figure
imagesc(gc_result)

% Construct all-in-focus image
aifImage = zeros(size(Ialigned{1}));
for i = 1:rows
    for j = 1:cols
        aifImage(i,j,:) = Ialigned{gc_result(i,j)}(i,j,:);
    end
end

figure
imshow(aifImage);
imwrite(aifImage, 'results/all_in_focus.png');
disp('Done generating all-in-focus-image...');


%% Generate map data for use in calibration optimization
%     Get difference map, blur map, and confidence map.

% Blur AIF image with disk point-spread functions of various sizes
testRadii = 0:0.25:6.5;
numPSFs = length(testRadii);
blurStack = cell(numPSFs,1);

aifGray = rgb2gray(aifImage);
blurStack{1} = aifGray;
for i = 2:numPSFs
    diskPSF = fspecial('disk',testRadii(i));
    blurStack{i} = conv2(aifGray, diskPSF, 'same');
end
    
% Compute all maps
mu_d = 11;                              % std dev of Gaussian neighborhood
psf_d = fspecial('gaussian',[2*mu_d-1 2*mu_d-1], mu_d); 
otf_d = psf2otf(psf_d,[rows,cols]);     % precompute Gaussian filter
        
R = cell(N,1);              % radius maps
tform = eye(3);             % initial affine transformation (img1)
for i = [1 N]
    
    % Difference map (D)
    Di = zeros([rows,cols,numPSFs]);
    for r = 1:numPSFs
        absDiff = abs(rgb2gray(Ialigned{i}) - blurStack{r});
        Di(:,:,r) = ifft2(fft2(absDiff) .* otf_d);
    end
    
    % Get delta (inverse magnification)
    if i == 1
        delta = 1.0;
    else
        % Propogate transform to current image
        tform = tform * affineTransforms{i};
        
        % Estimate scaling as average of sx and sy in transform
        delta = 1.0 / ((tform(1,1) + tform(2,2)) / 2);
    end
    
    % Get PSF radius corresonding to best match, i.e. radius map (R)
    [~,argminDi] = min(Di,[],3);
    R{i} = delta * testRadii(argminDi);
    
    
end

disp('Done generating maps...');

%% Get relative depth map estimate (kernel radii)
%     Use computed blur maps to generate a simple depth estimate. Values
%     in this estimate actually represent the blur kernel radius at each
%     pixel relative to the nearest focal plane. 

radiusMap = R{1} + (max(R{end}(:)) - R{end});

figure;
imagesc(radiusMap);

%% Perform refocusing

mkdir('results','refocus');
disp('Performing refocusing...');
for relativeDepth = 0:14
  
    refocusedImage = zeros(size(aifImage));
    for i = 15:rows-15
        for j=15:cols-15

            % Create blur kernel with appropriate radius
            blurRadius = abs(radiusMap(i,j) - relativeDepth);
            if blurRadius == 0
                refocusedImage(i,j,:) = aifImage(i,j,:);
                continue;
            end        
            blurKernel = fspecial('disk', blurRadius);

            % Perform spatially-variant blurring
            for c = 1:chans

                [kRows,kCols] = size(blurKernel);                
                aifROI = aifImage(i-(kRows-1)/2:i+(kRows-1)/2, ...
                                  j-(kCols-1)/2:j+(kCols-1)/2, c);

                refocusedImage(i,j,c) = sum(sum(aifROI .* blurKernel));

            end
        end
    end

    % Show/write results
    imwrite(refocusedImage, ...
        sprintf('results/refocus/refocus_%.2f.png', relativeDepth));
    
    disp(sprintf('Refocused at relative depth=%.2f', relativeDepth));

end
