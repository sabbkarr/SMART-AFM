% initialization
params.imgDim              = 512;      % Image dimensions (px)
params.imageSize           = 1000;     % Image dimensions (nm)

% preprocessing
params.maxDelta          = double(15e-9); % nm clipping range
params.upsampleFactor    = 4;             % upscaling
params.gaussianSigma     = 2;             % smoothing Gaussian sigma
params.fftHighPassCutoff = 4;             % high-pass cutoff (px)
params.sharpeningRadius  = 4;             % standard deviation of the Gaussian lowpass filter
params.sharpeningAmount  = 4;             % strength of the sharpening effect
params.medianFilterSize  = [7, 7];        % median filter window

params.pixelLength = params.imageSize/(params.imgDim * params.upsampleFactor);

% segmentation:
params.seg.window       = 1 + 10 * floor(30/params.pixelLength);
params.seg.neighborhood = [params.seg.window, params.seg.window]; % adaptive threshold window
params.seg.sensitivity  = 0.1;              % adaptthresh sensitivity
params.seg.minAreaFrac  = 0.0002;           % min area fraction

% refinement:
params.ref.minAreaFrac       = 0.0002;           % min area fraction
params.ref.intensityThreshLc = 30;               % 8-bit intensity cutoff
params.ref.intensityThreshGl = 25;               % 8-bit intensity cutoff
params.ref.localExpand       = 0.3;              % local bbox expansion
params.ref.seE1              = strel('disk', 4); % Morphological cleanup
params.ref.seD1              = strel('disk', 8); % Morphological cleanup
params.ref.se2               = strel('disk', 6); % Morphological cleanup

% grouping
params.group.edgeMargin          = 5;         % px
params.group.seS1                = strel('disk', 5);

% step 2
params.group.coarse.minARInit    = 1.5;       % initial screening of aspect ratio
params.group.coarse.minAR        = 1.9;       % screening of aspect ratio 
params.group.coarse.solidity     = 0.85;      % between 0 and 1
params.group.coarse.lengthSize   = [10, 300]; % nm
params.group.coarse.heightSize   = [10, 40];  % nm

% step 3
params.group.lateral.feretThresh = 10;        % px

% step 4
params.group.vertical.stepThresh = 0.33;      % relative height fraction
params.group.vertical.windowsize = 30;        % step detecting kernel window
params.group.vertical.margin     = 0.3;       % exclude margins from step-detectors

% measurement
params.groupList = {'isolated', 'vertical', 'lateral'};  % Choose any subset of: 'lateral','vertical','isolated'