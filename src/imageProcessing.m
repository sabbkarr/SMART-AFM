% Runs the full SMART‑AFM image processing pipeline:
%   1) (segmentation pre-processing)  preprocess raw height map
%   2) (segmentation)                 segment into a binary mask
%   3) (segmentation post-processing) refine the mask to remove false positives
%       3-a) Object‐wise re‐segmentation and morphological cleanup.
%       3-b) Global cleanup: intensity‐based removal, morphology, border‐clearing.
%
% Inputs:
%   imgRaw – 1×N cell array of AFM height values in nanometers (image) matrices (double)   
%   params – struct with fields
%   names  – image names
%
% Outputs:
%   imgPre      – pre‑processed, upscaled image (uint8 matrix)
%   imgSeg      – logical segmentation mask (boolean matrix)
%   imgRawUp    – upscaled image in nm (double matrix)
%   back_height – average background height of an image (double scalar)
%   noisee      – image noise (double scalar)
%   contrastt   – image contrast (double scalar)

function [imgPre, imgSeg, imgRawUp, back_height, noisee, contrastt] = imageProcessing(imgRaw, names, params)
    [imgPre, imgSeg, imgRawUp, back_height, noisee, contrastt] = deal({}, {}, {}, [], [], []); % Init

    outPre = fullfile(pwd, 'outputs/preprocessed_images'); % Directory to save preprocessed images
    if ~exist(outPre, 'dir'), mkdir(outPre); end       % Create directory

    outSeg = fullfile(pwd, 'outputs/segmented_images'); % Directory to save preprocessed images
    if ~exist(outSeg, 'dir'), mkdir(outSeg); end       % Create directory
    
    for i = 1:numel(imgRaw)
        [imgPre{end+1}, imgRawUp{end+1}] = preprocessImage(imgRaw{i}, params); % Preprocess
        imgBin = segmentImage(imgPre{i}, params.seg);                              % Segment
        imgSegRef = refineLocalSeg(imgPre{i}, imgBin, params.ref);                 % Postprocess
        imgSeg{end + 1} = refineGlobalSeg(imgSegRef, imgPre{i}, params.ref);

        back_height(end + 1) = mean(imgRawUp{i}(imgSeg{i} == 0));
        noisee(end + 1)      =  std(imgRawUp{i}(imgSeg{i} == 0));
        contrastt(end + 1)   = mean(imgRawUp{i}(imgSeg{i} == 1)) - mean(imgRawUp{i}(imgSeg{i} == 0));

        imwrite(uint8(imgPre{i}), fullfile(outPre, [names{i} '_preprocessed.png'])); % Save Preprocessed Image
        imwrite(uint8(imgSeg{i} .* 255), fullfile(outSeg, [names{i} '_segmented.png'])); % Save Segmented Image
    end 
end

%% Subfunction: preprocessImage
% Clip spikes, rescale, upscale & preprocess (smoothing, sharpening, FFT
%                               high-pass, adjusting and median-filtering)
%
% Steps:
%   a) Upsample by upsampleFactor
%   b) Clip values > min(imgRaw)+maxDelta to that threshold
%   b) Normalize to [0,1], convert to uint8 [0,255]
%   d) Preprocess image
%       i)   Gaussian blur (sigma = gaussianSigma)
%       ii)  Adjust image intensity values
%       iii) FFT high‑pass (cutoff = fftHighPassCutoff)
%       iv)  Sharpening - Unsharp Masking
%       iv)  Edge‑preserving median filter (size = medianFilterSize)

function [imgPre, imgRawUp] = preprocessImage(imgRaw, params)
    % a) Upsample
    imgUp = imresize(imgRaw, params.upsampleFactor);

    % b) Clip height spikes
    hMin = min(imgUp(:));
    imgRawUp = min(imgUp, hMin + params.maxDelta);

    % c) Normalize & convert--this is the *raw upscaled* image.
    img8 = uint8(255 * mat2gray(imgRawUp));
    
    % d) Preprocess image
    % d-i) Gaussian blur
    imgGauss = imgaussfilt(img8, params.gaussianSigma);

    % d-ii) Adjust image intensity values
    imgAdj = imadjust(imgGauss);

    % d-iii) FFT high‑pass
    imgFFT = log(abs(fftshift(fft2(imgAdj))));
    [M, N] = size(imgFFT);
    [X, Y] = meshgrid(1:N, 1:M);
    D = (X - N/2).^2 + (Y - M/2).^2;  
    newFFT = zeros(size(imgFFT));
    newFFT(D < params.fftHighPassCutoff) = imgFFT(D < params.fftHighPassCutoff);
    highPass = uint8(imgAdj - uint8(abs(ifft2(ifftshift(exp(newFFT))))));

    % d-iv) Image sharpening
    imgSharp = imsharpen(highPass,'radius', params.sharpeningRadius, 'Amount', params.sharpeningAmount);
    
    % d-v) Edge‑preserving median filter
    imgPre = medfilt2(imgSharp, params.medianFilterSize);
end

%% Subfunction: segmentImage
% Adaptive thresholding to get binary mask.
%
% Uses adaptthresh + imbinarize + small object removal:

function imgBin = segmentImage(imgPre, segParams)
    T = adaptthresh(imgPre, segParams.sensitivity, ...
                    'NeighborhoodSize', segParams.neighborhood);

    Seg = imbinarize(imgPre, T);
    imgBin = bwareaopen(Seg, floor(numel(Seg) * segParams.minAreaFrac));
end

%% Subfunction: Local re‐segmentation (Stage 1 Post-Processing)
% Object‐wise re‐segmentation and morphological cleanup.
%
%   Inputs:
%     imgPre     – preprocessed (grayscale) image (H×W), type uint8 or double
%     imgBin     – initial binary mask (H×W), logical or numeric
%     refParams  – struct with parameters
%
%   Output:
%     imgSegRef  – binary mask (H×W) after local re‐segmentation + cleanup
%
%   Description:
%     1. For each connected component in imgSegInit:
%         a. Expand its bounding box by refParams.localExpand (using expandBox).
%         b. Within that expanded region, replace all pixels NOT in the “main object”
%            with the mean background intensity.
%         c. Compute a 2‐level quantization via multithresh → imquantize:
%            keep only sub‐regions whose median intensity ≥ refParams.intensityThresh.
%         d. Paint all surviving pixels back into a global “newSeg” mask.
%     2. Perform morphological cleanup:
%         a. Erosion
%         b. Removal of small objects
%         c. Dilation

function imgSegRef = refineLocalSeg(imgPre, imgBin, refParams)
    % Initialization
    propsFirst = regionprops(imgBin, 'PixelIdxList', 'BoundingBox'); % Extract objects
    imgSegRef = zeros(size(imgBin));
    imgPre = double(imgPre);

    % Loop over objects
    for k = 1:numel(propsFirst)
        % Extract bounding box, crop patches and masks for kth object
        bbox = propsFirst(k).BoundingBox;
        exX = round(max(1,                     bbox(1) -     refParams.localExpand * bbox(3)));
        exY = round(max(1,                     bbox(2) -     refParams.localExpand * bbox(4)));
        exW = round(min(size(imgPre, 2) - exX, bbox(3) + 2 * refParams.localExpand * bbox(3)));
        exH = round(min(size(imgPre, 1) - exY, bbox(4) + 2 * refParams.localExpand * bbox(4)));
        
        regPre = imgPre(exY:(exY + exH - 1), exX:(exX + exW - 1));
        regMsk = imgBin(exY:(exY + exH - 1), exX:(exX + exW - 1));

        mainObjMask = zeros(size(imgBin));
        mainObjMask(propsFirst(k).PixelIdxList) = 1;
        regMainObjMsk = mainObjMask(exY:(exY + exH - 1), exX:(exX + exW - 1));

        meanBkg = mean(regPre(regMsk == 0)); % mean background intensity
        objPre = regPre; 
        objPre(regMainObjMsk == 0) = uint8(meanBkg); % mask out other objects

        % 2‐level quantization → objSeg ∈ {0,1,2}
        thresh = multithresh(objPre, 2);
        objSeg = imquantize(objPre, thresh) - 1;

        no2s = rem(objSeg, 2);
        all2s = zeros(size(objSeg));
        all2s(objSeg == 2) = 1;

        % Remove “outer” regions (where rem(objSeg, 2) == 1) if their median < intensityThresh
        outer = regionprops(no2s, 'PixelIdxList');
        for m = 1:numel(outer)
            if median(objPre(outer(m).PixelIdxList)) < refParams.intensityThreshLc
                objSeg(outer(m).PixelIdxList) = 0;
            end
        end

        % Remove “inner” regions (where objSeg==2) if their median < intensityThresh
        inner = regionprops(all2s, 'PixelIdxList');
        for m = 1:numel(inner)
            if median(double(objPre(inner(m).PixelIdxList))) < refParams.intensityThreshLc
                objSeg(inner(m).PixelIdxList) = 0;
            end
        end
        
        % Paint surviving pixels (objSeg > 0) back into imgSegRef
        [rr, cc] = find(objSeg > 0);
        for m = 1:numel(rr)
            imgSegRef(rr(m) + exY, cc(m) + exX) = 1;
        end
    end
    
    % Morphological cleanup
    imgSegRef = imerode(imgSegRef, refParams.seE1);
    imgSegRef = bwareaopen(imgSegRef, ceil(numel(imgSegRef)*refParams.minAreaFrac/10), 4);
    imgSegRef = imdilate(imgSegRef, refParams.seD1);
end


%% Subfunction: Global post‐processing (Stage 2 Post-Processing)
% Global cleanup: intensity‐based removal, morphology, border‐clearing.
%
%   Inputs:
%     newSeg      – binary mask after local re‐segmentation (H×W)
%     imgPre      – preprocessed (grayscale) image (H×W)
%     refParams   – struct with parameters
%
%   Output:
%     finalSeg    – binary mask after global post‐processing
%
%   Description:
%     1. Remove any connected component in newSeg whose median intensity in imgPre < intensityThresh.
%     2. Perform: erode(6px) → bwareaopen((minAreaFrac/10)×totalPixels, 4) → dilate(6px) → 
%                   bwareaopen((minAreaFrac)×totalPixels).
%     3. Remove any component touching the image border (imclearborder).

function imgSeg = refineGlobalSeg(newSeg, imgPre, refParams)
    % Intensity‐based removal
    imgSeg = newSeg;
    props = regionprops(imgSeg, 'PixelIdxList');
    for k = 1:numel(props)
        if median(imgPre(props(k).PixelIdxList)) < refParams.intensityThreshGl
            imgSeg(props(k).PixelIdxList) = false;
            continue;
        end
    end

    % Morphological cleanup
    imgSeg = imerode(imgSeg, refParams.se2);
    imgSeg = bwareaopen(imgSeg, ceil(numel(imgSeg)*(refParams.minAreaFrac/10)), 4);
    imgSeg = imdilate(imgSeg, refParams.se2);
    imgSeg = bwareaopen(imgSeg, floor(numel(imgSeg)*refParams.minAreaFrac));
end