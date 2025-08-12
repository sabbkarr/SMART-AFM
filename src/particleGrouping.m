%% ParticleGrouping.m
% Sequentially assign each connected component (CC) to one of:
%   1) border
%   2) coarse cluster
%   3) lateral bundle
%   4) vertical stack
%   5) isolated (remaining)
%
% For each CC in the input mask, we call step1→step4 functions. Once a CC
% is classified in step n, it is not tested in later steps.
%
% Inputs:
%   imgSeg     – logical 2D array: segmentation (foreground=1)
%   imgRawUp   – double 2D array: original AFM height data (same size)
%   params     – struct with fields
%   names      – image names
%
% Output:
%   groups – struct with fields:
%       .border    – CC indices assigned as border
%       .coarse    – CC indices assigned as coarse clusters
%       .lateral   – CC indices assigned as lateral bundles
%       .vertical  – CC indices assigned as vertical stacks
%       .isolated  – CC indices assigned as isolated CNCs
%   group_stat – struct with fields:
%       .pc_<group_name>   – <group_name> particle count
%       .paf_<group_name>  – <group_name> particle area fraction

function [groups, group_stat] = particleGrouping(imgSeg, imgRawUp, params, names)
    % 1) Find all connected components once
    CC = bwconncomp(imgSeg, 8);
    N  = CC.NumObjects;
    imgSize = size(imgSeg);
    
    % If no objects, return empty groups
    if N == 0
        groups = struct('border',[], 'coarse',[], 'lateral',[], 'vertical',[], 'isolated',[]);
        return;
    end
    
    % 2) Precompute regionprops needed by multiple steps
    statsAll = regionprops(CC, imgRawUp, ...
        'PixelIdxList', 'Solidity', 'MinorAxisLength', 'MajorAxisLength', ...
        'Orientation', 'BoundingBox', 'ConvexImage', 'Image', 'Area');
    
    % 3) Prepare empty index lists
    [borderIdx, coarseIdx, lateralIdx, verticalIdx, isolatedIdx] = deal([], [], [], [], []);
    
    % 4) Loop over each CC and test sequentially
    for k = 1:N
        props_k   = statsAll(k);
        compMask  = false(imgSize);
        compMask(props_k.PixelIdxList) = true;
        
        % Step 1: border?
        if step1(compMask)
            borderIdx(end+1) = k;
            
        % Step 2: coarse cluster?
        elseif step2(props_k, compMask, params) 
            coarseIdx(end+1) = k;

        % Step 3: lateral bundle?
        elseif step3(compMask, params.group.lateral)
            lateralIdx(end+1) = k;

        % Step 4: vertical stack?
        elseif step4(props_k, compMask, imgRawUp, params.group.vertical)
            verticalIdx(end+1) = k;

        % Else: isolated
        else
            isolatedIdx(end+1) = k;
        end
    end
    
    % 5) Return groups (as row vectors)
    groups = struct('border',   borderIdx(:)', ...
        'coarse',   coarseIdx(:)', 'lateral',  lateralIdx(:)', ...
        'vertical', verticalIdx(:)', 'isolated', isolatedIdx(:)');

    
    group_names = {'border','coarse','lateral', 'vertical', 'isolated'};
    group_stat = struct();
    for k = 1:numel(group_names)
        g = group_names{k};
        idxs = groups.(g);   
        group_stat.(['pc_' g]) = numel(idxs);
        if isempty(idxs)
            group_stat.(['paf_' g]) = 0;
        else
            group_stat.(['paf_' g]) = sum([statsAll(idxs).Area]/2048/2048*100);
        end
    end

    % 6) Save segmented images with grouping boundaries
    outBase = fullfile(pwd, 'outputs');
    if ~exist(outBase, 'dir'), mkdir(outBase); end

    outGrp = fullfile(outBase, 'grouped_images');
    if ~exist(outGrp, 'dir'), mkdir(outGrp); end
    saveSegmentedBoundaries(imgRawUp, imgSeg, groups, names, outGrp);
end

%% step1: is this CC touching the image border?
function yes = step1(compMask)
    interiorMask = imclearborder(compMask);
    % If something was removed, it was touching the border
    yes = ~isequal(interiorMask, compMask);
end

%% step2: is this CC a coarse cluster?
function yes = step2(props_k, compMask, params)
    p = params.group.coarse;
    % Compute aspect ratio and solidity
    AR = props_k.MajorAxisLength / props_k.MinorAxisLength;
    solidityPerc = props_k.Solidity;
    
    % Initial screening
    if solidityPerc < p.solidity || AR < p.minARInit
        yes = true;
        return;
    end

    % Measure precise major & minor axes
    % Preprocess object
    ccRot = imrotate(compMask, 90 - props_k.Orientation);
    ccRot = imerode(ccRot, params.group.seS1);
    ccRot = imdilate(ccRot, params.group.seS1);
    ccRot = imfill(ccRot,'holes');
    
    % Measure bounding box of the processed rotated object
    propsOnRot = regionprops(ccRot, 'BoundingBox');
    if isempty(propsOnRot)
        yes = true;
        return;
    end
    
    bb = propsOnRot.BoundingBox;
    [min_nm, maj_nm] = deal(bb(3) * params.pixelLength, bb(4) * params.pixelLength);
    
    % Check against thresholds
    
    if AR >= p.minAR && ...
       min_nm >= p.heightSize(1) && min_nm <= p.heightSize(2) && ...
       maj_nm >= p.lengthSize(1) && maj_nm <= p.lengthSize(2)
        yes = false;
    else
        yes = true;
    end
end

%% step3: is this CC a lateral bundle (pocket‐based)?
function yes = step3(compMask, params)
    % Compute convex hull of this component
    hullMask = bwconvhull(compMask);
    % “Pocket” = hull minus the object
    pocket = hullMask & ~compMask;
    
    % If no “pocket” pixels, cannot be lateral bundle
    if ~any(pocket(:))
        yes = false;
        return;
    end
    
    % Label each pocket region as a connected component
    pocketStats = regionprops(pocket, 'MinFeretProperties');
    
    % Check if any pocket’s MinFeretDiameter ≥ threshold
    ferets = [pocketStats.MinFeretDiameter];
    if any(ferets >= params.feretThresh)
        yes = true;
    else
        yes = false;
    end
end

%% step4: is this CC a vertical stack (height‐profile step)?
function yes = step4(props_k, compMask, imgRawUp, params)
    filterH = [1, -1];

    % 1) Rotate mask so that major axis ∥ x-axis & process object
    angle = props_k.Orientation;
    obj = imdilate(compMask, strel('line', 3, angle));
    obj = imfill(obj, 'holes');
    obj = imrotate(obj, -angle);
    objRot = imdilate(obj, strel('line', ceil(0.01 * (props_k.MajorAxisLength)), 0));
    objRot = imfill(objRot, 'holes');
    heiRawRot = imrotate(imgRawUp, -angle);
    
    % Find left and right boundary.
    filtered = imfilter(double(objRot), filterH, 'conv');

    [idx_lft, idx_rgt] = deal(find(filtered == 1), find(filtered == -1));

    [idlft_y, idlft_x] = ind2sub(size(objRot), idx_lft);
    [idrgt_y, idrgt_x] = ind2sub(size(objRot), idx_rgt);

    [newlfty, Ilfty] = sort(idlft_y);
    [~,       Irgty] = sort(idrgt_y);

    [newlftx, newrgtx] = deal(idlft_x(Ilfty), idrgt_x(Irgty));

    if numel(newrgtx) ~= numel(newlftx)
        yes = false;
        return
    end

    horizontals = newrgtx - newlftx; % Horizontal chords 

    % Average raw height values along the longest chords.
    lngH = find(horizontals == max(horizontals)); % horizontals all having max length
    [startIdx, endIdx] = deal(max(0, lngH(1) - 2), min(numel(horizontals), lngH(end) + 2));
    idh = startIdx:endIdx;                        % Chord indices
    strip = heiRawRot(newlfty(idh), max(1, min(newlftx(idh)) - 1): min(size(heiRawRot,2), max(newrgtx(idh)) + 1)); 
                                                  % Height values on the longest strip                                       
    strip = mean(strip, 1);                       % Averaged height values 

    % Convolve averaged strip with kernel    
    kernl = ones(1, floor(params.windowsize));
    kernl(1:floor(params.windowsize/2)) = -1;                % Kernel
    edgeSignal = conv(strip, kernl, 'same');                 % Convolution
    edgeSignal = rescale(abs(edgeSignal), 0, range(strip));  % Rescaling

    peaks = edgeSignal > range(strip) * params.margin;       % Limiting the selection to points far from 2 ends
    peaks = peaks + circshift(peaks, 1);                     % peaks is the threshold for 'edges'
    cross = (peaks == 1).* sign(gradient(edgeSignal));       % points that crossed the threshold

    [cross(1), cross(end)] = deal(0, 0);
    [lft, rgt] = deal(min(find(cross == -1)), max(find(cross == 1)));

    if sum(abs(cross(lft + 1:rgt - 1))) == 0
        yes = false;
        return
    else
        yes = true;
        return
    end
end

%% Helper: saveSegmentedBoundaries
% Overlay boundaries by group in colored PNGs under outSeg directory
function saveSegmentedBoundaries(imgRawUp, imgSeg, imgGrp, names, outGrp)
    % Define group to color mapping
    cmap = struct( ...
        'coarse',   [1, 0, 0], ...   % red
        'border',   [1, 0, 1], ...   % magenta
        'isolated', [0, 1, 0], ...   % green
        'vertical', [1, 1, 0], ...   % yellow
        'lateral',  [1, 0.5, 0]);    % orange

    segMask = imgSeg;
    CC = bwconncomp(segMask, 8);
    L = labelmatrix(CC);

    if size(imgRawUp,3) == 1
        bg = mat2gray(imgRawUp);       % normalize to [0,1]
        rgb = cat(3, bg, bg, bg);      % grayscale to RGB
    else
        rgb = im2double(imgRawUp);     % assume already RGB
    end

    figure('Visible','off'); imshow(rgb); hold on;
    groupnames = {'border', 'coarse', 'lateral', 'vertical', 'isolated'};

    % Overlay boundaries for each specified group
    for gi = 1:numel(groupnames)
        grp = groupnames{gi};
        if ~isfield(imgGrp, grp), continue; end
        idxs = imgGrp.(grp);
        if isempty(idxs), continue; end

        B = bwboundaries(ismember(L, idxs));
        for k = 1:numel(B)
            boundary = B{k};
            plot(boundary(:,2), boundary(:,1), 'Color', cmap.(grp), 'LineWidth', 1.5);
        end
    end
    hold off;

    % Save overlay image
    outFile = fullfile(outGrp, [names '_grouped.png']);
    frame = getframe(gca);
    imwrite(frame.cdata, outFile);
    close;
end