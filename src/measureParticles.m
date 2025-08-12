% Compute length, height & aspect ratio for desired particles.
%
% Inputs:
%   imgRawUp       – upscaled raw image (uint8 or double)
%   imgSeg         – logical mask of particles
%   imgGrp         – grouping of the particles
%   params         – struct with fields
%   names          – image names
%   back_height    – image background heights
%
% Output:
%   measTbl – table with columns with each row representing one particle:
%     particleID   – particle index
%     imageID      – image index
%     group        – particle group assignment
%     Length_nm    – length along major axis (nm)
%     Height_nm    – avg of top `topFrac` heights (nm)
%     AspectRatio  – Length / Height
%     Solidity     – Area / Convex area
%     Area 2D      – Object's projected area
%     Area 3D      – Object's 3D area
%   backgroundRoughness – image background roughness

function [tbl, backgroundRoughness] = measureParticles(imgRawUp, imgGrp, imgSeg, name, params, back_height)
    tbl = table();
    for g = 1:numel(params.groupList)
        % Get particle measurement table
        groupIdx  = imgGrp.(params.groupList{g}); % indices of that group
        if isempty(groupIdx), continue; end
        maskGroup = ismember(labelmatrix(bwconncomp(imgSeg, 8)), groupIdx);
        measTbl = measureImage(imgRawUp, maskGroup, params, back_height);

        % Add 'Group' and 'ImageName' columns to this batch
        nRows = height(measTbl);
        measTbl.Group = repmat({params.groupList{g}}, nRows, 1);
        measTbl.ImageName = repmat(string(name), nRows, 1);
        tbl = [tbl; measTbl];
    end
    A_proj = (size(imgRawUp, 1) - 1) * (size(imgRawUp, 2) - 1) * params.pixelLength^2;  % projected area in nm²
    backgroundRoughness = surfaceArea(imgRawUp .* (imgSeg == 0) * 1e9, params.pixelLength)/A_proj;
end

%% Measure particles in one image
function tbl = measureImage(imgRawUp, particleMask, params, back_height)
    CC = bwconncomp(particleMask);
    stats = regionprops(CC, 'PixelIdxList','Orientation', 'MajorAxisLength', 'Solidity', 'BoundingBox', 'Area');

    N = CC.NumObjects;
    [lengths, heights, AspRats, Solidity, Area3d, Area2d] = deal(zeros(N,1), zeros(N,1), zeros(N,1), zeros(N,1), zeros(N,1), zeros(N,1));

    for k = 1:N
        % 1) Preprocess object & its mask
        angle = stats(k).Orientation;
        obj = imdilate(ismember(labelmatrix(CC), k), strel('line', 5, angle));
        obj = imfill(obj, 'holes');
        obj = imdilate(obj, strel('disk', 1));
        obj_heights = imgRawUp(obj == 1);
        objRot = imrotate(obj, -angle);
        objRot = imdilate(objRot, strel('line', ceil(0.01 * (stats(k).MajorAxisLength)), 0));
        objRot = imfill(objRot, 'holes');
        nm_rotated = imrotate(imgRawUp, -angle);

        % 2) Apply convolution to find horizontal chords
        filtered = imfilter(double(objRot), [1, -1], 'conv');
        [idx_lft, idx_rgt] = deal(find(filtered == 1), find(filtered == -1));

        [idlft_y, idlft_x] = ind2sub(size(objRot), idx_lft);
        [idrgt_y, idrgt_x] = ind2sub(size(objRot), idx_rgt);

        [~, Ilfty] = sort(idlft_y);
        [~, Irgty] = sort(idrgt_y);

        [newlftx, newrgtx] = deal(idlft_x(Ilfty),  idrgt_x(Irgty));

        horizontals = newrgtx - newlftx;
        horizontals = horizontals(horizontals > 0.05 * range(horizontals));

        % 3) Measure
        heights(k)  = (mean(obj_heights(obj_heights > min(obj_heights) + 0.8 * range(obj_heights))) - back_height) * 1e9;
        lengths(k)  = params.pixelLength * max(horizontals);
        AspRats(k)  = lengths(k)/heights(k);
        Solidity(k) = stats(k).Solidity * 100;
        Area3d(k)   = surfaceArea(imcrop(imgRawUp .* obj * 1e9, stats(k).BoundingBox), params.pixelLength);
        Area2d(k)   = stats(k).Area * params.pixelLength * params.pixelLength;
    end

    tbl = table((1:N)', lengths, heights, AspRats, Solidity, Area3d, Area2d, ...
        'VariableNames', {'ID','Length_nm','Height_nm','AspectRatio', 'Solidity_%', 'Area 3D_nm2', 'Area 2D_nm2'});
end
%% Helper function-Particle 3D surface area
% Compute surface‐area of height map Z
% by summing two triangles per cell.
%
%   A = surfaceArea_tri_noloop(Z, dx, dy)
%
% Inputs:
%   Z   M×N matrix of heights--particle masked out
%   dx  grid spacing in the x–direction
%   dy  grid spacing in the y–direction
%
% Output:
%   A   total surface area

function A = surfaceArea(Z, dxy)
    [m, n] = size(Z);

    % build grid coordinates
    [I, J] = ndgrid(0:m-1, 0:n-1);
    [X, Y] = deal(I * dxy, J * dxy);

    % corner points of each cell
    P1 = cat(3, X(1:end-1,1:end-1), Y(1:end-1,1:end-1), Z(1:end-1,1:end-1));
    P2 = cat(3, X(2:end,  1:end-1), Y(2:end,  1:end-1), Z(2:end,  1:end-1));
    P3 = cat(3, X(1:end-1,2:end  ), Y(1:end-1,2:end  ), Z(1:end-1,2:end  ));
    P4 = cat(3, X(2:end,  2:end  ), Y(2:end,  2:end  ), Z(2:end,  2:end  ));

    % reshape to 3×N arrays for vectorized cross‐product
    P1v = reshape(P1, [], 3).';
    P2v = reshape(P2, [], 3).';
    P3v = reshape(P3, [], 3).';
    P4v = reshape(P4, [], 3).';

    C1  = cross(P2v - P1v, P3v - P1v, 1);  % first triangle: (P1,P2,P3)
    C2  = cross(P3v - P4v, P2v - P4v, 1);  % second triangle: (P4,P3,P2)
    
    % total area
    A = 0.5 * sum(sqrt(sum(C1.^2,1))) + 0.5 * sum(sqrt(sum(C2.^2,1)));
end
