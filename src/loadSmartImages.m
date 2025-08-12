%% loadSmartImages.m
% Loads all AFM images from a folder and extract names.
%
% Usage:
%   [imgRaw, fileNames] = loadSmartImages(imageDir)
%
% Inputs:
%   imageDir    – folder to cd into before listing files
%
% Outputs:
%   imgRaw      – 1xN cell array of image matrices (double)
%   fileNames   – 1xN cell array of original file name strings excluding file extension (e.g. 'ABC')

function [imgRaw, fileNames] = loadSmartImages(imageDir)
    files = dir(imageDir);
    [imgRaw, fileNames] = deal({}, {});

    for k = 1:numel(files)
        [~, name, ext] = fileparts(files(k).name);

        switch lower(ext)
            case '.mat'
                S = load(files(k).name);
                vars = fieldnames(S);
                imgRaw{end+1} = double(S.(vars{1}));
                fileNames{end + 1} = name;
            case {'.tif','.tiff','.jpg','.png'}
                imgRaw{end+1} = double(imread(files(k).name));
                fileNames{end + 1} = name;
            case '.txt'
                imgRaw{end+1} = load(files(k).name);
                fileNames{end + 1} = name;
        end
end