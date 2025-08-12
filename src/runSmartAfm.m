% Top-level pipeline: load, process, segment, group, measure & plot.
%
% Usage:
%   runSmartAfm; 

defineParams;                                     % read parameters
% imageDir = fullfile(pwd,'data');                  % directory for input images

outBase = fullfile(pwd, 'outputs');               % directory to save output 
if ~exist(outBase, 'dir'), mkdir(outBase); end    % create directory
% Initialization

[backgroundRoughness, particleMeasurements, grouping, group_stat] = deal([], table(), {}, {});
vars = {'Background Roughness','Noise (nm)','Contrast (nm)','Background Height (nm)', '# Border Particles',...
    '# Coarse Particles','# Lateral Particles','# Vertical Particles','# Isolated Particles', 'Area Border Particles %',...
    'Area Coarse Particles %','Area Lateral Particles %','Area Vertical Particles %','Area Isolated Particles %'};
% -------------------------------------- workflow starts
%% -------------------------------------- 1) Load images


fprintf('--- SMART‑AFM: Loading images from %s\n', imageDir);
tStartLoad = tic;
[imgRaw, names] = loadSmartImages(imageDir);
tLoad = toc(tStartLoad);
% imgRaw = imgRaw(6:10);
% names = names(6:10);
N = numel(imgRaw);
fprintf('    ✓ Loaded %d image(s) in %.2f s\n\n', N, tLoad);


%% -------------------------------------- 2) Preprocess & segment each image


fprintf('--- SMART‑AFM: Processing & segmenting images\n');
tStartLoad = tic;
[imgPre, imgSeg, imgRawUp, back_height, noisee, contrastt] = imageProcessing(imgRaw, names, params);
tLoad = toc(tStartLoad);
fprintf('    ✓ Processed & segmented %d image(s) in %.2f s\n\n', N, tLoad);


%% -------------------------------------- 3) Classify (group) & measure particles

[cumLen, cumHgt, imgIdx] = deal([], [], []);

% Setup live figure with two subplots
figure; fLive = figure('Name', 'Live Cumulative Averages', 'NumberTitle', 'off');

% Length plot
ax1 = subplot(2,1,1);
pLen = plot(ax1, NaN, NaN, '-o', 'LineWidth', 1.5); grid(ax1,'on');
xlabel('CNC Count'); ylabel('Mean Length (nm)');
title('Cumulative Mean Length');

% Height plot
ax2 = subplot(2,1,2);
pHgt = plot(ax2, NaN, NaN, '-or', 'LineWidth', 1.5); grid(ax2,'on');
xlabel('CNC Count'); ylabel('Mean Height (nm)');
title('Cumulative Mean Height');

fprintf('--- SMART‑AFM: Classify & measure particles...\n');
tStartLoad = tic;
for i = 1:N
    [grouping{i}, group_stat{i}] = particleGrouping(imgSeg{i}, imgRawUp{i}, params, names{i});
    [tbl, backgroundRoughness(end+1)] = ...
        measureParticles(imgRawUp{i}, grouping{i}, imgSeg{i}, names{i}, params, back_height(i));
    particleMeasurements = [particleMeasurements; tbl];

    cumLen(end+1) = mean(particleMeasurements.Length_nm);
    cumHgt(end+1) = mean(particleMeasurements.Height_nm);
    imgIdx(end+1) = numel(particleMeasurements.Length_nm);

    set(pLen, 'XData', imgIdx, 'YData', cumLen);
    set(pHgt, 'XData', imgIdx, 'YData', cumHgt);
    drawnow;

    % Progress bar
    pct = floor(100 * i / N);
    fprintf('\b\b\b\b%3d%%', pct);  % Overwrites last percentage
    if i == 1, fprintf('\n'); end   % new line after first percentage
end
tLoad = toc(tStartLoad);
fprintf('\n    ✓ Grouped and measured particles in %d image(s) in %.2f s\n\n', N, tLoad);

%% -------------------------------------- 4) Prepare Outputs

imageMetaData = table('Size', [N, numel(vars)], 'VariableTypes', ...
    repmat("double", 1, numel(vars)), 'VariableNames', vars, 'RowNames', names);

for i = 1:N
    imageMetaData(i, :) = { ...
        backgroundRoughness(i), ...
        noisee(i) * 1e9, ...
        contrastt(i) * 1e9, ...
        back_height(i) * 1e9, ...
        group_stat{i}.pc_border, ...
        group_stat{i}.pc_coarse, ...
        group_stat{i}.pc_lateral, ...
        group_stat{i}.pc_vertical, ...
        group_stat{i}.pc_isolated, ...
        group_stat{i}.paf_border, ...
        group_stat{i}.paf_coarse, ...
        group_stat{i}.paf_lateral, ...
        group_stat{i}.paf_vertical, ...
        group_stat{i}.paf_isolated };
end

prepareOutputs(particleMeasurements, params, imageMetaData);

% save the cumulative plots
fig = figure('Visible','off');
plot(imgIdx, cumLen, '-o','LineWidth',1.5);
xlabel('Particle Count'); ylabel('Cumulative Mean Length (nm)');
title('Cumulative Mean Length');
saveas(fig, fullfile(pwd,'outputs','measurements','all_cumulative_mean_length.png'));
close(fig);

fig = figure('Visible','off');
plot(imgIdx, cumHgt, '-o','LineWidth',1.5);
xlabel('Particle Count'); ylabel('Cumulative Mean Height (nm)');
title('Cumulative Mean Height');
saveas(fig, fullfile(pwd,'outputs','measurements','all_cumulative_mean_height.png'));
close(fig);
