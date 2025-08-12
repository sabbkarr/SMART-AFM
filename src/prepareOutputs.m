% Generate & save summary figures and tables.
%
% Usage:
%   prepareOutputs(imgRawUp, imgPre, names, imgSeg, imgGrp, particleMeasurements, params, imageMetaData)
%
% Inputs:
%   particleMeasurements - particle measurement table (length, height, etc.)
%   params               - parameters
%   imageMetaData        - image metadata table (noise, contrast, etc.)

function prepareOutputs(particleMeasurements, params, imageMetaData)
    % Main output folder in src/output/
    outBase = fullfile(pwd, 'outputs');
    if ~exist(outBase, 'dir'), mkdir(outBase); end

    % 1. Measurements by group (including plots)
    outMeas = fullfile(outBase, 'measurements');
    if ~exist(outMeas, 'dir'), mkdir(outMeas); end

    groupList = params.groupList(:)';
    groupTbl = cell(numel(groupList),1);

    for g = 1:numel(groupList)
        gName = groupList{g};
        groupFolder = fullfile(outMeas, gName);
        if ~exist(groupFolder, 'dir'), mkdir(groupFolder); end

        groupTbl = particleMeasurements(strcmp(particleMeasurements.Group, gName), :); % Filter rows by group
        
        % Round numeric columns to 2 significant digits
        varNames = groupTbl.Properties.VariableNames;
        for v = 1:numel(varNames)
            col = groupTbl.(varNames{v});
            if isnumeric(col) groupTbl.(varNames{v}) = round(col, 2); end
        end

        writetable(groupTbl, fullfile(groupFolder, [gName '_measurements.csv'])); % Save CSV
        savePlots(groupTbl, groupFolder, gName); % Plot and save distributions
    end

    % Aggregated measurements (across all groups)
    aggFolder = fullfile(outMeas, 'all');
    if ~exist(aggFolder, 'dir'), mkdir(aggFolder); end
    writetable(particleMeasurements, fullfile(aggFolder, 'all_measurements.csv'));
    savePlots(particleMeasurements, aggFolder, 'all');

    % Average statistics
    getaverages(particleMeasurements, params, outMeas);

    % 2. Image metadata
    outMeta = fullfile(outBase, 'image_metadata');
    if ~exist(outMeta, 'dir'), mkdir(outMeta); end
    writetable(imageMetaData, fullfile(outMeta, 'image_metadata.csv'));
end

%% Helper: average statistics: per‐group mean & std, plus overall
function getaverages(particleMeasurements, params, outMeas)
    % Identify the numeric variables to summarize
    allVars    = particleMeasurements.Properties.VariableNames;
    skipVars   = {'ID', 'Group', 'ImageName'};  % skip non‐numeric cols
    numericVars = setdiff(allVars, skipVars, 'stable');
    
    % Preallocate containers
    summaryRows = {};
    summaryData = [];
    summaryCount = [];  

    % Loop over each requested group + 'all'
    for nameCell = [params.groupList(:); {'all'}]'
        grp = nameCell{1};
        if strcmp(grp, 'all')
            T = particleMeasurements;
        else
            T = particleMeasurements(strcmp(particleMeasurements.Group, grp), :);
        end
        
        % Compute mean and std (returns table of 1×numel(numericVars))
        M = varfun(@mean, T, 'InputVariables', numericVars);
        S = varfun(@std,  T, 'InputVariables', numericVars);
        N = height(T);

        % Append to arrays
        summaryData  = [summaryData;  M{1,:};             S{1,:}];
        summaryCount = [summaryCount; N;                  N];  % one for avg‐row, one for std‐row
        summaryRows  = [summaryRows;  {[grp '_average']}; {[grp '_std']}]; 
    end
    
    summaryTbl = array2table(summaryData, 'RowNames', summaryRows, 'VariableNames', numericVars);
    summaryTbl.Count = summaryCount;
    writetable(summaryTbl, fullfile(outMeas, 'summary_measurements.csv'), 'WriteRowNames', true);
end

%% Helper: savePlots
% Plot Length, Height, Aspect Ratio distributions
function savePlots(tbl, outDir, label)
    if ~isempty(tbl)
        f = figure('Visible','off');
        histogram(tbl.Length_nm, 30);
        xlabel('Length (nm)'); ylabel('Count');
        title(['Length Distribution (' label ')']);
        saveas(f, fullfile(outDir, [label '_length_hist.png'])); close(f);

        f = figure('Visible','off');
        histogram(tbl.Height_nm, 30);
        xlabel('Height (nm)'); ylabel('Count');
        title(['Height Distribution (' label ')']);
        saveas(f, fullfile(outDir, [label '_height_hist.png'])); close(f);

        f = figure('Visible','off');
        histogram(tbl.AspectRatio, 30);
        xlabel('Aspect Ratio'); ylabel('Count');
        title(['Aspect Ratio Distribution (' label ')']);
        saveas(f, fullfile(outDir, [label '_aspect_ratio_hist.png'])); close(f);

        % 2D Heat Map: Length vs Height
        f = figure('Visible','off');
        histogram2(tbl.Length_nm, tbl.Height_nm, 'NumBins', [30, 30], 'DisplayStyle', 'tile', 'ShowEmptyBins', 'off');
        xlabel('Length (nm)'); ylabel('Height (nm)');
        title(['2D Histogram: Length vs Height (' label ')']);
        colorbar;
        saveas(f, fullfile(outDir, [label '_length_height_heatmap.png'])); close(f);
    end
end

