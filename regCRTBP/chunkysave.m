function chunkysave(atlas, savePath, saveName, varargin)
%CHUNKYSAVE - Custom save method for RegCRTBPAtlas class objects. Atlases over 1GB are broken into pieces and saved across
% several different .mat files. Use chunkyload to reassemble the atlas from these parts.
%
%   Inputs:
%       atlas - RegCRTBPAtlas object
%       savePath - Path to folder to save data
%       saveName - Name of file or subfolder where atlas data will be saved
%
%   Outputs:
%       output - Produces either a single .mat file with specified name, or a subfolder with specified name containing multiple chunks of the 
% atlas data

%   Author: Shane Kepley
%   email: s.kepley@vu.nl
%   Date: 06-Feb-2021; Last revision: 06-Feb-2021

%% ================================================== SECTION 1 ==================================================

% function input
% atlas = atlasP1B;
% savePath = [computerPath, 'Desktop/CRTBP_collision_manifolds/'];
% saveName = 'chunkysavetest';
% saveName = 'L4_stable_t7_CL4_4020_equalmass_default';
% saveName = ['P1_stable_', runID]



testChart = atlas.Chart(1);  % get a representative Chart to estimate the Atlas size.
testScalar = testChart.Coordinate;
scalarAttributes = whos('testScalar');
scalarSize = scalarAttributes.bytes;  % estimated size of a single Chart (in bytes)
atlasSize = scalarSize*atlas.Size*1e-9;  % estimated size of the atlas (in gigabytes)
nChunk = ceil(atlasSize);  % number of 1 GB chunks to chop the atlas in to

if nChunk == 1  % atlas is less than 1GB. Just use the regular MATLAB save command
    save([savePath, saveName, '.mat'], 'atlas')
    
else
    
    % Initialize a new empty atlas to store metadata and properties
    emptyAtlas = RegCRTBPAtlas();
    atlasProperties = fields(struct(atlas)); % extract all atlas properties including hidden properties
    for iProperty = 1:length(atlasProperties) % loop over all properties of the atlas except the Charts/Crashstack and assign them to copyObj
        thisProperty = atlasProperties{iProperty};
        if ~isequal(thisProperty, 'Chart') && ~isequal(thisProperty, 'CrashStack')
            emptyAtlas.(thisProperty) = atlas.(thisProperty); % assign this property to copyObj
        end
    end
    
    % prepare orphaned crash stack
    nCrash = length(atlas.CrashStack);
    crashParentIdx = arrayfun(@(j)parentindex(atlas, atlas.CrashStack(j)), 1:length(atlas.CrashStack));  % record parent handle index for Crash stack
    emptyAtlas.CrashStack = RegCRTBPChart();  % % initialize crash stack as array of empty Charts
    emptyAtlas.CrashStack(nCrash) = RegCRTBPChart(); % initialize correct size of crash stack array
    for iCrash = 1:nCrash
        emptyAtlas.CrashStack(iCrash) = killparent(atlas.CrashStack(iCrash));  % append an orphaned version of the Crash Chart
    end
    
    chartParentIdx = arrayfun(@(j)parentindex(atlas, atlas.Chart(j)), 1:atlas.Size);  % compute parent handle indices to save with the empty atlas
    
    % save empty atlas
    saveFolder = [savePath, saveName, '/'];  % save Atlas in chunks to a folder
    if ~isfolder(saveFolder)
        mkdir(saveFolder);
    end
    save([saveFolder, 'atlas'], 'emptyAtlas', 'crashParentIdx', 'chartParentIdx')  % save empty atlas with attributes
    
    % save the orphaned Charts in chunks
    orphanCharts(atlas.Size) = RegCRTBPChart();  % initialize array of orphaned charts
    for iChart = 1:atlas.Size
        orphanCharts(iChart) = killparent(atlas.Chart(iChart));
    end
    
    chunkIdx = floor(linspace(1, atlas.Size, 1 + nChunk)); % indices for charts in each chunk
    for iChunk = 1:nChunk
        chunk = orphanCharts(chunkIdx(iChunk):chunkIdx(1+iChunk));
        save([saveFolder, sprintf('chart%d',iChunk)], 'chunk')
    end
end



end % end chunkysave

function orphanChart = killparent(chart)
% Return a chart with its parent handle removed
orphanChart = RegCRTBPChart();
chartProperties = fields(struct(chart)); % extract all chart properties including hidden properties
for iProperty = 1:length(chartProperties) % loop over all properties of the atlas except the Charts/Crashstack and assign them to copyObj
    thisProperty = chartProperties{iProperty};
    if ~isequal(thisProperty, 'ParentHandle')
        orphanChart.(thisProperty) = chart.(thisProperty); % assign this property to the blank chart
    end
end
end

function idx = parentindex(atlas, chart)
% Return the index of the chart parent handle if it exists or 0 if it is a root node
idx = 0;  % default return value
if ~isempty(chart.ParentHandle)
    idx = find(atlas.Chart==chart.ParentHandle);
end
end

% Revision History:
%{

%}
