function atlas = chunkyload(savePath, saveName, varargin)
%CHUNKYLOAD - Custom load method for RegCRTBPAtlas class objects which are saved using chunkysave. Atlases saved across
% several different .mat files are reassembled and returned.
%
%   Inputs:
%       atlas - RegCRTBPAtlas object
%       savePath - Path to folder to save data
%       saveName - Name of file or subfolder where atlas data will be saved
%
%   Outputs:
%       output - Produces either a single .mat file with specified name, or a subfolder with specified name containing multiple chunks of the
% atlas data

if isfolder([savePath, saveName])
    subFolder = [savePath, saveName, '/'];
    load([subFolder, 'atlas'], 'emptyAtlas', 'crashParentIdx', 'chartParentIdx')
    
    % attach the charts to the emptyAtlas
    iChunk = 1;
    while isfile([subFolder, sprintf('chart%d.mat', iChunk)])
        load([subFolder, sprintf('chart%d', iChunk)], 'chunk')
        emptyAtlas.Chart = [emptyAtlas.Chart, chunk];
        iChunk = iChunk + 1;
    end
    
    % reattach parent handles to each crash chart
    nCrash = length(emptyAtlas.CrashStack);
    for iCrash = 1:nCrash
        parentID = crashParentIdx(iCrash);
        if parentID > 0
            emptyAtlas.CrashStack(iCrash).ParentHandle = emptyAtlas.Chart(parentID);
        end
    end
    
    % reattach parent handles to each crash chart
    for iChart = 1:emptyAtlas.Size
        parentID = chartParentIdx(iChart);
        if parentID > 0
            emptyAtlas.Chart(iChart).ParentHandle = emptyAtlas.Chart(parentID);
        end
    end
    atlas = emptyAtlas;
    
elseif isfile([savePath, saveName, '.mat'])
    load([savePath, saveName, '.mat'], 'atlas')
    
else
    error(['No file or folder named ', savePath, saveName])
end

end % end chunkyload

% Revision History:
%{

%}
