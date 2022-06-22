function sols = mine_intersections(stableAtlas, unstableAtlas)
%MINE_INTERSECTIONS - Mine a pair of RegCRTBP Atlases for intersections.
%
%   MINE_INTERSECTIONS() - A more detailed description of the function
%
%   Syntax:
%       output = MINE_INTERSECTIONS(input1, input2)
%       [output1, output2] = MINE_INTERSECTIONS(input1, input2, input3)
%
%   Inputs:
%       stableAtlas - A RegCRTBP Atlas which has been integrated backward in time. 
%       unstableAtlas - A RegCRTBP Atlas which has been integrated forward in time. 
%
%   Outputs:
%       sols - A collection of connecting orbits
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: s.kepley@vu.nl
%   Date: 26-May-2022;

% mine for connections by leapfrogging fundamental domains in each atlas
sols = {}; % initialize cell array for connections
genB = -1; % stable atlas generation index
genF = 0;  % unstable atlas generation index
genFCurrent = unstableAtlas.generation(genF); % initialize unstable fundamental domain as a vector of charts
nGenF = length(genFCurrent); % number of charts parameterizing the initial unstable fundamental domain

while genB < stableAtlas.LastGeneration || genF < unstableAtlas.LastGeneration
    
    % leapfrog to the next stable fundamental domain
    if genB < stableAtlas.LastGeneration
        genB = genB + 1;
        genBCurrent = stableAtlas.generation(genB); % parameterize stable fundamental domain as a vector of charts
        fprintf('Mining generations: %d-%d \n',[genB,genF])
        nGenB = length(genBCurrent); % number of charts parameterizing the current stable fundamental domain
        
        % double loop over charts in both current fundamental domains
        for iGenB = 1:nGenB
            for jGenF = 1:nGenF
                chart1 = genBCurrent(iGenB);  % stable chart
                chart2 = genFCurrent(jGenF);  % unstable chart
                [isTrue, ijSols] = check4intersection(chart1,chart2,5);
                if ~isempty(ijSols)
                    sols{end+1} = {chart1; chart2; ijSols; isTrue};
                    
                end
            end
        end
    end
    
    % leapfrog to the next unstable fundamental domain
    if genF < unstableAtlas.LastGeneration
        genF = genF + 1;
        genFCurrent = unstableAtlas.generation(genF); % parameterize unstable fundamental domain as a vector of charts
        fprintf('Mining generations: %d-%d \n',[genB,genF])
        nGenF = length(genFCurrent); % number of charts parameterizing the current unstable fundamental domain
        
        % double loop over charts in both current fundamental domains
        for iGenB = 1:nGenB
            for jGenF = 1:nGenF
                chart1 = genBCurrent(iGenB);  % stable chart
                chart2 = genFCurrent(jGenF);  % unstable chart
                [isTrue, ijSols] = check4intersection(chart1,chart2,5);
                if ~isempty(ijSols)
                    sols{end+1} = {chart1; chart2; ijSols; isTrue};
                end
            end
        end
    end
end
end % end mine_intersections

