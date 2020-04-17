%MINE_COLLISION_MANIFOLDS - load and mine the large and small collision manifolds for intersections
%
%   Other m-files required: none
%   MAT-files required: none
%
%   See also: OTHER_SCRIPT_NAME,  OTHER_FUNCTION_NAME

%   Author: Shane Kepley 
%   email: shane.kepley@rutgers.edu
%   Date: 21-May-2019; Last revision: 21-May-2019

% clear all
% clc
% load ('/Users/sk2011/Desktop/CRTBP_collision_manifolds/collisions_3_3dot7_equalmass') 


THIS SCRIPT IS REPLACED BY MINE_CONNECTIONS.M
%% ================================================== % MINE FOR INTERSECTIONS  ==================================================

sols = {};
genIdx1 = [atlas1.Chart.Generation];
genIdx2 = [atlas2.Chart.Generation];
gen1 = 0; % large collision atlas generation index
gen2 = 0;  % small collision atlas generation index
gen2Current = atlas2.generation(gen2);
nGen2 = length(gen2Current); 

while gen1 < atlas1.LastGeneration || gen2 < atlas2.LastGeneration
    gen1 = min(gen1 + 1, atlas1.LastGeneration);
    gen1Current = atlas1.generation(gen1);
    disp([gen1,gen2])
    nGen1 = length(gen1Current);
    
    % loop over current atlas1 generation
    for iGen1 = 1:nGen1
        for jGen2 = 1:nGen2
            chart1 = gen1Current(iGen1);
            chart2 = gen2Current(jGen2);
            if isequal(chart1.RegType, chart2.RegType)
                ijSols = check4intersection(chart1,chart2,5);
                if ~isempty(ijSols)
                    sols{end+1} = {gen1Current(iGen1); gen2Current(jGen2); ijSols};
                end
            end
        end
    end
    
    % leapfrog to a fundamental domain for the next generation of atlas 2
    gen2 = min(gen2 + 1, atlas2.LastGeneration);
    gen2Current = atlas2.generation(gen2);
    disp([gen1,gen2])
    nGen2 = length(gen2Current);
    
    % loop over backward gen
    for iGen1 = 1:nGen1
        for jGen2 = 1:nGen2
            chart1 = gen1Current(iGen1);
            chart2 = gen2Current(jGen2);
            if isequal(chart1.RegType, chart2.RegType)
                ijSols = check4intersection(chart1,chart2,5);
                if ~isempty(ijSols)
                    sols{end+1} = {gen1Current(iGen1); gen2Current(jGen2); ijSols};
                end
            end
        end
    end
end
    % save('self_collisions')
    
    
    
    
    
