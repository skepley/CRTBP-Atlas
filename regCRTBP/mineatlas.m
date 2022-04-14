function connectionData = mineatlas(atlas1, atlas2, atlasIndex, connectionData)
% mine atlases with specified index for connections and append them to connectionData

% Author: Shane Kepley
% Date: 3-Mar-2018


% INPUT:
% atlasIndex = [stableIdx,unstableIdx]: File identifiers for the loaded atlases
% stableAtlas, unstableAtlas: (un)stable atlases to mine
% connectionData: Existing connectionData to append new connections to


thisConnectionData = cell(0);
stableGen = 1;
unstableGen = 1; 

% compare initial boundaries pairwise
thisConnectionData = minethisgen(atlas1{stableGen},atlas2{unstableGen},atlasIndex,[stableGen,unstableGen],thisConnectionData);

% expand boundaries one generation at a time, comparing at each step
while stableGen < length(atlas1) || unstableGen < length(atlas2)
    
    if stableGen < length(atlas1) % advance stable domain and mine connections
        stableGen = stableGen + 1;
        thisConnectionData = minethisgen(atlas1{stableGen},atlas2{unstableGen},atlasIndex,[stableGen,unstableGen],thisConnectionData);
    end
    
    if unstableGen < length(atlas2)  % advance unstable domain and mine connections
        unstableGen = unstableGen + 1;
        thisConnectionData = minethisgen(atlas1{stableGen},atlas2{unstableGen},atlasIndex,[stableGen,unstableGen],thisConnectionData);
    end
end

connectionData = [connectionData,thisConnectionData]; % append new connections
end

function connectionData = minethisgen(stableGen,unstableGen,atlasIndex,genIndex,connectionData)
% Rules out intersections for pairwise charts in fundamental domain for each atlas.

% coarse connection check using ell^1 box in R^4.
for j = 1:length(stableGen) % loop through fundamental domain for stable manifold
     
    stableChart = stableGen(j);
    stableBox = stableChart.ellonebox();
    for k = 1:length(unstableGen) % loop through fundamental domain for unstable manifold
        unstableChart = unstableGen(k);
        unstableBox = unstableChart.ellonebox();
        thisDistance = abs(stableBox - unstableBox);
        
        % connection ruled out by ell^1 box or (un)stable chart already contains a connection
        ruledOut = any(inf(thisDistance) > 0) || chkblacklist([atlasIndex(1),genIndex(1),j],[atlasIndex(2),genIndex(2),k],connectionData);
        if ~ruledOut
            localSolution = checkconnection(stableChart,unstableChart,10); % does Newton to check for intersection
            if ~isempty(localSolution)
                connectionData{end+1}.StableIndex = [atlasIndex(1),genIndex(1),j];
                connectionData{end}.UnstableIndex = [atlasIndex(2),genIndex(2),k];
                connectionData{end}.LocalSolution = localSolution;
                connectionData = getglobalconnection(stableChart,unstableChart,localSolution,connectionData);
            end
        end
    end
end
end

function ruleOut = chkblacklist(stableIdx,unstableIdx,connectionData)
% rule out connection indices which have already produced connections
try
    if isempty(connectionData)
        ruleOut = false;
    else
        ruleOutStable = any(arrayfun(@(j)isequal(stableIdx,connectionData{j}.StableIndex),1:length(connectionData)));
        ruleOutUnstable = any(arrayfun(@(j)isequal(unstableIdx,connectionData{j}.UnstableIndex),1:length(connectionData)));
        ruleOut = ruleOutStable || ruleOutUnstable;
    end
catch ME
    disp('stophere')
end
end

function connectionData = getglobalconnection(stableChart,unstableChart,localSolution,connectionData)
% get global connection data and append it to connection structure

% local connection coordinates
stableLocalSpace = localSolution(1);
unstableLocalSpace = localSolution(3);
stableLocalTime = localSolution(2);

% global connection coordinates
stableGlobalSpace = .5*(sum(stableChart.MTCrange) + stableLocalSpace*diff(stableChart.MTCrange));
unstableGlobalSpace = .5*(sum(unstableChart.MTCrange) + unstableLocalSpace*diff(unstableChart.MTCrange));
stableGlobalTime = stableChart.TimeSpan(1) + stableLocalTime*stableChart.Tau;
unstableGlobalTime = unstableChart.TimeSpan(1);

% get full orbit from local unstable to local stable.
connectTime = unstableGlobalTime - stableGlobalTime;
fullTimeNode = linspace(0,connectTime,1000); % time vector for orbit between unstable and stable manifolds
unstableTimeNode = fullTimeNode(fullTimeNode <= unstableGlobalTime); % portion of flight in unstable atlas
stableTimeNode = fullTimeNode(length(unstableTimeNode)+1:end) - connectTime; % portion of flight in stable atlas
[stableOrbit,~] = taylororbit(stableGlobalSpace,stableChart,stableTimeNode);
[unstableOrbit,~] = taylororbit(unstableGlobalSpace,unstableChart,unstableTimeNode);

% append connection Data
connectionData{end}.ConnectTime = connectTime;
connectionData{end}.Orbit = [unstableOrbit;stableOrbit];
connectionData{end}.OrbitTime = fullTimeNode';

% fullConnectionTime = [unstableTime',stableTime'];
% intersectionPt = unstableChart.eval([unstableLocalSpace,unstableLocalTime])';
% stableTarget = stableAtlas{1}.eval([stableGlobalSpace,0]);
% unstableTarget = unstableAtlas{1}.eval([unstableGlobalSpace,0]);
end










