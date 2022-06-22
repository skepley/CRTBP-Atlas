function [output1,output2] = new_run_script(templateFilename, saveAsFilename, mu, energy, time, source, target)
% NEW_RUN_SCRIPT - Generate a new run script from a given template and a choice of mass, energy, source, target, and integration parameters

% NEW_RUN_SCRIPT() - A more detailed description of the function

% Syntax:
% NEW_RUN_SCRIPT(templateFilename, saveAsFilename, mu, energy, time, source, target) generates a script in the same folder
% named "saveAsFilename". Running this script will compute the stable/unstable manifolds associated with the target/source
% respectively at the given mass and energy value and then try to mine these atlases for intersections. 
% 
% Inputs:
% templateFilname - a plain text file containing the path information for the computer which will run the script created
% saveAsFilename - A string specifying the name which the new script should be called.
% mu - The mass ratio of the small primary
% energy - The energy level to search in for connections. If either the source or target is L4 then this must be the same as the L4 energy. 
% time - [t-, t+] specifying the integration time for the stable and unstable maniflds respectively.
% source - A string identifying the unstable object for the connection. This should be either 'P1', 'P2', or 'filename_where_L4_local_data_is_saved'
% target - A string identifying the stable object for the connection. This should be either 'P1', 'P2', or 'filename_where_L4_local_data_is_saved'
%
%   Subfunctions: none
%   Classes required: IMP library, RegCRTBPAtlas, RegCRTBPChart
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: s.kepley@vu.nl
%   Date: 30-May-2022;

% unpack backward/forward integration times for stable/unstable manifolds
timeBackward = time(1);
timeForward = time(2);

if ~any(strcmp(source, {'P1', 'P2'}))
    unstableLocalDataFrom = ['''', source, ''''];
    source = 'L4';
else
    unstableLocalDataFrom = ['''', source, ''''];
end

if ~any(strcmp(target, {'P1', 'P2'}))
    stableLocalDataFrom = ['''', target, ''''];
    target = 'L4';
else
    stableLocalDataFrom = ['''', target, ''''];
end

localPath = mfilename('fullpath');  % fetch local path of this function
templatePath = localPath(1:end-14);  % set location to save new script by stripping this script name from local path
new_script = [templatePath, templateFilename]; % set location and name of the template script

file = matlab.desktop.editor.newDocument();  % initialize new file
file.saveAs(getFileName(saveAsFilename));  % save file

replacementDictionary = struct('name', saveAsFilename,...
    'NAME', upper(saveAsFilename),...
    'ENERGY', num2str(energy),...
    'MU', num2str(mu),...
    'SOURCE', source,...
    'TARGET', target,...
    'TIME1', num2str(timeBackward),...
    'TIME2', num2str(timeForward),...
    'date', char(datetime(date)),...
    'UNSTABLELOCALDATAFROM', unstableLocalDataFrom,...
    'STABLELOCALDATAFROM', stableLocalDataFrom...
    );
fields = fieldnames(replacementDictionary);


fid = fopen(new_script, 'r');
line = fgetl(fid);
[~, saveAsFilename, ~] = fileparts(file.Filename);

try
    while ~feof(fid)  % test for end of file
        line = replace(line, fields, replacementDictionary);
        file.appendText([line 10]);
        line = fgetl(fid);  % grab next line
    end
    
    % catch exception
catch exception
    throw(exception)  % handle exceptions so the file closes properly on exceptions
end
fclose(fid);


    function line = replace(line, fields, replacementDictionary)
        % replace template text with parameter values for this run
        for j = 1:numel(fields)
            fname = fields{j};
            line = strrep(line, ['$', fname], getfield(replacementDictionary, fname));
        end
    end



    function filename = getFileName(shortFilename)
        
        dot = strfind(shortFilename, '.');
        if ~isempty(dot)
            filename = shortFilename(1 : dot(1) - 1);
        else
            filename = shortFilename;
        end
        filename = [pwd, filesep, matlab.lang.makeValidName(filename), '.m'];
    end

end % end new_run_script