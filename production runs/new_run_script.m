function [output1,output2] = new_run_script(templateFilename, saveAsFilename, mu, energy, time, source, target)
% NEW_RUN_SCRIPT - Generate a new run script from a given template and a choice of mass, energy, source, target, and integration parameters

% NEW_RUN_SCRIPT() - A more detailed description of the function

% Syntax:
% output = NEW_RUN_SCRIPT(input1, input2)
% [output1, output2] = NEW_RUN_SCRIPT(input1, input2, input3)
% 
% Inputs:
% input1 - Description
% input2 - Description
% input3 - Description
% 
% Outputs:
% output1 - Description
% output2 - Description
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: s.kepley@vu.nl
%   Date: 30-May-2022;

% INPUT
% saveAsFilename = 'test_script';
timeBackward = time(1);
timeForward = time(2);

% mu = 0.5;
% energy = 3;
% source = 'L4';
% target = 'P1';
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