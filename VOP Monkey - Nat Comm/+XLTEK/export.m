
function [data] = export(session,PATHFILE,varargin)
% Load and process the XLTEK
testType = {'DCS'};
contact = {'APB'};
setNum = []; % If empty, this will just load all of the available set data
% in order

assignopts(who, varargin);

data = [];
% Determine the number of files based on the conditions we identifed.
D = dir(fullfile(PATHFILE,session));
fileNames = {D.name};
posTest = contains(fileNames,testType);
posContact = contains(fileNames,contact);
if ~isempty(setNum)
    posSet = contains(fileNames,setNum);
else
    posSet = ones(1,length(D));
end

% Load the data
valLoc = find(posTest & posContact & posSet);
for n = 1:length(valLoc)
    temp = csvread(fullfile(D(valLoc(n)).folder,D(valLoc(n)).name));
    if n == 1
        data = temp(:,1);
    else
        data(:,n) = temp(:,1);
    end
end
end