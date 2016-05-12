function viewresults(dirName, str1, lo, up)
% myinit();
addpath(fullfile('results', dirName));
addpath('utils');
if nargin < 4
    up = 1; 
end 
if nargin < 3 
    lo = 0;
end 
if nargin < 2       
    str1 = '.mat';
end 

warning off;

% dirName = 'results_Dec28';
% dirName = 'results_Dec27_reduce';
ext = '.mat';

fileList = getAllFiles_ext(fullfile('results', dirName), ext);
numel(fileList)
cnt = 0;
for i =  numel(fileList):-1:1

    if ~numel(findstr(fileList{i}, str1))
        continue;
    end
    load(fileList{i}, 'acc', 'rt');
    [p, fn, ext] = fileparts(fileList{i});
    if exist('acc', 'var')
        [a, id] = max(acc);
        if a > lo && a < up
            fprintf('%80s: %5f, id: %3d', fn, a, id);
            if exist('rt', 'var')
                fprintf('|  rt = %5.1f (s)', rt);
                
            end 
            fprintf('\n');
            cnt = cnt + 1;
        end
    end
end 
fprintf('Total: %3d\n', cnt);
end 
