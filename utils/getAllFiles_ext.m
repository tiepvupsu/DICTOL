function fileList = getAllFiles_ext(dirName, ext)
  % get all files with extension EXT in the folder dirName
  % output: a cell of filenames. fileList{1} is the directory + name of the 1st file in dirName
  dirData = dir(strcat(dirName, '/*', ext));       % Get the data for the current directory
  dirIndex = [dirData.isdir]    ;   % Find the index for directories
  fileList = {dirData(~dirIndex).name}';  %' Get a list of the files in current folder
  % --------------- In subfolders -------------------------
  dirData = dir(dirName);                 % Get the data for the current directory
  dirIndex = [dirData.isdir];
  if ~isempty(fileList)
    fileList = cellfun(@(x) fullfile(dirName,x),...  % Prepend path to files
                       fileList,'UniformOutput',false);
  end
  subDirs = {dirData(dirIndex).name};           % Get a list of the subdirectories
  validIndex = ~ismember(subDirs,{'.','..'});  % Find index of subdirectories
                                                %   that are not '.' or '..'
  for iDir = find(validIndex)                   % Loop over valid subdirectories
    nextDir = fullfile(dirName,subDirs{iDir});    % Get the subdirectory path

    fileList = [fileList; getAllFiles_ext(nextDir, ext)];  % Recursively call getAllFiles
  end


% sfjjsafj%
end