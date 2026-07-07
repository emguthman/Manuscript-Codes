function updated_filedata = removeBogus(filedata, string)
%%%INPUTS:
%   filedata - file object containing name
%   string - string to target files to remove

%%%OUTPUTS:
%   updated_filedata - updated file objects

%inits
nFiles = length(filedata);

%ID bogus copies of file names
bogus_files = zeros(1, nFiles);
parfor versye = 1:nFiles
    if contains(filedata(versye).name, string)
        bogus_files(versye) = 1;
    end
end

%output only real filenames
updated_filedata = filedata(~bogus_files);