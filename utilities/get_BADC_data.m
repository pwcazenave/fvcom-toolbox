function tfiles = get_BADC_data(site, filepath, files, credentials)
% Child function to do the actual downloading from the BADC site via FTP.
% If the remote file doesn't exist, the function continues to the next file
% warning that the file couldn't be found. You may end up with gaps if the
% Met Office don't have those data (e.g. the 31st July 2010 doesn't have
% any model results, for some reason).
% 
% INPUTS:
%   site - FTP server name (e.g. ftp.ceda.ac.uk')
%   filepath - path to the files to download
%   files - cell array of a file or files to download
%   credentials - cell array of {'username', 'password'}.
% 
% OUTPUTS:
%   tfiles = cell array of files downloaded. NaN = failed to download.
%
% WARNING:
%   This function will indiscriminately overwrite files in the destination
%   directory (which is the system temporary directory appended with
%   'metum'). Really, the files in 

global ftbverbose

assert(iscell(files), 'Provide a cell array of files to download')

tdir = fullfile(tempdir, 'metum');

if exist(tdir, 'dir') ~= 7
    mkdir(tdir)
end

nf = length(files);

tfiles = cell(nf, 1);

for j = 1:nf
    % Open a remote connection to the FTP site. I found that the timeout on
    % the FTP connection to ftp.ceda.ac.uk is pretty low, so it's best to
    % explicitly re-open it each time we want a bunch of new files.
    remote = ftp(site, credentials(1), credentials(2));

    S = whos('remote');
    assert(strcmpi(S.class, 'ftp'), 'remote is not an FTP class. See HELP FTP.')
    clear S
    
    try
        cd(remote, filepath);
    catch
        warning('No such path %s (file %s)', filepath, files{j})
        continue
    end

    try
        if ftbverbose
            fprintf('Downloading %s to %s\n', fullfile(site, filepath, files{j}), tdir)
        end
        tfiles{j} = mget(remote, files{j}, tdir);
    catch
        tfiles{j} = nan;
        warning('Failed to fetch data from %s', fullfile(site, filepath, files{j}))
    end
    
    % Close the connection to the FTP server.
    close(remote)
end

