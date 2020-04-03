%% save_custom_ntf: 
function save_custom_ntf(NTF,filename,from_fdatool)

% First, export from fdatool into 'Num' and 'Den' variables

% filename = 'simple';
% filename = 'hyb4-v0';

% NTF = tf(Num,Den,fs);
if (from_fdatool)
	NTF = from_fdatool(Num,Den,fs);
end

libdir = '~/Documents/MATLAB/tools/dsm/';
outfile = [libdir 'custom_ntfs/' filename '.mat'];
if (exist(outfile,'file')==2)
	error('Don''t overwrite your old file!')
end

save(outfile, 'NTF')
