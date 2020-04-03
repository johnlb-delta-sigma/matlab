%% custom_ntf: loads custom-designed NTF from .mat file
function [NTF] = custom_ntf(filename)
	libdir = '~/Documents/MATLAB/tools/dsm/';
	load([libdir 'custom_ntfs/' filename '.mat']);
end