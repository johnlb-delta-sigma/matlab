%% Setup testing environment.


%% Make sure the matlab folder is on the path
thispath = mfilename('fullpath');
thispath_= split(thispath,filesep);
toolpath = ['/' fullfile(thispath_{1:end-1})];

pathCell = regexp(path, pathsep, 'split');
if ispc  % Windows is not case-sensitive
  onPath = any(strcmpi(toolpath, pathCell));
else
  onPath = any(strcmp(toolpath, pathCell));
end

if ~onPath
	disp('hi')
	addpath(toolpath);
end

clear all;