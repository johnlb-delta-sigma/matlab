# matlab
Some MATLAB Design Tools for Delta Sigma Modulators originally from John Bell's PhD thesis. More details on how these functions work can be found there.


# Installing
To install this toolbox:

1. Clone this into a new directory: '~/Documents/MATLAB/johnlb-delta-sigma' (or '%USERPROFILE%/MATLAB/johnlb-delta-sigma' on windows). 
2. If it doesn't already exist, add the file "startup.m" to your '~/Documents/MATLAB' folder with the following contents:

```matlab
tempPath = userpath;
if tempPath(end)=='/'
    tempPath(end) = [];
end
addpath([tempPath '/johnlb-delta-sigma/matlab'])

clear all;
```

The "startup.m" script runs every time MATLAB starts, and adds this repo's folder to MATLAB's path so you can call any of the functions like any other MATLAB function.
