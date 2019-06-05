function tensorlabinstalled = checkfortensorlab,
% Check if tensorlab is installed. If not, it opens a dialog window to select tensorlab path. 

tensorlabinstalled = 0; 

if exist('sdf_check') ~= 2, % this is a function from tensorlab
    warning('This code requires tensorlab. Please download tensorlab from "http://www.tensorlab.net/" (click <a href="http://www.tensorlab.net">here</a>).');

    % to manually add tensorlab path to the MATLAB path
    tlpath = uigetdir(pwd,'Select tensorlab path');
    if exist('tlpath'), addpath(tlpath); end
end

if exist('sdf_check') ~= 2,
    error('This code requires tensorlab.');
else
    tensorlabinstalled = 1;
end

end % function

