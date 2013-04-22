
function [result] = readFoamResultFile (fname)

	%Read one FLACS result file 
% % 	
% % 	fid = fopen(fname, 'rt');
% % 
% % 	%GRID
% % 	line = fgetl(fid);
% % 
% % 	%dimension
% % 	line = fgetl(fid);
% % 	line = fgetl(fid);

% % 
	result.data=load (fname);%fname;
	result.dcount=length(result.data);
    
% This information should be in another data file and read from it to 
% generalize and automate the process. Could be the same one with value for
% molecular viscosity.
% For now just writing it since the dimension is known.
    result.nx = 120;
    result.ny = 72;
    result.nz = 72;

% % 	fclose (fid);
