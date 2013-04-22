% Script for running modelReduction. 
% Main function:
% [ProdEstXY UXY VXY WXY BoundaryMatchingErrorXY RANSErrorXY]= modelReductionOpenFoamXY(locx, locy, flagSB, noBasis, Case, callBasis);
% Input variable description.
% i) locx, locy - % location of Center of back tile
% 
% ii) different basis for front and back tile.
%     flagSB = 0, same basis for all tiles.
%     flagSB = 1, separate basis for front and back tiles.
% 
% iii) noBasis = size of basis to be used for model reduction.
% 
% iv) planar movement of back tile. 
%     Case = 0, X/Transverse movement of back tile.
%     Case = 1, Y/Longitudinal movement of back tile.
%     Case = 2, Planar movement of back tile.
% 
% If one wants to compute basis(boundary basis), then use flag callBasis =1
% 
% Output from function (model reduced solution).
% ProdEstXY - Production Estimate
% UXY - x-component of velocity
% VXY - y-component of velocty
% WXY - z-component of velocity
% BoundaryMatchingErrorXY - Error in matching boundary at the interface of
% front and back tile.
% RANSErrorXY - Error in satisfying the RANS equations. [This needs to be corrected]
% 


clear all
clc

for i = 1:7
    a(i) = 0+(i-1)*3;
end

flagSB = 1; % Separate Basis flag
% location of Center of back tile
locx = 0;   % x-coordinate
locy = 0;   % y-coordinate

countBasis = 1;
Case = 2; % all snapshots
noBasis = 7;
callBasis = 1;
countLoc = 1;   

[ProdEstXY UXY VXY WXY BoundaryMatchingErrorXY RANSErrorXY] = modelReductionOpenFoamXY(locx, locy, flagSB, noBasis, Case, callBasis);
 