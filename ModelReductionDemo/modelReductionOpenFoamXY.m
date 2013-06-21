
function [ProdEstim U V W BoundaryError RANSError] = modelReductionOpenFoamXY(locx, locy, flagSB, noBasis, Case, callBasis)


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
% v) If one wants to compute basis(boundary basis), then use flag 
%    callBasis =1
% 
% Output from function (model reduced solution).
% ProdEstXY - Production Estimate
% UXY - x-component of velocity
% VXY - y-component of velocty
% WXY - z-component of velocity
% BoundaryMatchingErrorXY - Error in matching boundary at the interface of
% front and back tile.
% RANSErrorXY - Error in satisfying the RANS equations. [This needs to be corrected

    if (nargin == 5)
        callBasis = 0;        
    end
% 	CFD_ResultsPath = '..\VB-3tiles_OpenFOAM_data\'; %'
	CFD_ResultsPathH = '..\VBShearFlow\Horizontal\'; %'

	CFD_JobsStruct.hcount=5; % No of data files to be read in
	CFD_JobsStruct.listh(1,:) = 'dat_files_1';
	CFD_JobsStruct.listh(2,:) = 'dat_files_2';
	CFD_JobsStruct.listh(3,:) = 'dat_files_3';
	CFD_JobsStruct.listh(4,:) = 'dat_files_4';
	CFD_JobsStruct.listh(5,:) = 'dat_files_5';
    
    CFD_ResultsPathV = '..\VBShearFlow\Vertical\'; %'        
	CFD_JobsStruct.vcount=7; % No of data files to be read in
	CFD_JobsStruct.listv(1,:) = 'dat_files_1';
	CFD_JobsStruct.listv(2,:) = 'dat_files_2';
	CFD_JobsStruct.listv(3,:) = 'dat_files_3';
	CFD_JobsStruct.listv(4,:) = 'dat_files_4';
	CFD_JobsStruct.listv(5,:) = 'dat_files_5';    
	CFD_JobsStruct.listv(6,:) = 'dat_files_6';
	CFD_JobsStruct.listv(7,:) = 'dat_files_7';
    
    if (flagSB ==1 )        
        snapshotCollectionFileName = 'snapshotsCollection_4_openFoam_xy_SB.mat';
        frontTileFileName = 'front-tile_4_openFoam_horizontal_xy_SB.mat';
        if (Case == 0)        
            basisCollectionFileName = sprintf('%s%d%s','basisCollection_4_openFoam_horizontal_y_SB',noBasis,'.mat');
            RANS_OperatorsFileName = sprintf('%s%d%s','RANS_4_openFoam_horizontal_y_SB',noBasis,'.mat');
        elseif(Case == 1)
            basisCollectionFileName = sprintf('%s%d%s','basisCollection_4_openFoam_horizontal_x_SB',noBasis,'.mat');
            RANS_OperatorsFileName = sprintf('%s%d%s','RANS_4_openFoam_horizontal_x_SB',noBasis,'.mat');
        else
            basisCollectionFileName = sprintf('%s%d%s','basisCollection_4_openFoam_horizontal_xy_SB',noBasis,'.mat');
            RANS_OperatorsFileName = sprintf('%s%d%s','RANS_4_openFoam_horizontal_xy_SB',noBasis,'.mat');
        end
    else
        snapshotCollectionFileName = 'snapshotsCollection_4_openFoam_horizontal_xy.mat';
        frontTileFileName = 'front-tile_4_openFoam_horizontal_xy.mat';
        if (Case == 0)        
            basisCollectionFileName = sprintf('%s%d%s','basisCollection_4_openFoam_horizontal_y',noBasis,'.mat');
            RANS_OperatorsFileName = sprintf('%s%d%s','RANS_4_openFoam_horizontal_y',noBasis,'.mat');
        elseif(Case == 1)
            basisCollectionFileName = sprintf('%s%d%s','basisCollection_4_openFoam_horizontal_x',noBasis,'.mat');
            RANS_OperatorsFileName = sprintf('%s%d%s','RANS_4_openFoam_horizontal_x',noBasis,'.mat');
        else
            basisCollectionFileName = sprintf('%s%d%s','basisCollection_4_openFoam_horizontal_xy',noBasis,'.mat');
            RANS_OperatorsFileName = sprintf('%s%d%s','RANS_4_openFoam_horizontal_xy',noBasis,'.mat');
        end
    end
    
    disp('Please type (and not copy-paste) the following:')
    disp('makeSnapshotsCollecionOpenFoamXY(snapshotCollectionFileName, CFD_ResultsPathH, CFD_ResultsPathV, CFD_JobsStruct, flagSB);');
    keyboard;
%  	makeSnapshotsCollecionOpenFoamXY(snapshotCollectionFileName, CFD_ResultsPathH, CFD_ResultsPathV, CFD_JobsStruct, flagSB);
%     eval(sprintf('%s','makeSnapshotsCollecionOpenFoamXY(snapshotCollectionFileName, CFD_ResultsPathH, CFD_ResultsPathV, CFD_JobsStruct, flagSB)'))

	frontTileJob = 'dat_files_1';
 	readFrontTileOpenFoam(frontTileFileName, frontTileJob, CFD_ResultsPathH) % Either x or y direction front-tile should work. 
    
%   VB
    if (flagSB == 1)
        if (Case == 0 )
            snapshotList = [11 12 13 14 15 16 17 18 19 20 21 22 23 24 6 7 8 9 10 11 12]; % vertical/y-direction OpenFoam results            
        elseif (Case == 1)
            snapshotList = [1 2 3 4 5 6 7 8 9 10 1 2 3 4 5]; % horizontal/x-direction OpenFoam results            
        else
            snapshotList = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 1 2 3 4 5 6 7 8 9 10 11 12];
        end
    else
        snapshotList = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15];    
    end

    if (callBasis == 1) 
% %  Make basis from snapshots collected
%         fprintf('%s \n',basisCollectionFileName)
        makeBasisCollectionOpenFoam (snapshotCollectionFileName, basisCollectionFileName, snapshotList, CFD_ResultsPathH, flagSB, noBasis);
        RANSequationsOperators2(basisCollectionFileName,  CFD_ResultsPathH, RANS_OperatorsFileName, flagSB);        
    end

% %  RANS Operators:
%     RANSequationsOperators2(basisCollectionFileName,  CFD_ResultsPathH, RANS_OperatorsFileName, flagSB);

    alpha = 0; % significance of momentum error/boundary matching term.
%     keyboard;
    
		%Solve only for specified location
    [ProdEstim U V W BoundaryError RANSError] = modelReductionSolverOpenFoamXY(basisCollectionFileName, frontTileFileName, CFD_ResultsPathH, RANS_OperatorsFileName, alpha, locx, locy, flagSB, noBasis) ;%loc = 0...31
       
	fprintf(1,'Finished!\n');
	
	
	
