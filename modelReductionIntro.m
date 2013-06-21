%% InstantWind Project
% A step by step introduction to using InstantWind.
% 
% InstantWind was a project financed by Vindforsk, lead by Meventus.  
% The project group consisted of Meventus, CMR, Høgskolan på Gotland and
% OpenCFD
% For more information & contact info, see 
% InstantWind homepages <http://www.instantwind.com www.instantwind.com >
%
% Meventus AS <http://www.meventus.com www.meventus.com>
%
% CMR <http://www.cmr.com www.cmr.com>
%
% Høgskolan på Gotland <http://www.hgo.se www.hgo.se>
%
% OpenCFD <http://www.opencfd.com www.opencfd.com>
%    
%% Project Setup 
%
%   From <http://85.19.218.244/instantwind/matlab.html Matlab Code> download and unzip InstantWind.zip into your matlab folder    
%   The zip file contains the following folders and files  
%   InstantWind/                            - main folder 
%   InstantWind/VBShearFlow                 - should contain CFD results
%   InstantWind/VBShearFlow/bin             - store CFD results as .mat files  
%   InstantWind/ModelReductionDemo          - contains the matlab and c functions
%   InstantWind/ModelReductionDemo/figures  - results as figures
   
%% 
% It is suggested that the user goes through the theory section of the report <http://85.19.218.244/instantwind/documents.html InstantWind Report> to get an understanding of the tasks involved. Some of the keywords like 'snapshots', 'basis' are described in the report.
% For a summary of the algorithm please refer to chapter 2 in the above report. 
%
%   Download the CFD results at <http://85.19.218.244/instantwind/openfoams.html> (Download All files link) into the VBShearFlow folder. 
%   The .zip file include the subfolders 
%   
%   * Horizontal
%   * Vertical
%   * AdditionalCases
%   and each of the folders contains a /bin folder to store CFD results in
%   .mat files


%%    
% Before we can run the model, we need to compile the c functions in the ModelReductionDemo folder into mex files (for RANS operators) in matlab.  At the command prompt type the following: 
    mex A_mex_check2.c
    mex V_mex_check2.c
    mex P_mex_check3.c
    mex K_mex_check3.c
%%    
% In case your matlab is not associated with a compiler, you will be
% requested by matlab to do that.  We have tested that the mex files work
% on windows 32-bit versions of Matlab.  
% *NOTE:* The Windows 64-bit version of Matlab crashes when mexing.  This
% has been reported to Matlab, and will be solved later. 

%% Running InstantWind
% We now describe how to run the project. 
% It is suggested that the user goes through the theory section of the report <http://85.19.218.244/instantwind/documents.html InstantWind Report> to get an understanding of the tasks involved. Some of the keywords like 'snapshots', 'basis' are described in the report.
% For a summary of the algorithm please refer to chapter 2 in the above report. 

%%  
% 1) Start matlab and change current directory to InstantWind/ModelReductionDemo.
    cd \InstantWind/ModelReductionDemo
%% 
% 3) At the matlab command prompt type 'modelReductionDemo' and press
% the 'return' key.  This is the main function for the project
    modelReductionDemo

%% 
% 3) Clearly, the main script/function for the project is 
% _modelReductionDemo.m_ and the description of variables is described in
% it

%%
% 4) When the user is asked for an input to the keyboard, please type in
% (and not copy paste) 
%    makeSnapshotsCollectionOpenFoamXY(snapshotCollectionFileName, CFD_resultsPathH, CFD_ResultsPathV, CFD_JObsStruct, flagsB);
    
%%
% 5) The code goes through the following main functions: 
% * *makeSnapshotsCollectionOpenFoamXY* - Collects the CFD results
%   (snapshots)
%
%  
% * *readFrontTileOpenFoam* - Collects information of the frontTile, to be
%   used as a constraint (rhs of equations Ax=b)
%
% * *makeBasisCollectionOpenFoam* – Performs SVD on the snapshots (from the
% snapshotlist.)
%
% * *RANSEquationsOperators2* – Computers the operators (Advection, Viscous,
%   Pressure, Turbulent Kinetic Energy) for the RANS equations.
%
% * *modelReductionSolverOpenFoamXY* – computes the reduced solutions based on the
% information from the above functions.       