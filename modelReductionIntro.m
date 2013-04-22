%% InstantWind Project
% A step by step introduction to using InstantWind.
% 
% InstantWind was a project financed by Vindforsk, lead by Meventus.  
% The project group consisted of Meventus, CMR, Høgskolan på Gotland and
% OpenCFD
% For more information & contact info, see 
% InstantWind homepages www.instantwind.com 
% Meventus AS www.meventus.com 
% CMR www.cmr.com
% Høgskolan på Gotland www.hgo.se
% OpenCFD www.opencfd.com

%% Project setup 
% I) Project setup - file system/directory structure for running the project. 
% 
% 1) Create the following directories : 

    !mkdir InstantWind
    !mkdir InstantWind/VBShearFlow
    !mkdir InstantWind/VBShearFlow/bin
    !mkdir InstantWind/ModelReductionDemo
    !mkdir InstantWind/ModelReductinoDemo/figures
    
%% 
% * InstantWind is the main directory 
% * InstantWind/VBShearFlow directory must contain the CFD results, one CFD simulation or "snapshot" of the format dat_files_1, dat_files_2, ..., dat_files_7, .. etc.   
%   Download the CFD results at <http://85.19.218.244/instantwind/openfoams.html www.instantwind.com $\rightarrow$ CFD model $\rightarrow$ OpenFoams $\rightarrow$ Download all files> into the VBShearFlow folder. 
%   The .zip file include the subfolders 
%
%       Horizontal 
%       Horizontal/bin
%       Vertical
%       Vertical/bin
%       AdditionalCases
%       AdditionalCases/bin
%   
% * InstantWind/VBSshearFlow/bin shall be used to store files created during execution of the program, i.e. .mat files such basisCollection.mat, snapshotsCollection.mat, etc. These files are described later. 
% * InstantWind/ModelReductionDemo contains the matlab and c-files.  From <http://85.19.218.244/instantwind/matlab.html www.instantwind.com $\rightarrow$ code $\rightarrow$ matlab>
%   download & unzip InstantWindMatlab.zip into this folder
%   The folder structure should be 
%   (Klaus: type in final folder structure here) 
% * InstantWind/ModelReductionDemo/figures will also contain output figures from the model runs

%% Running InstantWind
% We now describe how to run the project. 
% It is suggested that the user goes through the theory section of the report <http://85.19.218.244/instantwind/documents.html InstantWind Report> to get an understanding of the tasks involved. Some of the keywords like 'snapshots', 'basis' are described in the report.
% For a summary of the algorithm please refer to chapter 2 in the above report. 

%%  
% 1) Start matlab and change current directory to InstantWind/ModelReductionDemo.
    cd \instantwind0/ModelReductionDemo
%%    
% 2) We need to create mex files for using c functions (for RANS operators) in matlab.  At the command prompt type the following: 
    mex A_mex_check2.c
    mex V_mex_check2.c
    mex P_mex_check3.c
    mex K_mex_check3.c
%%    
% In case your matlab is not associated with a compiler, you will be
% requested by matlab to do that.  We have tested that the mex files work
% on windows 32-bit versions of Matlab.  
% NOTE: The Windows 64-bit version of Matlab requires some additonal options, but crashes when mexing.  This
% has been reported to Matlab, and will be solved later. 
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