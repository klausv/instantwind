function [ProdEstim U V W BoundaryError RANSError] = modelReductionSolverOpenFoamXY(basisCollectionFileName, frontTileFileName, CFD_ResultsPath, RANS_OperatorsFileName, alpha, locx, locy, flagSB, noBasis) %loc = 0...31

%Solve for the flow field by model reduction when the back row turbine is positioned at grid position loc

   % Here both the tiles share the same basis of solution vectors
   % Hence we just have boundaryBasisVectors1_L and boundaryBasisVectors1_R
   % for both the front tile (1) and back tile (2).
% % Tile set-up
%       front (1)  back(2)
%     |---------|
%     |    A    |---------|
%     |_________|         | 
%     |         |    C    |
%     |    B    |---------|
%     |---------| 


% %     Orientation for basis :
% Y  ny ____N____ 
% ^    |         |
% |    W         E
% |  1 |____S____| ---> X axis
%      1        nx 

% Note all tiles have the same basis:
% B_A = B_B = B_C
%     |---------|
%     |  A(a1)  |---------|
%     |_______E1|W2       | 
%     |       E2|W1 C(a2) |
%     |  B(a3)  |---------|
%     |---------| 

	if nargin==7
		flagSB = 1;
        noBasis = 5;
    elseif nargin==8
        noBasis = 5;
	elseif nargin~=9
		fprintf(1, '\nFunction modelReductionSolver: Unknown usage')
		return;
    end
	
    if (alpha == 0) 
        fprintf(1, '\n Reducing Boundary Error Only')
    elseif (alpha == 1) 
        fprintf(1, '\n Reducing Momentum Error Only')
    else
        fprintf(1, '\n Adrresing both Boundary and Momentum Error')
    end
    
	fprintf(1, '\nStart modelReductionSolver, location = %d , %d\n', locx,locy)

	% Grid cell size
	delta = 3;

	tStart=tic;

% %     Boundary matching part:
% % 

	basisCollectionStruct = load(sprintf('%s%s%s', CFD_ResultsPath, 'bin\', basisCollectionFileName)); %'
	selssize=size(basisCollectionStruct.basisVectors);
    
    if (flagSB == 1)
        boundaryBasisCollectionFileName = 'boundaryBasisCollection_4_openFoam_SB.mat';
        boundaryBasisComputationXY(basisCollectionFileName, CFD_ResultsPath,  boundaryBasisCollectionFileName, locx, locy, flagSB);
        boundaryBasisCollectionStruct = load(sprintf('%s%s%s', CFD_ResultsPath, 'bin\', boundaryBasisCollectionFileName)); %'
    else
        boundaryBasisCollectionFileName = 'boundaryBasisCollection_4_openFoam_WS.mat';
        boundaryBasisComputationXY(basisCollectionFileName, CFD_ResultsPath,  boundaryBasisCollectionFileName, locx, locy, flagSB);
        boundaryBasisCollectionStruct = load(sprintf('%s%s%s', CFD_ResultsPath, 'bin\', boundaryBasisCollectionFileName)); %'
    end

	frontTileStruct = load(sprintf('%s%s%s', CFD_ResultsPath, 'bin\', frontTileFileName)); %'

	%Extract the ingoing boundary conditions
	%from the front-tile to the seccond tile
	%according to the position indicator loc
	%Center tile is starting from y-index 31
    % 	yStartIndex=31+loc; VB commented out
    yStartIndex = 1;
    bvalue = zeros(basisCollectionStruct.ny*basisCollectionStruct.nz*6,1);
    bvalue2 = zeros(basisCollectionStruct.ny*basisCollectionStruct.nz*6,1);
	for comp=1:6
		for i=1:basisCollectionStruct.nz
			for j=1:basisCollectionStruct.ny
				bvalue(j+(i-1)*basisCollectionStruct.ny+(comp-1)*basisCollectionStruct.ny*basisCollectionStruct.nz,1)=frontTileStruct.tileA_left_xsec(j+(yStartIndex-1),i,comp);
                bvalue2(j+(i-1)*basisCollectionStruct.ny+(comp-1)*basisCollectionStruct.ny*basisCollectionStruct.nz,1)=frontTileStruct.tileB_left_xsec(j+(yStartIndex-1),i,comp);
			end
		end
    end

%     VB

    tile1 = boundaryBasisCollectionStruct.boundaryBasisVectors_W'*bvalue;
    tile3 = boundaryBasisCollectionStruct.boundaryBasisVectors_W'*bvalue2;
    a1 =zeros(size(boundaryBasisCollectionStruct.boundaryBasisVectors_W,2),1);
    a2 =zeros(length(a1),1);%zeros(size(basisCollectionStruct.boundaryBasisVectors1_L,2),1);
    a3 =zeros(length(a1),1);

    rhs = zeros(length(a1)+length(a2)+length(a3),1);
    rhs(1:length(a1),1) = tile1;
    rhs(1+length(a1)+length(a2):end,1) = tile3;
    [aa_0] = inv(boundaryBasisCollectionStruct.A)*rhs;

   

 aa = aa_0;
 func = (aa'*boundaryBasisCollectionStruct.A*aa) - (2*aa'*rhs)+ ((bvalue'*bvalue)+(bvalue2'*bvalue2));
 disp(['Func:',num2str(func)]);

 [func, BoundaryError, RANSError, grad, hess] = combinedObjectiveFunctionOpenFoam ...
                     (aa, boundaryBasisCollectionStruct.A,  CFD_ResultsPath, RANS_OperatorsFileName, rhs, alpha,bvalue,bvalue2, flagSB);
RANSError = RANSError^2;                 
 disp(['Initial Objective function value:',num2str(func)])

    a1 = aa(1:length(a1),1);
    a2 = aa(1+length(a1):(length(a1)+length(a2)),1);
    a3 = aa(1+length(a1)+length(a2):end,1);
      
	% Start plotting
	close all
% 	iz=20;
 	iz=36; % 24*3 = 72 position of turbine in z direction 

    %Find the corresponding solution for all the tiles
    estimSol1=basisCollectionStruct.basisVectors*a1;
    if (flagSB == 1)
        estimSol2=basisCollectionStruct.basisVectors2*a2;
    else
        estimSol2=basisCollectionStruct.basisVectors*a2;
    end
    estimSol3=basisCollectionStruct.basisVectors*a3;
    
    
	for plotComponent=1:1
	

		%Horizontal section of component [1, 2, 3, 4]=[u, v, w, mu] from back-tile  at z-level iz
		low = 1+basisCollectionStruct.nx*basisCollectionStruct.ny*(iz-1)+(plotComponent-1)*basisCollectionStruct.nx*basisCollectionStruct.ny*basisCollectionStruct.nz;
		high= basisCollectionStruct.nx*basisCollectionStruct.ny*iz+(plotComponent-1)*basisCollectionStruct.nx*basisCollectionStruct.ny*basisCollectionStruct.nz;
		zvecEstimSol1=estimSol1(low:high);
        zvecEstimSol2=estimSol2(low:high);
		zvecEstimSol3=estimSol3(low:high);
		
		for i=1:basisCollectionStruct.ny
			zsecEstimSol1(:,i)=zvecEstimSol1(1+(i-1)*basisCollectionStruct.nx:i*basisCollectionStruct.nx);
            zsecEstimSol2(:,i)=zvecEstimSol2(1+(i-1)*basisCollectionStruct.nx:i*basisCollectionStruct.nx);
            zsecEstimSol3(:,i)=zvecEstimSol3(1+(i-1)*basisCollectionStruct.nx:i*basisCollectionStruct.nx);            
		end

%         if plotComponent==1
% 			cv=5.6:0.15:8.4;
% 		elseif plotComponent==2
% 			cv=-0.7:0.08:0.6;
% 		elseif plotComponent==3
% 			cv=-0.04:0.02:0.06;
% 		elseif plotComponent==4
% 			cv=15:0.4:20;
%         elseif plotComponent==5
% 			cv=-11:1.5:12;
%         else
%             cv=1.2:0.08:2.1;
%         end
        if plotComponent==1
			cv=6:0.15:10;
		elseif plotComponent==2
			cv=-0.9:0.08:0.6;
		elseif plotComponent==3
			cv=-1.4:0.02:0.04;
		elseif plotComponent==4
			cv=15:0.5:28;
        elseif plotComponent==5
			cv=-15:1.5:15;
        else
            cv=0.6:0.08:2.2;
        end

        
		figure;
        contourf((1:basisCollectionStruct.nx)*3,(1:basisCollectionStruct.ny)*3,zsecEstimSol1(:,:)',cv); % til1
        hold on
        contourf((1:basisCollectionStruct.nx)*3,(basisCollectionStruct.ny+1:2*basisCollectionStruct.ny)*3,zsecEstimSol3(:,:)',cv); % tile2
        hold on
%         contourf((basisCollectionStruct.nx:2*basisCollectionStruct.nx-1)*3,(1+loc+basisCollectionStruct.ny/2:loc+(basisCollectionStruct.ny/2+basisCollectionStruct.ny))*3,zsecEstimSol2(:,:)',cv); % tile2
        contourf((basisCollectionStruct.nx-locx:2*basisCollectionStruct.nx-1-locx)*3,(1+locy+basisCollectionStruct.ny/2:(locy+basisCollectionStruct.ny/2+basisCollectionStruct.ny))*3,zsecEstimSol2(:,:)',cv); % tile2
        axis([0 2*delta*basisCollectionStruct.nx 0 2*delta*basisCollectionStruct.ny])
        
		titletxt=sprintf('ModRedSol flow component %d, back row turbine production: %6.2f KW', plotComponent, computeTurbineEffectInTileUsingCpCurveFoam(estimSol2, basisCollectionStruct.nx, basisCollectionStruct.ny, basisCollectionStruct.nz));
		title(titletxt);
		axis('tight');
		cb=colorbar;
		if plotComponent==1
			set(get(cb,'ylabel'), 'String', 'Wind speed in x-direction (m/s)', 'Color', 'b');
		elseif plotComponent==2
			set(get(cb,'ylabel'), 'String', 'Wind speed in y-direction (m/s)', 'Color', 'b');
		elseif plotComponent==3
			set(get(cb,'ylabel'), 'String', 'Wind speed in z-direction (m/s)', 'Color', 'b');
		elseif plotComponent==4
			set(get(cb,'ylabel'), 'String', 'Total viscosity (Turbulent + molecular) (Pa\cdot s)', 'Color', 'b');
		elseif plotComponent==5
			set(get(cb,'ylabel'), 'String', 'Pressure (Pa\cdot s)', 'Color', 'b');
        else
            set(get(cb,'ylabel'), 'String', 'Turbulent Kinetic Energy (Pa\cdot s)', 'Color', 'b');            
		end
		
		xlabel('Length (x) [m]')
		ylabel('Width (y) [m]')
		hold off

		fname=sprintf('figures/movingTurbine-b%d-c%d-px%d-py%d-alpha-%d.png', noBasis, plotComponent,locx,locy,alpha);
		print(fname,'-dpng');
       
    end
    ProdEstim = computeTurbineEffectInTileUsingCpCurveFoam(estimSol2, basisCollectionStruct.nx, basisCollectionStruct.ny, basisCollectionStruct.nz);
%     Volumetric data for U
    plotComponent = 1;
    low = 1+(plotComponent-1)*basisCollectionStruct.nx*basisCollectionStruct.ny*basisCollectionStruct.nz;
    high= (plotComponent)*basisCollectionStruct.nx*basisCollectionStruct.ny*basisCollectionStruct.nz;    
    U = estimSol2(low:high);

%     Volumetric data for V
    plotComponent = 2;    
    low = 1+(plotComponent-1)*basisCollectionStruct.nx*basisCollectionStruct.ny*basisCollectionStruct.nz;
    high= (plotComponent)*basisCollectionStruct.nx*basisCollectionStruct.ny*basisCollectionStruct.nz;        
    V = estimSol2(low:high);
    
%     Volumetric data for V
    plotComponent = 3;    
    low = 1+(plotComponent-1)*basisCollectionStruct.nx*basisCollectionStruct.ny*basisCollectionStruct.nz;
    high= (plotComponent)*basisCollectionStruct.nx*basisCollectionStruct.ny*basisCollectionStruct.nz;        
    W = estimSol2(low:high);
  
	tElapsed = toc(tStart);
	fprintf(1, 'Time spent for reading data, calculations, and plotting is %f seconds\n\n',tElapsed);
    