
function readFrontTileOpenFoam(frontTileFileName, frontTileJob, CFD_ResultsPath)

	% Read a OpenFOAM result for the front-tile and save to binary file.
    % Left cross-section is constant for all positions of back tile.

	fprintf(1,'\nStart reading foam result for front-tile\n')
    fres = load (sprintf('%s%s%s', CFD_ResultsPath, frontTileJob, '\', frontTileJob));


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
    
    nx = fres.foamResult1.nx/2;
    ny = fres.foamResult1.ny/2;
    nz = fres.foamResult1.nz;
    ny1 = ny;

    tileA_right_xsec = zeros(ny1,nz,6);
    tileB_right_xsec = zeros(ny1,nz,6);
    tileA_left_xsec = zeros(ny1,nz,6);
    tileB_left_xsec = zeros(ny1,nz,6);               

	%Right boundary components of fres.data (column vector)
	%U-V-W-NU-components
	%Read into a matrix
    for comp=1:6
		for i=1:nz
			for j=1:ny1
% 				right_xsec(j,i,comp)=tmp_data(nx+(j-1)*nx+(i-1)*nx*ny1+(comp-1)*nx*ny1*nz);
                tileA_right_xsec(j,i,comp)=fres.foamResult1.data(nx+(j-1)*nx+(i-1)*nx*ny1+(comp-1)*nx*ny1*nz);
                tileB_right_xsec(j,i,comp)=fres.foamResult3.data(nx+(j-1)*nx+(i-1)*nx*ny1+(comp-1)*nx*ny1*nz);
			end
		end
    end
    

	%Left boundary components of fres.data (column vector)
	%U-V-W-NU-components
	%Read into a matrix
    for comp=1:6
		for i=1:nz
			for j=1:ny1
                tileA_left_xsec(j,i,comp)=fres.foamResult1.data(1+(j-1)*nx+(i-1)*nx*ny1+(comp-1)*nx*ny1*nz);
                tileB_left_xsec(j,i,comp)=fres.foamResult3.data(1+(j-1)*nx+(i-1)*nx*ny1+(comp-1)*nx*ny1*nz);
			end
		end
    end    

    
    ny = ny1;

	save(sprintf('%s%s%s', CFD_ResultsPath, 'bin\', frontTileFileName), ...
         'tileA_left_xsec', 'tileB_left_xsec', ...
         'tileA_right_xsec', 'tileB_right_xsec', 'nx', 'ny', 'nz');

	%figure;
	%contour(right_xsec(:,:,1));