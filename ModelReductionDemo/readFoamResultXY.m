function [foamResult1, foamResult2, foamResult3] = readFoamResultXY(basePath, dirn, jobNo, dirName)

	%Read components U(ux,ux,uz), p, k, epsilon from FOAM ascii result files
	
	if (nargin~=3 && nargin~=2)
		fprintf(1,'\n\nUsage 1 : readFlacsResult (basePath, jobNo, dirName\n');
		fprintf(1,'Usage 2 : readFlacsResult (basePath, jobNo\n\n');
		flacsResult=-1;
	else
		if nargin==3
			dirName=jobNo;
		elseif nargin == 4
			dirName = sprintf('%s%s%s', jobNo, '\', dirName); %'
		end

		fname_U = sprintf('%s%s%s%s%s',basePath, dirName,'\U_dat');
		fname_p = sprintf('%s%s%s%s%s',basePath, dirName,'\p_dat');
		fname_k = sprintf('%s%s%s%s%s',basePath, dirName,'\k_dat');
		fname_e = sprintf('%s%s%s%s%s',basePath, dirName,'\epsilon_dat');

		foamRes_U = readFoamResultFile (fname_U);
		foamRes_p = readFoamResultFile (fname_p);
		foamRes_k = readFoamResultFile (fname_k);
		foamRes_e = readFoamResultFile (fname_e);

  		foamResult1.dcount = foamRes_U.dcount;
  		foamResult1.nx = foamRes_U.nx;
        foamResult1.ny = foamRes_U.ny;
        foamResult1.nz = foamRes_U.nz;
        foamResult2 = foamResult1; % foamResult2 initialized
        foamResult3 = foamResult1; % foamResult3 initialized

%   Computing total effective viscosity 
%   Total effective viscosity = molecular viscosity + turbulent viscosity
%   Molecular viscosity must be supplied in a separate data file or 
%   maybe from one of the other data files - tke or epsilon?
%   Turbulent viscosity = c_mu*(k^2/epsilon)

        foamRes_mu_eff.data = zeros(size(foamRes_e.data));
        c_mu = 0.09;
        molecular_viscosity = (1.5e-5);
        foamRes_mu_eff.data = c_mu.*(foamRes_k.data.^2)./foamRes_e.data;
        foamRes_mu_eff.data = foamRes_mu_eff.data + molecular_viscosity;

        nx=foamResult1.nx;
        ny=foamResult1.ny;
        nz=foamResult1.nz;
%       Dimensions of the tiles - halved in x and y directions from the complete field        
        nx_by_2=nx/2;
        ny_by_2=ny/2;
        
% % Tile set-up
%       front (1)  back(2)
%     |---------|
%     |    A(1) |---------|  Y axis
%     |_________|         |   ^
%     |         |    C(2) |   |
%     |    B(3) |---------|   |
%     |---------|             |---> X axis



%     |---------|
%     |  A(a1)  |---------|
%     |_______E1|W2       | 
%     |       E2|W1 C(a2) |
%     |  B(a3)  |---------|
%     |---------| 
%       Taking the appropriate y-range for each of the three tiles. 
%       Will vary for tile3.

        digit=eval(jobNo(end));
        j3=0; % tile3
        j1=ny/2; % tile1

        if (dirn == 'h' || dirn == 'H')
            j2=(ny/4);%+(digit-1)*3; % tile2
        elseif (dirn == 'v' || dirn == 'V')
            j2=(ny/4)+(digit-1)*3; % tile2
        end
        
%       X-Location of front turbines(1 and 3) - 30
%       X-Location of back turbine (2)        - 90+(digit-5)*3
%       X-range for tile 1 and 3              - 1:60
%       To keep relative position of turbine wrt tiles fixed.
        if (dirn == 'h' || dirn == 'H')
            i3 = (nx/2) +(digit-5)*3;
        elseif (dirn == 'v' || dirn == 'V')
            i3 = nx_by_2;
        end
                      
        nnz=nx_by_2*ny_by_2*nz; % Total number of non-zeros in each tile
        
        tile1=zeros(nnz,6);
        tile2=zeros(nnz,6);
        tile3=zeros(nnz,6);
        
        
            for k=1:nz
                for j=1:ny_by_2
                    for i=1:nx_by_2
% %  For the three components of velocity                                                    
                        for ss = 1:3
%                 Velocity                            
                tile3(i+(j-1)*nx_by_2+(k-1)*nx_by_2*ny_by_2,ss)=foamRes_U.data(i + (j+j3-1)*nx + (k-1)*nx*ny,ss);
                tile1(i+(j-1)*nx_by_2+(k-1)*nx_by_2*ny_by_2,ss)=foamRes_U.data(i + (j+j1-1)*nx + (k-1)*nx*ny,ss);
                tile2(i+(j-1)*nx_by_2+(k-1)*nx_by_2*ny_by_2,ss)=foamRes_U.data((i3+i) + (j+j2-1)*nx + (k-1)*nx*ny,ss);
                        end
%                Pressure        
                tile3(i+(j-1)*nx_by_2+(k-1)*nx_by_2*ny_by_2,ss+2)=foamRes_p.data(i + (j+j3-1)*nx + (k-1)*nx*ny);
                tile1(i+(j-1)*nx_by_2+(k-1)*nx_by_2*ny_by_2,ss+2)=foamRes_p.data(i + (j+j1-1)*nx + (k-1)*nx*ny);
                tile2(i+(j-1)*nx_by_2+(k-1)*nx_by_2*ny_by_2,ss+2)=foamRes_p.data((i3+i) + (j+j2-1)*nx + (k-1)*nx*ny);
%                 TKE
                tile3(i+(j-1)*nx_by_2+(k-1)*nx_by_2*ny_by_2,ss+3)=foamRes_k.data(i + (j+j3-1)*nx + (k-1)*nx*ny);
                tile1(i+(j-1)*nx_by_2+(k-1)*nx_by_2*ny_by_2,ss+3)=foamRes_k.data(i + (j+j1-1)*nx + (k-1)*nx*ny);
                tile2(i+(j-1)*nx_by_2+(k-1)*nx_by_2*ny_by_2,ss+3)=foamRes_k.data((i3+i) + (j+j2-1)*nx + (k-1)*nx*ny);   
%                 nu_effective
                tile3(i+(j-1)*nx_by_2+(k-1)*nx_by_2*ny_by_2,ss+1)=foamRes_mu_eff.data(i + (j+j3-1)*nx + (k-1)*nx*ny);
                tile1(i+(j-1)*nx_by_2+(k-1)*nx_by_2*ny_by_2,ss+1)=foamRes_mu_eff.data(i + (j+j1-1)*nx + (k-1)*nx*ny);
                tile2(i+(j-1)*nx_by_2+(k-1)*nx_by_2*ny_by_2,ss+1)=foamRes_mu_eff.data((i3+i) + (j+j2-1)*nx + (k-1)*nx*ny);                                        
                    end
                end
            end
        end
        
		foamResult1.data=[tile1(:,1);tile1(:,2);tile1(:,3);tile1(:,4);tile1(:,5);tile1(:,6)];
		foamResult2.data=[tile2(:,1);tile2(:,2);tile2(:,3);tile2(:,4);tile2(:,5);tile2(:,6)];
		foamResult3.data=[tile3(:,1);tile3(:,2);tile3(:,3);tile3(:,4);tile3(:,5);tile3(:,6)];

        % save the field to be read in for the front tile
        if(digit == 1)
            save(sprintf('%s%s%s%s%s',basePath, dirName, '\',dirName), 'foamResult1', 'foamResult2', 'foamResult3')
        end

end