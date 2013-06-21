
function boundaryBasisComputationXY(basisCollectionFileName, CFD_ResultsPath,  boundaryBasisCollectionFileName, locx, locy, flagSB)

fprintf(1, '\nComputing Boundary Basis\n')   
basisCollectionStruct = load(sprintf('%s%s%s', CFD_ResultsPath, 'bin\', basisCollectionFileName)); %'

if (flagSB == 1)
   basisVectors = basisCollectionStruct.basisVectors;
   basisVectors2 = basisCollectionStruct.basisVectors2;
else
   basisVectors = basisCollectionStruct.basisVectors;
end
   nx= basisCollectionStruct.nx;
   ny= basisCollectionStruct.ny;
   nz= basisCollectionStruct.nz;
   su=size(basisVectors);
   
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



%    Because the tile 1 is split into two - tile A and tile B
%    and for the east and west boundary of tile 1 and 2 respectively  
% 
%    yloc=0
%     |---------|
%     |  A(a1)  |---------|
%     |_______E1|W2       | 
%     |       E2|W1 C(a2) |
%     |  B(a3)  |---------|
%     |---------| 
% 
%    yloc~=0
%     |---------|_ _ _ _ _
%     |  A(a1)  |W2       |
%     |_______E1|   C(a2) | 
%     |       E2|_W1_ _ _ |
%     |  B(a3)  |         
%     |---------| 



% To avoid indexing errors when mod(ny,2)~=0 we do:
   
   y1_begin = 1;
   if (mod(ny,2)~=0)
       y1_end = (floor(ny/2))+locy;   %For E1
       y1_Wend = y1_end+1-(2*locy);   %For W1
   else 
       y1_end =  (ny/2)+locy;         %For E1
       y1_Wend = y1_end-(2*locy);     %For W1
   end
   y2_begin = y1_end + 1;
   y2_Wbegin = y1_Wend + 1;
   y2_end = ny;
   y2_diff = y2_end - y2_begin + 1;
   y2_Wdiff = y2_end - y2_Wbegin + 1;
   
   boundaryBasisVectors_W = zeros(ny*nz*6,su(2));
   boundaryBasisVectors_E = zeros(ny*nz*6,su(2));
   boundaryBasisVectors_W_1 = zeros(y1_Wend*nz*6,su(2));
   boundaryBasisVectors_W_2 = zeros(y2_Wdiff*nz*6,su(2));
   boundaryBasisVectors_E_1 = zeros(y1_end*nz*6,su(2));
   boundaryBasisVectors_E_2 = zeros(y2_diff*nz*6,su(2));
   boundaryBasisVectors_N = zeros(nx*nz*6,su(2));
   boundaryBasisVectors_S = zeros(nx*nz*6,su(2));

if (flagSB == 1)
   %East/Right boundary components of basis vectors vectors
   %U-V-W-NU-P-TKE components
   %West/Left boundary components of basis vectors vectors
   %U-V-W-NU-P-TKE components
   for n=1:su(2)
     for comp=1:6
       for i=1:nz
		 for j=1:ny
             
		   boundaryBasisVectors_E(j+(i-1)*ny+(comp-1)*ny*nz,n)=basisVectors(nx+(j-1)*nx+(i-1)*nx*ny+(comp-1)*nx*ny*nz,n);             
		   boundaryBasisVectors_W(j+(i-1)*ny+(comp-1)*ny*nz,n)=basisVectors(1+(j-1)*nx+(i-1)*nx*ny+(comp-1)*nx*ny*nz,n);
% 		   boundaryBasisVectors_E(j+(i-1)*ny+(comp-1)*ny*nz,n)=basisVectors(nx+(j-1)*nx+(i-1)*nx*ny+(comp-1)*nx*ny*nz,n);           
           if (j <= y1_Wend) % To avoid indexing errors when mod(ny,2)~=0
           boundaryBasisVectors_W_1(j+(i-1)*y1_Wend+(comp-1)*y1_Wend*nz,n)=basisVectors2(1+(j-1)*nx+(i-1)*nx*ny+(comp-1)*nx*ny*nz,n);
%            boundaryBasisVectors_E_1(j+(i-1)*y1_end+(comp-1)*ny*nz,n)=basisVectors(nx+(j-1)*nx+(i-1)*nx*ny+(comp-1)*nx*ny*nz,n);
           end
           if (j >= y2_Wbegin)
           boundaryBasisVectors_W_2((j-y2_Wbegin+1)+(i-1)*y2_Wdiff+(comp-1)*y2_Wdiff*nz,n)=basisVectors2(1+(j-1)*nx+(i-1)*nx*ny+(comp-1)*nx*ny*nz,n);    
%            boundaryBasisVectors_E_2((j-y2_begin+1)+(i-1)*y2_diff+(comp-1)*ny*nz,n)=basisVectors(nx+(j-1)*nx+(i-1)*nx*ny+(comp-1)*nx*ny*nz,n);               
           end
           if (j <= y1_end)
           boundaryBasisVectors_E_1(j+(i-1)*y1_end+(comp-1)*y1_end*nz,n)=basisVectors(nx-locx+(j-1)*nx+(i-1)*nx*ny+(comp-1)*nx*ny*nz,n);
           end
           if (j >= y2_begin)
           boundaryBasisVectors_E_2((j-y2_begin+1)+(i-1)*y2_diff+(comp-1)*y2_diff*nz,n)=basisVectors(nx-locx+(j-1)*nx+(i-1)*nx*ny+(comp-1)*nx*ny*nz,n);    
           end 
           
		 end
       end
     end
   end
    
else   
   %East/Right boundary components of basis vectors vectors
   %U-V-W-NU-P-TKE components
   %West/Left boundary components of basis vectors vectors
   %U-V-W-NU-P-TKE components
   for n=1:su(2)
     for comp=1:6
       for i=1:nz
		 for j=1:ny
             
		   boundaryBasisVectors_E(j+(i-1)*ny+(comp-1)*ny*nz,n)=basisVectors(nx+(j-1)*nx+(i-1)*nx*ny+(comp-1)*nx*ny*nz,n);             
		   boundaryBasisVectors_W(j+(i-1)*ny+(comp-1)*ny*nz,n)=basisVectors(1+(j-1)*nx+(i-1)*nx*ny+(comp-1)*nx*ny*nz,n);
% 		   boundaryBasisVectors_E(j+(i-1)*ny+(comp-1)*ny*nz,n)=basisVectors(nx+(j-1)*nx+(i-1)*nx*ny+(comp-1)*nx*ny*nz,n);           
           if (j <= y1_Wend) % To avoid indexing errors when mod(ny,2)~=0
           boundaryBasisVectors_W_1(j+(i-1)*y1_Wend+(comp-1)*y1_Wend*nz,n)=basisVectors(1+(j-1)*nx+(i-1)*nx*ny+(comp-1)*nx*ny*nz,n);
%            boundaryBasisVectors_E_1(j+(i-1)*y1_end+(comp-1)*ny*nz,n)=basisVectors(nx+(j-1)*nx+(i-1)*nx*ny+(comp-1)*nx*ny*nz,n);
           end
           if (j >= y2_Wbegin)
           boundaryBasisVectors_W_2((j-y2_Wbegin+1)+(i-1)*y2_Wdiff+(comp-1)*y2_Wdiff*nz,n)=basisVectors(1+(j-1)*nx+(i-1)*nx*ny+(comp-1)*nx*ny*nz,n);    
%            boundaryBasisVectors_E_2((j-y2_begin+1)+(i-1)*y2_diff+(comp-1)*ny*nz,n)=basisVectors(nx+(j-1)*nx+(i-1)*nx*ny+(comp-1)*nx*ny*nz,n);               
           end
           if (j <= y1_end)
           boundaryBasisVectors_E_1(j+(i-1)*y1_end+(comp-1)*y1_end*nz,n)=basisVectors(nx-locx+(j-1)*nx+(i-1)*nx*ny+(comp-1)*nx*ny*nz,n);
           end
           if (j >= y2_begin)
           boundaryBasisVectors_E_2((j-y2_begin+1)+(i-1)*y2_diff+(comp-1)*y2_diff*nz,n)=basisVectors(nx-locx+(j-1)*nx+(i-1)*nx*ny+(comp-1)*nx*ny*nz,n);    
           end 
           
		 end
       end
     end
   end
   
end
%   VB   
   %North/South boundary components of basis vectors vectors
   %U-V-W-NU-components
   for n=1:su(2)
     for comp=1:6
       for i=1:nz
		 for j=1:nx
		   boundaryBasisVectors_N(j+(i-1)*nx+(comp-1)*nx*nz,n)=basisVectors(j+(nx*(ny-1))+(i-1)*nx*ny+(comp-1)*nx*ny*nz,n);
           boundaryBasisVectors_S(j+(i-1)*nx+(comp-1)*nx*nz,n)=basisVectors(j+(i-1)*nx*ny+(comp-1)*nx*ny*nz,n);
		 end
       end
     end
   end   

   
   A_11 = (boundaryBasisVectors_E_1'*boundaryBasisVectors_E_1) + ...
          (boundaryBasisVectors_S'*boundaryBasisVectors_S) + ...
          (boundaryBasisVectors_W'*boundaryBasisVectors_W);
   A_12 = -(boundaryBasisVectors_E_1'*boundaryBasisVectors_W_2);
   A_21 = -(boundaryBasisVectors_W_2'*boundaryBasisVectors_E_1);
   A_13 = -(boundaryBasisVectors_S'*boundaryBasisVectors_N);
   A_31 = -(boundaryBasisVectors_N'*boundaryBasisVectors_S);
   A_22 = (boundaryBasisVectors_W_2'*boundaryBasisVectors_W_2) + ...
          (boundaryBasisVectors_W_1'*boundaryBasisVectors_W_1);
   A_23 = -(boundaryBasisVectors_W_1'*boundaryBasisVectors_E_2);
   A_32 = -(boundaryBasisVectors_E_2'*boundaryBasisVectors_W_1);
   A_33 = (boundaryBasisVectors_E_2'*boundaryBasisVectors_E_2) + ...
          (boundaryBasisVectors_N'*boundaryBasisVectors_N) + ...
          (boundaryBasisVectors_W'*boundaryBasisVectors_W);

   A = [A_11 A_12 A_13; A_21 A_22 A_23; A_31 A_32 A_33];
   
   fprintf(1, 'Saving boundary basis\n')
   save(sprintf('%s%s%s', CFD_ResultsPath, 'bin\', boundaryBasisCollectionFileName), ...
        'A', 'boundaryBasisVectors_W', 'boundaryBasisVectors_E', ...
       'boundaryBasisVectors_E_1', 'boundaryBasisVectors_E_2');
    