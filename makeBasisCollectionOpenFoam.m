

 function makeBasisCollectionOpenFoam (snapshotCollectionFileName, basisCollectionFileName, snapshotList, CFD_ResultsPath, flagSB, noBasis)

	% Make a collection of basis vectors from selected snapshots and save to file
 
   fprintf(1, '\nConstructing basis\n')
 
   snapshotStruct = load(sprintf('%s%s%s', CFD_ResultsPath, 'bin\', snapshotCollectionFileName)); %'

   if (flagSB == 1)

%        1 and 3
       for i=1:2*length(snapshotList)/3
        selSnapshots(:, i) = snapshotStruct.snapshots(:,snapshotList(i));
       end

       for i=1:2*length(snapshotList)/3
        selBSnapshots(:, i) = snapshotStruct.bsnapshots(:,snapshotList(i));
       end
       
%        2
       for i=1+2*length(snapshotList)/3:length(snapshotList)
        selSnapshots2(:, i-(2*length(snapshotList)/3)) = snapshotStruct.snapshots2(:,snapshotList(i));
       end


       for i=1+2*length(snapshotList)/3:length(snapshotList)
        selBSnapshots2(:, i-2*length(snapshotList)/3) = snapshotStruct.bsnapshots2(:,snapshotList(i));
       end
       
   else
   
       for i=1:length(snapshotList)
        selSnapshots(:, i) = snapshotStruct.snapshots(:,snapshotList(i));
       end


       for i=1:length(snapshotList)
        selBSnapshots(:, i) = snapshotStruct.bsnapshots(:,snapshotList(i));
       end
   
   end
   
   ssize=size(snapshotStruct.snapshots);
   selssize=size(selSnapshots);
   
   fprintf(1, 'Start computing SVD\n');
   [basisVectors,S,V]=svd(selSnapshots,'econ');

   if (flagSB == 1)
       ssize=size(snapshotStruct.snapshots2);
       selssize=size(selSnapshots2);       
       [basisVectors2,S2,V2]=svd(selSnapshots2,'econ');       
       su2=size(basisVectors2);
       ss2=size(S2);
       sv2=size(V2);       
   end
   fprintf(1, 'Finished computing SVD\n');
   
   su=size(basisVectors);
   ss=size(S);
   sv=size(V);
   
   disp('Singular values:');
   
if (flagSB == 1)
   disp ('For tile 2');
   for i=1:su2(2)
       if (S2(i,i)/max(diag(S2)) <= 3e-4)
           break;
       end
   end
    
else
   (diag(S))'
   for i=1:su(2)
       if (S(i,i)/max(diag(S)) <= 3e-4)
           break;
       end
   end
end
   
   
   
   if (flagSB == 1)
       old=su2(2);
       basisVectors2 = basisVectors2(:,1:noBasis);
       su2=size(basisVectors2);
       fprintf(1,'Reducing the size of the basisVectors for tile2 based on their singular values from %d to %d \n',old, su2(2));
%        semilogy(1:length(diag(S2)),diag(S2))
%        xlabel('Singular Value/Basis Vector for tile2');
%        ylabel('Singular Values');
       figure();
       semilogy(1:length(diag(S2)),diag(S2), 'LineWidth', 2)
       hXLabel = xlabel('Singular Value/Basis Vector for tile2');
       hYLabel = ylabel('Singular Values');
       hTitle = title('Singular Values of tile 2 for OpenFoam');

       set([hXLabel, hYLabel], ...
           'FontName', 'AvantGarde', ...
           'FontSize', 10);
       set(hTitle, ...
           'FontName','AvantGarde',...
           'FontSize', 12, ...
           'FontWeight', 'bold')
       set(gca, ...
           'FontSize', 8);
       set(gca,'FontName','Helvetica');       
       
       basisVectors = basisVectors(:,1:noBasis);       
       figure();
%        semilogy(1:length(diag(S)),diag(S))
%        xlabel('Singular Value/Basis Vector for tile1');
%        ylabel('Singular Values');
       semilogy(1:length(diag(S)),diag(S), 'LineWidth', 2)
       hXLabel = xlabel('Singular Value/Basis Vector for tiles 1 and 3');
       hYLabel = ylabel('Singular Values');
       hTitle = title('Singular Values of tiles 1 and 3 (combined) for OpenFoam');
       set([hXLabel, hYLabel], ...
           'FontName', 'AvantGarde', ...
           'FontSize', 10);
       set(hTitle, ...
           'FontName','AvantGarde',...
           'FontSize', 12, ...
           'FontWeight', 'bold')
       set(gca, ...
           'FontSize', 8);
       set(gca,'FontName','Helvetica');

   else
       old=su(2);
       basisVectors = basisVectors(:,1:noBasis);
       su=size(basisVectors);
       fprintf(1,'Reducing the size of the basisVectors based on their singular values from %d to %d \n',old, su(2));
       semilogy(1:length(diag(S)),diag(S))
       xlabel('Singular Value/Basis Vector no');
       ylabel('Singular Values');
   end
%    keyboard;
 

   nx= snapshotStruct.nx;
   ny = snapshotStruct.ny;
   nz= snapshotStruct.nz;
   
   fprintf(1, 'Saving basis\n')
   if (flagSB == 1)
   save(sprintf('%s%s%s', CFD_ResultsPath, 'bin\', basisCollectionFileName), ...
       'selSnapshots', 'selBSnapshots', 'basisVectors', 'S',  ...
      'selSnapshots2', 'selBSnapshots2', 'basisVectors2', 'S2',  ...       
       'nx', 'ny', 'nz');       
   else   
   save(sprintf('%s%s%s', CFD_ResultsPath, 'bin\', basisCollectionFileName), ...
       'selSnapshots', 'selBSnapshots', 'basisVectors', 'S',  ...
       'nx', 'ny', 'nz');
   end
   
% %    fprintf(1, 'Saving basis\n')
% %    save(sprintf('%s%s%s', CFD_ResultsPath, 'bin\', basisCollectionFileName), ...
% %        'selSnapshots', 'selBSnapshots', 'basisVectors', 'S', 'A', ...
% %        'boundaryBasisVectors_W', 'boundaryBasisVectors_E', ...
% %        'boundaryBasisVectors_E_1', 'boundaryBasisVectors_E_2','nx', 'ny', 'nz');
    
   fprintf(1, 'Finished constructing basis\n')
   
