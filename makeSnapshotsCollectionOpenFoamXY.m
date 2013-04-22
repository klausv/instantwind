
function makeSnapshotsCollectionOpenFoamXY (snapshotCollectionFileName, CFD_ResultsPathH, CFD_ResultsPathV, CFD_JobsStruct, flagSB)

	% Read OpenFOAM results and save to file

	fprintf(1,'\n\nStart reading CFD results for snapshotsCollection\n')
    j=1;
    fprintf(1,'\n\nFirst reading Horizontal/x-direction CFD results for snapshotsCollection\n')
    dirn = 'h'; % Horizontal/x-direction
	for i=1:CFD_JobsStruct.hcount
		fprintf(1,'Reading job %s\n',CFD_JobsStruct.listh(i,:))
		[fres1, fres2, fres3] = readFoamResultXY(CFD_ResultsPathH, dirn, sprintf('%s',CFD_JobsStruct.listh(i,:)));
		fprintf(1,'Job %s has data count %d\n',sprintf('%s',CFD_JobsStruct.listh(i,:)),fres1.dcount)
		snapshots(:,j)=fres1.data; % Tile 1
        j=j+1;
        if (flagSB == 1)
        snapshots2(:,i)=fres2.data; % Tile 2
        else
        snapshots(:,j)=fres2.data; % Tile 2
        j=j+1;            
        end
        snapshots(:,j)=fres3.data; % Tile 3
        j=j+1;
		s=size(snapshots);
		fprintf(1,'Size snapshots: [%d %d]\n\n',s(1),s(2))
    end

    fprintf(1,'\n\nNext reading Vertical/y-direction CFD results for snapshotsCollection\n')    
    dirn = 'v'; % Vertical/y-direction    
	for i=1:CFD_JobsStruct.vcount
		fprintf(1,'Reading job %s\n',CFD_JobsStruct.listv(i,:))
		[fres1, fres2, fres3] = readFoamResultXY(CFD_ResultsPathV, dirn, sprintf('%s',CFD_JobsStruct.listv(i,:)));
		fprintf(1,'Job %s has data count %d\n',sprintf('%s',CFD_JobsStruct.listv(i,:)),fres1.dcount)
		snapshots(:,j)=fres1.data; % Tile 1
        j=j+1;
        if (flagSB == 1)
        snapshots2(:,i+CFD_JobsStruct.hcount)=fres2.data; % Tile 2
        else
        snapshots(:,j)=fres2.data; % Tile 2
        j=j+1;            
        end
        snapshots(:,j)=fres3.data; % Tile 3
        j=j+1;
		s=size(snapshots);
		fprintf(1,'Size snapshots: [%d %d]\n\n',s(1),s(2))
    end    
    nx = fres1.nx;
    ny = fres1.ny;
    nz = fres1.nz;
    nx1 = nx/2;
    ny_by_2 = ny/2;
       
    ssize=size(snapshots);
	%Left boundary components of snapshots
	%U-V-W-NU-components
	for n=1:ssize(2)
		for comp=1:6
			for i=1:nz
				for j=1:ny/2
					bsnapshots(j+(i-1)*ny_by_2+(comp-1)*ny_by_2*nz,n)=snapshots(1+(j-1)*nx1+(i-1)*nx1*ny_by_2+(comp-1)*nx1*ny_by_2*nz,n);
				end
			end
		end
    end

if (flagSB == 1)
    ssize2=size(snapshots2);

    for n=1:ssize2(2)
		for comp=1:6
			for i=1:nz
				for j=1:ny/2
					bsnapshots2(j+(i-1)*ny_by_2+(comp-1)*ny_by_2*nz,n)=snapshots2(1+(j-1)*nx1+(i-1)*nx1*ny_by_2+(comp-1)*nx1*ny_by_2*nz,n);
				end
			end
		end
    end
end
    
    nx = nx1;
    ny = ny/2;

    if (flagSB == 1)
        save(sprintf('%s%s%s', CFD_ResultsPathH, 'bin\', snapshotCollectionFileName), 'snapshots', 'bsnapshots','snapshots2', 'bsnapshots2', 'nx', 'ny', 'nz')
    else
    	save(sprintf('%s%s%s', CFD_ResultsPathH, 'bin\', snapshotCollectionFileName), 'snapshots', 'bsnapshots', 'nx', 'ny', 'nz')
    end