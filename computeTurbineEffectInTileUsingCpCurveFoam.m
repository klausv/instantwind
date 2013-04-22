
function [production, vxy]=computeTurbineEffectInTileUsingCpCurveFoam(solVec, nx, ny, nz)

%  FLACS Results
% % 	turbinePosXmeters=101.0;
% % 	turbinePosYmeters=75.0;
% % 	hubHeightmeters=50.0;
% % 	rotorDiameter=37.0;
% % 	gridDelta=2.5;
	turbinePosXmeters=90.0;
	turbinePosYmeters=54.0;
	hubHeightmeters=108.0;
	rotorDiameter=36.0;
	gridDelta=3.0;


	cp=[0.0 0.0;
		1.0 0.0;
		2.0 0.0;
		3.0 0.0;
		4.0 0.0;
		5.0 24.9;
		6.0 63.9;
		7.0 109.0;
		8.0 162.4;
		9.0 217.8;
		10.0 272.0;
		11.0 331.1;
		12.0 373.3;
		13.0 423.5;
		14.0 454.7;
		15.0 466.9;
		16.0 476.2;
		17.0 469.7;
		18.0 435.6;
		19.0 425.0;
		20.0 414.9;
		21.0 410.0;
		22.0 405.0;
		23.0 400.0;
		24.0 395.0;
		25.0 390.0];

	
    ix = round((turbinePosXmeters-rotorDiameter/2)/gridDelta);
	iy = round(turbinePosYmeters/gridDelta);
	iz = round(hubHeightmeters/gridDelta);
	
	comp=1;
	u = solVec(ix+(iy-1)*nx+(iz-1)*nx*ny+(comp-1)*nx*ny*nz);
	comp=2;
	v = solVec(ix+(iy-1)*nx+(iz-1)*nx*ny+(comp-1)*nx*ny*nz);
	
	vxy = sqrt(u^2+v^2);
	production = interp1(cp(:,1),cp(:,2),vxy);
	