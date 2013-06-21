function [ObjFuncVal, boundaryMatchingValue, reducedRANSValue, dF_dm, dF_dm_dn ]= combinedObjectiveFunctionOpenFoam(aa, A,  CFD_ResultsPath, RANS_OperatorsFileName, rhs, alpha, bvalue, bvalue2, flagSB)

    if nargin==8
		alpha=0;
	elseif nargin~=9
		fprintf(1, '\n Function modelReductionSolver: Unknown usage')
		return;
    end
    
RANS_MatrixStruct =  load(sprintf('%s%s%s', CFD_ResultsPath, 'bin\', RANS_OperatorsFileName)); 
% % clear basisCollectionStruct;
if (flagSB == 1)
T1 = RANS_MatrixStruct.T_ijk;
L1 = RANS_MatrixStruct.L_ij;
T2 = RANS_MatrixStruct.T_ijk2;
L2 = RANS_MatrixStruct.L_ij2;    
else
T = RANS_MatrixStruct.T_ijk;
L = RANS_MatrixStruct.L_ij;
end
N = RANS_MatrixStruct.N;
clear RANS_MatrixStruct basisCollectionStruct

% keyboard;

%  Function = 1/2[(1-alpha)*{||Aij*xj - bi||^2} + alpha*{||Tijk*xk*xj + Lij*xj||^2}]
%                            Boundary matching                   RANS
% Boundary matching : boundaryMatchingValue = ||Aij*xj - bi||^2
boundaryMatchingValue = tileBoundaryMatching_openFoam(A, aa, rhs, bvalue, bvalue2);

% disp(num2str(aa))
% RANS : reducedRANSValue = ||(Tijk*xk*xj + Lij*xj)|| 
if (flagSB == 1)
    reducedRANSValue = RANS_Reduction_OpenFoam(T1, L1, N, aa, flagSB,T2, L2);    
    M = size(T1,1);
else
    reducedRANSValue = RANS_Reduction_OpenFoam(T, L, N, aa, flagSB);
    M = size(T,1);
end
% boundaryMatchingValue = tileBoundaryMatching(basisCollectionFileName,  CFD_ResultsPath, rhs, aa);


x = zeros(M,1);
dF_dm = zeros(length(aa),1);

% Gradient 
% Expression: (1-alpha)*(Aij*xj-bi)Aim + 
%             alpha*(Tijk*xk*xj + Lij*xj)(Tipm*xp + Timp*xp + Lim)
for a = 1:3
    begin_x = 1 + (a-1)*M;
    end_x = begin_x + M -1;
    x=aa(begin_x:end_x);
    if (flagSB == 1)
        if (a == 2)
            T = T2;
            L = L2;
        else
            T = T1;
            L = L1;    
        end
    end
     
%     b=bvalue(begin_x:end_x);
    
    for m=1:M
        f1 = zeros(M,1);
        f2 = zeros(M,1);
        
        for i=1:M

            f2_1 = 0;
            f2_2 = 0;
            
            for j=1:M
                f10=0;

                for k=1:M
                    f10=f10+T(i,j,k)*x(k);         % Tijk*xk
                end

                f1(i) = f1(i) + (f10+L(i,j))*x(j); % ((Tijk*xk)+ Lij)xj
                
                f2_1 = f2_1 + T(i,j,m)*x(j);       %  Tipm*xp = Tijm*xj
                f2_2 = f2_2 + T(i,m,j)*x(j);       %  Timp*xp = Timj*xj

            end
            f2(i) = f2_1 + f2_2 + L(i,m);          % (Tipm*xp + Timp*xp + Lim)
        dF_dm((a-1)*M+m) = (dF_dm((a-1)*M+m) + f1(i)*f2(i)); % RANS Part        
        end
    end
    
    
end
% keyboard;
dF_dm = alpha*dF_dm;
dF_dm = dF_dm + ((1-alpha)*(A*aa-rhs));

% dF_dm_check = gradient(T,L,aa,M,noTile);

%     dF_dm = dF_dm';
        
% Hessian
% Expression: (1-alpha)*Ain*Aim + 
%             alpha*(Tijk*xk*xj + Lij*xj)(Tinm + Timn) +
%             alpha*(Tijn*xj + Tink*xk + Lin)(Tipm*xp + Timp*xp + Lim)

dF_dm_dn = zeros(length(aa),length(aa));

for a = 1:3
    begin_x = 1 + (a-1)*M;
    end_x = begin_x + M -1;
    x=aa(begin_x:end_x);
    if (flagSB == 1)
        if (a == 2)
            T = T2;
            L = L2;
        else
            T = T1;
            L = L1;    
        end
    end
    
    for n = 1:M
        for m=1:M
            f11 = zeros(M,1);
            f22 = zeros(M,1);
            
            for i=1:M

                f22_1 = 0;
                f22_2 = 0;
                f22_3 = 0;
                f22_4 = 0;
                for j=1:M
                    f10=0;

                    for k=1:M
                        f10=f10+T(i,j,k)*x(k);
                    end

                    f11(i) = f11(i) + (f10+L(i,j))*x(j);
                    f22_1 = f22_1 + T(i,j,m)*x(j);
                    f22_2 = f22_2 + T(i,m,j)*x(j);
                    f22_3 = f22_3 + T(i,n,j)*x(j);
                    f22_4 = f22_4 + T(i,j,n)*x(j);

                end
                f22(i) = (f22_1 + f22_2 + L(i,m))*(f22_3 + f22_4 + L(i,n));
                dF_dm_dn((a-1)*M+m,(a-1)*M+n) = dF_dm_dn((a-1)*M+m,(a-1)*M+n) + (f11(i)*(T(i,m,n)+T(i,n,m))) + f22(i);
            end
        end
    end
    
    

end
dF_dm_dn = alpha*dF_dm_dn;
dF_dm_dn = dF_dm_dn + (1-alpha)*A;
 
ObjFuncVal = 0.5*((1 - alpha)*(boundaryMatchingValue) + alpha*(reducedRANSValue^2));
end