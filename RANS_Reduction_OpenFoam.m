function Reduced_RANS = RANS_Reduction_OpenFoam(T_ijk1, L_ij1, N, aa_0, flagSB, T_ijk2, L_ij2)

if (nargin == 5)
    T_ijk2 = T_ijk1;
    L_ij2 = L_ij1;
elseif (nargin ~=7 )
fprintf(1, '\n Function RANS_Reduction_OpenFoam: Unknown usage')
return;
end
T_ijk = T_ijk1;
L_ij = L_ij1;
M = size(T_ijk,1);
% Reduced_RANS_func = zeros(size(M,M));
Reduced_RANS_Vector = zeros(size(3*M,1));

% size(aa_0)

for m = 1:3
    if (flagSB == 1)
        if (m == 2)
            T_ijk = T_ijk2;
            L_ij = L_ij2;
        else
            T_ijk = T_ijk1;
            L_ij = L_ij1;    
        end
    end
    
    for i=1:M
        Reduced_RANS_Vector(i+(m-1)*M) = 0;
        for j=1:M
            p=0;
            for k=1:M
                p=p+T_ijk(i,j,k)*aa_0(k+(m-1)*M,1);
            end
            Reduced_RANS_Vector(i+(m-1)*M) = Reduced_RANS_Vector(i+(m-1)*M)...
                                            + ((p + L_ij(i,j))*aa_0(j+(m-1)*M,1));
        end
    end
end

% keyboard;
% Reduced_RANS_Vector = Reduced_RANS_Vector + (Reduced_RANS_func*aa_0(1+M:2*M,1));% ...
    %+ (L_ij*aa_0(1+M:2*M,1));
% Reduced_RANS_Vector=Reduced_RANS_Vector;
% keyboard;
Reduced_RANS = norm(Reduced_RANS_Vector,'fro');
disp(['Frobenius Norm of Reduced_RANS:', num2str(Reduced_RANS)])                  
end
% keyboard;
