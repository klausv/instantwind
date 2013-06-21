// V_mex
//
//  Generation of sparse Viscosity matrix for the RANS equations.
//  Refer to RANSequations.m for the spesification of the matrix.
//  Compile with Matlab command ">>mex V_mex.c".
//
//  Note: The transposeed matrix is generated.
//
//  Use: V=V_mex(Phi,nx,ny,nz,hx,hy,hz);
//

#include "mex.h"

void 
mexFunction(
    int nlhs,       mxArray     *plhs[],        
    int nrhs, const mxArray     *prhs[]){
    
    // Input arguments
    //
    double  *Phi;       // Basis vector for fields ux,uy,uz,mu
    int     nx,ny,nz;   // # grid-points
    double  hx,hy,hz;   // grid spacing
    
    // Output sparse matrix.
    // Note that the data arrays are sorted so that the first entries are
    // the none-zero entries of column 1, then column 2, etc.
    //
    double  *A;     // Array of none-zeros elements of sparse matrix
    int     *Air,   // (row # of none-zero elements). Row #'s start at 0.
            *Ajc;   // Ajc[j]: index into A and Air of first none-zero element
                    // of column j (0,1,...).
    int     nnz;    // # none-zero elemnents in A
    
    int     idummy;
    int     i,j,k;
    int     index_data,index_col;
    int     offset_ux,offset_uy,offset_uz,offset_mu;
    int     index_xminus,index_xplus,index_yminus,index_yplus,index_zminus,index_zplus;
    int     index_xminus_yminus,index_xplus_yplus,index_xminus_zminus,index_xplus_zplus,index_yminus_zminus,index_yplus_zplus;
    int     ibegin, iend, jbegin, jend, kbegin, kend;
    
    
    // Check inputs
    //
//     if (nrhs == 13){
//         printf("location of turbines as input arguments provided.");};        
    if (nrhs != 13 && nrhs != 7){
        mexErrMsgTxt("7 input argument required; Use: V=V_mex(Phi,nx,ny,nz,hx,hy,hz)");};
    if (nlhs != 1) {
        mexErrMsgTxt("No output argument spesified; Use : V=V_mex(...)");};
    
    // Assign inputs
    //
    nx=(int)*mxGetPr(prhs[1]);ny=(int)*mxGetPr(prhs[2]);nz=(int)*mxGetPr(prhs[3]);
    hx=(double)*mxGetPr(prhs[4]);hy=(double)*mxGetPr(prhs[5]);hz=(double)*mxGetPr(prhs[6]);
    
    if (nrhs == 13){
        ibegin = (int)*mxGetPr(prhs[7]);iend = (int)*mxGetPr(prhs[8]);
        jbegin = (int)*mxGetPr(prhs[9]);jend = (int)*mxGetPr(prhs[10]);
        kbegin = (int)*mxGetPr(prhs[11]);kend = (int)*mxGetPr(prhs[12]);
    }
    else{
        ibegin=0;iend=0; jbegin=0; jend=0; kbegin=0; kend=0;
    };
    
    idummy=mxGetM(prhs[0]);
    if(idummy!=(4*nx*ny*nz))
        mexErrMsgTxt("Input argument 1 (Phi) must be a column vector with length 4*nx*ny*nz (Phi=[ux,uy,uz,mu])");
    Phi=mxGetPr(prhs[0]);
        
    // Create output sparse matrix A
    //
if (nrhs == 13){    
    nnz=(3*15*(nz-2)*(ny-2)*(nx-2))-(3*15*(kend-kbegin+1)*(jend-jbegin+1)*(iend-ibegin+1)); // + 3*((2*nx*ny) + (2*(nz-2)*ny)+(2*(nz-2)*(nx-2)));
}
else{
    nnz=(3*15*(nz-2)*(ny-2)*(nx-2));
};
    
// //     nnz=(3*15*(nz-2)*(ny-2)*(nx-2));//+(3*(2*(nx*ny) + 2*((nz-2)*ny)  + 2*((nz-2)*(nx-2))));
    plhs[0] = mxCreateSparse(3*nz*ny*nx,3*nz*ny*nx,nnz,mxREAL); // (#rows,#columns,nnz,mxREAL)
    A  = mxGetPr(plhs[0]);
    Air = mxGetIr(plhs[0]);
    Ajc = mxGetJc(plhs[0]);
    
    // Build the (transposed) viscosity V matrix
    //
    index_data=0;
    offset_ux=0;
    offset_uy=(nx*ny*nz);
    offset_uz=2*(nx*ny*nz);
    offset_mu=3*(nx*ny*nz);
    for(k=0;k<nz;k++)for(j=0;j<ny;j++)for(i=0;i<nx;i++){
        index_col=k*(nx*ny)+j*nx+i;
        if((k==0)||(k==nz-1)||(j==0)||(j==ny-1)||(i==0)||(i==nx-1))
        {
            Ajc[index_col]=index_data;
//             Air[index_data]=index_col;       A[index_data]=1.0;
//             index_data++;        
        }
        else if((k>=kbegin)&&(k<=kend)&&(j>=jbegin)&&(j<=jend)&&(i>=ibegin)&&(i<=iend)){
            Ajc[index_col]=index_data;
        }        
        else{
            index_xminus=   (k)*(nx*ny)+(j)*nx+i-1;
            index_xplus=    (k)*(nx*ny)+(j)*nx+i+1;
            index_yminus=   (k)*(nx*ny)+(j-1)*nx+i;
            index_yplus=    (k)*(nx*ny)+(j+1)*nx+i;
            index_zminus=   (k-1)*(nx*ny)+(j)*nx+i;
            index_zplus=    (k+1)*(nx*ny)+(j)*nx+i;
            index_xminus_yminus=    (k)*(nx*ny)+(j-1)*nx+i-1;
            index_xplus_yplus=      (k)*(nx*ny)+(j+1)*nx+i+1;
            index_xminus_zminus=    (k-1)*(nx*ny)+(j)*nx+i-1;
            index_xplus_zplus=      (k+1)*(nx*ny)+(j)*nx+i+1;
            index_yminus_zminus=    (k-1)*(nx*ny)+(j-1)*nx+i;
            index_yplus_zplus=      (k+1)*(nx*ny)+(j+1)*nx+i;
            
            // a=x
            //
            Ajc[index_col]=index_data;
            Air[index_data]=index_zminus;       A[index_data]=(1/(2*hz*hz))*Phi[offset_mu+index_zminus];
            index_data++;
            Air[index_data]=index_yminus;       A[index_data]=(1/(2*hy*hy))*Phi[offset_mu+index_yminus];
            index_data++;
            Air[index_data]=index_xminus;       A[index_data]=(1/(2*hx*hx))*Phi[offset_mu+index_xminus] + (1/(2*hx*hx))*   Phi[offset_mu+index_xminus];
            index_data++;
            Air[index_data]=index_col;          A[index_data]=-(1/(2*hx*hx))*(Phi[offset_mu+index_xminus]+Phi[offset_mu+index_xplus])
                                                            -(1/(2*hy*hy))*(Phi[offset_mu+index_yminus]+Phi[offset_mu+index_yplus])
                                                            -(1/(2*hz*hz))*(Phi[offset_mu+index_zminus]+Phi[offset_mu+index_zplus])
                                                            -(1/(2*hx*hx))*(Phi[offset_mu+index_xminus]+Phi[offset_mu+index_xplus]);
            index_data++;
            Air[index_data]=index_xplus;        A[index_data]=(1/(2*hx*hx))*Phi[offset_mu+index_xplus] + (1/(2*hx*hx))*Phi[offset_mu+index_xplus];
            index_data++;
            Air[index_data]=index_yplus;        A[index_data]=(1/(2*hy*hy))*Phi[offset_mu+index_yplus];
            index_data++;
            Air[index_data]=index_zplus;        A[index_data]=(1/(2*hz*hz))*Phi[offset_mu+index_zplus];
            index_data++;
            Air[index_data]=offset_uy+index_xminus_yminus;  A[index_data]=(1/(2*hy*hx))*Phi[offset_mu+index_yminus];
            index_data++;
            Air[index_data]=offset_uy+index_yminus;         A[index_data]=-(1/(2*hy*hx))*Phi[offset_mu+index_yminus];
            index_data++;
            Air[index_data]=offset_uy+index_yplus;          A[index_data]=-(1/(2*hy*hx))*Phi[offset_mu+index_yplus];
            index_data++;
            Air[index_data]=offset_uy+index_xplus_yplus;    A[index_data]=(1/(2*hy*hx))*Phi[offset_mu+index_yplus];
            index_data++;
            Air[index_data]=offset_uz+index_xminus_zminus;  A[index_data]=(1/(2*hz*hx))*Phi[offset_mu+index_zminus];
            index_data++;
            Air[index_data]=offset_uz+index_zminus;         A[index_data]=-(1/(2*hz*hx))*Phi[offset_mu+index_zminus];
            index_data++;
            Air[index_data]=offset_uz+index_zplus;          A[index_data]=-(1/(2*hz*hx))*Phi[offset_mu+index_zplus];
            index_data++;
            Air[index_data]=offset_uz+index_xplus_zplus;    A[index_data]=(1/(2*hz*hx))*Phi[offset_mu+index_zplus];
            index_data++;
        };
    };
    for(k=0;k<nz;k++)for(j=0;j<ny;j++)for(i=0;i<nx;i++){
        index_col=offset_uy+k*(nx*ny)+j*nx+i;
        if((k==0)||(k==nz-1)||(j==0)||(j==ny-1)||(i==0)||(i==nx-1))
        {
        //    Ajc[index_col]=index_data;
            Ajc[index_col]=index_data;
//             Air[index_data]=index_col;       A[index_data]=1.0;
//             index_data++;            
        }
        else if((k>=kbegin)&&(k<=kend)&&(j>=jbegin)&&(j<=jend)&&(i>=ibegin)&&(i<=iend)){
            Ajc[index_col]=index_data;
        }        
        else{
            index_xminus=   (k)*(nx*ny)+(j)*nx+i-1;
            index_xplus=    (k)*(nx*ny)+(j)*nx+i+1;
            index_yminus=   (k)*(nx*ny)+(j-1)*nx+i;
            index_yplus=    (k)*(nx*ny)+(j+1)*nx+i;
            index_zminus=   (k-1)*(nx*ny)+(j)*nx+i;
            index_zplus=    (k+1)*(nx*ny)+(j)*nx+i;
            index_xminus_yminus=    (k)*(nx*ny)+(j-1)*nx+i-1;
            index_xplus_yplus=      (k)*(nx*ny)+(j+1)*nx+i+1;
            index_xminus_zminus=    (k-1)*(nx*ny)+(j)*nx+i-1;
            index_xplus_zplus=      (k+1)*(nx*ny)+(j)*nx+i+1;
            index_yminus_zminus=    (k-1)*(nx*ny)+(j-1)*nx+i;
            index_yplus_zplus=      (k+1)*(nx*ny)+(j+1)*nx+i;
            
            // a=y
            //
            Ajc[index_col]=index_data;
            Air[index_data]=index_xminus_yminus;            A[index_data]=(1/(2*hx*hy))*Phi[offset_mu+index_xminus];
            index_data++;
            Air[index_data]=index_xminus;                   A[index_data]=-(1/(2*hx*hy))*Phi[offset_mu+index_xminus];
            index_data++;
            Air[index_data]=index_xplus;                    A[index_data]=-(1/(2*hx*hy))*Phi[offset_mu+index_xplus];
            index_data++;
            Air[index_data]=index_xplus_yplus;              A[index_data]=(1/(2*hx*hy))*Phi[offset_mu+index_xplus];
            index_data++;
            Air[index_data]=offset_uy+index_zminus;         A[index_data]=(1/(2*hz*hz))*Phi[offset_mu+index_zminus];
            index_data++;
            Air[index_data]=offset_uy+index_yminus;         A[index_data]=(1/(2*hy*hy))*Phi[offset_mu+index_yminus]+(1/(2*hy*hy))*Phi[offset_mu+index_yminus];
            index_data++;
            Air[index_data]=offset_uy+index_xminus;         A[index_data]=(1/(2*hx*hx))*Phi[offset_mu+index_xminus];
            index_data++;
            Air[index_data]=index_col;                      A[index_data]=-(1/(2*hx*hx))*(Phi[offset_mu+index_xminus]+Phi[offset_mu+index_xplus])
                                                                 -(1/(2*hy*hy))*(Phi[offset_mu+index_yminus]+Phi[offset_mu+index_yplus])
                                                                 -(1/(2*hz*hz))*(Phi[offset_mu+index_zminus]+Phi[offset_mu+index_zplus])
                                                                 -(1/(2*hy*hy))*(Phi[offset_mu+index_yminus]+Phi[offset_mu+index_yplus]);
            index_data++;
            Air[index_data]=offset_uy+index_xplus;          A[index_data]=(1/(2*hx*hx))*Phi[offset_mu+index_xplus];
            index_data++;
            Air[index_data]=offset_uy+index_yplus;          A[index_data]=(1/(2*hy*hy))*Phi[offset_mu+index_yplus]+ (1/(2*hy*hy))*Phi[offset_mu+index_yplus];
            index_data++;
            Air[index_data]=offset_uy+index_zplus;          A[index_data]=(1/(2*hz*hz))*Phi[offset_mu+index_zplus];
            index_data++;
            Air[index_data]=offset_uz+index_yminus_zminus;  A[index_data]=(1/(2*hz*hy))*Phi[offset_mu+index_zminus];
            index_data++;
            Air[index_data]=offset_uz+index_zminus;         A[index_data]=-(1/(2*hz*hy))*Phi[offset_mu+index_zminus];
            index_data++;
            Air[index_data]=offset_uz+index_zplus;         	A[index_data]=-(1/(2*hz*hy))*Phi[offset_mu+index_zplus];
            index_data++;
            Air[index_data]=offset_uz+index_yplus_zplus;    A[index_data]=(1/(2*hz*hy))*Phi[offset_mu+index_zplus];
            index_data++;      
        };
    };
    for(k=0;k<nz;k++)for(j=0;j<ny;j++)for(i=0;i<nx;i++){
        index_col=offset_uz+k*(nx*ny)+j*nx+i;
        if((k==0)||(k==nz-1)||(j==0)||(j==ny-1)||(i==0)||(i==nx-1))
        {
        //    Ajc[index_col]=index_data;
            Ajc[index_col]=index_data;
//             Air[index_data]=index_col;       A[index_data]=1.0;
//             index_data++;   
        }
        else if((k>=kbegin)&&(k<=kend)&&(j>=jbegin)&&(j<=jend)&&(i>=ibegin)&&(i<=iend)){
            Ajc[index_col]=index_data;
        }        
        else{
            index_xminus=   (k)*(nx*ny)+(j)*nx+i-1;
            index_xplus=    (k)*(nx*ny)+(j)*nx+i+1;
            index_yminus=   (k)*(nx*ny)+(j-1)*nx+i;
            index_yplus=    (k)*(nx*ny)+(j+1)*nx+i;
            index_zminus=   (k-1)*(nx*ny)+(j)*nx+i;
            index_zplus=    (k+1)*(nx*ny)+(j)*nx+i;
            index_xminus_yminus=    (k)*(nx*ny)+(j-1)*nx+i-1;
            index_xplus_yplus=      (k)*(nx*ny)+(j+1)*nx+i+1;
            index_xminus_zminus=    (k-1)*(nx*ny)+(j)*nx+i-1;
            index_xplus_zplus=      (k+1)*(nx*ny)+(j)*nx+i+1;
            index_yminus_zminus=    (k-1)*(nx*ny)+(j-1)*nx+i;
            index_yplus_zplus=      (k+1)*(nx*ny)+(j+1)*nx+i;
            
            // a=z
            //
            Ajc[index_col]=index_data;
            Air[index_data]=index_xminus_zminus;            A[index_data]=(1/(2*hx*hz))*Phi[offset_mu+index_xminus];
            index_data++;
            Air[index_data]=index_xminus;                   A[index_data]=-(1/(2*hx*hz))*Phi[offset_mu+index_xminus];
            index_data++;
            Air[index_data]=index_xplus;                    A[index_data]=-(1/(2*hx*hz))*Phi[offset_mu+index_xplus];
            index_data++;
            Air[index_data]=index_xplus_zplus;              A[index_data]=(1/(2*hx*hz))*Phi[offset_mu+index_xplus];
            index_data++;
            Air[index_data]=offset_uy+index_yminus_zminus;  A[index_data]=(1/(2*hy*hz))*Phi[offset_mu+index_yminus];
            index_data++;
            Air[index_data]=offset_uy+index_yminus;         A[index_data]=-(1/(2*hy*hz))*Phi[offset_mu+index_yminus];
            index_data++;
            Air[index_data]=offset_uy+index_yplus;          A[index_data]=-(1/(2*hy*hz))*Phi[offset_mu+index_yplus];
            index_data++;
            Air[index_data]=offset_uy+index_yplus_zplus;    A[index_data]=(1/(2*hy*hz))*Phi[offset_mu+index_yplus];
            index_data++;
            Air[index_data]=offset_uz+index_zminus;         A[index_data]=(1/(2*hz*hz))*Phi[offset_mu+index_zminus]+(1/(2*hz*hz))*Phi[offset_mu+index_zminus];
            index_data++;
            Air[index_data]=offset_uz+index_yminus;         A[index_data]=(1/(2*hy*hy))*Phi[offset_mu+index_yminus];
            index_data++;
            Air[index_data]=offset_uz+index_xminus;         A[index_data]=(1/(2*hx*hx))*Phi[offset_mu+index_xminus];
            index_data++;
            Air[index_data]=index_col;                      A[index_data]=-(1/(2*hx*hx))*(Phi[offset_mu+index_xminus]+Phi[offset_mu+index_xplus])
                                                                    -(1/(2*hy*hy))*(Phi[offset_mu+index_yminus]+Phi[offset_mu+index_yplus])
                                                                    -(1/(2*hz*hz))*(Phi[offset_mu+index_zminus]+Phi[offset_mu+index_zplus])
                                                                    -(1/(2*hz*hz))*(Phi[offset_mu+index_zminus]+Phi[offset_mu+index_zplus]);
            index_data++;
            Air[index_data]=offset_uz+index_xplus;          A[index_data]=(1/(2*hx*hx))*Phi[offset_mu+index_xplus];
            index_data++;
            Air[index_data]=offset_uz+index_yplus;          A[index_data]=(1/(2*hy*hy))*Phi[offset_mu+index_yplus];
            index_data++;
            Air[index_data]=offset_uz+index_zplus;          A[index_data]=(1/(2*hz*hz))*Phi[offset_mu+index_zplus]+(1/(2*hz*hz))*Phi[offset_mu+index_zplus];
            index_data++;
            
        };
    };
    
    //mexPrintf("\n>>> index_col=%d, index_data=%d, 3*nznynx=%d, nnz=%d\n",index_col,index_data,3*nz*ny*nx,nnz);
    
    if(index_data!=(nnz))
        mexErrMsgTxt("MISMATCH between nnz and index_data");
    
    Ajc[3*nz*ny*nx]=nnz;
    
    //for(i=0;i<nnz;i++) printf("ir(%2d): %d\n",i,Air[i]);
    
    }
