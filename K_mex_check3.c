// P_mex
//
//  Geneation of sparse Advection matrix for the RANS equations.
//  Refer to RANSequations.m for the spesification of the matrix.
//  Compile with Matlab command ">>mex P_mex_P_VB.c".
//
//  Note: The transposeed matrix is generated.
//
//  Use: P=P_mex(Phi,rho,nx,ny,nz,hx,hy,hz);
//

#include "mex.h"

void 
mexFunction(
    int nlhs,       mxArray     *plhs[],        
    int nrhs, const mxArray     *prhs[])
{
    
    // Input arguments
    //
    double  *Phi;       // Basis vector for fields ux,uy,uz.
    double  rho;
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
    int     offset_k,offset_ky,offset_kz;
    int     index_xminus,index_xplus,index_yminus,index_yplus,index_zminus,index_zplus;
    int     checksize;
    int     ibegin, iend, jbegin, jend, kbegin, kend;
    
    
    // Check inputs
    //
//     if (nrhs == 14){
//         printf("location of turbines as input arguments provided. \n");};                
    if (nrhs != 14 && nrhs != 8){
        mexErrMsgTxt("8 input argument required; Use: A=A_mex(Phi,rho,nx,ny,nz,hx,hy,hz)");};
    if (nlhs != 1) {
        mexErrMsgTxt("No output argument spesified; Use : A=A_mex(...)");};
    
    // Assign inputs
    //
    rho=(double)*mxGetPr(prhs[1]);
    nx=(int)*mxGetPr(prhs[2]);ny=(int)*mxGetPr(prhs[3]);nz=(int)*mxGetPr(prhs[4]);
    hx=(double)*mxGetPr(prhs[5]);hy=(double)*mxGetPr(prhs[6]);hz=(double)*mxGetPr(prhs[7]);

    if (nrhs == 14){
        ibegin = (int)*mxGetPr(prhs[8]);iend = (int)*mxGetPr(prhs[9]);
        jbegin = (int)*mxGetPr(prhs[10]);jend = (int)*mxGetPr(prhs[11]);
        kbegin = (int)*mxGetPr(prhs[12]);kend = (int)*mxGetPr(prhs[13]);
    }
    else{
        ibegin=0; iend=0; jbegin=0; jend=0; kbegin=0; kend=0;
    }    
    
    idummy=mxGetM(prhs[0]);
    if(idummy!=(nx*ny*nz))
        mexErrMsgTxt("Input argument 1 (Phi) must be a column vector with length nx*ny*nz");
    Phi=mxGetPr(prhs[0]);
        
    // Create output sparse matrix A
    //
if (nrhs == 14){    
    nnz=(3*2*(nz-2)*(ny-2)*(nx-2))-(3*2*(kend-kbegin+1)*(jend-jbegin+1)*(iend-ibegin+1)); // + 3*((2*nx*ny) + (2*(nz-2)*ny)+(2*(nz-2)*(nx-2)));
}
else{
    nnz=(3*2*(nz-2)*(ny-2)*(nx-2));
};
    
// //     nnz=(3*2*(nz-2)*(ny-2)*(nx-2));// + 3*((2*nx*ny) + (2*(nz-2)*ny)+(2*(nz-2)*(nx-2)));

//     printf("nnz is %d",nnz);
    plhs[0] = mxCreateSparse(3*nz*ny*nx,3*nz*ny*nx,nnz,mxREAL); // (#rows,#columns,nnz,mxREAL)
    A  = mxGetPr(plhs[0]);
    Air = mxGetIr(plhs[0]);
    Ajc = mxGetJc(plhs[0]);
    
    checksize = sizeof(A)/sizeof(double *);
//     printf("size of A is %d",checksize);
    checksize = sizeof(Air)/sizeof(int *);
//     printf("size of Air is %d",checksize);
    checksize = sizeof(Ajc)/sizeof(int *);
//     printf("size of Ajc is %d",checksize);
    // Build the (transposed) advection A matrix
    //
    index_data=0;
    offset_k=0;
//     printf("Check 1");
    for(k=0;k<nz;k++)for(j=0;j<ny;j++)for(i=0;i<nx;i++){
        index_col=k*(nx*ny)+j*nx+i;
        if((k==0)||(k==nz-1)||(j==0)||(j==ny-1)||(i==0)||(i==nx-1)){
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
//             index_yminus=   (k)*(nx*ny)+(j-1)*nx+i;
//             index_yplus=    (k)*(nx*ny)+(j+1)*nx+i;
//             index_zminus=   (k-1)*(nx*ny)+(j)*nx+i;
//             index_zplus=    (k+1)*(nx*ny)+(j)*nx+i;
            
            Ajc[index_col]=index_data;
//             Air[index_data]=index_zminus;   A[index_data]=-(rho/(2*hz))*Phi[offset_k+index_zminus];    index_data++;
//             Air[index_data]=index_yminus;   A[index_data]=-(rho/(2*hy))*Phi[offset_k+index_yminus];    index_data++;
            Air[index_data]=index_xminus;   A[index_data]=-(rho/(2*hx));    index_data++;
            Air[index_data]=index_xplus;    A[index_data]= (rho/(2*hx));     index_data++;            
//             Air[index_data]=index_xminus;   A[index_data]=-(rho/(2*hx))*Phi[offset_k+index_xminus];    index_data++;
//             Air[index_data]=index_xplus;    A[index_data]= (rho/(2*hx))*Phi[offset_k+index_xplus];     index_data++;
//             Air[index_data]=index_yplus;    A[index_data]= (rho/(2*hy))*Phi[offset_k+index_yplus];     index_data++;
//             Air[index_data]=index_zplus;    A[index_data]= (rho/(2*hz))*Phi[offset_k+index_zplus];     index_data++;
            
        };
    };
    offset_ky=(nx*ny*nz);
//     printf("Check 2");    
    for(k=0;k<nz;k++)for(j=0;j<ny;j++)for(i=0;i<nx;i++){
        index_col=offset_ky+k*(nx*ny)+j*nx+i;
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
//             index_xminus=   (k)*(nx*ny)+(j)*nx+i-1;
//             index_xplus=    (k)*(nx*ny)+(j)*nx+i+1;
            index_yminus=   (k)*(nx*ny)+(j-1)*nx+i;
            index_yplus=    (k)*(nx*ny)+(j+1)*nx+i;
//             index_zminus=   (k-1)*(nx*ny)+(j)*nx+i;
//             index_zplus=    (k+1)*(nx*ny)+(j)*nx+i;
            
            Ajc[index_col]=index_data;
//             Air[index_data]=index_zminus;   A[index_data]=-(rho/(2*hz))*Phi[offset_k+index_zminus];    index_data++;
//             Air[index_data]=index_yminus;   A[index_data]=-(rho/(2*hy))*Phi[offset_k+index_yminus];    index_data++;
            Air[index_data]=index_yminus;   A[index_data]=-(rho/(2*hy));    index_data++;            
//             Air[index_data]=index_xminus;   A[index_data]=-(rho/(2*hx))*Phi[offset_k+index_xminus];    index_data++;
//             Air[index_data]=index_xplus;    A[index_data]= (rho/(2*hx))*Phi[offset_k+index_xplus];     index_data++;
//             Air[index_data]=index_yplus;    A[index_data]= (rho/(2*hy))*Phi[offset_k+index_yplus];     index_data++;
            Air[index_data]=index_yplus;    A[index_data]= (rho/(2*hy));     index_data++;            
//             Air[index_data]=index_zplus;    A[index_data]= (rho/(2*hz))*Phi[offset_k+index_zplus];     index_data++;
            
        };
    };
    offset_kz=2*(nx*ny*nz);
//     printf("Check 3");    
    for(k=0;k<nz;k++)for(j=0;j<ny;j++)for(i=0;i<nx;i++){
        index_col=offset_kz+k*(nx*ny)+j*nx+i;
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
//             index_xminus=   (k)*(nx*ny)+(j)*nx+i-1;
//             index_xplus=    (k)*(nx*ny)+(j)*nx+i+1;
//             index_yminus=   (k)*(nx*ny)+(j-1)*nx+i;
//             index_yplus=    (k)*(nx*ny)+(j+1)*nx+i;
            index_zminus=   (k-1)*(nx*ny)+(j)*nx+i;
            index_zplus=    (k+1)*(nx*ny)+(j)*nx+i;
            
            Ajc[index_col]=index_data;
//             Air[index_data]=index_zminus;   A[index_data]=-(rho/(2*hz))*Phi[offset_k+index_zminus];    index_data++;
            Air[index_data]=index_zminus;   A[index_data]=-(rho/(2*hz));    index_data++;            
//             Air[index_data]=index_yminus;   A[index_data]=-(rho/(2*hy))*Phi[offset_k+index_yminus];    index_data++;
//             Air[index_data]=index_xminus;   A[index_data]=-(rho/(2*hx))*Phi[offset_k+index_xminus];    index_data++;
//             Air[index_data]=index_xplus;    A[index_data]= (rho/(2*hx))*Phi[offset_k+index_xplus];     index_data++;
//             Air[index_data]=index_yplus;    A[index_data]= (rho/(2*hy))*Phi[offset_k+index_yplus];     index_data++;
//             Air[index_data]=index_zplus;    A[index_data]= (rho/(2*hz))*Phi[offset_k+index_zplus];     index_data++;
            Air[index_data]=index_zplus;    A[index_data]= (rho/(2*hz));     index_data++;            
            
        };
    };
//     printf("Check 4");   
    printf("index_data, nnz %d %d", index_data, nnz);
    if(index_data!=(nnz))
        mexErrMsgTxt("MISMATCH between nnz and index_data");
    
    Ajc[3*nz*ny*nx]=nnz;
//     printf("Check 5");    
    
    //mexPrintf("\n>>> index_col=%d, index_data=%d, 3*nznynx=%d, nnz=%d\n",index_col,index_data,3*nz*ny*nx,nnz);
    
    
    }