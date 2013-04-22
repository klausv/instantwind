function boundaryMatchingValue = tileBoundaryMatching_openFoam(A, aa, rhs, bvalue, bvalue2)

    boundaryMatchingVector = (aa'*A*aa) - (2*aa'*rhs)+ ((bvalue'*bvalue)+(bvalue2'*bvalue2));
    boundaryMatchingValue = boundaryMatchingVector;
    disp(['Frobenius Norm of Boundary matching term', num2str(boundaryMatchingValue)])