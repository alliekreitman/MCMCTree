
typedef struct matrixsolution{
    gsl_matrix *m;
    double t;
}matrixsolution;

typedef struct matArray{
    matrixsolution *transMat;
}matArray;

matArray mats[MUNUM];

matrixsolution *findsolution(int nmax,int num,int index)
{
    matrixsolution *transMat=malloc(sizeof(matrixsolution)*nmax);
    gsl_matrix *onestep=makematrix(num,omega[index]);
    transMat[0].m=onestep;
    transMat[0].t=0.0;
    int i;
    for(i=1;i<nmax;i++){
        transMat[i].m=gsl_matrix_alloc(num,num);
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                        1.0, transMat[i-1].m, onestep,
                        0.0, transMat[i].m);
        transMat[i].t=i*odeh;
    }
    return transMat;
}
void freesolution(){
    int i,j;
    for(i=0;i<MUNUM;i++){
        for(j=0;j<maxn;j++){
            gsl_matrix_free(mats[i].transMat[j].m);
        }
        free(mats[i].transMat);
    }
}
gsl_matrix *ProbMat(double t,int index){
    
    int k=t/odeh;
    //printf("%lf %d\n",t,k);
    return mats[index].transMat[k].m;
}
