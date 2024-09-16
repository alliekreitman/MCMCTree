#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>


typedef struct Nodename{
    char name;
    int str[3];
}Nodename;
Nodename mnode[61];

typedef struct TheMatrix{
    gsl_matrix *onestep;
    gsl_matrix *evec;
    gsl_matrix *diag;
    gsl_matrix *inv;
}TheMatrix;

TheMatrix Matrixes[MUNUM];

int diffstr(int *str1,int *str2){
    int i;
    int ct=0;
    for(i=0;i<3;i++)
        if (str1[i]!=str2[i])
            ct++;
    return ct;
}
// 10 11 14
int findindex(char *str){
    
    int n;
    int trans[3];
    for (n=0;n<3;n++){
        if (str[n]=='T') trans[n]=0;
        else if (str[n]=='C') trans[n]=1;
        else if (str[n]=='A') trans[n]=2;
        else if (str[n]=='G') trans[n]=3;
    }
    n=trans[0]*16+trans[1]*4+trans[2];
    if (n<10) return n;
    if (n==10 || n==11 || n==14) return -1;
    if (n==12) return 10;
    if (n==13) return 11;
    
    return n-3;
}

char acidname(int *str){
    
    char *tmp;
    int i;
    tmp=calloc(3,sizeof(char));
    for (i=0;i<3;i++){
        if (str[i]==0) strcat(tmp,"T");
        else if (str[i]==1) strcat(tmp,"C");
        else if (str[i]==2) strcat(tmp,"A");
        else if (str[i]==3) strcat(tmp,"G");
    }
    //printf("%s\n",tmp);
    if (strcmp(tmp,"GCT")==0||strcmp(tmp,"GCC")==0 || strcmp(tmp,"GCA")==0 || strcmp(tmp,"GCG")==0)
    {free(tmp);return 'A';}
    else if (strcmp(tmp,"CGT")==0||strcmp(tmp,"CGC")==0||strcmp(tmp,"CGA")==0||strcmp(tmp,"CGG")==0||strcmp(tmp,"AGG")==0||strcmp(tmp,"AGA")==0)
    {free(tmp);return 'R';}
    else if (strcmp(tmp,"AAT")==0||strcmp(tmp,"AAC")==0)
    {free(tmp);return 'N';}
    else if (strcmp(tmp,"GAT")==0||strcmp(tmp,"GAC")==0)
    {free(tmp);return 'D';}
    else if (strcmp(tmp,"TGC")==0||strcmp(tmp,"TGT")==0)
    {free(tmp);return 'C';}
    else if (strcmp(tmp,"CAA")==0||strcmp(tmp,"CAG")==0)
    {free(tmp);return 'Q';}
    else if (strcmp(tmp,"GAA")==0||strcmp(tmp,"GAG")==0)
    {free(tmp);return 'E';}
    else if (strcmp(tmp,"GGT")==0||strcmp(tmp,"GGA")==0||strcmp(tmp,"GGC")==0||strcmp(tmp,"GGG")==0)
    {free(tmp);return 'G';}
    else if (strcmp(tmp,"CAT")==0||strcmp(tmp,"CAC")==0)
    {free(tmp);return 'H';}
    else if (strcmp(tmp,"ATT")==0||strcmp(tmp,"ATC")==0||strcmp(tmp,"ATA")==0)
    {free(tmp);return 'I';}
    else if (strcmp(tmp,"CTA")==0||strcmp(tmp,"CTT")==0||strcmp(tmp,"CTC")==0||strcmp(tmp,"CTG")==0||strcmp(tmp,"TTG")==0||strcmp(tmp,"TTA")==0)
    {free(tmp);return 'L';}
    else if (strcmp(tmp,"AAG")==0||strcmp(tmp,"AAA")==0)
    {free(tmp);return 'K';}
    else if (strcmp(tmp,"ATG")==0)
    {free(tmp);return 'M';}
    else if (strcmp(tmp,"TTC")==0||strcmp(tmp,"TTT")==0)
    {free(tmp);return 'F';}
    else if (strcmp(tmp,"CCT")==0||strcmp(tmp,"CCA")==0||strcmp(tmp,"CCC")==0||strcmp(tmp,"CCG")==0)
    {free(tmp);return 'P';}
    else if (strcmp(tmp,"AGC")==0||strcmp(tmp,"AGT")==0||strcmp(tmp,"TCG")==0||strcmp(tmp,"TCC")==0||strcmp(tmp,"TCA")==0||strcmp(tmp,"TCT")==0)
    {free(tmp);return 'S';}
    else if (strcmp(tmp,"ACA")==0||strcmp(tmp,"ACT")==0||strcmp(tmp,"ACG")==0||strcmp(tmp,"ACC")==0)
    {free(tmp);return 'T';}
    else if (strcmp(tmp,"TGG")==0)
    {free(tmp);return 'W';}
    else if (strcmp(tmp,"TAC")==0||strcmp(tmp,"TAT")==0)
    {free(tmp);return 'Y';}
    else if (strcmp(tmp,"GTG")==0||strcmp(tmp,"GTA")==0||strcmp(tmp,"GTT")==0||strcmp(tmp,"GTC")==0)
    {free(tmp);return 'V';}
    else
    {free(tmp);return '*';}
}


void inverse(gsl_matrix *m,int n,gsl_matrix *invm){
    
    int i,j,s;
    gsl_matrix *tmp=gsl_matrix_alloc(n,n);
    for(i=0;i<n;i++)
        for(j=0;j<n;j++){
            gsl_matrix_set(tmp,i,j,gsl_matrix_get(m,i,j));
        }
    gsl_permutation * p = gsl_permutation_alloc (n);
    gsl_linalg_LU_decomp (tmp, p, &s);
    gsl_linalg_LU_invert (tmp, p, invm);
    gsl_permutation_free(p);
    gsl_matrix_free(tmp);
    
}

void eigen(gsl_matrix *m,int n,gsl_matrix *diagm, gsl_matrix *evecm){
    
    gsl_vector *eval = gsl_vector_alloc (n);
    gsl_eigen_symmv_workspace * w =gsl_eigen_symmv_alloc (n);
    gsl_eigen_symmv (m, eval, evecm, w);
    gsl_eigen_symmv_free (w);
    int i;
    for(i=0;i<n;i++) gsl_matrix_set(diagm,i,i,gsl_vector_get (eval, i));
    
    gsl_vector_free (eval);
    
}
gsl_matrix *TransitMatrix(int n,double t,int Mindex)
{
    gsl_matrix *result,*tmpmat;
    tmpmat = gsl_matrix_alloc(n,n);
    result = gsl_matrix_alloc(n,n);
    int i,j;
    double tmp;
    for(i=0;i<n;i++){
        tmp=gsl_matrix_get(Matrixes[Mindex].diag,i,i);
        tmp=exp(tmp*t);
        for(j=0;j<n;j++)
            gsl_matrix_set(tmpmat,j,i,tmp*gsl_matrix_get(Matrixes[Mindex].evec,j,i));
    }
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                    1.0, tmpmat, Matrixes[Mindex].inv,
                    0.0, result);
    gsl_matrix_free(tmpmat);

    return result;
}

void MakeMatrixNode(){
    int nodecount=0;
    int i,j,k;
    for (i=0;i<4;i++)
        for(j=0;j<4;j++)
            for(k=0;k<4;k++)
            {
                mnode[nodecount].str[0]=i;
                mnode[nodecount].str[1]=j;
                mnode[nodecount].str[2]=k;
                mnode[nodecount].name=acidname(mnode[nodecount].str);
                if (mnode[nodecount].name!='*') nodecount++;
                // else {printf("%d\n",count);}
            }
}
gsl_matrix *makematrix(int n,double omeg){

    gsl_matrix *oneStep;
    oneStep=gsl_matrix_alloc(n,n);
    int i,j;
    for (i=0;i<n;i++)
        for(j=0;j<n;j++){
            if(diffstr(mnode[i].str,mnode[j].str)==1){
                if(mnode[i].name==mnode[j].name) gsl_matrix_set(oneStep,i,j,lambda);
                else gsl_matrix_set(oneStep,i,j,lambda*omeg);
            }else{
                gsl_matrix_set(oneStep,i,j,0.0);
            }
        }
    for (i=0;i<n;i++){
        double sum=0.0;
        for(j=0;j<n;j++)
            sum-=gsl_matrix_get(oneStep,i,j);
        gsl_matrix_set(oneStep,i,i,sum);
    }
    return oneStep;
}

