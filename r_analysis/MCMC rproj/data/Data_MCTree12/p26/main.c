#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#define MAXDATASIZE 500
#define MAXMCTIMES 50000
#define INF 1e300
#define LOWBODUND -1e10
#define MIN(x,y)  (((x)<(y)) ?  (x) : (y))
#define MAX(x,y)  (((x)>(y)) ?  (x) : (y))
#define MUNUM 2

//global par
double THEH;
double start_time;
double lambda=0.001;  //mutation rate
double omega1=0.1;
double omega2=0.1;
double Nbegin=10000.0;
int seqlength=-1;
int flagcount=0;
double drug_time=14.0;
#include "header/fasta.h"
#include "header/tree.h"
#include "init/init_p26.h"  //change here
#include "header/randTree.h"
#include "header/matrix.h"
#include "header/Coleselikelihood.h"
#include "header/Treelikelihood.h"


void MCMC(){

    int ncount=0;
    int tcount=0;
    double ALPHA,r;
    double llhood;
    double pllhood=LOWBODUND;
    BiNode *TREE,*TEMPTREE;
    int id=getpid();
    char fname1[30],fname2[30];
    double tlambda,to1,to2,Ntmp;

    srand(time(0)*id);
    TREE=init();
    MakeMatrixNode();
    datatime[0]=-200;
    
    sprintf(fname1,"result/p26_%d_mu.txt",id); // change here
    sprintf(fname2,"result/p26_%d_omega.txt",id); //change here
    
    FILE *fp1,*fp2;
    fp1=fopen(fname1,"w+");
    fp2=fopen(fname2,"w+");
    
    while (ncount<MAXMCTIMES){

        tlambda=lambda;
        r=rand()/(RAND_MAX+0.0)-0.5;
        lambda*=exp(r);

        to1=omega1;
        r=rand()/(RAND_MAX+0.0)-0.5;
        omega1*=exp(r);
        
        to2=omega2;
        r=rand()/(RAND_MAX+0.0)-0.5;
        omega2*=exp(r);
        
        Ntmp=Nbegin;
        r=rand()/(RAND_MAX+0.0)-0.5;
        Nbegin*=exp(0.5*r);

        // random change tree
        TEMPTREE=cloneTree(TREE,NULL);
        RandTime(TEMPTREE);
        TEMPTREE=SPRtheTree(TEMPTREE);

        //get out the time value
        flagcount=0;
        TimeListN=0;
        PrintTime(TEMPTREE);
        sortTimeN=sortTime();
        start_time=TEMPTREE->NodeTime;
        THEH=drug_time-start_time;
        
        if(THEH<0){printf("drug time too small\n"); exit(1);}

        //calculate likelihood
        double tt1=LogTreeLikelihood(TEMPTREE);
        double tt2=coallikelihood(sortedTime,sortTimeN);
        llhood=tt1+tt2;
        //printf("%lf %lf %d\n",tt,llhood-tt,ncount);
        ALPHA=MIN(0,llhood-pllhood);
        r=log(rand()/(RAND_MAX+0.0));

        if (r<ALPHA){
             pllhood=llhood;
             DestroyTree(TREE);
             TREE=cloneTree(TEMPTREE,NULL);
             DestroyTree(TEMPTREE);
             tcount++;
             printf("%lf %lf count=%d\n",tt1,tt2,tcount);
             fprintf(fp1,"%lf  %lf\n",lambda,Nbegin);
             fprintf(fp2,"%lf  %lf\n",omega1,omega2);
        }
        else{
            DestroyTree(TEMPTREE);
            lambda=tlambda;
            omega1=to1;
            omega2=to2;
            Nbegin=Ntmp;
        }

        ncount++;

    }
     fclose(fp1);
     fclose(fp2);
     printTree(TREE);

}


int main(){

    MCMC();
    return 0;

}
