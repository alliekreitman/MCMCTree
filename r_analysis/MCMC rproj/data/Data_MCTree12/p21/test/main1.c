#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#define MAXDATASIZE 500
#define MAXMCTIMES 10000
#define INF 1e300
#define LOWBODUND -1e10
#define MIN(x,y)  (((x)<(y)) ?  (x) : (y))
#define MAX(x,y)  (((x)>(y)) ?  (x) : (y))
#define MUNUM 5

//global par
double THEH;
double start_time;
double lambda=0.001;  //mutation rate
double birth[MUNUM]={0.01,0.001,0.001,0.001,0.001};
double omega[MUNUM]={0.01,0.01,0.01,0.01,0.01};
double Nbegin=10000.0;
int seqlength=-1;
int flagcount=0;
#include "header/fasta.h"
#include "header/tree.h"
#include "init/init_p1n.h"
#include "header/randTree.h"
#include "header/matrix.h"
#include "header/Coleselikelihood.h"
#include "header/Treelikelihood.h"


void MCMC(){

    int i;
    int ncount=0;
    int tcount=0;
    double ALPHA,r;
    double llhood;
    double pllhood=LOWBODUND;
    BiNode *TREE,*TEMPTREE;
    int id=getpid();
    char fname1[20],fname2[20];
    char ft[10],ft2[20];
    double tlambda,tbirth[MUNUM],tomega[MUNUM],tlam,Ntmp;

    srand(time(0)*id);
    TREE=init();
    MakeMatrixNode();
    datatime[0]=-200;

    strcpy(ft2,"result/p1n_");
    strcpy(fname1,ft2);
    strcpy(fname2,ft2);
    sprintf(ft,"%d",id);
    strcat(fname1,ft);
    strcat(fname2,ft);
    strcat(fname1,"mu.txt");
    strcat(fname2,"lam.txt");
    FILE *fp1,*fp2;
    fp1=fopen(fname1,"w+");
    fp2=fopen(fname2,"w+");
    int printflag=0;

    while (ncount<MAXMCTIMES){

        tlambda=lambda;
        r=rand()/(RAND_MAX+0.0)-0.5;
        lambda*=exp(r);

        for(i=0;i<MUNUM;i++){
            tbirth[i]=birth[i];
            r=rand()/(RAND_MAX+0.0)-0.5;
            birth[i]*=exp(r);
        }

        for(i=0;i<MUNUM;i++){
            tomega[i]=omega[i];
            r=rand()/(RAND_MAX+0.0)-0.5;
            omega[i]*=exp(r);
        }

        Ntmp=Nbegin;
        r=rand()/(RAND_MAX+0.0)-0.5;
        Nbegin*=exp(r);

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
        THEH=(BiggestTime()-start_time+1.0)/MUNUM;

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
             //printflag=1;
        }
        else{
             DestroyTree(TEMPTREE);
             lambda=tlambda;
             for (i=0;i<MUNUM;i++) {
                birth[i]=tbirth[i];
                omega[i]=tomega[i];
             }
            Nbegin=Ntmp;
            //printflag=0;
        }

        ncount++;

    if (printflag){
            fprintf(fp1,"%lf %lf %lf %lf %lf\n",birth[0],birth[1],birth[2],birth[3],birth[4]);
            fprintf(fp2,"%lf %lf %lf %lf %lf\n",omega[0],omega[1],omega[2],omega[3],omega[4]);

        }

    }
     fclose(fp1);
     fclose(fp2);
     printTree(TREE);

}


int main(){

    MCMC();

}
