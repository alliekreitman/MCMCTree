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
#define MUNUM 2

//global par
int seed = 76807;  // define seed
double THEH;
double start_time;
double lambda=0.01;  //mutation rate init
double omega1=0.1;    // pre ratio init
double omega2=0.1;    // post ratio init
double Nbegin=10000.0; // effective pop size init
int seqlength=-1;
int flagcount=0;
double drug_time=14;
#include "header/fasta.h"
#include "header/tree.h"
#include "init/init_p1.h"  //change here
#include "header/randTree.h"
#include "header/matrix.h"
#include "header/Coleselikelihood.h"
#include "header/Treelikelihood.h"

double gen_normal(double mean, double sd){
    double u1 = 1.0 - rand() / (RAND_MAX+1.0);
    double u2 = 1.0 - rand() / (RAND_MAX+1.0);
    double z0 = sqrt(-2.0 * log(u1))*cos(2.0*M_PI*u2);
    return mean+z0*sd;
}    

double gen_lognormal(double mu, double sigma){
    return exp(gen_normal(mu,sigma));
}    

double lognormal_density(double x, double mux, double sigmax) {
    if (x <= 0.000001) {
        return 0;
    }
    return (1 / (x * sigmax * sqrt(2 * M_PI))) * exp(-pow(log(x) - mux, 2) / (2 * sigmax * sigmax));
}

double log_lognormal_density(double x, double mux, double sigmax) {
     return -log(x)-log(sigmax)-0.5*log(2 * M_PI)-(pow(log(x) - mux, 2) / (2 * sigmax * sigmax));
}

void MCMC(){
    int ncount=0;
    int tcount=0;
    double propmean = 0;
    double propse = 0.3;
    double ALPHA,r;
    double lnr,lnrd1,lnrd2,lnrd3,lnrd4,lnrd5,lnrd6,lnrd7,lnrd8;
    double llhood; 
    double pllhood=LOWBODUND;
    BiNode *TREE,*TEMPTREE;
    int id=getpid();
    char fname1[30],fname2[30];
    double tlambda,to1,to2,Ntmp;

//   srand(time(0)*id);
    srand(seed);    // set seed
    TREE=init();
    MakeMatrixNode();
    datatime[0]=-200;
    
    sprintf(fname1,"result/p1_%d_%d_mu.txt",seed,id); // change here
    sprintf(fname2,"result/p1_%d_%d_omega.txt",seed,id); //change here
    
    FILE *fp1,*fp2;
    fp1=fopen(fname1,"w+");
    fp2=fopen(fname2,"w+");
    
    while (ncount<MAXMCTIMES){

// proposals
        tlambda=lambda;
        r=rand()/(RAND_MAX+0.0)-0.5;
	lnr = gen_lognormal(propmean,propse);
	lnrd1 = log_lognormal_density(lambda,propmean,propse);
        lambda*=lnr;
	lnrd2 = log_lognormal_density(lambda,propmean,propse);
	
        to1=omega1;
        r=rand()/(RAND_MAX+0.0)-0.5;
	lnr = gen_lognormal(propmean,propse);
	lnrd3 = log_lognormal_density(omega1,propmean,propse);
        omega1*=lnr;
	lnrd4 = log_lognormal_density(omega1,propmean,propse);
        
        to2=omega2;
        r=rand()/(RAND_MAX+0.0)-0.5;
	lnr = gen_lognormal(propmean,propse);
	lnrd5 = log_lognormal_density(omega2,propmean,propse);
        omega2*=lnr;
	lnrd6 = log_lognormal_density(omega2,propmean,propse);
        
        Ntmp=Nbegin;
        r=rand()/(RAND_MAX+0.0)-0.5;
	lnr = gen_lognormal(propmean,propse);
	lnrd7 = log_lognormal_density(Nbegin,propmean,propse);
        Nbegin*=lnr;
	lnrd8 = log_lognormal_density(Nbegin,propmean,propse);

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
        ALPHA=MIN(0,llhood-pllhood+(lnrd1-lnrd2+lnrd3-lnrd4+lnrd5-lnrd6+lnrd7-lnrd8));
        r=log(rand()/(RAND_MAX+0.0));

        if (r<ALPHA){
             pllhood=llhood;
             DestroyTree(TREE);
             TREE=cloneTree(TEMPTREE,NULL);
             DestroyTree(TEMPTREE);
             tcount++;
        }
        else{
            DestroyTree(TEMPTREE);
            lambda=tlambda;
            omega1=to1;
            omega2=to2;
            Nbegin=Ntmp;
        }
        
        // save output
            ncount++;
             printf("%lf %lf Unique_Tree_count=%d\n Iter_count=%d\n",tt1,tt2,tcount,ncount);
             fprintf(fp1,"%lf  %lf\n Unique_Tree_count=%d\n Iter_count=%d\n",lambda,Nbegin,tcount,ncount);
             fprintf(fp2,"%lf  %lf\n Unique_Tree_count=%d\n Iter_count=%d\n",omega1,omega2,tcount,ncount);
 
 
    }
     fclose(fp1);
     fclose(fp2);
     printTree(TREE);

}


int main(){

    MCMC();
    return 0;

}
