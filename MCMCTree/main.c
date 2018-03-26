#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#define MAXDATASIZE 200
#define MAXMCTIMES 1000
#define MAXDATASIZE 200
#define INF 99999.0
#define LOWBODUND -99999.0
#define ODENUM 5000  // longest day * ode.n^-1
#define MIN(x,y)  (((x)<(y)) ?  (x) : (y))
#define MAX(x,y)  (((x)>(y)) ?  (x) : (y))
#define DATADAYS 8
double mu=0.001,lambda=0.001;
int seqlength=-1;
#include "fasta.h"
#include "tree.h"
#include "init.h"
#include "randTree.h"
#include "ode.h"
#include "likelihood.h"

void MCMC(){

    //initial data
    double datatime[DATADAYS+1]={0,4,28,56,112,140,252,392,480};
    int datanum[DATADAYS+1]={0,6,1,8,7,5,8,8,6};
    double par[9]={1000000,0,2000,1000,10000,30,1,10,200}; // S0 I0 v0 lambda alpha beta ds di dv

    //initial par
    int selectnum=4;
    int parnum=2; //the number of parameters that changing in ode
    int selectparnum=parnum*(selectnum+1);
    double par2[10]={10000,30,22000,38,23000,40,28000,40,31000,38}; //changing parameters
    double SelectTime[4]={60,90,130,160}; //select time

    int count=0,i;
    double ALPHA,r;
    double llhood,pllhood=LOWBODUND;
    double TEMPar[selectparnum],TEMPsTime[selectnum];
    double a,b;
    BiNode *TREE,*TEMPTREE;
    //VirusSolution virus;

    srand(time(0));
    TREE=init();
    datatime[0]=TREE->NodeTime-50;

    while (count<MAXMCTIMES){
        for (i=0;i<selectparnum;i++){
            if (i%2==0)
               TEMPar[i]=par2[i]+rand()%100*0.001;
            else
               TEMPar[i]=par2[i]+rand()%10000*0.001;
        }
        for (i=0;i<selectnum;i++){
            if (i==0){
                a=SelectTime[i]/2;
                b=(SelectTime[i+1]-SelectTime[i])/2;
            }
            else if (i==selectnum){
                a=(SelectTime[i]-SelectTime[i-1])/2;
                b=(datatime[DATADAYS]-SelectTime[i])/2;
            }
            else{
                a=(SelectTime[i]-SelectTime[i-1])/2;
                b=(SelectTime[i+1]-SelectTime[i])/2;
            }

            TEMPsTime[i]=SelectTime[i]+rand()%100*0.01*(a+b)-a;
        }
        //TEMPTREE=cloneTree(TREE);
        TEMPTREE=RearrangeTree(TREE);
        RandTime(TEMPTREE,datatime[0]);
        VirusSolution virus=virusNumber(par,TEMPsTime,selectnum,TEMPar,parnum,TEMPTREE->NodeTime);
        TimeListN=0;
        PrintTime(TEMPTREE);
        llhood=LogTreeLikelihood(TREE);
       // printf("%lf\n",llhood);
        llhood+=Logselectllhood(virus,TEMPsTime,selectnum,TEMPTREE->NodeTime);
       // printf("%lf\n",llhood);
        llhood+=logcoalikelihood(TimeList,TimeListN,datatime,DATADAYS+1,datanum,DATADAYS+1,virus);
        printf("%lf\n",llhood);
        ALPHA=MIN(0,llhood-pllhood);
        r=log(rand()/(RAND_MAX+0.0))*100; //change to simulated annealing
        if (r<ALPHA){
             for (i=0;i<selectparnum;i++) par2[i]=TEMPar[i];
             for (i=0;i<selectnum;i++) SelectTime[i]=TEMPsTime[i];
             pllhood=llhood;
             //DestroyTree(&TREE);
             TREE=TEMPTREE;
             count++;
             printf("count=%d\n",count);
        }
       // else
           // DestroyTree(&TEMPTREE);

    }

}


int main(){

    MCMC();

}
