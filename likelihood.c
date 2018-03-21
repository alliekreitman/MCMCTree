#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#define MAXDATASIZE 200
#define INF 99999.0


typedef struct BiNode{

    struct BiNode *lchi, *rchi;
    double NodeTime;
    char *sequence;

}BiNode;

double timelist[MAXDATASIZE];
int listlen=0;

typedef struct SIVSolution{

    double S0,I0,v0; // initial point
    double h; // time interval
    double lam,alpha,beta,ds,di,dv; // parameters
    int n; // Number of Node
    double S,I,v; //solution

}SIVSolution;

SIVSolution findSolution(SIVSolution ode){ // one step

    ode.S=ode.lam/(ode.beta*ode.v0+ode.ds)+(ode.S0-ode.lam/(ode.beta*ode.v0+ode.ds))*exp(-(ode.beta*ode.v0+ode.ds)*ode.h);
    ode.I=ode.beta*ode.S*ode.v0/ode.di+(ode.I0-ode.beta*ode.S*ode.v0/ode.di)*exp(-ode.di*ode.h);
    ode.v=ode.alpha*ode.di*ode.I/ode.dv+(ode.v0-ode.alpha*ode.di*ode.I/ode.dv)*exp(-ode.dv*ode.h);
    ode.S0=ode.S;
    ode.I0=ode.I;
    ode.v0=ode.v;

    return ode;

}

int PrintTime(BiNode *T){

    if ((T->lchi==NULL)&&(T->rchi==NULL)){

          timelist[listlen]=T->NodeTime;
          listlen++;
          return 1;
    }

    timelist[listlen]=T->NodeTime;
    listlen++;
    PrintTime(T->lchi);
    PrintTime(T->rchi);
    return 1;
}

int DoMatchTime(double *t,int N,
                double NodeT,
                int Tend)
{

    if (t[N-1] <= NodeT ) return N-1;
    int i;
    for (i=Tend;i>=0;i--)
        if ((t[i-1]<NodeT && t[i]>NodeT)|| fabs(t[i]-NodeT)<0.001)
           return i;
	return 0;
}

double coallikelihood(double *timelist,int n,
	                  double *datatime, int m, //begin with -1
			          int *datanum, int k,//begin with 0
			          SIVSolution ode)
{
	n--;m--;
	double M=1.0;
	int num=datanum[--k];

	int i;
    double v[ode.n],t[ode.n];
    double N=1.0;
    int matchT1=ode.n-1,matchT2=ode.n-1;

    v[0]=ode.v0;t[0]=0.0;
    for (i=1;i<ode.n;i++){
        ode=findSolution(ode);
        v[i]=ode.v;
        t[i]=(i+1)*ode.h;
		 //printf("v[%d]=%lf\n",i,v[i]);
      }
  //  matchT1=DoMatchTime(t,ode.n,T,matchT2);
  //  matchT2=DoMatchTime(t,ode.n,T-timelist[n-1],matchT1);


	while(n>0){
        matchT1=DoMatchTime(t,ode.n,timelist[n],matchT2);
        matchT2=DoMatchTime(t,ode.n,timelist[n-1],matchT1);
		 //printf("%d,%d\n",matchT2,matchT1);

		if (timelist[n-1]>datatime[m-1]){

            //for (i=matchT2;i<=matchT1;i++) N+=1/v[i];
            //M=num*(num-1)*M*exp(-num*(num-1)/2*N)/v[matchT2]/2;
			for (i=matchT2;i<matchT1;i++) N=N*(1-1/v[i]*num*(num-1)/2);
			N=N*1/v[i]*num*(num-1)/2;
			M=M*N;
            //printf("calculate coalescence at time %lf, number is %d, M=%lf, N=%lf\n",
			//timelist[n-1],num-1,log(M),N);
			num--;
		}
		else{
			m--;
            if (num>1) {
                //for (i=matchT2;i<=matchT1;i++) N+=1/v[i];
                //M=M*exp(-num*(num-1)/2*timelist[n]*N);
				for (i=matchT2;i<=matchT1;i++) N=N*(1-1/v[i]*num*(num-1)/2);
				M=M*N;
            }
            //printf("calculate no coalescence at time %lf, number is %d, and M=%lf, N=%lf\n",
			//timelist[n-1],num-1,log(M),N);
			num+=datanum[--k];
		}
     N=1;
	 n--;
	}

	return M;
}



void sort(double *TimeList, int len){

    int i,j;
    double temp;
    int judge[len];

    for (i=0;i<len;i++) judge[i]=1;

    //for (i=0;i<len;i++) printf("%lf ",TimeList[i]);
    //printf("\n");

    double min=INF,pmin=INF-1;
    int n=-1,m=0,count=0;
    double result[100];

    while(count<len){

        for (i=0;i<len;i++)
           if (TimeList[i]<min && judge[i]){
             min=TimeList[i];
             n=i;
           }

        if (min!=pmin){
            result[m]=min;
            pmin=min;
            m++;
            judge[n]=0;
           // printf("%lf\n",min);
            min=INF;
        }
        else{
            judge[n]=0;
            min=INF;
        }


        count++;

    }
   // printf("%d\n",m);
   // for (i=0;i<len;i++) printf("[%d]=%d ",i,judge[i]);
   // printf("\n");
   // for (i=0;i<m;i++) printf("%lf ",result[i]);
   // printf("\n");
\
    double datatime[3]={-1,3,6};
    int datanum[3]={0,2,3};
    SIVSolution ode;
    ode.S0=10;ode.I0=0;ode.v0=20;
    ode.lam=10;ode.alpha=1000;ode.beta=3;ode.ds=0.1;ode.di=1;ode.dv=20;
    ode.h=0.01;ode.n=800;

   // printf("\n");
    printf("likelihood=%lf \n",log(coallikelihood(result,m,datatime,3,datanum,3,ode)));
    //coallikelihood(result1,m-1,datatime,3,datanum,3,ode);

}

int main(){

   int NN=9;
   BiNode *T[NN];
   int i=0,j;
   for (i=0;i<NN;i++){
       T[i]=malloc(sizeof(BiNode));
       T[i]->lchi=NULL; T[i]->rchi=NULL;
   }
   T[0]->lchi=T[1];T[0]->rchi=T[5];T[0]->NodeTime=0.0;
   T[1]->lchi=T[2];T[1]->rchi=T[3];T[1]->NodeTime=1.0;
   T[2]->lchi=T[4];T[2]->rchi=T[6];T[2]->NodeTime=2.0;
   T[3]->lchi=NULL;T[3]->rchi=NULL;T[3]->NodeTime=3.0;
   T[4]->lchi=NULL;T[4]->rchi=NULL;T[4]->NodeTime=3.0;
   T[5]->lchi=T[7];T[5]->rchi=T[8];T[5]->NodeTime=4.0;
   T[6]->lchi=NULL;T[6]->rchi=NULL;T[6]->NodeTime=6.0;
   T[7]->lchi=NULL;T[7]->rchi=NULL;T[7]->NodeTime=6.0;
   T[8]->lchi=NULL;T[8]->rchi=NULL;T[8]->NodeTime=6.0;

    for (i=0;i<MAXDATASIZE;i++) timelist[i]=0;

    listlen=0;
    PrintTime(T[0]);

    double TimeList[listlen];

    for (i=0;i<listlen;i++) {
        TimeList[i]=timelist[i];
    //    printf("%lf ",timelist[i]);
    }
    //printf("\n");

    sort(TimeList,listlen);


}
