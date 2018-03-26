
typedef struct SIVSolution{

    double S0,I0,v0; // initial point
    double h; // time interval
    double lam,alpha,beta,ds,di,dv; // parameters
    int n; // Number of Node
    double S,I,v; //solution

}SIVSolution;

typedef struct VirusSolution{

    double *v;
    double *t;
    int *sTime;

}VirusSolution;


SIVSolution findSolution(SIVSolution ode){ // one step

    ode.S=(ode.lam+ode.S0/ode.h)/(1/ode.h+ode.beta*ode.v0+ode.ds);
    ode.I=(ode.I0/ode.h+ode.beta*ode.S*ode.v0)/(ode.di+1/ode.h);
    ode.v=(ode.alpha*ode.di*ode.I+ode.v0/ode.h)/(ode.dv+1/ode.h);
    ode.S0=ode.S;
    ode.I0=ode.I;
    ode.v0=ode.v;

    return ode;
}


VirusSolution virusNumber(double *par,
                          double *SelecTime,int timenum,
                          double *par2,int parnum,
                          double startTime)
{
    SIVSolution ode;
    ode.S0=par[0];ode.I0=par[1];ode.v0=par[2];
    ode.lam=par[3];ode.alpha=par2[0];ode.beta=par2[1];ode.ds=par[6];ode.di=par[7];ode.dv=par[8];
    ode.h=0.1;ode.n=ODENUM;


    int i,tn=0,flag=1,j,k=0;
    //int sflag=0,count=0;

    double prob,fichange;

    VirusSolution virus;
    virus.v=malloc(sizeof(double)*(ode.n+1));
    virus.t=malloc(sizeof(double)*(ode.n+1));
    virus.sTime=malloc(sizeof(int)*(timenum+1));

    virus.v[0]=ode.v0;
    virus.t[0]=0.0;
    for (i=1;i<ode.n;i++){

        ode=findSolution(ode);
        virus.v[i]=ode.v;
        virus.t[i]=(i+1)*ode.h+startTime;

       // if( virus.v[i]<0 ) printf("wrong at %d t=%lf\n",i,virus.t[i]);

       if (flag && fabs(virus.t[i]-SelecTime[tn])<0.05){

            if (flag){
                 if (tn==0) prob=mu*(1-exp(-lambda*(SelecTime[tn]-startTime)));
                 else  prob=mu*(1-exp(-lambda*(SelecTime[tn]-SelecTime[tn-1])));
                 ode.I0=(pow(prob,ode.I0)*(ode.I0*log(prob)-1)+1)/pow(log(prob),2);
                 ode.v0=(pow(prob,ode.v0)*(ode.v0*log(prob)-1)+1)/pow(log(prob),2);
                // printf("p=%lf,I0=%lf,v0=%lf\n",lambda,ode.I0,ode.v0);

                 ode.I0=MAX(800,ode.I0);
                 ode.v0=MAX(800,ode.v0);
                 ode.alpha=par2[(tn+1)*parnum];
                 ode.beta=par2[(tn+1)*parnum+1];
                 virus.sTime[k]=i;
                 k++;
                 //sflag=1;
                 //fitchange=;
                 //printf("%d\n",i);
            }
           // printf("%d\n",tn);
            tn++;
            if (tn==timenum) flag=0;
        }

    }

    return virus;

}
