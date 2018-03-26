

typedef struct LNode{

  double TCAG[4];

}LNode;

double Logselectllhood(VirusSolution virus, double *SelectTime, int num, double startTime){

    double sum;
    int i;
    sum=log(1-pow(1-mu*(1-exp(-lambda*(SelectTime[0]-startTime))),virus.v[virus.sTime[0]]));

    for (i=1;i<num;i++){
        sum+=log(1-pow(1-mu*(1-exp(-lambda*(SelectTime[i]-SelectTime[i-1]))),virus.v[virus.sTime[i]]));
       // printf("%.20lf\n",sum);
    }

    return sum;
}

double ProbJC69(int i, int j, double t){
    if (i==j) return 0.25+0.75*exp(-4*lambda*t);
    else  return 0.25-0.25*exp(-4*lambda*t);
}

LNode TreeLikelihood(BiNode *T, int n){ //the time here is the molecular clock

    LNode like;

    like.TCAG[0]=0.0; like.TCAG[1]=0.0; like.TCAG[2]=0.0; like.TCAG[3]=0.0;

	if (isTerminal(T)){
		if (T->sequence[n]=='T') {like.TCAG[0]=1.0; return like;}
		if (T->sequence[n]=='C') {like.TCAG[1]=1.0; return like;}
		if (T->sequence[n]=='A') {like.TCAG[2]=1.0; return like;}
		if (T->sequence[n]=='G') {like.TCAG[3]=1.0; return like;}

		printf("WRONG!");
	}

	LNode plikeLeft,plikeRight;
	int i,j;
	double tLeft=(T->lchi->NodeTime - T->NodeTime)*lambda;
	double tRight=(T->rchi->NodeTime - T->NodeTime)*lambda;

	plikeLeft=TreeLikelihood(T->lchi,n);
	plikeRight=TreeLikelihood(T->rchi,n);

    for (i=0;i<4;i++)
        for (j=0;j<4;j++)
            like.TCAG[i]+=ProbJC69(i,j,tLeft)*plikeLeft.TCAG[j];
    double temp=0.0;

    for (i=0;i<4;i++){
        for (j=0;j<4;j++)
        temp+=ProbJC69(i,j,tRight)*plikeRight.TCAG[j];

        like.TCAG[i]=like.TCAG[i]*temp;
        temp=0.0;
    }

	return like;
}
double LogTreeLikelihood(BiNode *root){// p(D|G,t)

	double m=0,m2=0;
	int i,j;
	LNode temp;
	for(i=0;i<seqlength;i++){
		temp=TreeLikelihood(root,i);
		m2=temp.TCAG[0]+temp.TCAG[1]+temp.TCAG[2]+temp.TCAG[3];
		//m2=m2*0.25;
		//printf("%.10lf,%.10lf,%.10lf,%.10lf\n",temp.TCAG[0],temp.TCAG[1],temp.TCAG[2],temp.TCAG[3]);
        m+=log(m2);
	}
	return m;

}

int DoMatchTime(double *t,int N,
                double NodeT,
                int Tend)
{

    if (t[N-1] <= NodeT ) return N-1;
    int i;
    for (i=Tend;i>=0;i--)
        if ((t[i-1]<NodeT && t[i]>NodeT)|| fabs(t[i]-NodeT)<0.05)
           return i;
	return 0;
}

double coallikelihood(double *timelist,int n,
	                  double *datatime, int m, //begin with data[0]-1
			          int *datanum, int k,//begin with 0
			          VirusSolution virus) // double *v, double *t
{
	n--;m--;k--;
	double M=0.0;
	int num=datanum[k];
    double MK=0.0;
	int i;
    double N=0.0;
    int matchT1=ODENUM-1,matchT2=ODENUM-1;


	while(n>0){
        matchT1=DoMatchTime(virus.t,ODENUM,timelist[n],matchT2);
        matchT2=DoMatchTime(virus.t,ODENUM,timelist[n-1],matchT1);
		// printf("%d,%d\n",matchT2,matchT1);

		if (timelist[n-1]>datatime[m-1]){

            //for (i=matchT2;i<=matchT1;i++) N+=1/v[i];
            //M=num*(num-1)*M*exp(-num*(num-1)/2*N)/v[matchT2]/2;
			for (i=matchT2;i<matchT1;i++)
            {
                MK=log(1-1/virus.v[i]*num*(num-1)/2);
               // printf("%lf, v=%lf, a=%lf num=%d\n",MK,virus.v[i],1/virus.v[i]*num*(num-1)/2,num);
                N+=MK;
            }
			N+=log(1/virus.v[i]*num*(num-1)/2);
			M+=N;
            //printf("calculate coalescence at time %lf, number is %d, M=%lf, N=%lf\n",timelist[n-1],num-1,M,N);
			num--;
		}
		else{
			m--;
            if (num>1) {
                //for (i=matchT2;i<=matchT1;i++) N+=1/v[i];
                //M=M*exp(-num*(num-1)/2*timelist[n]*N);
				for (i=matchT2;i<=matchT1;i++)
                {
                    MK=log(1-1/virus.v[i]*num*(num-1)/2);
                   // printf("%lf, v=%lf, a=%lf num=%d\n",MK,virus.v[i],1/virus.v[i]*num*(num-1)/2,num);
                    N+=MK;

                }
				M+=N;
            }
            //printf("calculate no coalescence at time %lf, number is %d, and M=%lf, N=%lf\n",timelist[n-1],num-1,M,N);
			k--;
			num+=datanum[k];
		}

     N=0.0;
	 n--;
	// printf("num=%d,n=%d\n",num,n);
	}

	return M;
}



double logcoalikelihood(double *timelist, int len,
                        double *datatime, int n1,
                        int *datanum, int n2,
                        VirusSolution virus)
{

    int i,j;
    double temp;
    int judge[len];
    for (i=0;i<len;i++) judge[i]=1;

    //printf("len=%d\n",len);
    //for (i=0;i<len;i++) printf("%lf\n",timelist[i]);
    //printf("\n");

    double rmin=INF,pmin=INF-1;
    int n=-1,m=0,rcount=0;
    double result[100];

    while(rcount<len){

        for (i=0;i<len;i++)
           if (timelist[i]<rmin && judge[i]){
             rmin=timelist[i];
             n=i;
           }

        if (rmin!=pmin){
            result[m]=rmin;
            pmin=rmin;
            m++;
            judge[n]=0;
            //printf("%lf\n",rmin);
            rmin=INF;
        }
        else{
            judge[n]=0;
            rmin=INF;
        }
        rcount++;
    }

    //for (i=0;i<m;i++) printf("%lf\n",result[i]);
    //printf("\n");
    //return 3.1415926;
     return coallikelihood(result,m,datatime,n1,datanum,n2,virus);

}
