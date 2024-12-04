double coallikelihood(double *timelist,int n)
{
	n--;
    int m=DATADAYS;
    int k=DATADAYS;
	double M=0.0;
	int num=datanum[k];
    double MK=0.0;
	int i;
    double N=0.0;
    double t1,t2;
    double EXPmu[MUNUM]={0};
    double Totalmu=0.0;
    for(i=0;i<MUNUM;i++){
        Totalmu+=mu[i]*THEH;
        EXPmu[i]=exp(Totalmu)*Nbegin;
    }
	while(n>0){
        t1=timelist[n-1]-start_time;
        t2=timelist[n]-start_time;
        
        if(t1<0 || t2<0) {printf("Coalesecnce Time is Wrong!\n"); exit(1);}
        
        int k1=t1/THEH;
        int k2=t2/THEH;
        
        if(k1==k2){
            N=lam/(mu[k1]-lam)*(exp((mu[k1]-lam)*(t2-t1))-1);
        }
        else{
            N=lam/(mu[k1]-lam)*(exp((mu[k1]-lam)*((k1+1)*THEH-t1))-1);
            int THEINTER=k2-k1;
            for (i=1;i<THEINTER;i++){
                N+=lam/(mu[k1+i]-lam)*(exp((mu[k1+i]-lam)*THEH)-1);
            }
            N+=lam/(mu[k2]-lam)*(exp((mu[k2]-lam)*(t2-THEH*k2))-1);
        }
        
        if (k1!=0) N*=EXPmu[k1-1];
        N*=exp(mu[k1]*(t1-k1*THEH));
        N*=num*(num-1)/2.0;
        
		if (timelist[n-1]>datatime[m-1]){
            MK=EXPmu[k2]*exp(-lam*t2)*exp(mu[k2]*(t2-k2))*num*(num-1)/2.0;
            M+=log(MK)-N;
            num--;
		}
		else{
            M-=N;
            m--;
			k--;
			num+=datanum[k];
		}

     N=0.0;
	 n--;
	}

	return M;
}
