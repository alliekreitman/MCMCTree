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
    double EXPb[MUNUM]={0};
    double Totalb=0.0;
    double end_time=timelist[n];

    for(i=0;i<MUNUM;i++){
        Totalb+=birth[i]*THEH;
        EXPb[i]=exp(Totalb);
    }
	while(n>0){
        t1=end_time-timelist[n];
        t2=end_time-timelist[n-1];

        if(t1<0 || t2<0) {printf("Coalesecnce Time is Wrong!\n"); exit(1);}

        int k1=t1/THEH;
        int k2=t2/THEH;

        if(k1==k2){
            N=(exp(birth[k1]*(t2-t1))-1)/birth[k1];
        }
        else{
            N=(exp(birth[k1]*((k1+1)*THEH-t1))-1)/birth[k1];
            int THEINTER=k2-k1;
            for (i=1;i<THEINTER;i++){
                N+=(exp(birth[k1+i]*THEH)-1)/birth[k1+i];
            }
            N+=(exp(birth[k2]*(t2-THEH*k2))-1)/birth[k2];
        }

        if (k1!=0) N*=EXPb[k1-1];
        N*=exp(birth[k1]*(t1-k1*THEH));
        N*=num*(num-1)/2.0/Nbegin;
        

        if (timelist[n-1]>datatime[m-1]){
            MK=EXPb[k2]*exp(birth[k2]*(t2-k2*THEH))*num*(num-1)/2.0/Nbegin;
            M+=log(MK)-N-100;
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
