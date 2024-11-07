double coallikelihood(double *timelist,int n)
{
	n--;
    int m=DATADAYS;
    int k=DATADAYS;
	double M=0.0;
	int num=datanum[k];
    double MK=0.0;
    double N=0.0;
    double t1,t2;
    double end_time=timelist[n];

	while(n>0){
        t1=end_time-timelist[n];
        t2=end_time-timelist[n-1];

        if(t1<0 || t2<0 || t1>t2) {printf("Coalesecnce Time is Wrong!\n"); exit(1);}
        
        N=num*(num-1)/2.0/Nbegin*(t2-t1);
        
        if (timelist[n-1]>datatime[m-1] && num>1){
            MK=num*(num-1)/2.0/Nbegin;
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
