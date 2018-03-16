#include <stdio.h>
#include <math.h>

double coallikelihood(double *timelist,int n
	           double *datatime, int m //begin with 0
			   int *datanum, int k)
{
	double T=datatime[--m];
	double M=1.0;
	int num=datanum[--k];
	T=T-timelist[--n];
	while(1){
		
		if (T>datatime[m-1]){
			M=M*exp(-num*(num-1)/2*timelist[n]/N)/N;
			T=T-timelist[--n];
			num--;
		}
		else{
			m--;
			M=M*exp(-num*(num-1)/2*timelist[n]/N);
			T=T-timelist[--n];
			num+=datanum[--k];
		}
		
		if (num==1 && m==1) break;
		
	}
	return M;
}


int main(){
	
}
