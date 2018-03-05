#include <stdio.h>

long double *TreeLikelihood(BiNode *T, int n){

	if (isLeave(T)){
		if (T->sequence[n]=='A') return {1,0,0,0};
		if (T->sequence[n]=='T') return {0,1,0,0};
		if (T->sequence[n]=='C') return {0,0,1,0};
		if (T->sequence[n]=='G') return {0,0,0,1};
	}
	
	long double like1[4]={0,0,0,0},like2[4]={0,0,0,0};
	long double plikeLeft[4],plikeRight[4];
	int i,j;
	plikeLeft=Treelikelihood(T->lchi,n);
	plikeRight=Treelikelihood(T->rchi,n);
	
	for (i=0;i<4;i++)
		for (j=0;j<4;j++)
			like1[i]+=Prob(i,j,T->t)*plikeLeft[j];
	for (i=0;i<4;i++)
		for (j=0;j<4;j++)
			like2[i]+=Prob(i,j,T->t)*plikeRight[j];
	
	for (i=0;i<4;i++)
		like1[i]=like1[i]*like2[i];	
		
	return like1;
}


long double likelihood(BiNode *root){
	
	double m=1;
	int i,j;
	long double temp[4];
	for(i=0;i<seqlength;i++){
		temp[]=Treelikelihood(root,i);
    m=m*(temp[0]+temp[1]+temp[2]+temp[3])*0.25;
	}
	return m;
	
}



int main(){
	
	printf("%lf",likelihood(T));
	
}
