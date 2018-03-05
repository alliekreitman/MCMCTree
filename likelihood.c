#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#define lambda 0.2
extern int seqlength;

typedef struct BiNode{

    struct BiNode *lchi, *rchi;
    int data;
    double NodeTime;
    char *sequence;

}BiNode;

typedef struct LNode{

    double T;
    double C;
    double A;
    double G;

}LNode;

int isTerminal(BiNode* T){

    if (T->lchi==NULL && T->rchi==NULL)
        return 1;
    else
        return 0;
}


double ProbJC69(int i, int j, double t){
    if (i==j) return 0.25+0.75*exp(-4*lambda*t);
    else  return 0.25-0.25*exp(-4*lambda*t);
}

LNode TreeLikelihood(BiNode *T, int n){

    LNode like;

    like.T=0.0; like.C=0.0; like.A=0.0; like.G=0.0;

	if (isTerminal(T)){
		if (*(T->sequence+n)=='T') {like.T=1.0; return like;}
		if (*(T->sequence+n)=='C') {like.C=1.0; return like;}
		if (*(T->sequence+n)=='A') {like.A=1.0; return like;}
		if (*(T->sequence+n)=='G') {like.G=1.0; return like;}

		printf("WRONG!");
	}

	LNode plikeLeft,plikeRight;
	int i,j;
	double tLeft=T->lchi->NodeTime - T->NodeTime;
	double tRight=T->rchi->NodeTime - T->NodeTime;

	plikeLeft=TreeLikelihood(T->lchi,n);
	plikeRight=TreeLikelihood(T->rchi,n);


like.T=(ProbJC69(0,0,tLeft)*plikeLeft.T+ProbJC69(0,1,tLeft)*plikeLeft.C+ProbJC69(0,2,tLeft)*plikeLeft.A+ProbJC69(0,3,tLeft)*plikeLeft.G);
like.C=(ProbJC69(1,0,tLeft)*plikeLeft.T+ProbJC69(1,1,tLeft)*plikeLeft.C+ProbJC69(1,2,tLeft)*plikeLeft.A+ProbJC69(1,3,tLeft)*plikeLeft.G);
like.A=(ProbJC69(2,0,tLeft)*plikeLeft.T+ProbJC69(2,1,tLeft)*plikeLeft.C+ProbJC69(2,2,tLeft)*plikeLeft.A+ProbJC69(2,3,tLeft)*plikeLeft.G);
like.G=(ProbJC69(3,0,tLeft)*plikeLeft.T+ProbJC69(3,1,tLeft)*plikeLeft.C+ProbJC69(3,2,tLeft)*plikeLeft.A+ProbJC69(3,3,tLeft)*plikeLeft.G);

like.T=like.T*(ProbJC69(0,0,tRight)*plikeRight.T+ProbJC69(0,1,tRight)*plikeRight.C+ProbJC69(0,2,tRight)*plikeRight.A+ProbJC69(0,3,tRight)*plikeRight.G);
like.C=like.C*(ProbJC69(1,0,tRight)*plikeRight.T+ProbJC69(1,1,tRight)*plikeRight.C+ProbJC69(1,2,tRight)*plikeRight.A+ProbJC69(1,3,tRight)*plikeRight.G);
like.A=like.A*(ProbJC69(2,0,tRight)*plikeRight.T+ProbJC69(2,1,tRight)*plikeRight.C+ProbJC69(2,2,tRight)*plikeRight.A+ProbJC69(2,3,tRight)*plikeRight.G);
like.G=like.G*(ProbJC69(3,0,tRight)*plikeRight.T+ProbJC69(3,1,tRight)*plikeRight.C+ProbJC69(3,2,tRight)*plikeRight.A+ProbJC69(3,3,tRight)*plikeRight.G);

    //printf("%lf  %lf  %lf  %lf\n",like.T,like.C,like.A,like.G);

	return like;
}

double CoalLikelihood(BiNode *T){ //p(G,t|v(\theta))

    odesolution();


}



double TheTreeLikelihood(BiNode *root){// p(D|G,t)

	double m=1,m2=0;
	int i,j;
	LNode temp;
	for(i=0;i<seqlength;i++){
		temp=TreeLikelihood(root,i);
		m2=temp.T+temp.C+temp.A+temp.G;
		m2=m2*0.25;
		printf("%.20lf\n",m2);
        m=m*m2;
	}
	return m;

}


int main(){

   BiNode *T0, *T1, *T2, *T3, *T4, *T5;

   T0=malloc(sizeof(BiNode));
   T1=malloc(sizeof(BiNode));
   T2=malloc(sizeof(BiNode));
   T3=malloc(sizeof(BiNode));
   T4=malloc(sizeof(BiNode));
   T5=malloc(sizeof(BiNode));

    T0->lchi=T1; T0->rchi=T4; T0->NodeTime=1; T0->sequence=NULL;
    T1->lchi=T2; T1->rchi=T3; T1->NodeTime=2; T1->sequence=NULL;
    T2->lchi=NULL; T2->rchi=NULL; T2->NodeTime=3; T2->sequence="AAAA";
    T3->lchi=NULL; T3->rchi=NULL; T3->NodeTime=3; T3->sequence="AAAT";
    T4->lchi=NULL; T4->rchi=NULL; T4->NodeTime=3; T4->sequence="AGTA";

	seqlength=strlen(T2->sequence);

	printf("%.20lf",TheTreeLikelihood(T0));

}
