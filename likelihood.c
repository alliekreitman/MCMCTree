#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define MAXDATASIZE 100

typedef struct BiNode{
	struct BiNode *lchi,*rchi;
	int data;
}BiNode;

int printTree(BiNode* T){
	if (T->lchi==NULL && T->rchi==NULL) {
		printf("%d",T->data);
		return 1;
	}
	printf("(");
	//printf("%d(",T->data);
	printTree(T->lchi);
	printf(",");
	printTree(T->rchi);
	printf(")");
	return 1;
}


BiNode *init(){
	
	int n=0;
	BiNode* S[MAXDATASIZE];
	int top=-1;
	while(n<8){
        S[n]=malloc(sizeof(BiNode));
		S[n]->lchi=NULL;
		S[n]->rchi=NULL;
		S[n]->data=n;
		n++;
	}
	top=n-1;
	int i,j,len=n;
	int judge[MAXDATASIZE];
	//for (i=0;i<len;i++) printf("S[%d]=%d\n",i,S[i]->data);
	
	for (i=0;i<MAXDATASIZE;i++) judge[i]=i;
	srand(time(0));	
	int r;
	while(len>1){ //shuffle
		top++;
		S[top]=malloc(sizeof(BiNode));
		S[top]->data=top;
		r=rand()%len;
		j=judge[r];
		S[top]->lchi=S[j];
		len--;
		judge[r]=judge[len];
		//printf("j=%d\n",j);
	
		r=rand()%len;
		j=judge[r];
		S[top]->rchi=S[j];
		len--;
		judge[r]=judge[len];
		judge[len]=top;
		len++;
		//printf("j=%d\n\n",j);

	}
	
	//for (i=0;i<top+1;i++) printf("S[%d]=%d\n",i,S[i]->data);
		
	return S[top];
}

int main(){
	
	printTree(init());
	printf("\n");
}
