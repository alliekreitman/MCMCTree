typedef struct BiNode{
	struct BiNode *lchi,*rchi;
	struct BiNode *parent;
	char *name;
	char *sequence;
	double NodeTime;
}BiNode;

double TimeList[800]={0.0};
int TimeListN=0;
double sortedTime[500]={0.0};
int sortTimeN=0;

double logplus(double loga,double logb){
    
    if(loga>logb){
        return logb+log(1+exp(loga-logb));
    }else{
        return loga+log(1+exp(logb-loga));
    }
    
}

int printTree(BiNode* T){
	if (T->lchi==NULL && T->rchi==NULL) {
		printf("%s",T->name);
		//printf(" %lf",T->NodeTime-T->parent->NodeTime);
		//printf("%c ",T->sequence[0]);
		return 1;
	}
	printf("(");
	//printf("%d(",T->data);
	printTree(T->lchi);
   // printf(" %lf",T->lchi->NodeTime-T->NodeTime);
	printf(",");
	printTree(T->rchi);
   // printf(" %lf",T->lchi->NodeTime-T->NodeTime);
	printf(")");
	return 1;
}
int PrintTime(BiNode *T){

    if (T->lchi==NULL && T->rchi==NULL){

          TimeList[TimeListN]=T->NodeTime;
          TimeListN++;
          return 1;

    }

    TimeList[TimeListN]=T->NodeTime;
    TimeListN++;
    PrintTime(T->lchi);
    PrintTime(T->rchi);
    return 1;
}
int sortTime(){
    
    int i;
    int len=TimeListN;
    int judge[len];
    for (i=0;i<len;i++) judge[i]=1;
    
    double rmin=INF,pmin=INF-1;
    int n=-1,m=0,rcount=0;
    
    while(rcount<len){
        
        for (i=0;i<len;i++)
            if (TimeList[i]<rmin && judge[i]){
                rmin=TimeList[i];
                n=i;
            }
        
        if (rmin!=pmin){
            sortedTime[m]=rmin;
            pmin=rmin;
            m++;
            judge[n]=0;
            rmin=INF;
        }
        else{
            judge[n]=0;
            rmin=INF;
        }
        rcount++;
    }
    return m;
}

double BiggestTime(){
    return sortedTime[sortTimeN-1];
}
void DestroyTree(BiNode *T)
{
    if (T==NULL) return;

    DestroyTree(T->lchi);
    DestroyTree(T->rchi);

    if (T->sequence!=NULL) free(T->sequence);
    free(T->name);
    free(T);
}

int isTerminal(BiNode* T){

    if (T->lchi==NULL && T->rchi==NULL)
        return 1;
    else
        return 0;
}
int isRoot(BiNode *T){

    if (T->parent==NULL)
        return 1;
    else
        return 0;
}


int TreeSize(BiNode* T){

    if (T->lchi==NULL && T->rchi==NULL)
        return 1;
    if (T->lchi==NULL)
        return 1+TreeSize(T->rchi);
    if (T->rchi==NULL)
        return 1+TreeSize(T->lchi);
    return 1+TreeSize(T->lchi)+TreeSize(T->rchi);

}
int TreeLeaveSize(BiNode *T){

    if (T->lchi==NULL && T->rchi==NULL)
        return 1;
    if (T->lchi==NULL)
        return TreeLeaveSize(T->rchi);
    if (T->rchi==NULL)
        return TreeLeaveSize(T->lchi);
    return TreeLeaveSize(T->lchi)+TreeLeaveSize(T->rchi);

}
BiNode *cloneTree(BiNode *T,BiNode *Tp){


    if (T==NULL) return NULL;

    BiNode *tmp;
    tmp=calloc(1,sizeof(BiNode));
    tmp->NodeTime=T->NodeTime;
    tmp->name=calloc(10,sizeof(char));
    strcpy(tmp->name,T->name);
    if(T->sequence==NULL){
        tmp->sequence=NULL;
    }
    else{
        tmp->sequence=calloc(1000,sizeof(char));
        strcpy(tmp->sequence,T->sequence);
    }
    tmp->parent=Tp;
    tmp->lchi=cloneTree(T->lchi,tmp);
    tmp->rchi=cloneTree(T->rchi,tmp);


    return tmp;

}

