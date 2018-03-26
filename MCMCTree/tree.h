typedef struct BiNode{
	struct BiNode *lchi,*rchi;
	char *name;
	char *sequence;
	double NodeTime;
}BiNode;

double TimeList[200]={0};
int TimeListN=0;


int printTree(BiNode* T){
	if (T->lchi==NULL && T->rchi==NULL) {
		printf("%s",T->name);
		//printf("%c ",T->sequence[0]);
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
void DestroyTree(BiNode **T)
{
    if(*T)
    {
        if((*T)->lchi)
            DestroyTree(&(*T)->lchi);
        if((*T)->rchi)
            DestroyTree(&(*T)->rchi);
        free(*T);
        *T=NULL;
    }
}

int isTerminal(BiNode* T){

    if (T->lchi==NULL && T->rchi==NULL)
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
BiNode *cloneTree(BiNode *T){


    if (T==NULL) return NULL;

    BiNode *tmp;
    tmp=malloc(sizeof(BiNode));
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
    tmp->lchi=cloneTree(T->lchi);
    tmp->rchi=cloneTree(T->rchi);

    return tmp;

}

