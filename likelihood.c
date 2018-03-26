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
