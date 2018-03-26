
int RandTime(BiNode *T, double fatherTime){

    if (isTerminal(T)) return 1;

    double temp=T->NodeTime;
    double childTime=MIN(T->lchi->NodeTime,T->rchi->NodeTime);
    double a,b;
    double r=rand()/(RAND_MAX+0.0);

    a=(temp-fatherTime)/2;
    b=(childTime-temp)/2;

   // printf("a=%lf,b=%lf\n",a,b);
    r=r*(a+b)-a;
    T->NodeTime=temp+r;

    RandTime(T->lchi,temp);
    RandTime(T->rchi,temp);

    return 1;

}

int RescaleTime(BiNode *T,double E){

    if (isTerminal(T)) return 1;

    T->NodeTime=T->NodeTime*E;

    RescaleTime(T->lchi,E);
    RescaleTime(T->rchi,E);

    return 1;

}

BiNode* NewTree(BiNode* Left, BiNode* Right, BiNode *cT){

    BiNode* T;
    T=malloc(sizeof(BiNode));
    T->name=calloc(10,sizeof(char));
    strcpy(T->name,cT->name);
    T->sequence=NULL;
    T->lchi=Left;
    T->rchi=Right;
    T->NodeTime=cT->NodeTime;

    return T;
}
BiNode* DoRerrange(BiNode* T){

    BiNode *gchi1,*gchi2;
    int r;
    if (T->lchi->NodeTime>T->rchi->NodeTime){

        if (isTerminal(T->rchi)) return T;

        gchi1=T->rchi->lchi;
        gchi2=T->rchi->rchi;
        r=rand()%2;
        if (r==0)
           return NewTree(gchi2,NewTree(gchi1,T->lchi,T->rchi),T);
        else
           return NewTree(gchi1,NewTree(T->lchi,gchi2,T->rchi),T);
    }

    if (T->lchi->NodeTime<T->rchi->NodeTime){

        if (isTerminal(T->lchi)) return T;

        gchi1=T->lchi->lchi;
        gchi2=T->lchi->rchi;
        r=rand()%2;
        if (r==0)
           return NewTree(NewTree(T->rchi,gchi2,T->lchi),gchi1,T);
        else
           return NewTree(NewTree(gchi1,T->rchi,T->lchi),gchi2,T);

    }

    return T;

}

BiNode* RearrangeTree(BiNode* T){

    int r=0;
    int LleaveSize=TreeLeaveSize(T->lchi);
    int RleaveSize=TreeLeaveSize(T->rchi);

    if (LleaveSize<3 && RleaveSize<3)
       return DoRerrange(T);

    if (LleaveSize < 3){
        r=rand()% (RleaveSize- 1);
        if (r==0)
            return DoRerrange(T);
        else
            return NewTree(T->lchi,RearrangeTree(T->rchi),T);
    }
    if (RleaveSize < 3){
        r=rand()%(LleaveSize-1);
        if (r==0)
            return DoRerrange(T);
        else
            return NewTree(RearrangeTree(T->lchi),T->rchi,T);
    }
    r=rand()% (RleaveSize+LleaveSize-3);
    if (r==0)
        return DoRerrange(T);
    else if (r<(LleaveSize-1))
        return NewTree(RearrangeTree(T->lchi),T->rchi,T);
    else
        return NewTree(T->lchi,RearrangeTree(T->rchi),T);


}
