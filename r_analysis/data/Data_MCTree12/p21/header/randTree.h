int RandTime(BiNode *T){

    if (isTerminal(T)) return 1;

    double temp=T->NodeTime;
    double childTime=MIN(T->lchi->NodeTime,T->rchi->NodeTime);
    double fatherTime;
    double a,b;
    double r=rand()/(RAND_MAX+0.0);
    if (isRoot(T)){
        fatherTime=temp+log(1-rand()/(RAND_MAX+1.0))*10.0;
    }
    else{
        fatherTime=T->parent->NodeTime;
    }
    a=(temp-fatherTime)/2.0;
    b=(childTime-temp)/2.0;
    r=r*(a+b)-a;
    T->NodeTime=temp+r;
    RandTime(T->lchi);
    RandTime(T->rchi);

    return 1;

}

int RescaleTime(BiNode *T,double E){

    if (isTerminal(T)){
        T->NodeTime=T->NodeTime*E;
        return 1;
    }

    T->NodeTime=T->NodeTime*E;

    RescaleTime(T->lchi,E);
    RescaleTime(T->rchi,E);

    return 1;

}

BiNode* NewTree(BiNode* Left, BiNode* Right, BiNode *cT){

    BiNode* T;
    T=calloc(1,sizeof(BiNode));
    T->name=calloc(10,sizeof(char));
    strcpy(T->name,cT->name);
    T->sequence=NULL;
    T->lchi=Left;
    T->rchi=Right;
    T->NodeTime=cT->NodeTime;
    T->parent=cT->parent;
    free(cT->name);
    free(cT);
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

BiNode* NNItheTree(BiNode* T){

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
            return NewTree(T->lchi,NNItheTree(T->rchi),T);
    }
    if (RleaveSize < 3){
        r=rand()%(LleaveSize-1);
        if (r==0)
            return DoRerrange(T);
        else
            return NewTree(NNItheTree(T->lchi),T->rchi,T);
    }
    r=rand()% (RleaveSize+LleaveSize-3);
    if (r==0)
        return DoRerrange(T);
    else if (r<(LleaveSize-1))
        return NewTree(NNItheTree(T->lchi),T->rchi,T);
    else
        return NewTree(T->lchi,NNItheTree(T->rchi),T);


}

BiNode* searchNode(BiNode* T){

    int tsize=TreeSize(T);
    if (tsize==3) return T;
    int R=rand()%tsize;

    if (isRoot(T)){
    if (R>=1 && R<=TreeSize(T->lchi) && T->lchi!=NULL)
        return searchNode(T->lchi);
    else
        return searchNode(T->rchi);
    }

    if (R==0)
        return T;
    else if (R>=1 && R<=TreeSize(T->lchi) && T->lchi!=NULL)
        return searchNode(T->lchi);
    else
        return searchNode(T->rchi);

}

BiNode *searchNodeTime(BiNode *T,double tp){

    int tsize=TreeSize(T);
    int R=rand()%tsize;

    if(isRoot(T)){
       if (tsize==3) return T;
       if (R==0)
            return T;
       else if (R>=1 && R<=TreeSize(T->lchi) && T->lchi!=NULL)
            return searchNodeTime(T->lchi,tp);
       else
            return searchNodeTime(T->rchi,tp);

    }

    if (T->parent->NodeTime<tp){

        if (tsize==3) return T;
        if (R==0)
              return T;
        else if (R>=1 && R<=TreeSize(T->lchi) && T->lchi!=NULL)
            return searchNodeTime(T->lchi,tp);
        else
            return searchNodeTime(T->rchi,tp);

        }
        else{
          return T->parent;
        }

}


BiNode *SPRtheTree(BiNode *T)
{
    BiNode *tmp1=searchNode(T); // tmp1 can't be the root
    BiNode *p=tmp1->parent;

    //cut off
    if (isRoot(p)){
       // printf("ok1\n");
        if (p->lchi==tmp1){
             p->rchi->parent=NULL;
             T=p->rchi;
        }
        else{
             p->lchi->parent=NULL;
             T=p->lchi;
         }


    }
    else{
        //printf("ok2\n");
        BiNode  *pp=tmp1->parent->parent;
        if (p->lchi==tmp1 && pp->lchi==p)
        {
            pp->lchi=p->rchi;
            p->rchi->parent=pp;
        }
        else if (p->lchi==tmp1 && pp->rchi==p)
        {
            pp->rchi=p->rchi;
            p->rchi->parent=pp;
        }
        else if (p->rchi==tmp1 && pp->lchi==p)
        {
            pp->lchi=p->lchi;
            p->lchi->parent=pp;
        }
        else if (p->rchi==tmp1 && pp->rchi==p)
        {
            pp->rchi=p->lchi;
            p->lchi->parent=pp;
        }


    }
   BiNode *tmp2=searchNodeTime(T,tmp1->NodeTime);
   if (isRoot(tmp2)){

      if (tmp1==p->lchi)
      {
         p->rchi=tmp2;
         tmp2->parent=p;
      }
      else
      {
         p->lchi=tmp2;
         tmp2->parent=p;
      }
      p->NodeTime=MIN(tmp1->NodeTime,tmp2->NodeTime)+log(1-rand()/(RAND_MAX+1.0))*10.0;
      p->parent=NULL;

      return p;
   }

   BiNode *p2=tmp2->parent;
   if (tmp2==p2->lchi){

        p2->lchi=p;
        p->parent=p2;

   }
   else {

        p2->rchi=p;
        p->parent=p2;

   }

   if (p->lchi==tmp1){

         p->rchi=tmp2;
         tmp2->parent=p;
    }
   else{

        p->lchi=tmp2;
        tmp2->parent=p;
    }
    double a,b;
    a=MIN(tmp1->NodeTime,tmp2->NodeTime);
    b=p2->NodeTime;
    //printf("a=%lf,b=%lf\n",a,b);
    p->NodeTime=rand()/(RAND_MAX+1.0)*(a-b)+b;

  return T;
}





