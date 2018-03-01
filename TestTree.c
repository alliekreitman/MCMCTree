#include <stdio.h>
#include <stdlib.h>
#include <time.h>

typedef struct BiNode{

    struct BiNode *lchi, *rchi;
    int data;
    double NodeTime;

}BiNode;

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
void PrintBTree(BiNode* BT)
{
    if (BT != NULL)
    {
        if (isTerminal(BT))
            printf("%d", BT->data);

        if (BT->lchi != NULL || BT->rchi != NULL)
        {
            printf("(");
            PrintBTree(BT->lchi);
            if (BT->rchi != NULL)
                printf(",");
            PrintBTree(BT->rchi);
            printf(")");
        }
    }
}

BiNode* NewTree(BiNode* Left, BiNode* Right){

    BiNode* T;
    T=malloc(sizeof(BiNode));
    T->data=0;
    T->lchi=Left;
    T->rchi=Right;

    return T;
}
BiNode* DoRerrange(BiNode* T){

    BiNode *gchi1,*gchi2,*gchi3,*gchi4;
    int r=0;

    if (isTerminal(T->lchi)&&isTerminal(T->rchi))
        return T;

    if (isTerminal(T->lchi)){

        gchi1=T->rchi->lchi;
        gchi2=T->rchi->rchi;
        r=rand()%2;
        if (r==0)
           return NewTree(NewTree(T->lchi,gchi1),gchi2);
        else
           return NewTree(NewTree(T->lchi,gchi2),gchi1);

    }
    if (isTerminal(T->rchi)){

        gchi1=T->lchi->lchi;
        gchi2=T->lchi->rchi;
        r=rand()%2;
        if (r==0)
           return NewTree(gchi1,NewTree(gchi2,T->rchi));
        else
           return NewTree(gchi2,NewTree(gchi1,T->rchi));
    }
        gchi1=T->lchi->lchi;
        gchi2=T->lchi->rchi;
        gchi3=T->rchi->lchi;
        gchi4=T->rchi->rchi;
        r=rand()%14;
        if (r==0)
          return NewTree(NewTree(gchi1,gchi3),NewTree(gchi2,gchi4));
        if (r==1)
          return NewTree(NewTree(gchi1,gchi4),NewTree(gchi2,gchi3));
        if (r==2)
          return NewTree(NewTree(NewTree(gchi1,gchi2),gchi3),gchi4);
        if (r==3)
          return NewTree(NewTree(NewTree(gchi1,gchi3),gchi2),gchi4);
        if (r==4)
          return NewTree(NewTree(NewTree(gchi2,gchi3),gchi1),gchi4);
        if (r==5)
          return NewTree(NewTree(NewTree(gchi1,gchi2),gchi4),gchi3);
        if (r==6)
          return NewTree(NewTree(NewTree(gchi1,gchi4),gchi2),gchi3);
        if (r==7)
          return NewTree(NewTree(NewTree(gchi4,gchi2),gchi1),gchi3);
        if (r==8)
          return NewTree(NewTree(NewTree(gchi1,gchi3),gchi4),gchi2);
        if (r==9)
          return NewTree(NewTree(NewTree(gchi1,gchi4),gchi3),gchi2);
        if (r==10)
          return NewTree(NewTree(NewTree(gchi3,gchi4),gchi1),gchi2);
        if (r==11)
          return NewTree(NewTree(NewTree(gchi4,gchi2),gchi3),gchi1);
        if (r==12)
          return NewTree(NewTree(NewTree(gchi3,gchi4),gchi2),gchi1);
        if (r==13)
          return NewTree(NewTree(NewTree(gchi2,gchi3),gchi4),gchi1);

}

BiNode* RearrangeTree(BiNode* T){

    int r=0;
    int lsize=TreeSize(T->lchi);
    int rsize=TreeSize(T->rchi);
    int Tsize=TreeSize(T);
    if (Tsize==5 || Tsize==6)
       return DoRerrange(T);

    if (lsize < 4){
        r=rand()% rsize;
        if (r==0)
            return DoRerrange(T);
        else
            return NewTree(T->lchi,RearrangeTree(T->rchi));
    }
    if (rsize < 4){
        r=rand()%lsize;
        if (r==0)
            return DoRerrange(T);
        else
            return NewTree(RearrangeTree(T->lchi),T->rchi);
    }
    r=rand()% Tsize;
    printf("3   %d\n",r);
    if (r==0)
        return DoRerrange(T);
    else if (r>=1 && r<=lsize)
        return NewTree(RearrangeTree(T->lchi),T->rchi);
    else
        return NewTree(T->lchi,RearrangeTree(T->rchi));


}



BiNode* RandSelecNode(BiNode* T){

    int R=rand()%TreeSize(T);
    if (R==0)
        return T;
    else if (R>=1 && R<=TreeSize(T->lchi) && T->lchi!=NULL)
        return RandSelecNode(T->lchi);
    else
        return RandSelecNode(T->rchi);

}





int main(){

   BiNode *T[15];
   int i=0,n=1;
   for (i=0;i<15;i++){
       BiNode *tmp;
       tmp=malloc(sizeof(BiNode));
       tmp->data=i;
       tmp->lchi=NULL; tmp->rchi=NULL;
       T[i]=tmp;
   }
   T[0]->lchi=T[1];T[0]->rchi=T[2];
   T[1]->lchi=T[3];T[1]->rchi=T[4];
   T[2]->lchi=T[5];T[2]->rchi=T[6];
   T[3]->lchi=T[7];T[3]->rchi=T[8];
   T[4]->lchi=T[9];T[4]->rchi=T[10];
   T[5]->lchi=T[11];T[5]->rchi=T[12];
   T[6]->lchi=T[13];T[6]->rchi=T[14];


   PrintBTree(T[0]);
   printf("\n");
   srand(time(0));
   PrintBTree(RearrangeTree(T[0]));


}

