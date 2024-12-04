typedef struct LNode{
    double codon[61];
}LNode;
typedef struct resultlist{
    gsl_matrix *m;
}resultlist;

resultlist leftlist[500];
int leftN=0;
int leftcount=0;
resultlist rightlist[500];
int rightN=0;
int rightcount=0;

LNode TreeLikelihood(BiNode *T, int n,int TheN){ //the time here is the molecular clock

    LNode like;
    int i,j;
    for (i=0;i<TheN;i++) like.codon[i]=0.0;

	if (isTerminal(T)){
        char *tmp=calloc(3,sizeof(char));
        strncpy(tmp,T->sequence+n,3);
        i=findindex(tmp);
        free(tmp);
        if (i==-1) {
            printf("The codon is the stop codon!");
            exit(1);
        }
        like.codon[i]=1.0;
        return like;
	}

	LNode plikeLeft,plikeRight;
	plikeLeft=TreeLikelihood(T->lchi,n,TheN);
	plikeRight=TreeLikelihood(T->rchi,n,TheN);

    double tNode=T->NodeTime-start_time;
    double tL=T->lchi->NodeTime-start_time;
    double tR=T->rchi->NodeTime-start_time;

    if(tNode>tL || tNode>tR){
        printf("Time WRONG!\n");
        exit(1);
    }

    double tmp;
    gsl_matrix *leftresult;
    //caculate the left
    if (tNode<THEH && tL<THEH){
        leftresult=TransitMatrix(TheN,tL-tNode,0);
    }
    else if (tNode<THEH && tL>THEH){
        gsl_matrix *mm1=TransitMatrix(TheN,THEH-tNode,0);
        gsl_matrix *mm2=TransitMatrix(TheN,tL-THEH,1);
        leftresult=gsl_matrix_alloc(TheN,TheN);
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                        1.0, mm1, mm2,
                        0.0, leftresult);
        gsl_matrix_free(mm1);
        gsl_matrix_free(mm2);
    }
    else{
        leftresult=TransitMatrix(TheN,tL-tNode,1);
    }
    
    for (i=0;i<TheN;i++)
        for (j=0;j<TheN;j++)
            like.codon[i]+=gsl_matrix_get(leftresult,i,j)*plikeLeft.codon[j];
    leftlist[leftN].m=leftresult;
    leftN++;
   // gsl_matrix_free(leftresult);
    gsl_matrix *rightresult;
    //caculate the right
    if (tNode<THEH && tR<THEH){
        rightresult=TransitMatrix(TheN,tR-tNode,0);
    }
    else if (tNode<THEH && tR>THEH){
        gsl_matrix *mm1=TransitMatrix(TheN,THEH-tNode,0);
        gsl_matrix *mm2=TransitMatrix(TheN,tR-THEH,1);
        rightresult=gsl_matrix_alloc(TheN,TheN);
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                        1.0, mm1, mm2,
                        0.0, rightresult);
        gsl_matrix_free(mm1);
        gsl_matrix_free(mm2);
    }
    else{
        rightresult=TransitMatrix(TheN,tR-tNode,1);
    }
    
    for (i=0;i<TheN;i++){
        tmp=0.0;
        for (j=0;j<TheN;j++)
            tmp+=gsl_matrix_get(rightresult,i,j)*plikeRight.codon[j];

        like.codon[i]*=tmp;
    }
    rightlist[rightN].m=rightresult;
    rightN++;
    //gsl_matrix_free(rightresult);

	return like;
}

LNode TreeLikelihood2(BiNode *T, int n,int TheN){ //the time here is the molecular clock

    LNode like;
    int i,j;
    for (i=0;i<TheN;i++) like.codon[i]=0.0;

    if (isTerminal(T)){
        char *tmp=calloc(3,sizeof(char));
        strncpy(tmp,T->sequence+n,3);
        i=findindex(tmp);
        if (i==-1) {
            printf("The codon is the stop codon at site %d!\n",n);
            printf("Its codon is %s and name is %s\n",tmp,T->name);
            exit(1);
        }
        free(tmp);
        like.codon[i]=1.0;
        return like;
    }

    LNode plikeLeft,plikeRight;
    plikeLeft=TreeLikelihood2(T->lchi,n,TheN);
    plikeRight=TreeLikelihood2(T->rchi,n,TheN);

    for (i=0;i<TheN;i++)
        for (j=0;j<TheN;j++)
            like.codon[i]+=gsl_matrix_get(leftlist[leftcount].m,i,j)*plikeLeft.codon[j];
    leftcount++;
    double tmp;
    for (i=0;i<TheN;i++){
        tmp=0.0;
        for (j=0;j<TheN;j++)
            tmp+=gsl_matrix_get(rightlist[rightcount].m,i,j)*plikeRight.codon[j];

        like.codon[i]*=tmp;
    }
    rightcount++;
    return like;
}

double LogTreeLikelihood(BiNode *root){// p(D|G,t)

    double m=0.0,m2=0.0;
    int n=61;
    int i,j;
    LNode temp;

    for(i=0;i<MUNUM;i++){
       // Matrixes[i].onestep = makematrix(n,omega[i]);
        Matrixes[i].evec = gsl_matrix_alloc(n,n);
        Matrixes[i].diag = gsl_matrix_alloc(n,n);
        Matrixes[i].inv = gsl_matrix_alloc(n,n);
    }
    Matrixes[0].onestep = makematrix(n,omega1);
    Matrixes[1].onestep = makematrix(n,omega2);
    
    for(i=0;i<MUNUM;i++){
        eigen(Matrixes[i].onestep,n,Matrixes[i].diag,Matrixes[i].evec);
        inverse(Matrixes[i].evec,n,Matrixes[i].inv);
    }
    leftN=0;
    rightN=0;
    temp=TreeLikelihood(root,0,n);
    for (j=0; j<n; j++) m2+=temp.codon[j];
    if(m2>1){
        printf("The Tree likelihood wrong! The codon likelihood:\n");
        for (j=0; j<n; j++)  printf("%e\n",temp.codon[j]);
        printf("%e %e",gsl_matrix_get(rightlist[3].m,3,3),gsl_matrix_get(rightlist[3].m,6,3));
        exit(1);
    }
    m+=log(m2);
    for(i=3;i<seqlength;i+=3){
        m2=0.0;
        leftcount=0;
        rightcount=0;
        if (i+3>seqlength) break;
        temp=TreeLikelihood2(root,i,n);
        for (j=0; j<n; j++) m2+=temp.codon[j];
        if(m2>1){
            printf("The Tree likelihood wrong! The codon likelihood:\n");
            for (j=0; j<n; j++)  printf("%e\n",temp.codon[j]);
            printf("%e %e",gsl_matrix_get(rightlist[3].m,3,3),gsl_matrix_get(rightlist[3].m,6,3));
            exit(1);}
        m+=log(m2);
        //printf("%lf %e\n",m,m2);
	}
    for(i=0;i<MUNUM;i++){
        gsl_matrix_free(Matrixes[i].onestep);
        gsl_matrix_free(Matrixes[i].evec);
        gsl_matrix_free(Matrixes[i].diag);
        gsl_matrix_free(Matrixes[i].inv);
    }
    for(i=0;i<leftN;i++) gsl_matrix_free(leftlist[i].m);
    for(i=0;i<rightN;i++) gsl_matrix_free(rightlist[i].m);

	return m;

}
