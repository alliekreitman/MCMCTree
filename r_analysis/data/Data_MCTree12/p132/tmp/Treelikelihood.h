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
//5 5 7
int findindex(char *str){
    
    int n;
    int trans[3];
    for (n=0;n<3;n++){ // U A G C
        if (str[n]=='U') trans[n]=0;
        else if (str[n]=='A') trans[n]=1;
        else if (str[n]=='G') trans[n]=2;
        else if (str[n]=='C') trans[n]=3;
    }
    n=trans[0]*16+trans[1]*4+trans[2];
    if (n<5) return n;
    if (n==5 || n==6 || n==9) return -1;
    if (n==7) return 5;
    if (n==8) return 6;
    
    return n-3;
}

LNode TreeLikelihood(BiNode *T, int n,int TheN){ //the time here is the molecular clock

    LNode like;
    int i,j;
    for (i=0;i<TheN;i++) like.codon[i]=0.0;
    
	if (isTerminal(T)){
        char *tmp=calloc(3,sizeof(char));
        strncpy(tmp,T->sequence+n,3);
        i=findindex(tmp);
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
        
    int kNode=tNode/THEH;
    int kL=tL/THEH;
    int kR=tR/THEH;
    double tmp;
    gsl_matrix *leftresult;
    //caculate the left
    if (kNode==kL){
        leftresult=TransitMatrix(TheN,tL-tNode,kNode);
    }
    else{
        gsl_matrix *mm1=TransitMatrix(TheN,(kNode+1)*THEH-tNode,kNode);
        gsl_matrix *mm2;
        int interval=kL-kNode;
        for(i=1;i<interval;i++){
            mm2=TransitMatrix(TheN,THEH,kNode+i);
            leftresult=gsl_matrix_alloc(TheN,TheN);
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                            1.0, mm1, mm2,
                            0.0, leftresult);
            gsl_matrix_free(mm1);
            gsl_matrix_free(mm2);
            mm1=leftresult;
        }
        mm2=TransitMatrix(TheN,tL-kL*THEH,kL);
        leftresult=gsl_matrix_alloc(TheN,TheN);
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                        1.0, mm1, mm2,
                        0.0, leftresult);
        free(mm1);
        free(mm2);
    }
    for (i=0;i<TheN;i++)
        for (j=0;j<TheN;j++)
            like.codon[i]+=gsl_matrix_get(leftresult,i,j)*plikeLeft.codon[j];
    leftlist[leftN].m=leftresult;
    leftN++;
   // gsl_matrix_free(leftresult);
    gsl_matrix *rightresult;
    //caculate the right
    if(kNode==kR){
        rightresult=TransitMatrix(TheN,tR-tNode,kNode);
    }
    else{
        gsl_matrix *mm1=TransitMatrix(TheN,(kNode+1)*THEH-tNode,kNode);
        gsl_matrix *mm2;
        int interval=kR-kNode;
        for(i=1;i<interval;i++){
            mm2=TransitMatrix(TheN,THEH,kNode+i);
            rightresult=gsl_matrix_alloc(TheN,TheN);
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                            1.0, mm1, mm2,
                            0.0, rightresult);
            gsl_matrix_free(mm1);
            gsl_matrix_free(mm2);
            mm1=rightresult;
        }
        mm2=TransitMatrix(TheN,tR-kR*THEH,kR);
        rightresult=gsl_matrix_alloc(TheN,TheN);
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                        1.0, mm1, mm2,
                        0.0, rightresult);
        free(mm1);
        free(mm2);
    }
    for (i=0;i<TheN;i++){
        tmp=0.0;
        for (j=0;j<TheN;j++)
            tmp+=gsl_matrix_get(rightresult,i,j)*plikeRight.codon[j];
        
        like.codon[i]*=tmp;
    }
    rightlist[leftN].m=rightresult;
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
            printf("The codon is the stop codon!");
            exit(1);
        }
        like.codon[i]=1.0;
        return like;
    }
    
    LNode plikeLeft,plikeRight;
    plikeLeft=TreeLikelihood(T->lchi,n,TheN);
    plikeRight=TreeLikelihood(T->rchi,n,TheN);
    
    gsl_matrix *leftresult=leftlist[leftcount].m;
    leftcount++;
    for (i=0;i<TheN;i++)
        for (j=0;j<TheN;j++)
            like.codon[i]+=gsl_matrix_get(leftresult,i,j)*plikeLeft.codon[j];

    gsl_matrix *rightresult=rightlist[rightcount].m;
    rightcount++;
    for (i=0;i<TheN;i++){
        tmp=0.0;
        for (j=0;j<TheN;j++)
            tmp+=gsl_matrix_get(rightresult,i,j)*plikeRight.codon[j];
        
        like.codon[i]*=tmp;
    }
    return like;
}

double LogTreeLikelihood(BiNode *root){// p(D|G,t)

	double m=0.0,m2=0.0;
    int n=61;
    int i,j;
    
    for(i=0;i<MUNUM;i++){
        Matrixes[i].onestep = makematrix(n,omega[i]);
        Matrixes[i].evec = gsl_matrix_alloc(n,n);
        Matrixes[i].diag = gsl_matrix_alloc(n,n);
        Matrixes[i].inv = gsl_matrix_alloc(n,n);
    }
    
    for(i=0;i<MUNUM;i++){
        eigen(Matrixes[i].onestep,n,Matrixes[i].diag,Matrixes[i].evec);
        inverse(Matrixes[i].evec,n,Matrixes[i].inv);
    }
    leftN=0;
    rightN=0;
    TreeLikelihood(root,0,n);
    
    LNode temp;
    int tmplength=seqlength-3;
    for(i=3;i<tmplength;i+=3){
        m2=0.0;
        leftcount=0;
        rightcount=0;
        if (i+3>seqlength) break;
        temp=TreeLikelihood2(root,i,n);
        for (j=0; j<n; j++) m2+=temp.codon[j];
        
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
