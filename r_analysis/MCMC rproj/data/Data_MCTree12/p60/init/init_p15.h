int IsInclude(char *a,char b[20][10],int n){

    int i;
    for (i=0;i<n;i++)
        if (strcmp(a,b[i])==0)
            return 1;

    return 0;

}
int checkseq(char *seqen,int n){
    int i;
    char *tmp=calloc(3,sizeof(char));
    for(i=0;i<n;i+=3){
        strncpy(tmp,seqen+i,3);
        if (strcmp(tmp,"TAA")==0 || strcmp(tmp,"TGA")==0 || strcmp(tmp,"TAG")==0)
            return 0;
    }
    free(tmp);
    return 1;
}
//
#define DATADAYS 8
double datatime[DATADAYS+1]={-1,0,13,31,58,85,113,142,175};
int datanum[DATADAYS+1]={0,7,8,7,2,1,2,6,6};

BiNode *init()
{

    FASTAFILE *ffp;
    char *seq;
    char *name;
    int   L;
    char *filename="fasta/p15sequence.fasta";  //

	int n=0;
	BiNode* S[MAXDATASIZE];
	int top=-1;
	ffp = OpenFASTA(filename);

    if (ffp==NULL) {

        printf("wrong input!\n");
        exit(1);
    }


	// data name
	char b0[7][10]={"AY000282","AY000283","AY000284","AY000285","AY000288","AY000293","AY000315"};
	char b1[8][10]={"AY000286","AY000287","AY000289","AY000291","AY000294","AY000299","AY000314","AY000318"};
	char b2[7][10]={"AY000296","AY000300","AY000302","AY000304","AY000306","AY000316","AY000320"};
	char b3[2][10]={"AY000307","AY000313"};
    char b4[1][10]={"AY000297"};
    char b5[2][10]={"AY000290","AY000292"};
    char b6[6][10]={"AY000295","AY000308","AY000309","AY000310","AY000311","AY000312"};
    char b7[6][10]={"AY000298","AY000301","AY000303","AY000305","AY000317","AY000319"};
    //char b8[7][10]={"AY000544","AY000546","AY000549","AY000551","AY000552","AY000561","AY000585"};



	while(ReadFASTA(ffp, &seq, &name, &L)){
        S[n]=calloc(1,sizeof(BiNode));
		S[n]->lchi=NULL;
		S[n]->rchi=NULL;
		S[n]->parent=NULL;
        S[n]->sequence=calloc(1000,sizeof(char));
        S[n]->name=calloc(10,sizeof(char));
	    strcpy(S[n]->sequence,seq);
	    strncpy(S[n]->name,name,8);

        //
	    if (IsInclude(S[n]->name,b0,datanum[1])) {S[n]->NodeTime=datatime[1];}
	    else if (IsInclude(S[n]->name,b1,datanum[2])) {S[n]->NodeTime=datatime[2];}
	    else if (IsInclude(S[n]->name,b2,datanum[3])) {S[n]->NodeTime=datatime[3];}
	    else if (IsInclude(S[n]->name,b3,datanum[4])) {S[n]->NodeTime=datatime[4];}
	    else if (IsInclude(S[n]->name,b4,datanum[5])) {S[n]->NodeTime=datatime[5];}
	    else if (IsInclude(S[n]->name,b5,datanum[6])) {S[n]->NodeTime=datatime[6];}
	    else if (IsInclude(S[n]->name,b6,datanum[7])) {S[n]->NodeTime=datatime[7];}
	    else if (IsInclude(S[n]->name,b7,datanum[8])) {S[n]->NodeTime=datatime[8];}
	    //else if (IsInclude(S[n]->name,b8,datanum[9])) {S[n]->NodeTime=datatime[9];}
	    else {
                printf("%s have no time\n",S[n]->name);
                S[n]->NodeTime=252;
	    }
	    free(name);
	    free(seq);

        if(checkseq(S[n]->sequence,strlen(S[n]->sequence))) n++;
        else{
            free(S[n]->sequence);
            free(S[n]->name);
            free(S[n]);
        }
	}
	seqlength=strlen(S[1]->sequence);
	//printf("number of node=%d\n",n);
	top=n-1;
	//printf("length %d\n",strlen(S[top]->sequence));
	int i,j,len=n;
	int judge[MAXDATASIZE];
	//for (i=0;i<len;i++) printf("S[%d]=%d\n",i,S[i]->data);

	for (i=0;i<MAXDATASIZE;i++) judge[i]=i;
	//srand(time(0));
	int r;
	while(len>1){ //shuffle algorithm to build the Tree
		top++;
		S[top]=calloc(1,sizeof(BiNode));
		S[top]->name=calloc(10,sizeof(char));
		strcpy(S[top]->name,"AYNULL");
		S[top]->sequence=NULL;
		r=rand()%len;
		j=judge[r];
		S[top]->lchi=S[j];
        S[j]->parent=S[top];
		len--;
		judge[r]=judge[len];
		//printf("j=%d\n",j);

		r=rand()%len;
		j=judge[r];
		S[top]->rchi=S[j];
		S[j]->parent=S[top];
		len--;
		judge[r]=judge[len];
		judge[len]=top;
		len++;

		S[top]->NodeTime=MIN(S[top]->lchi->NodeTime,S[top]->rchi->NodeTime)-1;
		S[top]->parent=NULL;
		//printf("j=%d\n\n",j);

	}
    //printf("top=%d\n",top);
	//for (i=0;i<top+1;i++) printf("S[%d]=%d\n",i,S[i]->data);
	return S[top];
}
