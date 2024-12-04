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
#define DATADAYS 7
double datatime[DATADAYS+1]={-1,0,28,56,84,159,236,413};
int datanum[DATADAYS+1]={0,8,1,4,6,8,8,7};

BiNode *init()
{

    FASTAFILE *ffp;
    char *seq;
    char *name;
    int   L;
    char *filename="fasta/p126sequence.fasta";  //

	int n=0;
	BiNode* S[MAXDATASIZE];
	int top=-1;
	ffp = OpenFASTA(filename);

    if (ffp==NULL) {

        printf("wrong input!\n");
        return NULL;
    }


	// data name
	char b0[][10]={"AY001621","AY001627","AY001631","AY001639","AY001649","AY001656","AY001658","AY001632"};
	char b1[][10]={"AY001641"};
    char b2[][10]={"AY001626","AY001635","AY001654","AY001660"};
	char b3[][10]={"AY001619","AY001623","AY001625","AY001628","AY001651","AY001655"};
    char b4[][10]={"AY001622","AY001629","AY001636","AY001638","AY001644","AY001652","AY001653","AY001659"};
    char b5[][10]={"AY001620","AY001624","AY001630","AY001633","AY001634","AY001637","AY001640","AY001657"};
    char b6[][10]={"AY001642","AY001643","AY001645","AY001646","AY001647","AY001648","AY001650"};
    //char b7[6][10]={"AY000572","AY000574","AY000575","AY000576","AY000577","AY000582"};
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
	    //else if (IsInclude(S[n]->name,b7,datanum[8])) {S[n]->NodeTime=datatime[8];}
	    //else if (IsInclude(S[n]->name,b8,datanum[9])) {S[n]->NodeTime=datatime[9];}

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
