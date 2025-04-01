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

// change these 3 lines for each patient
#define DATADAYS 8
double datatime[DATADAYS+1]={-1,0,28,56,86,139,171,258,392};
int datanum[DATADAYS+1]={0,6,6,4,6,4,6,1,3};

BiNode *init()
{

    FASTAFILE *ffp;
    char *seq;
    char *name;
    int   L;
    char *filename="fasta/complete_data_patients/p100sequence.fasta";  // change for each patient

	int n=0;
	BiNode* S[MAXDATASIZE];
	int top=-1;
	ffp = OpenFASTA(filename);

    if (ffp==NULL) {

        printf("wrong input!\n");
        return NULL;
    }


	// data name - change for each patient
	char b0[][10]={"AY001671","AY001672","AY001674","AY001675","AY001676","AY001679"};
    	char b1[][10]={"AY001662","AY001677","AY001678","AY001681","AY001682","AY001686"};
	char b2[][10]={"AY001664","AY001670","AY001680","AY001684"};
	char b3[][10]={"AY001665","AY001683","AY001687","AY001690","AY001695","AY001699"};
  	char b4[][10]={"AY001666","AY001685","AY001689","AY001694"};
   	char b5[][10]={"AY001668","AY001669","AY001688","AY001692","AY001693","AY001698"};
   	char b6[][10]={"AY001673"};
    	char b7[][10]={"AY001661","AY001663","AY001667"};



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
