
#include <ctype.h>

#define FASTA_MAXLINE 1024	/* Requires FASTA file lines to be <1024 characters */

typedef struct fastafile_s {
    FILE *fp;
    char  buffer[FASTA_MAXLINE];
} FASTAFILE;

extern FASTAFILE *OpenFASTA(char *seqfile);
extern int        ReadFASTA(FASTAFILE *fp, char **ret_seq, char **ret_name, int *ret_L);
extern void       CloseFASTA(FASTAFILE *ffp);


FASTAFILE *OpenFASTA(char *seqfile){
  FASTAFILE *ffp;

  ffp = malloc(sizeof(FASTAFILE));
  ffp->fp = fopen(seqfile, "r");
  if (ffp->fp == NULL) { free(ffp); return NULL; }
  if ((fgets(ffp->buffer, FASTA_MAXLINE, ffp->fp)) == NULL)
    { free(ffp); return NULL; }
  return ffp;
}

int ReadFASTA(FASTAFILE *ffp, char **ret_seq, char **ret_name, int *ret_L){
  char *s;
  char *name;
  char *seq;
  int   n;
  int   nalloc;

  if (ffp->buffer[0] != '>') return 0;

  s  = strtok(ffp->buffer+1, "\n");
  name = malloc(sizeof(char) * (strlen(s)+1));
  strcpy(name, s);

  seq = malloc(sizeof(char) * 128);
  nalloc = 128;
  n = 0;
  while (fgets(ffp->buffer, FASTA_MAXLINE, ffp->fp)){

      if (ffp->buffer[0] == '>') break;

      for (s = ffp->buffer; *s != '\0'; s++){

	      if (! isalpha(*s)) continue;

	      seq[n] = *s;
	      n++;
	      if (nalloc == n){

              nalloc += 128;
              seq = realloc(seq, sizeof(char) * nalloc);
	    }
	}
    }
     seq[n] = '\0';

    *ret_name = name;
    *ret_seq  = seq;
    *ret_L    = n;
    return 1;
}

void CloseFASTA(FASTAFILE *ffp){

    fclose(ffp->fp);
    free(ffp);
}
