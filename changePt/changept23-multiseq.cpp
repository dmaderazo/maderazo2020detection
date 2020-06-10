//From changept23.cpp
//Allow for different alphabets in the different input sequences
//Multiple chains and "pair" functionality removed - left as comments so can be made to work if needed
//Printing pis and alphas to screen removed - takes a long time


#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<math.h>
#include "MTrandom.h"
#include "keytree3.h"
#include "lfunc.h"

#define SCREEN 1
#define LOG 1
#define MAXGRPS 100
#define MAXSEQS 10
#define MAXALPHSIZE 40
#define MAXCP 4000000 //Could make this dynamic
#define FEPS 1e-7  //Check smallest possible value that can be subtracted from 1


char **inputfilename;
char *outputfilename,logfilename[500];
char *segfilename;
char *treefilename;
FILE *logfile;
unsigned long NUMBURN,NUMSAMP,SAMPBLOCK,sampsize,burnt,cnt,element;
unsigned long PF,NF;
char **alphabet;            //CHANGE: add additional layer of pointers
unsigned short *alphsize;   //CHANGE: make as pointer for vector
unsigned short effalphsize;
unsigned short sumalphsize; //CHANGE: add sumalphsize to replace nseq*alphsize
double p;
double pa,pb;
unsigned short ng;
unsigned long seed;
double phi;
char pair;

//TO DO: Add to parameters?
unsigned short numchains;
double heat;
double **chains;
double temperatures[10];

#define MAXLEN 300000000 //Max sequence length
unsigned short nseq;
unsigned short **seq;
unsigned long len;
unsigned short ESIZE[MAXSEQS],NUMBITS[MAXSEQS],CHPERE[MAXSEQS],MASK[MAXSEQS];           //CHANGE: initialise as vectors with size MAXSEQS

typedef struct cp* pcp;  //cutpoint
struct cp{
  char fixed;
  unsigned long c;        //Cut point
  pcp next; 
};

#define cElsPerBlock 100000 
// spcpAllocBlock's are allocated with this many cp elements
// and linked together as new spcpAllocBlock's may be needed
typedef struct TAG_spcpAllocBlock {
	struct TAG_spcpAllocBlock *pPrev; // pts to previously alloc'd spcpAllocBlock
	struct cp cps[cElsPerBlock];	// to hold cps
} spcpAllocBlock;

class parameters{
  private:
    pcp ppcpFree;				// pts to next cp avail for allocation, NULL if none free
    spcpAllocBlock *ppcpAllocBlocks;	// pts to MRA block
    FILE *out;
    spcpAllocBlock *AllocNewCppBlock(void);
    pcp GetCP();
    unsigned long numindex;
    unsigned long *counts;
    unsigned long ***totals;        //CHANGE: add layer to pointer
    double *part;
    double ***gamrat;               //CHANGE: add layer to pointer
    double ***sumgamrat;            //CHANGE: add layer to pointer
    unsigned long maxseglen;
    unsigned long initindex;
  public:
    pcp list;
    unsigned long numcp,numfixed;
    double ***alpha;
    double pi[MAXGRPS]; //Actually logs of pi's
    parameters();
    ~parameters();
    void init();
    void output();
    pcp insert(pcp cur);
    void remove(pcp cur);
    double *curpi;      //also logs of pis
    void setcountsandtotals();
    void setgamrat(unsigned short m,unsigned short j,double *alpha);
    void setcurpiandgamrat(double *params);
    double lnlike();
    void setpart(unsigned short g);
    double lnlike(unsigned short j);
} *seg;


typedef struct ratnode* pratnode;
struct ratnode{
  double *weights;
  double *ratio;
  pratnode *child;
};


class probabilityvector{
  private:
    unsigned long iAllocdScores;
    unsigned long numind;
    long numnodes;
    pratnode nodearray;
    pratnode root;
    void resettree(pratnode tree);
  public:
    double *probs; //Conditional Probabilities for Gibbs options
    probabilityvector();
    ~probabilityvector();
    void reset();
    void getscores(unsigned long numscores,unsigned long L,unsigned long ncp,unsigned long ncpfree);
    unsigned long choose(unsigned long numscores);
} *pr;


char *index(char *str,const char ch)
{
  while(*str != 0){
    if(ch == *str) return str;
    str++;
  }

  return NULL;
}


void getmtot(unsigned long *mtot,unsigned long L,unsigned long n,unsigned short m)      //CHANGE: use m for seq num instead of seq[m] as input argument
{
  unsigned short i;
  unsigned long j,k;
  unsigned short numch,term,ch;
  
  for(i=0;i<alphsize[m];i++)
    mtot[i]=0;
  k=L/CHPERE[m];
  numch=(unsigned short)(L%CHPERE[m]);
  term=seq[m][k]>>(numch*NUMBITS[m]);
  for(j=0;j<n;j=j+1){
    ch=term & MASK[m];
    term=term>>NUMBITS[m];
	numch++;
    if(numch==CHPERE[m]){term=seq[m][++k];numch=0;}
    mtot[ch]++; 
  }
}


unsigned long getsegmentation(FILE *in,unsigned short numgrps)
{
  unsigned long i;
  unsigned long num,numcp,idummy;
  double ddummy;

  if(fscanf(in,"%lu. %lu ",&num,&numcp)!=2){numcp=0;return 0;}
  
  if(numcp>MAXCP){
    printf("Too many cut-points in segmentation.\n");
    fflush(stdout);
    exit(1);
  }
  
  for(i=0;i<(unsigned long)(numgrps*(1+sumalphsize));i++){  //CHANGE: use sumalphsize
    if(fscanf(in,"%lf ",&ddummy)!=1){numcp=0;return 0;}
  }

  for(i=0;i<numcp+2;i++){
    if(fscanf(in,"%lu ",&idummy)!=1){numcp=0;return 0;}
  }

  return num;
}

// AllocNewCppBlock
// spcpAllocBlock allocator, called if new block required as when there are no pcp's free
spcpAllocBlock *parameters::AllocNewCppBlock(void)
{
  pcp pFree;
  spcpAllocBlock *pBlock = new spcpAllocBlock;
  int Idx;
  if(pBlock == NULL)
  {
	printf("Unable to allocate memory for pcp block\n");
    fflush(stdout);
	exit(1);
  }
  pBlock->pPrev = NULL;
  pFree = &pBlock->cps[0];
  for(Idx = 0; Idx < cElsPerBlock-1; Idx++,pFree++)
    pFree->next = pFree+1;
  pFree->next = NULL;
  return(pBlock);
}


pcp parameters::GetCP()
{
  pcp newcp;
  spcpAllocBlock *pNew;

  if(ppcpFree == NULL)		// alloc new cpp blocks if no pcp free
  {
	pNew = AllocNewCppBlock();
	pNew->pPrev = ppcpAllocBlocks;
	ppcpAllocBlocks = pNew;
	ppcpFree = &pNew->cps[0];
  }
  newcp = ppcpFree;
  ppcpFree = ppcpFree->next;

  newcp->fixed=0;
  newcp->next=NULL;
  return newcp;
}


parameters::parameters()
{
  unsigned long i;
  char outfn[102];
  unsigned short j,k;
  
  //Allocate alphas
  alpha=new double**[nseq];
  if(!alpha){printf("Could not allocate alpha in parameters::parameters.\n");fflush(stdout);exit(1);}
  for(j=0;j<nseq;j++){
    alpha[j]= new double*[ng];
    if(!alpha[j]){printf("Could not allocate alpha in parameters::parameters.\n");fflush(stdout);exit(1);}
    for(k=0;k<ng;k++){
      alpha[j][k]=new double[alphsize[j]+1];        //CHANGE: alphsize from vector
      if(!alpha[j][k]){printf("Could not allocate alpha in parameters::parameters.\n");fflush(stdout);exit(1);}
    }
  }

  ppcpAllocBlocks = AllocNewCppBlock();
  ppcpFree = &ppcpAllocBlocks->cps[0];

  //Initialise model list
  list=GetCP();
  list->c=0;
  list->fixed=1;
  list->next=GetCP();
  numcp=0;
  numfixed=0;

  //Note MAXCP should be much larger than numindex, but still feasible
  //Shouldn't need to do this for just one group, so can economise here... if(ng>1)...

  counts=new unsigned long[MAXCP+1];
  if(!counts){printf("Memory allocation problem in parameters.\n");fflush(stdout);exit(1);}

  totals=new unsigned long**[(int)nseq];             //CHANGE: vector for each seq
  for(unsigned short m=0;m<nseq;m++){                          //CHANGE: loop over seqs
      totals[m]=new unsigned long*[MAXCP+1];
      if(!totals[m]){printf("Memory allocation problem in parameters.\n");fflush(stdout);exit(1);}
  }

  part=new double[MAXCP+1];
  if(!part){printf("Memory allocation problem in parameters.\n");fflush(stdout);exit(1);}

  for(unsigned short m=0;m<nseq;m++){  //CHANGE: loop over sequences for totals
      for(i=0;i<=MAXCP;i++){
        totals[m][i]=new unsigned long[alphsize[m]+1];
        if(!totals[m][i]){printf("Memory allocation problem in parameters.\n");fflush(stdout);exit(1);}
      }
  }

  if(NF>1) sprintf(outfn,"%s.1",outputfilename);
  else strcpy(outfn,outputfilename);
  if(!(out=fopen(outfn,"w"))){
    printf("Could not open file %s.\n",outfn);
    fflush(stdout);
    exit(1);
  }
}


parameters::~parameters()
{
  spcpAllocBlock *pOld;
  unsigned long i;
  unsigned short j,k,m;
  
  fclose(out);

  delete [] part;

  for(k=0;k<nseq;k++){              //CHANGE: delete for each seq
      for(i=0;i<=MAXCP;i++)
        delete [] totals[k][i];
      delete [] totals[k];
  }
  delete [] counts;
  delete [] totals;

  for(m=0;m<nseq;m++){              //CHANGE: delete for each seq

      for(k=0;k<=alphsize[m];k++){
        if(gamrat[m][k]) delete [] gamrat[m][k];
      }
      if(gamrat[m]) delete [] gamrat[m];

      for(k=0;k<ng;k++){
        if(sumgamrat[m][k]) delete [] sumgamrat[m][k];
      }
      if(sumgamrat[m]) delete [] sumgamrat[m];
  }
  delete [] gamrat;
  delete [] sumgamrat;

  while((pOld = ppcpAllocBlocks)!=NULL)
  {
	ppcpAllocBlocks = pOld->pPrev;
	delete pOld;
  }

  for(j=0;j<nseq;j++){
    for(k=0;k<ng;k++)
      delete [] alpha[j][k];
    delete [] alpha[j];
  }
  delete [] alpha;
}


void parameters::init()
{
  unsigned long i;
  unsigned short j,k,m;
  pcp cur,newcp;
  double temp[MAXGRPS+MAXALPHSIZE],dpi[MAXGRPS];
  unsigned short chn;
  double sum;
  FILE *in;
  unsigned long seglen,pos,cpcnt;
  unsigned long segcnt;
  char line[10000];
  unsigned short numgrps;
  unsigned short nzero;
  
  //Initialise changepoint list
  if(segfilename==NULL){
    cur=list;
    for(i=1;i<len;i++){
      if(i==cur->next->c) cur=cur->next;
      else if(uniform()<p){
        newcp=insert(cur);
        newcp->c=i;
        cur=newcp;
      }
    }
    initindex=0;
  }
  else{
    //if(nseq>1){printf("Can't support -sf option if num seqs > 1\n"); exit(1);} //CHANGE: add no support for initial segmentation if nseq>1
    in=fopen(segfilename,"r");
    if(!in){printf("Could not open file %s in init.\n",segfilename);fflush(stdout);exit(1);}

    //First work out the number of groups
    if(!fgets(line,10000,in)){printf("Invalid segmentation file in init.\n");fflush(stdout);exit(1);}
    if(!fgets(line,10000,in)){printf("Invalid segmentation file in init.\n");fflush(stdout);exit(1);}
    numgrps=0;
    for(j=0;j<strlen(line);j++){
      if(line[j]==' ') numgrps++;
    }
    if(numgrps>ng){printf("Number of groups in segmentation file greater than ng.\n");fflush(stdout);exit(1);}
    else if(numgrps<ng){printf("Warning: The number of groups in the segmentation file is less than ng.\n");}
    //NB: I assume nseq is the same
    
    //Count segmentations and discard all but last
    if(SCREEN) printf("Counting segmentations.\n");
    fflush(stdout);
    if(LOG) fprintf(logfile,"Counting segmentations.\n");
    segcnt=0;
    rewind(in);    
    while(getsegmentation(in,numgrps)) segcnt++;
    if(segcnt==0){printf("Invalid segfile format in init.\n");fflush(stdout);exit(1);}
    if(SCREEN){
      printf("There are %lu valid segmentations in segfile.\n",segcnt);
      fflush(stdout);
      printf("Reading last segmentation.\n");
      fflush(stdout);
    }

    if(LOG){
      fprintf(logfile,"There are %lu valid segmentations in segfile.\n",segcnt);
      fprintf(logfile,"Reading last segmentation.\n");
    }
    rewind(in);
    for(i=1;i<segcnt;i++) getsegmentation(in,numgrps);
    
    //Read in numcp, pi and alpha
    if(fscanf(in,"%lu. %lu",&initindex,&cpcnt)!=2){printf("Invalid segfile format in init.\n");fflush(stdout);exit(1);}
    if(cpcnt>MAXCP){
      printf("Too many cut-points in segmentation.\n");
      exit(1);
    }
    if(SCREEN) printf("Initial segmentation index = %lu.\n",initindex);
    fflush(stdout);
    if(LOG) fprintf(logfile,"Initial segmentation index = %lu.\n",initindex);
  
    //First read pi
    for(j=0;j<numgrps;j++){
      if(fscanf(in,"%lf",&chains[0][j])!=1){printf("Invalid segfile format in init. Could not read distribution parameters.\n");fflush(stdout);exit(1);}
    }
    for(;j<ng;j++) chains[0][j]=0;
    //Next read alphas
    for(m=0;m<nseq;m++){
      for(j=0;j<numgrps*alphsize[m];j++){       //CHANGE: use alphsize vector
        if(fscanf(in,"%lf",&chains[0][ng+m*ng*alphsize[m]+j])!=1){printf("Invalid segfile format in init. Could not read distribution parameters.\n");fflush(stdout);exit(1);}    //CHANGE: use alphsize vector
      }
      for(k=0;k<alphsize[m];k++)    //CHANGE: use alphsize vector
        temp[k]=1.0;
      for(j=numgrps;j<ng;j++){
        sum=uniform(); 
        sum=1/(sum*sum)-1;
        dirichlet(chains[0]+ng+(m*ng+j)*alphsize[m],temp,alphsize[m]);  //CHANGE: use alphsize vector
        for(k=0;k<alphsize[m];k++)                                      //CHANGE: use alphsize vector
          chains[0][ng+(m*ng+j)*alphsize[m]+k]=exp(chains[0][ng+(m*ng+j)*alphsize[m]+k])*sum;   //CHANGE: use alphsize vector
      }
    }
    sum=0;
    nzero=0;
    for(j=0;j<ng;j++){
      if(chains[0][j]>=FEPS) sum+=chains[0][j];
      else nzero++;
    }
    if(nzero==ng){printf("All pi values zero in init!\n");fflush(stdout);exit(1);}
    sum/=1-nzero*FEPS;
    for(j=0;j<ng;j++){
      if(chains[0][j]<FEPS) chains[0][j]=log(FEPS);
      else chains[0][j]=log(chains[0][j]/sum);
    }
    for(j=0;j<ng;j++){
      pi[j]=chains[0][j];
      for(m=0;m<nseq;m++){
        sum=0;
        for(k=0;k<alphsize[m];k++){                 //CHANGE: use alphsize vector
          sum+=chains[0][ng+(m*ng+j)*alphsize[m]+k];    //CHANGE: use alphsize vector
        }
        if(sqrt(1/(1+sum))<FEPS) sum=1/(FEPS*FEPS)-1;
        else if(sqrt(1/(1+sum))>1-FEPS) sum=1/((1-FEPS)*(1-FEPS))-1;
        alpha[m][j][alphsize[m]]=sum;   //CHANGE: use alphsize vector
        sum=0;
        nzero=0;
        for(k=0;k<alphsize[m];k++){         //CHANGE: use alphsize vector
          if(chains[0][ng+(m*ng+j)*alphsize[m]+k] >= FEPS*alpha[m][j][alphsize[m]]) sum+=chains[0][ng+(m*ng+j)*alphsize[m]+k];  //CHANGE: use alphsize vector
          else nzero++;
        }
        for(k=0;k<alphsize[m];k++){        //CHANGE: use alphsize vector
          if(chains[0][ng+(m*ng+j)*alphsize[m]+k] >= FEPS*alpha[m][j][alphsize[m]])
            chains[0][ng+(m*ng+j)*alphsize[m]+k]=chains[0][ng+(m*ng+j)*alphsize[m]+k]*(1-nzero*FEPS)*alpha[m][j][alphsize[m]]/sum;  //CHANGE: use alphsize vector
          else 
            chains[0][ng+(m*ng+j)*alphsize[m]+k]=FEPS*alpha[m][j][alphsize[m]];     //CHANGE: use alphsize vector
          alpha[m][j][k]=chains[0][ng+(m*ng+j)*alphsize[m]+k];  //CHANGE: use alphsize vector
        }
      }
    }

    //Read in change-points
    if(fscanf(in,"%lu",&pos)!=1){printf("Invalid segfile format in init. Could not read change-point 1.\n");fflush(stdout);exit(1);}
    cur=list;
    for(i=0;i<cpcnt;i++){
      if(fscanf(in,"%lu",&seglen)!=1){printf("Invalid segfile format in init. Could not read change-points.\n");fflush(stdout);exit(1);}
      pos+=seglen;
      if(pos==cur->next->c) cur=cur->next;
      else{
        newcp=insert(cur);
        newcp->c=pos;
        cur=newcp;
      }
    }
    if(fscanf(in,"%lu",&seglen)!=1){printf("Invalid segfile format in init. Could not read last change-point.\n");fflush(stdout);exit(1);}
    pos+=seglen;
    if(pos!=len){printf("Invalid segfile format in init. pos!=len.\n");fflush(stdout);exit(1);}
    fclose(in);
  }
  
  for(chn=0;chn<numchains;chn++){
    if(segfilename!=NULL && chn==0) continue; //skips this section when using an input seg file
    
    for(j=0;j<ng;j++)
      temp[j]=1.0;
    dirichlet(dpi,temp,ng);
    sum=0;
    nzero=0;
    for(j=0;j<ng;j++){
      dpi[j]=exp(dpi[j]);
      if(dpi[j]>=FEPS) sum+=dpi[j];
      else nzero++;
    }
    if(nzero==ng){printf("All pi values zero in init!\n");fflush(stdout);exit(1);}
    sum/=1-nzero*FEPS;
    for(j=0;j<ng;j++){
      if(dpi[j]<FEPS) dpi[j]=log(FEPS);
      else dpi[j]=log(dpi[j]/sum);
    }
    for(j=0;j<ng;j++){
      chains[chn][j]=dpi[j];
      if(chn==0) pi[j]=dpi[j]; 
    }
    
    unsigned short curchainspos=0;  //CHANGE: new variable to keep track of position in chains vector
                                    // necessary because of differnt alphsizes
    for(m=0;m<nseq;m++){
      for(k=0;k<alphsize[m];k++) temp[k]=1.0; //CHANGE: move this line inside loop over sequences
      for(j=0;j<ng;j++){
        //OLD LINE: dirichlet(chains[chn]+ng+(m*ng+j)*alphsize[m],temp,alphsize[m]); 
        dirichlet(chains[chn]+ng+curchainspos,temp,alphsize[m]);   //CHANGE: use alphsize vector and curchainspos
        sum=0;
        nzero=0;
        for(k=0;k<alphsize[m];k++){
          //OLD LINE: chains[chn][ng+(m*ng+j)*alphsize+k]=exp(chains[chn][ng+(m*ng+j)*alphsize+k]);
          chains[chn][ng+curchainspos+k]=exp(chains[chn][ng+curchainspos+k]);
          //OLD LINE: if(chains[chn][ng+(m*ng+j)*alphsize+k]>=FEPS) sum+=chains[chn][ng+(m*ng+j)*alphsize+k];
          if(chains[chn][ng+curchainspos+k]>=FEPS) sum+=chains[chn][ng+curchainspos+k];
          else nzero++;
        }
        for(k=0;k<alphsize[m];k++){
          //OLD LINE: if(chains[chn][ng+(m*ng+j)*alphsize+k]>=FEPS) chains[chn][ng+(m*ng+j)*alphsize+k]*=(1-nzero*FEPS)/sum;
          if(chains[chn][ng+curchainspos+k]>=FEPS) chains[chn][ng+curchainspos+k]*=(1-nzero*FEPS)/sum;
          //OLD LINE: else chains[chn][ng+(m*ng+j)*alphsize+k]=FEPS;
          else chains[chn][ng+curchainspos+k]=FEPS;
        }
        sum=uniform(); 
        if(sum<FEPS) sum=FEPS;
        else if(sum>1-FEPS) sum=1-FEPS;
        sum=1/(sum*sum)-1;
        for(k=0;k<alphsize[m];k++)
          //OLD LINE: chains[chn][ng+(m*ng+j)*alphsize+k]=chains[chn][ng+(m*ng+j)*alphsize+k]*sum;
          chains[chn][ng+curchainspos+k]=chains[chn][ng+curchainspos+k]*sum;
        if(chn==0) {
          alpha[m][j][alphsize[m]]=sum;
          for(k=0;k<alphsize[m];k++)
            //OLD LINE: alpha[m][j][k]=chains[chn][ng+(m*ng+j)*alphsize+k];
            alpha[m][j][k]=chains[chn][ng+curchainspos+k];
        }
        
        curchainspos+=alphsize[m];      //CHANGE: update curchainspos by adding the current alphsize
      }
    }
  }

  numindex=0;

  gamrat=new double**[nseq];
  if(!gamrat){printf("Memory allocation problem in init.\n");fflush(stdout);exit(1);}
  sumgamrat=new double**[nseq];
  if(!sumgamrat){printf("Memory allocation problem in init.\n");fflush(stdout);exit(1);}

  //CHANGE: loop over seqs for gamrat
  for(m=0;m<nseq;m++){
      gamrat[m]=new double*[alphsize[m]+1];
      if(!gamrat[m]){printf("Memory allocation problem in init.\n");fflush(stdout);exit(1);}
      for(k=0;k<=alphsize[m];k++) gamrat[m][k]=NULL;

      sumgamrat[m]=new double*[ng];
      if(!sumgamrat[m]){printf("Memory allocation problem in init.\n");fflush(stdout);exit(1);}
      for(k=0;k<ng;k++) sumgamrat[m][k]=NULL;
  }


}


void parameters::output()
{
  pcp cur;
  unsigned long prev;
  unsigned short i,j,k;
  double LL;
  static unsigned long outfilecount=1;
  char outfn[104];
  

  if(burnt<NUMBURN){
    if(cnt==SAMPBLOCK-1){cnt=0; burnt++;}
    else cnt++;
  }
  else{
    if(cnt==element && sampsize<NUMSAMP) {
      sampsize++;
      fprintf(out,"%lu. %lu\n",sampsize+initindex,numcp);
      for(j=0;j<ng;j++)
        fprintf(out,"%lf ",exp(pi[j]));
      fprintf(out,"\n");
      for(k=0;k<nseq;k++){
        for(j=0;j<ng;j++){
          for(i=0;i<alphsize[k];i++)                //CHANGE: use alphsize vector
	        fprintf(out,"%lf ",alpha[k][j][i]);
          fprintf(out,"\n");
        }
      }

      cur=list;
      prev=0;
      while(cur){
        fprintf(out,"%lu\n",cur->c-prev);
	    prev=cur->c;
        cur=cur->next;
      }
      fprintf(out,"\n");
      fflush(out);
  
      if(sampsize%PF==0 && sampsize!=NUMSAMP){
        fclose(out);
        if(outfilecount==NF) outfilecount=1; //Re-write over old files if too many segmentations
        else outfilecount++;
        if(NF>1) sprintf(outfn,"%s.%lu",outputfilename,outfilecount);
        else strcpy(outfn,outputfilename);
        out=fopen(outfn,"w");
      }

      if(LOG || SCREEN){
        setcountsandtotals();
        setcurpiandgamrat(chains[0]);
        LL=lnlike(); 
      }

      if(SCREEN){
        printf("Sample %lu numcp=%lu ",sampsize+initindex,numcp);
        fflush(stdout);

        /* CHANGE: printing pis and alphas takes a lot of time
        for(j=0;j<ng;j++)
          printf("pi[%hu]=%lf ",j,exp(pi[j]));
        printf("\n");
        for(k=0;k<nseq;k++){
          for(j=0;j<ng;j++){
            for(i=0;i<alphsize[k];i++)                  //CHANGE: use alphsize vector
	          printf("alpha[%hu][%hu][%hu]=%lf ",k,j,i,alpha[k][j][i]);
	        printf("\n");
          }
        }
        */

        printf("Ln-likelihood=%lf\n",LL);
        fflush(stdout);
      }
      if(LOG){
        fprintf(logfile,"%lu. %lu ",sampsize+initindex,numcp);
        for(j=0;j<ng;j++)
          fprintf(logfile,"%lf ",exp(pi[j]));
        for(k=0;k<nseq;k++){
          for(j=0;j<ng;j++){
            for(i=0;i<alphsize[k];i++)                  //CHANGE: use alphsize vector
	          fprintf(logfile,"%lf ",alpha[k][j][i]);
          }
        }
        fprintf(logfile,"%lf ",LL);
        fprintf(logfile,"\n");
        fflush(logfile);
      }
    }
    if(cnt==SAMPBLOCK-1){cnt=0;element=(unsigned long)floor(uniform()*SAMPBLOCK);}
    else cnt++;
  }
}


pcp parameters::insert(pcp cur)
//Inserts change-point after cur
{
  pcp newcp;
  
  newcp=GetCP();
  newcp->next=cur->next;
  cur->next=newcp;
  numcp++;
  
  return newcp;
}


void parameters::remove(pcp cur)
//Removes changepoint AFTER cur
{
  pcp next;
  
  next=cur->next->next;
  cur->next->next = ppcpFree;
  ppcpFree = cur->next;
  cur->next=next;
  numcp--;
}


void parameters::setcountsandtotals()
{
  pcp cur;
  unsigned long seglen;
  unsigned long mtot[sumalphsize];        //CHANGE: mtot for sumalphsize
  pwtree tottree;
  unsigned long i,m;
  unsigned short k;
  long ind;

  //Set totals

  numindex=0;
  tottree=new wtree(sumalphsize*sizeof(unsigned long));     //CHANGE: use sumalphsize for word length

  cur=list;
  for(i=0;i<=numcp;i++){
    seglen=cur->next->c - cur->c;
    
    unsigned short curpos=0;                //CHANGE: use curpos to track position in mtot
    for(m=0;m<nseq;m++){
      getmtot(mtot+curpos,cur->c,seglen,m);             //CHANGE: use m as input argument
      curpos+=alphsize[m];
    }

    ind=tottree->getindex((void*)mtot);
    if(ind==-1){
      curpos=0;                             //CHANGE: use curpos to track position in mtot
      for(m=0;m<nseq;m++){
        for(k=0;k<alphsize[m];k++)
          totals[m][numindex][k]=mtot[curpos+k];
        totals[m][numindex][alphsize[m]]=seglen;
        curpos+=alphsize[m];
      }
      counts[numindex]=1;
      tottree->insertword((void*)mtot,(long)numindex++);
    }
    else
      counts[ind]++;


    /*OLD CODE
    ind=tottree->getindex((void*)mtot);
    if(ind==-1){
      for(k=0;k<nseq*alphsize;k++)
        totals[numindex][k]=mtot[k];
      totals[numindex][nseq*alphsize]=seglen;
      counts[numindex]=1;
      tottree->insertword((void*)mtot,(long)numindex++);
    }
    else 
      counts[ind]++;
    */
    
    cur=cur->next;
  }

  delete tottree;
}


void parameters::setgamrat(unsigned short m,unsigned short j,double *alpha)//CHANGE: this function is changed for new gamrat and sumgamrat
{
  unsigned long i;
  unsigned short k;
  double rat;
  // unsigned long *tot;

  for(k=0;k<=alphsize[m];k++){
    gamrat[m][k][0]=rat=0;
    for(i=1;i<=maxseglen;i++){
      rat+=log(i-1+alpha[k]);
      gamrat[m][k][i]=rat;
    }
  }

  //OLD CODE: j=j*nseq+m;
  for(i=0;i<numindex;i++){
    sumgamrat[m][j][i]=-gamrat[m][alphsize[m]][totals[m][i][alphsize[m]]];
    //OLD CODE: tot=totals[i]+m*alphsize;
    for(k=0;k<alphsize[m];k++)
      sumgamrat[m][j][i]+=gamrat[m][k][totals[m][i][k]];
  }
}


void parameters::setcurpiandgamrat(double *params)
{
  unsigned long i;
  unsigned short j,k,m; //,n
  double *curalpha;
  double sum[MAXGRPS];
  double alpha;
  double temp;
  double rat;
  unsigned short nzero;
  //unsigned long *tot;

  curpi=params;

  //renormalise curpi
  temp=0;
  nzero=0;
  for(j=0;j<ng;j++){
    curpi[j]=exp(curpi[j]);
    if(curpi[j]>=FEPS) temp+=curpi[j];
    else nzero++;
  }
  if(nzero==ng){printf("All pi values zero in init!\n");fflush(stdout);exit(1);}
  temp/=1-nzero*FEPS;
  for(j=0;j<ng;j++){
    if(curpi[j]<FEPS) curpi[j]=log(FEPS);
    else curpi[j]=log(curpi[j]/temp);
  }

  //update maxseglen
  maxseglen=0;
  for(i=0;i<numindex;i++){
    if(totals[0][i][alphsize[0]]>maxseglen) maxseglen=totals[0][i][alphsize[0]]; //CHANGE: seglengths the same for each seq, just check the first seq
  }

  //Reallocate gamrat and sumgamrat
  for(m=0;m<nseq;m++){      //CHANGE: loop over seqs for gamrat and sumgamrat
      for(k=0;k<=alphsize[m];k++){
        if(gamrat[m][k]) delete [] gamrat[m][k];
        gamrat[m][k]=new double[maxseglen+1];
        if(!gamrat[m][k]){printf("Allocation problem in in setcurpiandgamrat.\n");fflush(stdout);exit(1);}
      }
      for(k=0;k<ng;k++){
        if(sumgamrat[m][k]) delete [] sumgamrat[m][k];
        sumgamrat[m][k]=new double[numindex];
        if(!sumgamrat[m][k]){printf("Allocation problem in setcurpiandgamrat.\n");fflush(stdout);exit(1);}
      }
  }

  unsigned short curparamspos = 0;      //CHANGE: use curparamspos to get position in params vector with multiseq
  for(m=0;m<nseq;m++){
    //OLD CODE: curalpha=params+ng+m*ng*alphsize;
    curalpha=params+ng+curparamspos;
    //sum alphas  
    for(j=0;j<ng;j++){
      sum[j]=0;
      for(k=0;k<alphsize[m];k++)
        sum[j]+=curalpha[j*alphsize[m]+k];      //CHANGE: use alphsize vector
    }

    //Set all gamrats
    for(j=0;j<ng;j++){
      for(k=0;k<=alphsize[m];k++){
        if(k<alphsize[m]) alpha=curalpha[j*alphsize[m]+k];
        else alpha=sum[j];
        gamrat[m][k][0]=rat=0;
        for(i=1;i<=maxseglen;i++){
          rat+=log(i-1+alpha);
          gamrat[m][k][i]=rat;
        }
      }
      //OLD CODE: n=j*nseq+m;
      for(i=0;i<numindex;i++){
        sumgamrat[m][j][i]=-gamrat[m][alphsize[m]][totals[m][i][alphsize[m]]];
        //tot=totals[i]+m*alphsize;
        for(k=0;k<alphsize[m];k++)
          sumgamrat[m][j][i]+=gamrat[m][k][totals[m][i][k]];
      }
    }

    curparamspos+=ng*alphsize[m];       //CHANGE: update curparamspos
  }
}


double parameters::lnlike()
{
  unsigned long i;
  unsigned short j,m;
  double temp1,temp2,temp3;
  double retval;
  
  retval=0;
  temp2=0;
  for(i=0;i<numindex;i++){
      temp2=curpi[0];
      for(m=0;m<nseq;m++){
          temp2+=sumgamrat[m][0][i];
      }

      for(j=1;j<ng;j++){
          temp1=curpi[j];
          for(m=0;m<nseq;m++){
              temp1+=sumgamrat[m][j][i];
          }

          temp3=temp2-temp1;
          if(temp3<0){temp3=-temp3;temp2=temp1;}
          temp2=temp2+lfunc2(temp3);
      }

      retval+=counts[i]*temp2;
  }

  /* OLD CODE
  retval=0;
  for(i=0;i<numindex;i++){
    temp2=curpi[0];
    for(m=0;m<nseq;m++) temp2+=sumgamrat[m][i];
    for(j=1;j<ng;j++){
      temp1=curpi[j];
      for(;m<(j+1)*nseq;m++) temp1+=sumgamrat[m][i];
        
      //temp2=lnsum2(temp2,temp1);
      temp3=temp2-temp1;
      if(temp3<0){temp3=-temp3;temp2=temp1;}
      temp2=temp2+lfunc2(temp3);
    }
    retval+=counts[i]*temp2;
  } */

  return retval;
}

void parameters::setpart(unsigned short g)
{
  unsigned long i;
  unsigned short j,m;
  double temp1,temp2,temp3;

  for(i=0;i<numindex;i++){
      if(g==0) {j=1;}
      else {j=0;}
      temp2=curpi[j];
      for(m=0;m<nseq;m++){
          temp2+=sumgamrat[m][j][i];
      }
      for(j++;j<ng;j++){
          if(j==g){continue;}
          temp1=curpi[j];
          for(m=0;m<nseq;m++){
              temp1+=sumgamrat[m][j][i];
          }
          temp3=temp2-temp1;
          if(temp3<0){temp3=-temp3;temp2=temp1;}
          temp2=temp2+lfunc2(temp3);
      }
      part[i]=temp2;
  }

  /* OLD CODE
  for(i=0;i<numindex;i++){
    if(g==0) {j=1;m=nseq;}
    else {j=0;m=0;}
    temp2=curpi[j];
    for(;m<(j+1)*nseq;m++) temp2+=sumgamrat[m][i];
    for(j++;j<ng;j++){
      if(j==g) {m+=nseq; continue;}
      temp1=curpi[j];
      for(;m<(j+1)*nseq;m++) temp1+=sumgamrat[m][i];
        
      //temp2=lnsum2(temp2,temp1);
      temp3=temp2-temp1;
      if(temp3<0){temp3=-temp3;temp2=temp1;}
      temp2=temp2+lfunc2(temp3);
    }
    part[i]=temp2;
  } 
  */
}


double parameters::lnlike(unsigned short j)
{
  unsigned long i;
  unsigned short m;
  double temp1,temp2,temp3;
  double retval;
  
  retval=0;
  for(i=0;i<numindex;i++){
      temp1=curpi[j];
      for(m=0;m<nseq;m++){
          temp1+=sumgamrat[m][j][i];
      }
      if(ng>1){
          temp2=part[i];
          temp3=temp2-temp1;
          if(temp3<0){temp3=-temp3;temp2=temp1;}
          temp2=temp2+lfunc2(temp3);
          retval+=counts[i]*temp2;
      }
      else{
          retval+=counts[i]*temp1;
      }
  }

  /* OLD CODE
  for(i=0;i<numindex;i++){
    temp1=curpi[j];
    for(m=j*nseq;m<(j+1)*nseq;m++) temp1+=sumgamrat[m][i];

    if(ng>1){
      temp2=part[i];
  
      //temp2=lnsum2(temp2,temp1);
      temp3=temp2-temp1;
      if(temp3<0){temp3=-temp3;temp2=temp1;}
      temp2=temp2+lfunc2(temp3);
  
      retval+=counts[i]*temp2;
    }
    else retval+=counts[i]*temp1;
  } 
  */

  return retval;
}


probabilityvector::probabilityvector()
{
  FILE *in;
  long nodeID;
  short i;
  char leaf;
  long *childID;

  if(treefilename){
    childID=new long[effalphsize];
    if(!childID){printf("Could not allocate childID in probabilityvector::probabilityvector.\n");fflush(stdout);exit(1);}

    in=fopen(treefilename,"r");
    if(!in){printf("Could not open input file tree.\n");fflush(stdout); exit(1);}

    if(fscanf(in,"%ld",&numnodes)!=1){printf("Invalid tree file.\n");fflush(stdout);exit(1);}

    nodearray=new ratnode[numnodes];
    if(!nodearray){printf("Could not allocate nodearray in probabilityvector::probabilityvector");fflush(stdout);exit(1);}

    nodeID=0;
    while(fscanf(in,"%ld ",childID)==1){
      if(childID[0]) leaf=0;
      else leaf=1;
      for(i=1;i<effalphsize;i++){
        if(fscanf(in,"%ld ",childID+i)!=1){printf("Invalid input file.\n");fflush(stdout);exit(1);}
        if(childID[i]) leaf=0;
      }
      if(leaf) nodearray[nodeID].child=NULL;
      else{
        nodearray[nodeID].child=new pratnode[effalphsize];
        if(!nodearray[nodeID].child){printf("Could not allocate in probabilityvector::probabilityvector.\n");fflush(stdout);exit(1);}
        for(i=0;i<effalphsize;i++){
          if(childID[i])
            nodearray[nodeID].child[i]=nodearray+(childID[i]-1);
          else
            nodearray[nodeID].child[i]=NULL;
        }
      }
      nodearray[nodeID].ratio=new double[effalphsize];
      if(!nodearray[nodeID].ratio){printf("Could not allocate in probabilityvector::probabilityvector.\n");fflush(stdout);exit(1);}
      for(i=0;i<effalphsize;i++)
        nodearray[nodeID].ratio[i]=0;
      //Delete these two lines if you want to allocate them dynamically
      nodearray[nodeID].weights=new double[ng];
      if(!nodearray[nodeID].weights){printf("Could not allocate in probabilityvector::probabilityvector.\n");exit(1);}
      //nodearray[nodeID].weights=NULL;
      nodeID++;
    }
    fclose(in);
    delete [] childID;
  }
  else{
    numnodes=1;
    nodearray=new ratnode[1];
    if(!nodearray){printf("Could not allocate nodearray in probabilityvector::probabilityvector");fflush(stdout);exit(1);}
    nodearray[0].child=NULL;
    nodearray[0].ratio=new double[effalphsize]; //This could be a problem if effalphsize gets large
    if(!nodearray[0].ratio){printf("Could not allocate in probabilityvector::probabilityvector.\n");fflush(stdout);exit(1);}
    nodearray[0].weights=new double[ng];
    if(!nodearray[0].weights){printf("Could not allocate in probabilityvector::probabilityvector.\n");fflush(stdout);exit(1);}
  }
  root=nodearray;
  
  probs=NULL;
  iAllocdScores=0;
}


probabilityvector::~probabilityvector()
{
  long i;

  for(i=0;i<numnodes;i++){
    if(nodearray[i].child) delete [] nodearray[i].child;
    if(nodearray[i].ratio) delete [] nodearray[i].ratio;
    if(nodearray[i].weights) delete [] nodearray[i].weights;
  }
  delete [] nodearray;

  if(probs)
	delete [] probs;
}


void probabilityvector::resettree(pratnode tree)
{
  short i;

  if(!tree) return;
  for(i=0;i<effalphsize;i++){
    if(tree->ratio[i]>0) break;
  }
  if(i==effalphsize) return;  //Already reset

  //Could this be done more economically?

  for(i=0;i<effalphsize;i++)
    tree->ratio[i]=0;
  if(tree->child){
    for(i=0;i<effalphsize;i++)
      resettree(tree->child[i]);
  }
}


void probabilityvector::reset()
{
  short i;

  for(i=0;i<ng;i=i+1)
    root->weights[i]=exp(seg->pi[i]); //Normalise?

  if(root->child){
    for(i=0;i<effalphsize;i++) resettree(root->child[i]);
  }
}


void probabilityvector::getscores(unsigned long numscores,unsigned long L,unsigned long ncp,unsigned long ncpfree)
{
  unsigned short i,m;
  unsigned long k[MAXSEQS];                             //CHANGE: have a k for each seq
  unsigned long l;
  unsigned short ch[MAXSEQS],term[MAXSEQS],numch[MAXSEQS];      //CHANGE: have a numch for each seq
  unsigned long nummoves;
  unsigned long mtot[MAXSEQS][MAXALPHSIZE];
  double weights[MAXGRPS],ratio;
  unsigned long R;
  pratnode cur,parent;
  char newcalc;
  unsigned short index;
  double *alf;
  unsigned long cnt[MAXSEQS];


  //Allocate work space
  // dynamically grow workspace as needed...
  if(probs == NULL || numscores > iAllocdScores)
  {
	if(probs != NULL)
		delete [] probs;
	iAllocdScores = (numscores*4)/3;	
	probs = new double[iAllocdScores];
  }	
  if(probs == NULL)
  {
    printf("Memory allocation problem in getscores.\n");fflush(stdout);
    exit(1);
  }
  
  //Get scores for positions
  //OLD LINE: nummoves=3*ncp+1 + numchains*(ng*(alphsize + 1)) + numchains-1;
  //CHECK: was this correct for multiple seqs?
  nummoves=3*ncp+1 + numchains*(ng*(sumalphsize + 1)) + numchains-1; //CHANGE: use sumalphsize
  if(ng>1){
    if(pair) nummoves+=numchains*ng*ng; //NOT sure pifn2 is correct for multiple sequences, so don't use pair for that case
    else nummoves+=numchains*ng;    //CHECK: why add ng again?
  }
  if(phi>0)
    probs[0]=((1.0-phi)/phi) * (1.0 + 3.0/nummoves);
  else
    probs[0]=((len-2.0-ncp+pb)/(ncpfree+pa))*(1.0 + 3.0/nummoves);
  //Note: includes prior and GGS adjustment

  //Calculate ratios for left part

  for(m=0;m<nseq;m++){
      k[m]=L/CHPERE[m];                            //CHANGE: Bring inside loop over seqs and use CHPERE and k vector
      numch[m]=(unsigned short)(L%CHPERE[m]);      //CHANGE: Bring inside loop over seqs and use CHPERE and numch vector
      for(i=0;i<alphsize[m];i=i+1)            //CHANGE: use alhpsize vector
      mtot[m][i]=0;

    term[m]=seq[m][k[m]]>>(numch[m]*NUMBITS[m]);
  }

  //Do left segment ratio calculations using stored results
  parent=root;
  newcalc=1;
  for(l=0;l<numscores-1;l=l+1){
    index=0;
    for(m=0;m<nseq;m++){
      ch[m]=term[m] & MASK[m];
      index=index*alphsize[m]+ch[m]; //CHECK       //CHANGE: use alphsize vector
      term[m]=term[m]>>NUMBITS[m];

      numch[m]++;
      if(numch[m]==CHPERE[m]){
        k[m]++;
        term[m]=seq[m][k[m]];
        numch[m]=0;
      }
    }

    /* OLD CODE - CHANGE: MOVED INSIDE ABOVE LOOP
    numch++;
    if(numch==CHPERE){
      k++;
      for(m=0;m<nseq;m++)
        term[m]=seq[m][k];
      numch=0;
    }
    */

    if(!parent->child || !parent->child[index]) break; //End of path in tree

    cur=parent->child[index];
    if(cur->ratio[index]>0){ //Is stored
      probs[l+1]=cur->ratio[index];
      newcalc=1;
    }
    else{
      if(newcalc){ //Initialise weights
        for(i=0;i<ng;i=i+1)
          weights[i]=parent->weights[i];
        newcalc=0;
      }
      ratio=0;
      for(m=0;m<nseq;m=m+1) cnt[m]=mtot[m][ch[m]];
      for(i=0;i<ng;i=i+1){
        for(m=0;m<nseq;m=m+1){
          alf=seg->alpha[m][i];
          weights[i]=weights[i]*(cnt[m]+alf[ch[m]])/(l+alf[alphsize[m]]); //Most time spent here //CHANGE: use alpfsize vector
        }
        ratio+=weights[i];
      }
      for(i=0;i<ng;i=i+1){
        weights[i]=weights[i]/ratio;
        cur->weights[i]=weights[i]; //Store them
      }
      probs[l+1]=cur->ratio[index]=ratio;
    }
    for(m=0;m<nseq;m++)
	  mtot[m][ch[m]]++; 
    parent=cur;
  }

  //Reinitialise weights
  if(newcalc && l<numscores-1){
    for(i=0;i<ng;i=i+1)
      weights[i]=parent->weights[i];
  }

  //Do rest of left segment ratio calculation without storage
  for(;l<numscores-1;l=l+1){
    ratio=0;
    for(m=0;m<nseq;m=m+1) cnt[m]=mtot[m][ch[m]];
    for(i=0;i<ng;i=i+1){
      for(m=0;m<nseq;m=m+1){
        alf=seg->alpha[m][i];
        weights[i]=weights[i]*(cnt[m]+alf[ch[m]])/(l+alf[alphsize[m]]);     //CHANGE: use alphsize vector
      }
      ratio+=weights[i];
    }
    for(i=0;i<ng;i=i+1)
      weights[i]=weights[i]/ratio;
    probs[l+1]=ratio;
    for(m=0;m<nseq;m++)
      mtot[m][ch[m]]++; 

    for(m=0;m<nseq;m++){
      ch[m]=term[m] & MASK[m];
      term[m]=term[m]>>NUMBITS[m];

      numch[m]++;
      if(numch[m]==CHPERE[m]){
        k[m]++;
        term[m]=seq[m][k[m]];
        numch[m]=0;
      }
    }

    /* OLD CODE CHANGE: MOVE IN ABOVE LOOP
    numch++;
    if(numch==CHPERE){
      k++;
      for(m=0;m<nseq;m++)
        term[m]=seq[m][k];
      numch=0;
    }
    */
  }

  //Calculate ratios for right part
  R=L+numscores-1;

  for(m=0;m<nseq;m++){
      k[m]=R/CHPERE[m];                 //CHANGE: Bring inside loop over seqs and use CHPERE and k vector
      numch[m]=(unsigned short)(CHPERE[m]-R%CHPERE[m]-1);

      for(i=0;i<alphsize[m];i=i+1)              //CHANGE: use alphsize vector
      mtot[m][i]=0;

      term[m]=seq[m][k[m]]<<(numch[m]*NUMBITS[m]);
  }

  //Do right segment ratio calculations using stored results
  parent=root;
  newcalc=1;
  for(l=0;l<numscores;l=l+1){
    index=0;
    for(m=0;m<nseq;m++){
      ch[m]=term[m] >> (ESIZE[m]-NUMBITS[m]);   //CHANGE: use ESIZE and NUMBITS vector
      index=index*alphsize[m]+ch[m];            //CHANGE: use alphsize vector
      term[m]=term[m]<<NUMBITS[m];

      numch[m]++;
      if(numch[m]==CHPERE[m]){
        k[m]--;
        term[m]=seq[m][k[m]];
        numch[m]=0;
      }
    }

    /* OLD CODE - CHANGE: MOVE INTO ABOVE LOOP
    numch++;
    if(numch==CHPERE){
      k--;
      for(m=0;m<nseq;m++)
        term[m]=seq[m][k];
      numch=0;
    }
    */

    if(!parent->child || !parent->child[index]) break; //End of path in tree

    cur=parent->child[index];
    if(cur->ratio[index]>0){ //Is stored
      if(l!=0) probs[numscores-l]/=cur->ratio[index]; //First one doesn't get used
      newcalc=1;
    }
    else{
      if(newcalc){//Initialise weights
        for(i=0;i<ng;i=i+1)
          weights[i]=parent->weights[i];
        newcalc=0;
      }
      ratio=0;
      for(m=0;m<nseq;m=m+1) cnt[m]=mtot[m][ch[m]];
      for(i=0;i<ng;i=i+1){
        for(m=0;m<nseq;m=m+1){
          alf=seg->alpha[m][i];
          weights[i]=weights[i]*(cnt[m]+alf[ch[m]])/(l+alf[alphsize[m]]);   //CHANGE: use alphsize vector
        }
        ratio+=weights[i];
      }
      for(i=0;i<ng;i=i+1){
        weights[i]=weights[i]/ratio;
        cur->weights[i]=weights[i]; //Store them
      }
      cur->ratio[index]=ratio;
      if(l!=0) probs[numscores-l]/=ratio; //First one doesn't get used
    }
    for(m=0;m<nseq;m++)
      mtot[m][ch[m]]++; 
    parent=cur;
  }

  //Reinitialise weights
  if(newcalc && l<numscores){
    for(i=0;i<ng;i=i+1)
      weights[i]=parent->weights[i];
  }

  //Do rest of right segment ratio calculation without storage
  for(;l<numscores;l=l+1){
    ratio=0;
    for(m=0;m<nseq;m=m+1) cnt[m]=mtot[m][ch[m]];
    for(i=0;i<ng;i=i+1){
      for(m=0;m<nseq;m=m+1){
        alf=seg->alpha[m][i];
        weights[i]=weights[i]*(cnt[m]+alf[ch[m]])/(l+alf[alphsize[m]]);     //CHANGE: use alphsize vector
      }
      ratio+=weights[i];
    }
    for(i=0;i<ng;i=i+1)
      weights[i]=weights[i]/ratio;
    if(l!=0) probs[numscores-l]/=ratio; //First one doesn't get used
    for(m=0;m<nseq;m++)
      mtot[m][ch[m]]++; 

    for(m=0;m<nseq;m++){
      ch[m]=term[m] >> (ESIZE[m]-NUMBITS[m]);
      term[m]=term[m]<<NUMBITS[m];

      numch[m]++;
      if(numch[m]==CHPERE[m]){
        k[m]--;
        term[m]=seq[m][k[m]];
        numch[m]=0;
      }
    }

    /* OLD CODE  - CHANGE: MOVE INTO ABOVE LOOP
    numch++;
    if(numch==CHPERE){
      k--;
      for(m=0;m<nseq;m++)
        term[m]=seq[m][k];
      numch=0;
    }
    */
  }

  //Accumulate
  for(l=2;l<numscores;l=l+1)
    probs[l]*=probs[l-1];

  for(l=2;l<numscores;l=l+1)
    probs[l]+=probs[l-1];
}


unsigned long probabilityvector::choose(unsigned long numscores)
{
  double x;
  unsigned long a,b,c;

  //Choose option - binary search
  a=1;b=numscores-1;
  x=probs[b]*uniform();
  if(probs[a]>=x) return a;
  while(b>a+1){
    c=a+(b-a)/2;
    if(probs[c]<x) a=c;
    else b=c;
  }

  return b;
}


int TryInsertion(pcp cur)
//Assumes uniform prior
{
  pcp newcp;
  pcp next;
  unsigned long L;
  unsigned long i;
  unsigned long numscores;
  
  //Get left and right ends
  next=cur->next; 
  L=cur->c;
  numscores=next->c-L;
  if(numscores==1) return 0;

  pr->getscores(numscores,L,seg->numcp,seg->numcp-seg->numfixed);
  if(pr->probs[numscores-1]<pr->probs[0] && uniform()>pr->probs[numscores-1]/pr->probs[0]){ 
    seg->output();
    return 0;
  }
  else {
    if(seg->numcp==MAXCP){
      printf("Too many change-points. Try increasing MAXCP.\n");fflush(stdout);
      exit(1);
    }
    //Initialise new segment
    newcp=seg->insert(cur);
    i=pr->choose(numscores);
    newcp->c=L+i;
    seg->output();
    i=pr->choose(numscores);
    newcp->c=L+i;
    seg->output();
    return 1;
  }
}


int TryDeletion(pcp cur)
{
  pcp next,mid;
  unsigned long L;
  unsigned long i;
  unsigned long numscores;

  //Get left and right ends
  mid=cur->next;
  next=mid->next; 
  L=cur->c;
  numscores=next->c-L;

  //Choose to delete or not
  if(mid->fixed==0) pr->getscores(numscores,L,seg->numcp-1,seg->numcp-1-seg->numfixed);
  if(mid->fixed==0 && (pr->probs[numscores-1]<pr->probs[0] 
                       || uniform()<pr->probs[0]/pr->probs[numscores-1])){
    seg->remove(cur);
    seg->output();
    return 1;
  }
  else{
    seg->output();
    if(mid->fixed==0){
      i=pr->choose(numscores);
      mid->c=L+i;
    }
    seg->output();
    return 0;
  }
}


unsigned short g1,g2;
unsigned short s1;
double scale;
double gam[MAXGRPS+MAXALPHSIZE]; //max{MAXGRPS,MAXALPHSIZE} would be large enough
double maxx;
double alfsum;
double pifn1(double x)
{
  seg->curpi[g1]=log(x);
  seg->curpi[g2]=log(maxx-x);
  return seg->lnlike()/scale;
}


double tc,K,C1[MAXSEQS],C2[MAXSEQS];
double A[MAXSEQS][MAXALPHSIZE],B[MAXSEQS][MAXALPHSIZE];
double pifn2(double x)
{
  /*
  double y,t1,t2,prod1,prod2,z1[MAXSEQS],z2[MAXSEQS];
  double mu1[MAXSEQS][MAXALPHSIZE],mu2[MAXSEQS][MAXALPHSIZE];
  unsigned short k,m;
  double retval;
  double alpha1[MAXALPHSIZE+1],alpha2[MAXALPHSIZE+1];
  
  y=maxx-x;
  t1=tc-sqrt(y*K/x);
  t2=tc+sqrt(x*K/y);
  retval=1;
  for(m=0;m<nseq;m++){
    prod1=prod2=1;
    for(k=0;k<alphsize;k++){
      mu1[m][k]=A[m][k]*t1+B[m][k];
      mu2[m][k]=A[m][k]*t2+B[m][k];
      prod1*=mu1[m][k];
      prod2*=mu2[m][k];
    }
    z1[m]=pow(prod1/(C1[m]*x*x),1.0/(alphsize-1.0)); //Need -1?
    z2[m]=pow(prod2/(C2[m]*y*y),1.0/(alphsize-1.0)); //Need -1?
    retval*=1/((1+z1[m])*(1+z2[m]));
  }
  seg->curpi[g1]=log(x);
  seg->curpi[g2]=log(y);
  for(m=0;m<nseq;m++){
    for(k=0;k<alphsize;k++){
      alpha1[k]=mu1[m][k]*z1[m];
      alpha2[k]=mu2[m][k]*z2[m];
    }
    alpha1[alphsize]=z1[m];
    alpha2[alphsize]=z2[m];
    seg->setgamrat(m,g1,alpha1);
    seg->setgamrat(m,g2,alpha2);
  }
  return 0.5*log(retval)-(nseq*alphsize-1)*log(x*y)/2.0+seg->lnlike()/scale;  //Probably not correct for nseq>1
  //Need to check first two terms. Should second term have alphsize instead of nseq*alphsize?
  */
}


double pifn3(double x)
{
  double temp;
  unsigned short k;
  
  temp=log(1-x);
  for(k=0;k<ng;k++)
    seg->curpi[k]=gam[k]+temp;
  seg->curpi[g1]=log(x);
  return (ng-2)*temp+seg->lnlike()/scale; //Note: assumes ng>=2
}


double zfn(double x)
{
  unsigned short k;
  double z;
  double alpha[MAXALPHSIZE+1];

  z=1/(x*x)-1;
  for(k=0;k<alphsize[s1];k++)               //CHANGE: use alphsize vector
    alpha[k]=gam[k]*z;
  alpha[alphsize[s1]]=z;                    //CHANGE: use alphsize vector
  seg->setgamrat(s1,g1,alpha);  
  return /*(ac-1)*log(z)-ad*z-x*x*x+*/seg->lnlike(g1)/scale;  
}


//Use g1 as group number, g2 as char number
double mufn(double x)
{
  unsigned short k;
  double y;
  double alpha[MAXALPHSIZE+1];
  
  y=1-x;
  alpha[g2]=x*alfsum;
  for(k=0;k<alphsize[s1];k++){                  //CHANGE: use alphsize vector
    if(k!=g2) alpha[k]=gam[k]*y*alfsum;
  }
  alpha[alphsize[s1]]=alfsum;                   //CHANGE: use alphsize vector
  seg->setgamrat(s1,g1,alpha);                  //CHANGE: use alphsize vector
  return (alphsize[s1]-2)*log(y)+seg->lnlike(g1)/scale;
}


void slicebnd(double (*f)(double x),double L,double R,double *xnew,double *fxnew)
//Slice sampler for distribution (*f) on (L,R)
//At input *xnew is current position and *fxnew=(*f)(*xnew). Assumes *xnew>0
//Assumes (*f)(x) finite for all x>0 and is the log of a distribution up to normalisation constant
//At input, xnew should be strictly between L and R
{
  int i;
  double y,x1,fx1;
  
  y=-rand_Exponential()+*fxnew; //Assumes *f returns the log of the distribution and fx=(*f)(x)
                  //Only place log is assumed
  
  //Sample uniformly
  for(i=0;i<250;i++){
    x1=L+uniform()*(R-L);
    fx1=(*f)(x1);
    if(y<=fx1){*xnew=x1;*fxnew=fx1;return;}
    if(x1<*xnew) L=x1;
    else R=x1;
  }
  printf("L=%lf R=%lf x1=%lf fx1=%lf\n xnew=%lf fxnew=%lf y=%lf\n",L,R,x1,fx1,*xnew,*fxnew,y);
  printf("Maximum iterations exceeded in slice.\n");
  fflush(stdout);
  exit(1);
}


void GibbsIter()
{
  unsigned short j,k,l,m;
  pcp cur;
  double x,y,gg;
  double *curalpha[MAXSEQS];
  double sum[MAXSEQS][MAXGRPS];
  unsigned short chn;
  double max;

  /* unused variables when not using chains
  double *pdarc;
  double oldlnlike,newlnlike;
  double mu1[MAXSEQS][MAXALPHSIZE],mu2[MAXSEQS][MAXALPHSIZE];
  double prod1,prod2,L,R;
  double t1,t2;
  double temp;
  */


  pr->reset();
  
  //Do insertions/deletions/slides
  cur=seg->list;
  do {
	//Try insertion
	while(TryInsertion(cur))
      cur=cur->next;
	
	//Try deletion/substitution
	while(cur->next->next && TryDeletion(cur));
	cur=cur->next;
  }
  while(cur->next!=NULL);
  
  //Update hyperparameters
  seg->setcountsandtotals();
  for(chn=0;chn<numchains;chn++){
    seg->setcurpiandgamrat(chains[chn]);

    unsigned short curparamspos=0;      //CHANGE: use position marker for alpha
    for(m=0;m<nseq;m++){
      //OLD CODE: curalpha[m]=chains[chn]+ng+m*ng*alphsize;
      curalpha[m]=chains[chn]+ng+curparamspos;
      curparamspos+=ng*alphsize[m];
    }

    //sum alphas  
    for(m=0;m<nseq;m++){
      for(j=0;j<ng;j++){
        sum[m][j]=0;
        for(k=0;k<alphsize[m];k++)      //CHANGE: use alphsize vector
          sum[m][j]+=curalpha[m][j*alphsize[m]+k];
      }
    }

    scale=temperatures[chn];

    if(pair && ng>1){
      printf("No support for \"pair\" in multi-seq\n");
      fflush(stdout);
      exit(1);
        /*
        //Update pis
      for(j=0;j<ng-1;j++){
        g1=j;
        x=exp(seg->curpi[j]);
        for(k=j+1;k<ng;k++){
          g2=k;
          y=exp(seg->curpi[k]);
          maxx=x+y; 
          gg=pifn1(x);
          slicebnd(&pifn1,FEPS,maxx-FEPS,&x,&gg); 
          if(chn==0) {
            seg->pi[j]=seg->curpi[j];
            seg->pi[k]=seg->curpi[k];
          }
          seg->output();

          //Not currently using this move, but I need to modify it so that the alphas stay within bounds
          //pis should stay in bounds
          //I also need to check the Jacobian for more than one sequence
          y=maxx-x;
          tc=y/maxx;
          K=(x/maxx)*(y/maxx);
          for(m=0;m<nseq;m++){
            prod1=prod2=1;
            for(l=0;l<alphsize;l++){
              B[m][l]=curalpha[m][j*alphsize+l]/sum[m][j];
              A[m][l]=curalpha[m][k*alphsize+l]/sum[m][k];
              prod1*=B[m][l];
  		      prod2*=A[m][l];
              A[m][l]-=B[m][l];
            }
  		    C1[m]=prod1/(x*x*pow(sum[m][j],alphsize-1)); //Need +1?
  		    C2[m]=prod2/(y*y*pow(sum[m][k],alphsize-1)); //Need +1?  
          }
          gg=pifn2(x);
          //Compute upper and lower limits 
          L=0;R=1;
          for(m=0;m<nseq;m++){
            for(l=0;l<alphsize;l++){
              if(A[m][l]>0){
                temp=-B[m][l]/A[m][l];
                if(L==0 || temp>L) L=temp;
                temp=(1-B[m][l])/A[m][l];
                if(R==1 || temp<R) R=temp;
              }
              else if(A[m][l]<0){
                temp=(1-B[m][l])/A[m][l];
                if(L==0 || temp>L) L=temp;
                temp=-B[m][l]/A[m][l];
                if(R==1 || temp<R) R=temp;
              }
            }
          }
          L=maxx*K/(K+(tc-L)*(tc-L));
          if(L<FEPS) L=FEPS;
          R=maxx*(R-tc)*(R-tc)/(K+(R-tc)*(R-tc));
          if(R>maxx-FEPS) R=maxx-FEPS;
          if(x<L || x>R){printf("x out of bounds in GibbsIter.\n");exit(1);} //Impossible?
          slicebnd(&pifn2,L,R,&x,&gg); 
          y=maxx-x;
          t1=tc-sqrt(y*K/x);
          t2=tc+sqrt(x*K/y);
          for(m=0;m<nseq;m++){
            prod1=prod2=1;
            for(l=0;l<alphsize;l++){
              mu1[m][l]=A[m][l]*t1+B[m][l];
              mu2[m][l]=A[m][l]*t2+B[m][l];
              prod1*=mu1[m][l];
              prod2*=mu2[m][l];
            }
            sum[m][j]=pow(prod1/(C1[m]*x*x),1.0/(alphsize-1.0)); //Need -1?
            sum[m][k]=pow(prod2/(C2[m]*y*y),1.0/(alphsize-1.0)); //Need -1?
            for(l=0;l<alphsize;l++){
              curalpha[m][j*alphsize+l]=mu1[m][l]*sum[m][j];
              curalpha[m][k*alphsize+l]=mu2[m][l]*sum[m][k];
            } 
          }
          if(chn==0) {
            seg->pi[j]=seg->curpi[j];
            seg->pi[k]=seg->curpi[k];
            for(m=0;m<nseq;m++){
              for(l=0;l<alphsize;l++){
                seg->alpha[m][j][l]=curalpha[m][j*alphsize+l];
                seg->alpha[m][k][l]=curalpha[m][k*alphsize+l];
              }
              seg->alpha[m][j][alphsize]=sum[m][j];
              seg->alpha[m][k][alphsize]=sum[m][k];
            }
          }
          seg->output();
        }
      }
        */
    }
    
    for(j=0;j<ng;j++){
      g1=j;

      if(ng>1){
        //Update pis
        x=exp(seg->curpi[j]);
        max=1;
        for(k=0;k<ng;k++){
          if(k!=j && 1-FEPS*(1-x)/exp(seg->curpi[k])<max) max=1-FEPS*(1-x)/exp(seg->curpi[k]);
        }
        y=log(1-x);
        for(k=0;k<ng;k++)
          gam[k]=seg->curpi[k]-y;
        gg=pifn3(x);
        //printf("pifn3: chn=%hu j=%hu x=%lf gg=%lf\n",chn,j,x,gg);
        slicebnd(&pifn3,FEPS,max,&x,&gg); 
        if(chn==0) {
          for(k=0;k<ng;k++)
            seg->pi[k]=seg->curpi[k];
        }
        seg->output();

        seg->setpart(g1);
      }

      for(m=0;m<nseq;m++){
        s1=m;

        //Update sum[m][j]
        if(sum[m][j]==0){printf("Zero sum in GibbsIter.\n");fflush(stdout);exit(1);}
        for(k=0;k<alphsize[m];k++)                                                 //CHANGE: use alphsize vector
          gam[k]=curalpha[m][j*alphsize[m]+k]/sum[m][j];                            //CHANGE: use alphsize vector
        x=sqrt(1/(1+sum[m][j]));
        gg=zfn(x);
        //printf("zfn: chn=%hu j=%hu m=%hu x=%lf gg=%lf\n",chn,j,m,x,gg);
        //for(k=0;k<alphsize;k++)
        //  printf("%lf ",gam[k]);
        //printf("\n");
        slicebnd(&zfn,FEPS,1-FEPS,&x,&gg); 
        sum[m][j]=1/(x*x)-1;
        for(k=0;k<alphsize[m];k++)                                             //CHANGE: use alphsize vector
          curalpha[m][j*alphsize[m]+k]=gam[k]*sum[m][j];                        //CHANGE: use alphsize vector
        if(chn==0){
          for(k=0;k<alphsize[m];k++)                                           //CHANGE: use alphsize vector
            seg->alpha[m][j][k]=curalpha[m][j*alphsize[m]+k];                   //CHANGE: use alphsize vector
          seg->alpha[m][j][alphsize[m]]=sum[m][j];                              //CHANGE: use alphsize vector
        }
        seg->output();

        //Update alphas
        alfsum=sum[m][j];
        if(alfsum==0){printf("Zero alfsum in GibbsIter.\n");fflush(stdout);exit(1);}
        for(k=0;k<alphsize[m];k++){                                     //CHANGE: use alphsize vector
          g2=k;
          x=curalpha[m][j*alphsize[m]+k]/alfsum;                        //CHANGE: use alphsize vector
          y=1-x;

          max=1;
          for(l=0;l<alphsize[m];l++){                                   //CHANGE: use alphsize vector
            if(l!=k && 1-FEPS*y*alfsum/curalpha[m][j*alphsize[m]+l]<max) max=1-FEPS*y*alfsum/curalpha[m][j*alphsize[m]+l];  //CHANGE: use alphsize vector
          }

          for(l=0;l<alphsize[m];l++)                                   //CHANGE: use alphsize vector
            gam[l]=(curalpha[m][j*alphsize[m]+l]/alfsum)/y;            //CHANGE: use alphsize vector
          gg=mufn(x);
          //printf("mufn: chn=%hu j=%hu m=%hu k=%hu x=%lf gg=%lf\n",chn,j,m,k,x,gg);
          slicebnd(&mufn,FEPS,max,&x,&gg); 
          y=1-x;
          for(l=0;l<alphsize[m];l++)                                //CHANGE: use alphsize vector
            curalpha[m][j*alphsize[m]+l]=(gam[l]*y)*alfsum;         //CHANGE: use alphsize vector
          curalpha[m][j*alphsize[m]+k]=x*alfsum;                    //CHANGE: use alphsize vector
          if(chn==0){
            for(l=0;l<alphsize[m];l++)
              seg->alpha[m][j][l]=curalpha[m][j*alphsize[m]+l];     //CHANGE: use alphsize vector
          }
          seg->output();
        }
      }
    }
  }


  if(numchains>1){
      printf("Multiple chains not supported for multi-seq\n");
      fflush(stdout);
      exit(1);
  }
  //Do swaps
  /* COMMENT OUT FOR CHAINS
  seg->setcurpiandgamrat(chains[0]);
  oldlnlike=seg->lnlike();  
  for(chn=0;chn<numchains-1;chn++){
    seg->setcurpiandgamrat(chains[chn+1]);
    newlnlike=seg->lnlike();
    temp=(1.0/temperatures[chn]-1.0/temperatures[chn+1])*(newlnlike-oldlnlike);
    oldlnlike=newlnlike;
    if(temp>=1 || -rand_Exponential()<temp){
      pdarc=chains[chn];chains[chn]=chains[chn+1];chains[chn+1]=pdarc;
      if(SCREEN) printf("Swapping %hu - %hu: %lf - %lf.\n",chn,chn+1,oldlnlike,newlnlike);
      if(LOG) fprintf(logfile,"Swapping %hu - %hu: %lf - %lf.\n",chn,chn+1,oldlnlike,newlnlike);
    }
    if(chn==0){
      for(m=0;m<nseq;m++)
        curalpha[m]=chains[0]+ng+m*ng*alphsize;
      for(j=0;j<ng;j++){
        seg->pi[j]=chains[0][j];
        for(m=0;m<nseq;m++){
          sum[m][j]=0;
          for(k=0;k<alphsize;k++){
            seg->alpha[m][j][k]=curalpha[m][j*alphsize+k];
            sum[m][j]+=curalpha[m][j*alphsize+k];
          }
          seg->alpha[m][j][alphsize]=sum[m][j];
        }
      }
    }
    seg->output();
  }
  */
}


int compchar(const void*e1,const void *e2)
{
  char *p1,*p2;

  p1=(char*)e1;
  p2=(char*)e2;
  return (int)(*p1-*p2);
}


void getalphabet(char **filename,unsigned long *len) //CHANGE: add layer to filename pointer and remove alphabet and alphsize as arguments                                                        
{

  for(unsigned short seqnum=0;seqnum<nseq;seqnum++){           //CHANGE: add loop over multiple seqs with counter seqnum
      FILE *in;
      char ch;
      char tempalphabet[MAXALPHSIZE+1];


      in=fopen(filename[seqnum],"r");               //CHANGE: get the file for the given seq
      if(!in){printf("Could not open file %s in getalphabet.\n",filename[seqnum]);fflush(stdout);exit(1);}

      alphsize[seqnum]=0;                           //CHANGE: changed the whole function using new alphsize pointer
      tempalphabet[0]=0;
      *len=0;
      while(fscanf(in,"%c",&ch)==1){
        if((ch>='0' && ch<='9')||(ch>='a' && ch<='z')||(ch>='A' && ch<='Z' && ch!='I' && ch!='J' && ch!='K' && ch!='L' && ch!='M' && ch!='N' && ch!='O')){
          if(*len==MAXLEN){
            printf("Sequence too long in getalphabet. Truncating.\n");
            fflush(stdout);
            break;
          }
          (*len)++;
          if(!index(tempalphabet,ch)){
            if(alphsize[seqnum]==MAXALPHSIZE){printf("Alphabet too large in getalphabet.\n");fflush(stdout);exit(1);}
            tempalphabet[alphsize[seqnum]]=ch;
            alphsize[seqnum]++;
            tempalphabet[alphsize[seqnum]]=0;
          }
        }
      }
      fclose(in);

      qsort(tempalphabet,alphsize[seqnum],sizeof(char),&compchar);

      if(*len==0){printf("Zero length segment in getalphabet.\n");fflush(stdout);exit(1);}
      if(alphsize[seqnum]<2){printf("Less than two characters in alphabet in getalphabet.\n");fflush(stdout);exit(1);}

      alphabet[seqnum]=new char[alphsize[seqnum]+1];
      strcpy(alphabet[seqnum],tempalphabet);

      unsigned short tempalphsize = alphsize[seqnum];
      NUMBITS[seqnum]=(unsigned short)ceil(log(log((double)tempalphsize)/log(2.0))/log(2.0));  //CHANGE: use vectors of variables
      NUMBITS[seqnum]=(unsigned short)pow(2.0,NUMBITS[seqnum]);
      MASK[seqnum]=(unsigned short)pow(2.0,NUMBITS[seqnum])-1;
      ESIZE[seqnum]=8*sizeof(unsigned short);
      CHPERE[seqnum]=ESIZE[seqnum]/NUMBITS[seqnum];
  }
}


void readseq(char **inputfilename,unsigned short *nseq,unsigned short ***seq)
//seg must already exist
{
  FILE *in;
  char ch;
  unsigned short shch,term,bits;
  unsigned long i;
  pcp cur,newcp;
  unsigned long templen;
  unsigned short j;


  //Allocate
  *seq=new unsigned short*[*nseq];
  if(!(*seq)){printf("Could not allocate seq in main.\n");fflush(stdout);exit(1);}
  for(j=0;j<*nseq;j++){
    if(!((*seq)[j]=new unsigned short[len/CHPERE[j]+1])){               //CHANGE: use jth reference in vector CHPERE
      printf("Memory allocation problem in readseq.\n");
      fflush(stdout);
      exit(1);
    }
  }

  //Read first sequence and get fixed changepoints
  //Open file
  if(!(in=fopen(inputfilename[0],"r"))){
    printf("Could not open file %s.\n",inputfilename[0]);
    exit(1);
  }

  //Read in
  len=0;
  i=0;bits=0;term=0;
  cur=seg->list;
  while(fscanf(in,"%c",&ch)==1){
    if(ch=='#'){
      newcp=seg->insert(cur);
      newcp->c=len;
      newcp->fixed=1;
      cur=newcp;
      seg->numfixed++;
    }
    else if((ch>='0' && ch<='9')||(ch>='a' && ch<='z')||(ch>='A' && ch<='Z' && ch!='I' && ch!='J' && ch!='K' && ch!='L' && ch!='M' && ch!='N' && ch!='O')){
      if(len==MAXLEN){
        printf("Sequence too long in readseq. Truncating.\n");
        fflush(stdout);
        break;
      }
      for(shch=0;shch<alphsize[0];shch++){      //CHANGE: using alphsize vector ref
        if(alphabet[0][shch]==ch) break;        //CHANGE: using alphabet matrix
      }
      if(shch==alphsize[0]){printf("Unrecognized character in readseq.\n");fflush(stdout);exit(1);}        //CHANGE: using alphsize vector
      term=term | (shch<<bits);
      bits+=NUMBITS[0];                                                                 //CHANGE: using numbits vector
      if(bits==ESIZE[0]) {(*seq)[0][i++]=term;bits=0;term=0;}                           //CHANGE: using esize vector
      len++;
    }
  }
  fclose(in);

  if(bits>0) (*seq)[0][i]=term;
  
  cur->next->c=len;
  cur->next->fixed=1;
  
  if(SCREEN) printf("Sequence length=%lu.\n",len);
  fflush(stdout);
  if(LOG) fprintf(logfile,"Sequence length=%lu.\n",len);

  //Check there are no zero length segments
  cur=seg->list;
  for(i=0;i<seg->numcp;i++){
    if(cur->next->c == cur->c){
      printf("Zero segment length in readseq.\n");
      fflush(stdout);
      exit(1);
    }
    cur=cur->next;
  }

  //Read in remaining sequences
  for(j=1;j<*nseq;j++){
    //Open file
    if(!(in=fopen(inputfilename[j],"r"))){
      printf("Could not open file %s.\n",inputfilename[j]);
      exit(1);
    }

    //Read in
    templen=0;
    i=0;bits=0;term=0;
    while(fscanf(in,"%c",&ch)==1){
      if((ch>='0' && ch<='9')||(ch>='a' && ch<='z')||(ch>='A' && ch<='Z' && ch!='I' && ch!='J' && ch!='K' && ch!='L' && ch!='M' && ch!='N' && ch!='O')){
        if(templen==MAXLEN){
          printf("Sequence too long in readseq. Truncating.\n");
          fflush(stdout);
          break;
        }
        for(shch=0;shch<alphsize[j];shch++){            //CHANGE: using alphsize vector
          if(alphabet[j][shch]==ch) break;                 //CHANGE: using alphabet matrix
        }
        if(shch==alphsize[j]){printf("Unrecognized character in readseq.\n");fflush(stdout);exit(1);}      //CHANGE: using alphsize vector
        term=term | (shch<<bits);
        bits+=NUMBITS[j];                   //CHANGE: using numbits vector
        if(bits==ESIZE[j]) {(*seq)[j][i++]=term;bits=0;term=0;}     //CHANGE: using esize vector
        (templen)++;
      }
    }
    fclose(in);

    if(bits>0) (*seq)[j][i]=term;
    if(templen!=len){printf("Unequal sequence lengths in readseq.\n");fflush(stdout);exit(1);}
  }
}


void usage()
{
  printf("\n");
  printf("To correctly use this program, you must have the following parameters:\n");
  printf("  -i input file names (may segment up to %d sequences simultaneously)\n",MAXSEQS);
  printf("  -sf initial segmentation file (optional)\n");
  printf("  -o output_file_name\n");
  printf("  -n num_samples (integer, >=1)\n");
  printf("  -b num_burn (optional, integer, <num_iterations, default 0)\n");
  printf("  -s sampling_block_size (optional, integer, >=1, default 1)\n");
  printf("  -p prop_changepts_init (optional, double, 0<p<1, default 0.01)\n");
  printf("  -pa hyperparameter for phi (optional, double, pa>0, default 1.0)\n");
  printf("  -pb hyperparameter for phi (optional, double, pb>0, default 1.0)\n");
  printf("  -phi fixed value for phi (optional, double, 0<phi<1, overrides pa and pb settings)\n");
  printf("  -ng number of groups (optional, integer, 0<ng<=100, default 1)\n");
  // printf("  -nc number of chains (optional, integer, 1<=nc<=10, default 1)\n");
  printf("  -hp heating parameter for AP (optional, double, hp>=0, default 1.0)\n");
  printf("  -r random number seed (optional, integer, default=time)\n");
  printf("  -pf samples_per_file (optional, integer, >=1, default num_samples)\n");
  printf("  -nf number_of_output_files (optional, integer, >=1, default num_samples/samples_per_file)\n");
  // printf("  -pair use paired group moves (optional, flag, default: don't use paired group moves)\n");
  printf("  -t tree storage structure file name (optional, <99 char, default: don't use tree storage)\n"); 
  fflush(stdout);
  exit(1);
}


void handleargs(int argc,char **argv)
{
  int i,j;

  nseq=0;
  outputfilename=NULL;
  segfilename=NULL;
  treefilename=NULL;
  NUMSAMP=0;
  NUMBURN=0;
  SAMPBLOCK=1;
  PF=0;
  NF=0;
  p=0.01;
  pa=1.0;
  pb=1.0;
  phi=0;
  ng=1;
  pair=0;
  numchains=1;
  heat=1.0;
  seed=(unsigned long)time(NULL);
  if(argc==1) usage();
  for(i=1;i<argc;i+=2){
    if(!strcmp(argv[i],"-i")){
      while(argv[i+nseq+1][0]!='-'){
        if(nseq==MAXSEQS){printf("Too many input sequences.\n");fflush(stdout);exit(1);}
        nseq++;
      }
      inputfilename=new char*[nseq];
      if(!inputfilename){printf("Could not allocate inputfilename in handleargs.\n");fflush(stdout);exit(1);}
      for(j=0;j<nseq;j++) inputfilename[j]=argv[i+j+1];
      i+=nseq-1;
    } 
    else if(!strcmp(argv[i],"-sf")) 
      segfilename=argv[i+1]; 
    else if(!strcmp(argv[i],"-o")) {
      outputfilename=argv[i+1];
      sprintf(logfilename,"%s.log",outputfilename);
    }
    else if(!strcmp(argv[i],"-n")){
      if(sscanf(argv[i+1],"%lu",&NUMSAMP)!=1) usage();
    }
    else if(!strcmp(argv[i],"-b")){
      if(sscanf(argv[i+1],"%lu",&NUMBURN)!=1) usage();
    }
    else if(!strcmp(argv[i],"-s")){
      if(sscanf(argv[i+1],"%lu",&SAMPBLOCK)!=1) usage();
      if(SAMPBLOCK<1) usage();
    }
    else if(!strcmp(argv[i],"-p")){
      if(sscanf(argv[i+1],"%lf",&p)!=1) usage();
	  if(p<=0 || p>=1) usage();
    }
    else if(!strcmp(argv[i],"-pa")){
      if(sscanf(argv[i+1],"%lf",&pa)!=1) usage();
	  if(pa<=0) usage();
    }
    else if(!strcmp(argv[i],"-pb")){
      if(sscanf(argv[i+1],"%lf",&pb)!=1) usage();
	  if(pb<=0) usage();
    }
    else if(!strcmp(argv[i],"-phi")){
      if(sscanf(argv[i+1],"%lf",&phi)!=1) usage();
	  if(phi<=0 || phi>=1) usage();
    }
    else if(!strcmp(argv[i],"-ng")){
      if(sscanf(argv[i+1],"%hu",&ng)!=1) usage();
	  if(ng<1 || ng>MAXGRPS) usage();
    }
    /*
    else if(!strcmp(argv[i],"-nc")){
      if(sscanf(argv[i+1],"%hu",&numchains)!=1) usage();
	  if(numchains<1 || numchains>10) usage();
    }
    */
    else if(!strcmp(argv[i],"-hp")){
      if(sscanf(argv[i+1],"%lf",&heat)!=1) usage();
	  if(heat<0) usage();
    }
    else if(!strcmp(argv[i],"-r")){
      if(sscanf(argv[i+1],"%lu",&seed)!=1) usage();
    }
    else if(!strcmp(argv[i],"-pf")){
      if(sscanf(argv[i+1],"%lu",&PF)!=1) usage();
      if(PF<1) usage();
    }
    else if(!strcmp(argv[i],"-nf")){
      if(sscanf(argv[i+1],"%lu",&NF)!=1) usage();
      if(NF<1) usage();
    }
    /*
    else if(!strcmp(argv[i],"-pair")){
      pair=1;
      i--;
    }
    */
    else if(!strcmp(argv[i],"-t"))
      treefilename=argv[i+1]; 
    else usage();
  }
  
  if(nseq==0) usage();
  if(outputfilename==NULL) usage();
  if(NUMSAMP==0) usage();
  if(PF==0) PF=NUMSAMP;
  if(NF==0) NF=(unsigned long)ceil(1.0*NUMSAMP/PF);
}


int main(int argc,char **argv)
{
  clock_t t1,t2;
  unsigned short i;

  if(SCREEN){
    for(i=0;i<argc;i++)
      printf("%s ",argv[i]);
    printf("\n");
    fflush(stdout);
  }
  
  handleargs(argc,argv);

  if(LOG){
    if(!(logfile=fopen(logfilename,"a"))){
      printf("Could not open file %s.\n",logfilename);
      fflush(stdout);
      exit(1);
    }

    for(i=0;i<argc;i++)
      fprintf(logfile,"%s ",argv[i]);
    fprintf(logfile,"\n");
  }
  
  initlfunc2();

  //Seed random number generator
  init_genrand(seed);
  if(SCREEN) printf("Random number seed=%lu\n",seed);
  fflush(stdout);
  if(LOG) fprintf(logfile,"Random number seed=%lu\n",seed);

  
  //CHANGE: allocate alphabet and alphsize
  alphsize = new unsigned short[nseq];
  if(!alphsize){printf("Could not allocate alphsize in main.\n");fflush(stdout);exit(1);}
  alphabet = new char*[nseq];
  if(!alphabet){printf("Could not allocate alphabet in main.\n");fflush(stdout);exit(1);}
  
  //First determine alphabet and sequence length
  getalphabet(inputfilename,&len); //CHANGE: initial argument is pointer to all file names
                                                                //remove alphabet and alphsize as arguments
  
  //CHANGE: printing alphabet for multiple seqs
  printf("alphabet=");
  if(LOG){
    fprintf(logfile,"alphabet=");
  }
  for(int seqnum=0;seqnum<nseq;seqnum++){
      if(seqnum==0){
          printf("%s",alphabet[seqnum]);
          if(LOG){
            fprintf(logfile,"%s",alphabet[seqnum]);
          }
      }
      else{
          printf(" %s",alphabet[seqnum]);
          if(LOG){
            fprintf(logfile," %s",alphabet[seqnum]);
          }
      }
  }
  printf("\n");
  fflush(stdout);
  if(LOG){
      fprintf(logfile,"\n");
  }
         
  printf("alphsize=");
  if(LOG){
    fprintf(logfile,"alphsize=");
  }
  for(int seqnum=0;seqnum<nseq;seqnum++){
      if(seqnum==0){
          printf("%hu",alphsize[seqnum]);
          if(LOG){
            fprintf(logfile,"%hu",alphsize[seqnum]);
          }
      }
      else{
          printf(" %hu",alphsize[seqnum]);
          if(LOG){
            fprintf(logfile," %hu",alphsize[seqnum]);
          }
      }
  }
  printf("\n");
  fflush(stdout);
  if(LOG){
      fprintf(logfile,"\n");
  }
  
  //CHANGE: add sumalphsize calculation
  sumalphsize=0;
  for(int seqnum=0;seqnum<nseq;seqnum++){
      sumalphsize+=alphsize[seqnum];
  }

  seg=new parameters;
  
  //Input data
  if(SCREEN) printf("Reading sequence.\n");
  fflush(stdout);
  if(LOG) fprintf(logfile,"Reading sequence.\n");
  readseq(inputfilename,&nseq,&seq);
  for(int i=0;i<nseq;i++) 
      delete [] alphabet[i];        //CHANGE: delete rows in alphabet matrix
  delete [] alphabet;
  if(SCREEN) printf("numfixed=%lu\n",seg->numfixed);
  fflush(stdout);
  if(LOG) fprintf(logfile,"numfixed=%lu\n",seg->numfixed);
  
  if(pair && nseq>1) printf("WARNING: I'm not sure that pifn2 is correct for nseq>1. Needs checking.");

  effalphsize=1;
  for(i=0;i<nseq;i++) effalphsize*=alphsize[i];     //CHANGE: effalphsize multiplied by the actual alfsize for each seq using the vector

  //Allocate
  chains=new double*[numchains];
  if(!chains){printf("Could not allocate chains in main.\n");fflush(stdout);exit(1);}
  for(i=0;i<numchains;i++){
    chains[i]=new double[ng*(sumalphsize+1)]; //CHANGE: use sumalphsize
    if(!chains[i]){printf("Could not allocate chains in main.\n");fflush(stdout);exit(1);}
  }
  for(i=0;i<numchains;i++)
    temperatures[i]=1.0+i*heat;
  
  seg->init();
  
  pr=new probabilityvector;
 
  //Markov chain Monte Carlo
  if(SCREEN) printf("Beginning MCMC.\n");
  fflush(stdout);
  if(LOG) fprintf(logfile,"Beginning MCMC.\n");
  t1=clock();
  burnt=0;
  cnt=0;
  element=(unsigned long)floor(uniform()*SAMPBLOCK);
  sampsize=0;
  while(sampsize<NUMSAMP)
    GibbsIter();
  t2=clock();
  if(SCREEN) printf("total time=%.6lf\n",(double)(t2-t1)/CLOCKS_PER_SEC);
  fflush(stdout);
  if(LOG) fprintf(logfile,"total time=%.6lf\n",(double)(t2-t1)/CLOCKS_PER_SEC);

  //Clean Up
  if(LOG) fclose(logfile);
  for(i=0;i<nseq;i++)
    delete [] seq[i];
  delete [] seq;
  delete seg;
  delete pr;
  for(i=0;i<numchains;i++)
    delete [] chains[i];
  delete [] chains;
  cleanuplfunc2();
  delete [] inputfilename;

  return 0;
}
