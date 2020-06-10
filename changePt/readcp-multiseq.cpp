//From readcp4 - modified to allow for output from changept-multiseq
//Reads output of change point and produces profile of sequence
//Output histogram of theta values
//From readcp3 - uses averages for pi and theta instead of throwing values
//Recently changed so that no longer assumes alphsize==2. However, no longer sorts groups by mean.
//Recently modified to allow larger character set

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<time.h>
#include "keytree3.h"
#include "MTrandom.h"
#include "lfunc.h"


#define MAXLEN 300000000 //Max sequence length
#define MAXCP 3000000
#define MAXSLEN 500000
#define MAXGRPS 100
#define MAXALPHSIZE 40

unsigned short ESIZE[10],NUMBITS[10],CHPERE[10],MASK[10];

char seqfilename[10][500],segfilename[150][500],pfilename[500],thetafilename[500],cpfilename[500],lenfilename[500]; //CHANGE: added up to 10 seqfiles
unsigned long NUMBURN,NUMSKIP;
char **alphabet;    //CHANGE: added level of pointer
unsigned short ng,alphsize[10];    //CHANGE: make alphsize an array
unsigned short sumalphsize;     //CHANGE: added numseqs
unsigned short numsegfiles;
unsigned short numseqs;         //CHANGE: added numseqs
unsigned short pg[MAXGRPS],numpg;
char dotheta,dop;

unsigned short **seq; //CHANGE: added level of pointer
unsigned long *len;   //CHANGE: added level of pointer
unsigned long numcp;
double pi[MAXGRPS],***alpha; //CHANGE: added level of pointer tp alpha


char *index(char *str,const char ch)
{
  while(*str != 0){
    if(ch == *str) return str;
    str++;
  }

  return NULL;
}


unsigned long getsegmentation(FILE *in,unsigned long *numcp,unsigned long *cp,double *pi,double ***alpha)
{
  unsigned long i,j;
  unsigned long num;

  if(fscanf(in,"%lu. %lu ",&num,numcp)!=2){*numcp=0;return 0;}
  
  if(*numcp>MAXCP){
    printf("Too many change-points in segmentation.\n");
    fflush(stdout);
    exit(1);
  }
  
  for(i=0;i<ng;i++){
    if(fscanf(in,"%lf ",pi+i)!=1){*numcp=0;return 0;}
    pi[i]=log(pi[i]);
  }

  for(int seqnum=0;seqnum<numseqs;seqnum++){
      for(i=0;i<ng;i++){
        for(j=0;j<alphsize[seqnum];j++){
          if(fscanf(in,"%lf ",&(alpha[seqnum][i][j]))!=1){*numcp=0;return 0;}
        }
      }
  }

  for(i=0;i<*numcp+2;i++){
    if(fscanf(in,"%lu ",cp+i)!=1){
      printf("Invalid file format in getsegmentation.\n");
      printf("numcp=%lu i=%lu\n",*numcp,i);
      if(i!=0) printf("prevcp=%lu\n",cp[i-1]);
      fflush(stdout);
      exit(1);
    }
  }

  for(i=1;i<*numcp+2;i++)
    cp[i]+=cp[i-1];

  return num;
}


int compchar(const void*e1,const void *e2)
{
  char *p1,*p2;

  p1=(char*)e1;
  p2=(char*)e2;
  return (int)(*p1-*p2);
}


void getalphabet(int seqnum)
{
  FILE *in;
  char ch;
  char tempalphabet[MAXALPHSIZE+1];

  in=fopen(seqfilename[seqnum],"r");
  if(!in){printf("Could not open file %s in getalphabet.\n",seqfilename[seqnum]);fflush(stdout);exit(1);}

  alphsize[seqnum]=0;
  tempalphabet[0]=0;
  len[seqnum]=0;
  while(fscanf(in,"%c",&ch)==1){
    if((ch>='0' && ch<='9')||(ch>='a' && ch<='z')||(ch>='A' && ch<='Z' && ch!='I' && ch!='J' && ch!='K' && ch!='L' && ch!='M' && ch!='N' && ch!='O' )){
      if(len[seqnum]==MAXLEN){
        printf("Sequence too long in getalphabet. Truncating.\n");
        fflush(stdout);
        break;
      }
      len[seqnum]++;
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

  if(len[seqnum]==0){printf("Zero length segment in getalphabet.\n");fflush(stdout);exit(1);}
  if(*alphsize<2){printf("Less than two characters in alphabet in getalphabet.\n");fflush(stdout);exit(1);}

  alphabet[seqnum]=new char[alphsize[seqnum]+1];
  strcpy(alphabet[seqnum],tempalphabet);

  NUMBITS[seqnum]=(unsigned short)ceil(log(log((double)alphsize[seqnum])/log(2.0))/log(2.0));
  NUMBITS[seqnum]=(unsigned short)pow(2.0,NUMBITS[seqnum]);
  MASK[seqnum]=(unsigned short)pow(2.0,NUMBITS[seqnum])-1;
  ESIZE[seqnum]=8*sizeof(unsigned short);
  CHPERE[seqnum]=ESIZE[seqnum]/NUMBITS[seqnum];
}


//void readseq(char *inputfilename,unsigned short **seq,unsigned long *len)
void readseq(int seqnum)
{
  FILE *in;
  char ch;
  unsigned short shch,term,bits;
  unsigned long i;

  //Open file
  if(!(in=fopen(seqfilename[seqnum],"r"))){
    printf("Could not open file %s.\n",seqfilename[seqnum]);
    fflush(stdout);
    exit(1);
  }

  //Allocate
  if(!(seq[seqnum]=(unsigned short*)malloc(((len[seqnum])/CHPERE[seqnum]+1)*sizeof(unsigned short)))){
    printf("Memory allocation problem in readline.\n");
    fflush(stdout);
    exit(1);
  }

  //Read in
  len[seqnum]=0;
  i=0;bits=0;term=0;
  while(fscanf(in,"%c",&ch)==1){
    if((ch>='0' && ch<='9')||(ch>='a' && ch<='z')||(ch>='A' && ch<='Z' && ch!='I' && ch!='J' && ch!='K' && ch!='L' && ch!='M' && ch!='N' && ch!='O')){
      if(len[seqnum]==MAXLEN){
        printf("Sequence too long in readseq. Truncating.\n");
        fflush(stdout);
        break;
      }
      for(shch=0;shch<alphsize[seqnum];shch++){
        if(alphabet[seqnum][shch]==ch) break;
      }
      if(shch==alphsize[seqnum]){printf("Unrecognized character in readseq.\n");fflush(stdout);(1);}
      term=term | (shch<<bits);
      bits+=NUMBITS[seqnum];
      if(bits==ESIZE[seqnum]) {seq[seqnum][i++]=term;bits=0;term=0;}
      len[seqnum]++;
    }
  }
  if(bits>0) (*seq)[i]=term;
  
  fclose(in);
}


void usage()
{
  printf("\n");
  printf("To correctly use this program, you must have the following parameters:\n");
  printf("  -i change-point_input_file_name\n");
  printf("  -c change_point_output_file_name1 ... (up to 150 file names)\n");
  printf("  -b num_burn (optional, default 0)\n");
  printf("  -s num_skip (optional, default 0)\n");
  printf("  -ng num_groups (optional, default 2)\n");
  printf("  -pg list of profiled group numbers (optional, default ng-1)\n");
  printf("  -notheta suppress theta profile (optional)\n");
  printf("  -nop suppress p profile (optional)\n");
  fflush(stdout);
  exit(1);
}


void handleargs(int argc,char **argv)
{
  int i;

  NUMBURN=0;
  NUMSKIP=0;
  ng=2;
  numpg=0;
  dotheta=dop=1;
  numseqs=0;
  if(argc==1) usage();
  for(i=1;i<argc;i+=2){
    if(!strcmp(argv[i],"-i")) {
      //strcpy(seqfilename,argv[i+1]); //Check length, same for outputfilename
      while(argv[i+numseqs+1][0]!='-'){
        numseqs++;
      }
      for(int j=0;j<numseqs;j++) strcpy(seqfilename[j],argv[i+j+1]);
      i+=numseqs-1;
	}
    else if(!strcmp(argv[i],"-c")){
      numsegfiles=0;
      while(i+1<argc && argv[i+1][0]!='-'){
        if(numsegfiles==150) usage();
        strcpy(segfilename[numsegfiles++],argv[i+1]);
        i++;
      }
      if(numsegfiles==0) usage();
      i--;
    }
    else if(!strcmp(argv[i],"-b")){
      if(sscanf(argv[i+1],"%lu",&NUMBURN)!=1) usage();
    }
    else if(!strcmp(argv[i],"-s")){
      if(sscanf(argv[i+1],"%lu",&NUMSKIP)!=1) usage();
    }
    else if(!strcmp(argv[i],"-ng")){
      if(sscanf(argv[i+1],"%hu",&ng)!=1) usage();
	  if(ng<1 || ng>MAXGRPS) usage();
      if(numpg==0) {pg[0]=ng-1;numpg=1;}
    }
    else if(!strcmp(argv[i],"-pg")){
      numpg=0;
      while(i+1<argc && argv[i+1][0]!='-'){
        if(numpg==100) usage();
        if(sscanf(argv[i+1],"%hu",&pg[numpg++])!=1) usage();
        i++;
      }
      if(numpg==0) usage();
      i--;
    }
    else if(!strcmp(argv[i],"-notheta")){
      dotheta=0;
      i--;
    }
    else if(!strcmp(argv[i],"-nop")){
      dop=0;
      i--;
    }
    else usage();
  }

  for(i=0;i<numpg;i++){
    if(pg[i]>=ng) usage();
  }


  sprintf(pfilename,"%s.%luburnin.p",segfilename[0],NUMBURN);

  for(i=0;i<numpg;i++)
    sprintf(pfilename,"%s%hu",pfilename,pg[i]);
  sprintf(thetafilename,"%s.%luburnin.th",segfilename[0],NUMBURN);
  sprintf(cpfilename,"%s.%luburnin.cp",segfilename[0],NUMBURN);
  sprintf(lenfilename,"%s.%luburnin.len",segfilename[0],NUMBURN);
  for(i=0;i<numpg;i++)
    sprintf(lenfilename,"%s%hu",lenfilename,pg[i]);
}


//void getmtot(unsigned short *seq,unsigned long *mtot,unsigned long L,unsigned long n)
void getmtot(unsigned short seqnum,unsigned long *mtot,unsigned long L,unsigned long n)
{
  unsigned short i;
  unsigned long j,k;
  unsigned short numch,term,ch;
  
  for(i=0;i<alphsize[seqnum];i++)
    mtot[i]=0;
  k=L/CHPERE[seqnum];
  numch=(unsigned short)(L%CHPERE[seqnum]);
  term=seq[seqnum][k]>>(numch*NUMBITS[seqnum]);
  for(j=0;j<n;j=j+1){
    ch=term & MASK[seqnum];
    term=term>>NUMBITS[seqnum];
	numch++;
    if(numch==CHPERE[seqnum]){term=seq[seqnum][++k];numch=0;}
    mtot[ch]++; 
  }
}


unsigned long maxseglen;
unsigned long numindex;
unsigned long *indices;
unsigned long *counts;
unsigned long ***totals; //CHANGE: add pointer to totals
float ***gamrat;         //CHANGE: add pointer to totals
float **probs;

//void setcountsandtotals(unsigned short *seq,unsigned long numcp,unsigned long *cp)
void setcountsandtotals(unsigned long numcp,unsigned long *cp)
{
  unsigned long seglen;
  unsigned long mtot[MAXALPHSIZE+1];
  pwtree tottree;
  unsigned long i;
  unsigned short j,k;
  long ind;

  //Set totals
  maxseglen=0;
  numindex=0;
  tottree=new wtree(sumalphsize*sizeof(unsigned long)); //CHANGE: use sumalphsize
  for(i=0;i<numcp+1;i++){
    seglen=cp[i+1] - cp[i];
    if(seglen>maxseglen) maxseglen=seglen;
    
    unsigned short curpos=0;
    for(unsigned short seqnum=0;seqnum<numseqs;seqnum++){
        getmtot(seqnum,mtot+curpos,cp[i],seglen);
        curpos+=alphsize[seqnum];
    }
    ind=tottree->getindex((void*)mtot);
    if(ind==-1){
      curpos=0;
      for(unsigned short seqnum=0;seqnum<numseqs;seqnum++){
          for(k=0;k<alphsize[seqnum];k++)
            totals[seqnum][numindex][k]=mtot[curpos+k];
          totals[seqnum][numindex][alphsize[seqnum]]=seglen;
          curpos+=alphsize[seqnum];
      }
      counts[numindex]=1;
      indices[i]=numindex;
      for(j=0;j<ng;j++){
          probs[numindex][j]=0;
      }
      tottree->insertword((void*)mtot,(long)numindex++);
    }
    else{ 
      counts[ind]++;
      indices[i]=(unsigned long)ind;
    }
  }
  if(maxseglen>MAXSLEN){printf("Segment too long in setcountsandtotals.\n");fflush(stdout);exit(1);}
  delete tottree;
}


void setprobs(double *curpi,double ***curalpha)
{
  unsigned long i;
  unsigned short j,k,seqnum;
  float *gr;
  double sum;
  double alpha;
  double temp1,temp2,temp3;
  unsigned long *tot;
  double rat;

  //renormalise curpi
  temp2=curpi[0];
  for(j=1;j<ng;j++){
    temp1=curpi[j];
    //temp2=logsum(temp2,temp1);
    temp3=temp2-temp1;
    if(temp3<0){temp3=-temp3;temp2=temp1;}
    temp2=temp2+lfunc2(temp3);
  }
  for(j=0;j<ng;j++)
    curpi[j]-=temp2;

  for(seqnum=0;seqnum<numseqs;seqnum++){
      for(j=0;j<ng;j++){
        sum=0;
        for(k=0;k<alphsize[seqnum];k++)
          sum+=curalpha[seqnum][j][k];

        for(k=0;k<=alphsize[seqnum];k++){
          gr=gamrat[seqnum][k];
          if(k<alphsize[seqnum]) alpha=curalpha[seqnum][j][k];
          else alpha=sum;
          gr[0]=0;
          rat=0;
          for(i=1;i<=maxseglen;i++){
            rat+=log(i-1+alpha);
            gr[i]=(float)rat;
          }
        }

        for(i=0;i<numindex;i++){
          tot=totals[seqnum][i];
          temp1=curpi[j]-gamrat[seqnum][alphsize[seqnum]][tot[alphsize[seqnum]]];
          for(k=0;k<alphsize[seqnum];k++)
            temp1+=gamrat[seqnum][k][tot[k]];
          probs[i][j]+=(float)temp1;
        }
      }
  }

  for(i=0;i<numindex;i++){
    temp2=probs[i][0]; 
    for(j=1;j<ng;j++){
      temp1=probs[i][j];

      //temp2=logsum(temp2,temp1);
      temp3=temp2-temp1;
      if(temp3<0){temp3=-temp3;temp2=temp1;}
      temp2=temp2+lfunc2(temp3);
    }
    probs[i][ng]=(float)temp2;

    for(j=0;j<ng;j++)
      probs[i][j]=(float)exp(probs[i][j]-temp2);
  } 
}


double lnlike()
{
  unsigned long i;
  double retval;
  
  retval=0;
  for(i=0;i<numindex;i++)
    retval+=counts[i]*probs[i][ng];

  return retval;
}


void setalpha(double ***alpha)
{
  unsigned long i,j,seqnum;
  for(seqnum=0;seqnum<numseqs;seqnum++){
      for(i=0;i<ng;i++) {
        alpha[seqnum][i][alphsize[seqnum]]=0;
        for(j=0;j<alphsize[seqnum];j++)
          alpha[seqnum][i][alphsize[seqnum]]+=alpha[seqnum][i][j];
      }
  }
}


int main(int argc,char **argv)
{
  FILE *in,*out,*seqfile;   //,*freq;
  unsigned long i,j,k;
  /*
  unsigned short **seq; //CHANGE: added level of pointer
  unsigned long *len;   //CHANGE: added level of pointer
  unsigned long numcp;
  double pi[MAXGRPS],***alpha; //CHANGE: added level of pointer tp alpha
  */
  unsigned long *cp;
  unsigned long sampsize,redsize;
  float *pprof;
  float *theta;
  float *thetaprof;
  unsigned short *cpcnts;
  unsigned long blocknum;
  char ch;
  double tempalpha,tempbeta;
  double clen;
  unsigned long num;
  double lnl;
  float p;
  unsigned long seed;
  unsigned short grp;
  double sum;
  double th[MAXALPHSIZE];
  double alf[MAXALPHSIZE];
  unsigned long ind;
  unsigned short l;
  double x;


  initlfunc2();

  //Seed random number generator
  seed=(unsigned long)time(NULL);
  init_genrand(seed);
  printf("Random number seed=%lu\n",seed);
  fflush(stdout);

  cp=(unsigned long*)malloc((MAXCP+2)*sizeof(unsigned long));
  if(!cp){printf("Could not allocate cp in main.\n");fflush(stdout);exit(1);}
  
  handleargs(argc,argv);

  //First determine alphabet and sequence length
  sumalphsize=0;
  alphabet=new char*[numseqs];
  len=new unsigned long[numseqs];
  for(int seqnum=0;seqnum<numseqs;seqnum++){
      getalphabet(seqnum);
      printf("alphabet%hu=%s\n",seqnum+1,alphabet[seqnum]);
      printf("alphsize=%hu\n",alphsize[seqnum]);
      fflush(stdout);
      sumalphsize+=alphsize[seqnum];
  }
  if(numseqs>1){dotheta=0;}
  else{
    if(alphsize[0]!=2) dotheta=0;
  }

  indices=(unsigned long*)malloc((MAXCP+1)*sizeof(unsigned long*));
  if(!indices){printf("Memory allocation problem in main.\n");fflush(stdout);exit(1);}
  counts=(unsigned long*)malloc((MAXCP+1)*sizeof(unsigned long*));
  if(!counts){printf("Memory allocation problem in main.\n");fflush(stdout);exit(1);}

  totals=new unsigned long**[(int)numseqs];
  for(int seqnum=0;seqnum<numseqs;seqnum++){
      totals[seqnum]=(unsigned long**)malloc((MAXCP+1)*sizeof(unsigned long*));
      if(!totals){printf("Memory allocation problem in main.\n");fflush(stdout);exit(1);}
      for(i=0;i<=MAXCP;i++){
        totals[seqnum][i]=(unsigned long*)malloc((alphsize[seqnum]+1)*sizeof(unsigned long));
        if(!totals[seqnum][i]){printf("Memory allocation problem in main.\n");fflush(stdout);exit(1);}
      }
  }

  gamrat= new float**[numseqs];
  for(int seqnum=0;seqnum<numseqs;seqnum++){
      gamrat[seqnum]=(float**)malloc((alphsize[seqnum]+1)*sizeof(float*));
      if(!gamrat){printf("Memory allocation problem in main.\n");fflush(stdout);exit(1);}
      for(k=0;k<=alphsize[seqnum];k++){
        gamrat[seqnum][k]=(float*)malloc((MAXSLEN+1)*sizeof(float));
        if(!gamrat[seqnum][k]){printf("Memory allocation problem in main.\n");fflush(stdout);exit(1);}
      }
  }

  probs=(float**)malloc((MAXCP+1)*sizeof(float*));
  if(!probs){printf("Memory allocation problem in main.\n");fflush(stdout);exit(1);}
  for(i=0;i<=MAXCP;i++){
    probs[i]=(float*)malloc((ng+1)*sizeof(unsigned long));
    if(!probs[i]){printf("Memory allocation problem in main.\n");fflush(stdout);exit(1);}
  }

  //Input sequence
  seq=new unsigned short*[numseqs];
  printf("Reading sequence.\n");
  fflush(stdout);
  //readseq(seqfilename,&seq,&len);
  for(int seqnum=0;seqnum<numseqs;seqnum++){
      readseq(seqnum);
  }
  if(len[0]>MAXLEN){
    printf("Input sequence too long. Try increasing MAXLEN.\n");
    fflush(stdout);
    exit(1);
  }
  if(numseqs>1){
      for(int seqnum=1;seqnum<numseqs;seqnum++){
          if(len[seqnum]!=len[0]){
              printf("Unequal sequence lengths.\n");
              fflush(stdout);
              exit(1);
          }
      }
  }

  printf("sequence length=%lu\n",len[0]);
  fflush(stdout);
  for(int seqnum=0;seqnum<numseqs;seqnum++){
      delete [] alphabet[seqnum];
  }
  delete [] alphabet;
  
  //Allocate alphas
  alpha = new double**[numseqs];
  for(int seqnum=0;seqnum<numseqs;seqnum++){
      alpha[seqnum]=(double**)malloc((ng)*sizeof(double*));
      for(i=0;i<ng;i++){
        alpha[seqnum][i]=(double*)malloc((alphsize[seqnum]+1)*sizeof(double));
        if(!alpha[seqnum][i]){printf("Memory allocation problem in main.\n");fflush(stdout);exit(1);}
      }
  }

  if(dop){
    out=fopen(lenfilename,"w");
    if(!out){
      printf("Could not open %s in main.\n",lenfilename);
      fflush(stdout);
	  exit(1);
    }

    printf("Producing profile.\n");
    fflush(stdout);

    //Allocate conservation profile
    pprof=(float*)malloc(len[0]*sizeof(float));
    if(!pprof){printf("Memory allocation problem in main.\n");fflush(stdout);exit(1);}
    for(i=0;i<len[0];i++) pprof[i]=0;

    redsize=0;
    for(k=0;k<numsegfiles;k++){ 
      printf("Processing file %s.\n",segfilename[k]);
      fflush(stdout);
      //Open segfile
      in=fopen(segfilename[k],"r");
      if(!in){
        printf("Could not open file %s.\n",segfilename[k]);
        fflush(stdout);
        exit(1);
      }

      /*
      freq=fopen("freq.txt","w");
	  if(!freq){printf("Could not open freq.txt.\n");exit(1);}
      */

      sampsize=0;
      while(num=getsegmentation(in,&numcp,cp,pi,alpha)){
        if(num>NUMBURN && (num-1-NUMBURN)%(NUMSKIP+1)==0){
          redsize++;
          setalpha(alpha);
      
          setcountsandtotals(numcp,cp);
          setprobs(pi,alpha);
		  
          /*
		  for(i=0;i<=numcp;i++){
		    fprintf(freq,"%6lu\t%6lu\t",cp[i],cp[i+1]-1);
		    for(j=0;j<alphsize;j++)
			  fprintf(freq,"%.6lf\t",1.0*totals[indices[i]][j]/totals[indices[i]][alphsize]);
			fprintf(freq,"\n");
		  }
          */
		  
          lnl=lnlike();
	      printf("Processing sample %lu. numcp=%lu ln-likelihood=%lf\n",num,numcp,lnl);
          fflush(stdout);
       
          //Get profile for groups in pg
          for(i=0;i<numcp+1;i++){
            p=probs[indices[i]][pg[0]];
            for(l=1;l<numpg;l++)
              p+=probs[indices[i]][pg[l]];
            for(j=cp[i];j<cp[i+1];j++) pprof[j]+=p; //Pretty sure this is the slow part
	      }

          //Get length in group pg
          clen=0;
          for(i=0;i<numindex;i++){
            p=probs[i][pg[0]];
            for(l=1;l<numpg;l++)
              p+=probs[i][pg[l]];
            for(int seqnum=0;seqnum<numseqs;seqnum++)
                clen+=counts[i]*p*totals[seqnum][i][alphsize[seqnum]];
          }
          fprintf(out,"%.3lf\n",clen);
        }
        sampsize++;
      }
    
      //fclose(freq);
      fclose(in);
    }
    fclose(out);

    printf("Writing to p files.\n");
    fflush(stdout);
    out=fopen(pfilename,"w");
    if(!out){
      printf("Could not open %s in main.\n",pfilename);
      fflush(stdout);
	  exit(1);
    }
    seqfile=fopen(seqfilename[0],"r");
    if(!seqfile){
      printf("Could not open %s in main.\n",seqfilename[0]);
	  exit(1);
    }
    blocknum=1;i=0;
    fprintf(out,"1");
    while(fscanf(seqfile,"%c",&ch)==1){
      if(ch=='#') fprintf(out,"\n%ld",++blocknum);
      else if(ch=='I') fprintf(out,",-1");
      else if(ch=='J') fprintf(out,",-2");
      else if(ch=='K') fprintf(out,",-3");
      else if(ch=='L') fprintf(out,",-4");
      else if(ch=='M') fprintf(out,",-5");
      else if(ch=='N') fprintf(out,",-6");
      else if(ch=='O') fprintf(out,",-7");
      else if((ch>='0' && ch<='9')||(ch>='a' && ch<='z')||(ch>='A' && ch<='Z')) 
        fprintf(out,",%.3f",1.0f*pprof[i++]/redsize);

    }
    fprintf(out,"\n");
    fclose(seqfile);
    fclose(out);
    free(pprof);
  }

  if(dotheta){
    printf("Producing theta profile.\n");
    fflush(stdout);

    //Allocate theta profile
    thetaprof=(float*)malloc(len[0]*sizeof(float));
    if(!thetaprof){printf("Memory allocation problem in main.\n");fflush(stdout);exit(1);}
    for(i=0;i<len[0];i++) thetaprof[i]=0;
    theta=(float*)malloc((MAXCP+1)*sizeof(float));
    if(!theta){printf("Memory allocation problem in main.\n");fflush(stdout);exit(1);}

    redsize=0;
    for(k=0;k<numsegfiles;k++){ 
      printf("Processing file %s.\n",segfilename[k]);
    
      //Open segfile
      in=fopen(segfilename[k],"r");
      if(!in){
        printf("Could not open file %s.\n",segfilename[k]);
        exit(1);
      }

      sampsize=0;
      while(num=getsegmentation(in,&numcp,cp,pi,alpha)){
        if(num>NUMBURN && (num-1-NUMBURN)%(NUMSKIP+1)==0){
          redsize++;
          setalpha(alpha);
        
          setcountsandtotals(numcp,cp);
          setprobs(pi,alpha);
        
          //Assumes alphsize=2
          for(i=0;i<numindex;i++){
            theta[i]=0;
            for(j=0;j<ng;j++){ 
              tempalpha=totals[0][i][1]+alpha[0][j][1];
              tempbeta=totals[0][i][0]+alpha[0][j][0];
              theta[i]+=probs[i][j]*(float)(tempalpha/(tempalpha+tempbeta));
            }
          }
        
          lnl=0;  
          for(i=0;i<numcp+1;i++){
            ind=indices[i];
            p=theta[ind];
            for(j=cp[i];j<cp[i+1];j++) thetaprof[j]+=p; //Pretty sure this is the slow part
            
            //Choose grp
            x=genrand_real3();
            sum=0;
            for(grp=0;grp<ng-1;grp++){
              sum+=probs[ind][grp];
              if(sum>x) break;
            }
            
            //Choose theta
            for(l=0;l<alphsize[0];l++)
              alf[l]=totals[0][ind][l]+alpha[0][grp][l];
            dirichlet(th,alf,alphsize[0]);
            
            for(l=0;l<alphsize[0];l++)
              lnl+=totals[0][ind][l]*th[l];
          }
          printf("Processed sample %lu. lnl=%lf\n",num,lnl);
          fflush(stdout);
        }
        sampsize++;
      }

      fclose(in);
    }
    free(theta);
  
  
    printf("Writing theta profile.\n");
    fflush(stdout);
    out=fopen(thetafilename,"w");
    if(!out){
      printf("Could not open %s in main.\n",thetafilename);
      fflush(stdout);
      exit(1);
    }
    seqfile=fopen(seqfilename[0],"r");
    if(!seqfile){
      printf("Could not open %s in main.\n",seqfilename[0]);
      fflush(stdout);
	  exit(1);
    }
    blocknum=1;i=0;
    fprintf(out,"1");
    while(fscanf(seqfile,"%c",&ch)==1){
      if(ch=='#') fprintf(out,"\n%ld",++blocknum);
      else if(ch=='I') fprintf(out,",-1");
      else if(ch=='J') fprintf(out,",-2") ;
        else if(ch=='K') fprintf(out,",-3");
        else if(ch=='L') fprintf(out,",-4");
        else if(ch=='M') fprintf(out,",-5");
        else if(ch=='N') fprintf(out,",-6");
        else if(ch=='O') fprintf(out,",-7");
      else if((ch>='0' && ch<='9')||(ch>='a' && ch<='z')||(ch>='A' && ch<='Z')) {
        thetaprof[i]/=redsize;
        fprintf(out,",%.3f",thetaprof[i++]);
      }
    }
    fprintf(out,"\n");
    fclose(seqfile);
    fclose(out);
    free(thetaprof);
  }

  printf("Producing change-point profile.\n");
  fflush(stdout);

  //Allocate cut-point profile
  cpcnts=(unsigned short*)malloc(len[0]*sizeof(unsigned short));
  if(!cpcnts){printf("Memory allocation problem in main.\n");fflush(stdout);exit(1);}
  for(i=0;i<len[0];i++) cpcnts[i]=0;

  redsize=0;
  for(k=0;k<numsegfiles;k++){ 
    //Open segfile
    in=fopen(segfilename[k],"r");
    if(!in){
      printf("Could not open file %s.\n",segfilename[k]);
      fflush(stdout);
      exit(1);
    }

    sampsize=0;
    while(num=getsegmentation(in,&numcp,cp,pi,alpha)){
      if(num>NUMBURN && (num-NUMBURN)%(NUMSKIP+1)==0){
        printf("Processing sample %lu.\n",num);
        fflush(stdout);
        redsize++;
        for(i=1;i<=numcp;i++)
          cpcnts[cp[i]]++;
      }
      sampsize++;
    }

    fclose(in);
  }
  
  printf("Writing change-point profile.\n");
  fflush(stdout);
  out=fopen(cpfilename,"w");
  if(!out){
    printf("Could not open %s in main.\n",cpfilename);
    fflush(stdout);
	exit(1);
  }
  seqfile=fopen(seqfilename[0],"r");
  if(!seqfile){
    printf("Could not open %s in main.\n",seqfilename[0]);
    fflush(stdout);
	exit(1);
  }
  blocknum=1;i=0;
  fprintf(out,"1");
  while(fscanf(seqfile,"%c",&ch)==1){
    if(ch=='#') fprintf(out,"\n%ld",++blocknum);
    else if(ch=='I') fprintf(out,",-1");
    else if(ch=='J')fprintf(out,",-2") ;
    else if(ch=='K') fprintf(out,",-3");
    else if(ch=='L') fprintf(out,",-4");
    else if(ch=='M') fprintf(out,",-5");
    else if(ch=='N') fprintf(out,",-6");
    else if(ch=='O') fprintf(out,",-7");
    else if((ch>='0' && ch<='9')||(ch>='a' && ch<='z')||(ch>='A' && ch<='Z')) 
      fprintf(out,",%.3f",1.0f*cpcnts[i++]/redsize);
  }
  fprintf(out,"\n");
  fclose(seqfile);
  fclose(out);
  free(cpcnts);


  //clean up
  for(j=0;j<numseqs;j++){
    free(seq[j]);
  }
  delete [] seq;
  delete [] len;
  
  for(j=0;j<numseqs;j++){
      for(i=0;i<ng;i++)
        free(alpha[j][i]);
      free(alpha[j]);
  }
  delete [] alpha;

  free(cp);
  
  free(indices);
  free(counts);
  for(j=0;j<numseqs;j++){
      for(i=0;i<=MAXCP;i++)
        free(totals[j][i]);
      free(totals[j]);
      for(k=0;k<=alphsize[j];k++)
        free(gamrat[j][k]);
      free(gamrat[j]);
  }
  delete [] totals;
  delete [] gamrat;

  for(i=0;i<=MAXCP;i++)
    free(probs[i]);
  free(probs);

  return 0;
}
