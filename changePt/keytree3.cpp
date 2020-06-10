#include <stdio.h>
#include <stdlib.h>
#include "keytree3.h"

// AllocNewwnodeBlock
// wnodeAllocBlock allocator, called if new block required as when there are no wnode's free
wnodeAllocBlock *wtree::AllocNewwnodeBlock()
{
  wnodeAllocBlock *pBlock = (wnodeAllocBlock *)malloc(sizeof(wnodeAllocBlock));
  if(pBlock == NULL)
  {
	printf("Unable to malloc memory for pwnode block.\n");
	exit(1);
  }
  pBlock->pPrev = NULL;
  return(pBlock);
}


pwnode wtree::makenode()
{
  pwnode newwnode;
  wnodeAllocBlock *pNew;

  if(pwnodeFree == NULL)		// alloc new wnode blocks if no pwnode free
  {
	pNew = AllocNewwnodeBlock();
	pNew->pPrev = pwnodeAllocBlocks;
	pwnodeAllocBlocks = pNew;
	pwnodeFree = &pNew->wnodes[0];
  }
  newwnode = pwnodeFree;

  numalloc++;
  if(numalloc%wnodesPerBlock>0) pwnodeFree++;
  else pwnodeFree=NULL;

  newwnode->child[0]=newwnode->child[1]=NULL;

  return newwnode;
}


wtree::wtree(unsigned short wl)
{
  wordlen=wl;
  numalloc=0;

  pwnodeAllocBlocks = AllocNewwnodeBlock();
  pwnodeFree = &pwnodeAllocBlocks->wnodes[0];
  root=makenode();
}


wtree::~wtree()
{
  wnodeAllocBlock *pOld;

  while((pOld = pwnodeAllocBlocks)!=NULL)
  {
	pwnodeAllocBlocks = pOld->pPrev;
	free(pOld);
  }
}


void wtree::insertword(void *word,long index)
{
  pwnode cur;
  unsigned short i,j;
  int bit;
  unsigned char *pterm;
  unsigned char term;
  

  cur=root;
  pterm=(unsigned char*)word;
  for(i=0;i<wordlen;i++){
    term=pterm[i];
    for(j=0;j<8;j++){
      bit = (int)(term & 1);
      term=term >> 1;
      
      if(cur->child[bit]==NULL) 
        cur->child[bit]=makenode();
        
      cur=cur->child[bit];
    }
  }
  
  cur->index=index;
}


long wtree::getindex(void *word)
{
  unsigned short i,j;
  pwnode cur;
  int bit;
  unsigned char *pterm;
  unsigned char term;

  cur=root;
  pterm=(unsigned char*)word;
  for(i=0;i<wordlen;i++){
    term=pterm[i];
    for(j=0;j<8;j++){
      bit = (int)(term & 1);
      term=term >> 1;
      if(cur->child[bit]==NULL) return -1; //word not in tree
      cur=cur->child[bit];
    }
  }
  
  return cur->index;
}


void wtree::opnum()
{
	printf("number of nodes=%lu\n",numalloc);
}
