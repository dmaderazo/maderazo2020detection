//Code to construct binary keyword trees

//insertword inserts a word in the tree and sets its index to user-supplied value
//if word is already in the tree its index is changed to the new value

//getindex returns index if word is in tree, otherwise -1
//Clearly, -1 should never be used as an index

//makenode is used internally to allocate nodes

//opnum used for debugging purposes


typedef struct wnode* pwnode;
struct wnode{
  long index; //Used only at leaf node
  pwnode child[2];
};

#define wnodesPerBlock 100000 
// wnodeAllocBlock's are allocated with this many wnode elements
// and linked together as new wnodeAllocBlock's may be needed
typedef struct TAG_wnodeAllocBlock {
	struct TAG_wnodeAllocBlock *pPrev; // pts to previously alloc'd spwnodeAllocBlock
	struct wnode wnodes[wnodesPerBlock];	// to hold wnodes
} wnodeAllocBlock;

typedef class wtree* pwtree;
class wtree{
private:
  unsigned short wordlen;
  pwnode root;
  unsigned long numalloc;
  pwnode pwnodeFree;				// pts to next wnode avail for allocation, NULL if none free
  wnodeAllocBlock *pwnodeAllocBlocks;	// pts to MRA block
  wnodeAllocBlock *AllocNewwnodeBlock();
  pwnode makenode();  
public:
  wtree(unsigned short wl);
  ~wtree();
  void insertword(void *word,long index);
  long getindex(void *word);
  void opnum();
};
