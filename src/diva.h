/*JN 03/14/2007 02:35:59 PM CET*/
/*Bayes-DIVA*/


/*#pragma once*/

/*** define magic numbers  ****/

#define MAX 250    		/*** unsigned char. indicating impossibility  ***/
#define KEEPMAX 32767   /**** max. no. alternatives to keep at each node ***/ 
#define MAXLINE 10000	/**** max input line ***/
#define MAXTAXA 180		/*** max. no. taxa ***/
#define	DIVABLOCK 0		/*** block markers ***/
#define	DATABLOCK 1
#define	TREEBLOCK 2
#define FOREIGNBLOCK 3
#define OUTOFBLOCK 4
#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

/***  function prototypes ***/

char help(unsigned long,FILE *);
char instring (FILE *, char *);
char rarefy(unsigned short);
char stringtotree (char *);
short distribution (char *);
unsigned short disttono (char *);
void printanc(FILE *, unsigned short *d[]);
void printrec (FILE *);
char command (char *);
char option (char *);
char eprint (char *);
void statprint (char *);
char optimize (void);
unsigned short getno (unsigned char *, unsigned char);
void equis (void);
char area (unsigned short);
char sum(void);
void moveadd(char, unsigned short, unsigned short, unsigned short);
void transfer (void);
void printsum (FILE *,char);
void unambiguous (void);
void notodist (unsigned short, char *);
char pause (FILE *);
void vicariance (char, unsigned short, unsigned short, unsigned short);
void intoa (int n, char *s);
void OKExit(int sig);
int OKStop(void);
void sigHandler (int sig);
unsigned short atovno (unsigned short b, unsigned short a);
unsigned short vnotoa (unsigned short n, unsigned short a);
void agetoclass (void);
char interrupt (void);
char yes (void);
char matrix (char *s);
char charlabels (char *s);
short taxlabels (char *s);
short translate (char *s);
void showareas (FILE *fp);
void showtaxa (FILE *fp);
void ignore (int sig);
void setup();
char halt (void);

/***  global variables  ***/

FILE *infile, *outfile,*rarein,*rareout;		/*** file pointers ***/
char echo,status,tree,dist,isheuristic,ambiguous;   /** status variables **/
char input[MAXLINE+1]; 								 /** input variable **/
char nareas,maxareas,nareashigh,printrecs,proc;			 /** status variables **/
short ntax;					    		 /** number of taxa in last tree or distribution **/
long amax,amaxhigh;  				 /*** max distribution number ***/
short *l, *r, *root,*smal,*large;   /*** tree globals ***/
float *nodeage;						/*** tree global **/
char *nodeclass;					/*** tree global **/
short min;  /**** optimal reconstruction length  ****/
unsigned short *di;    /*** observed & reconstructed distributions (one poss.) ****/
short setbound;   /*** set reconstruction length ***/
unsigned short keep;     /**** number of alternative reconstructions to keep ***/
float *move[16][16], *tempmove[16][16];		/*** dispersals ***/
float terms, ancs, disps,ntrees;     /*** counters for summary statistics ***/
float weight;			/*** weight to be applied to each optimization ***/
float age; 				/*** relative age of a clade ***/
float *freq[6];			/*** frequency of distributions in time intervals ***/	
short *tempfreq[6];				/*** temporary summary variable ***/
float *tempvicfreq[5][256];		/*** temporary summary variable ***/
unsigned short nposs;		/*** number of reconstructions **/
char absolute,classes,sumareas,integer;	/*** status variables **/
short summax;			 /***max distribution number +1 corr. to sumareas ***/
float interval;					/*** segment bounds variable ***/
char taxlabel[MAXTAXA+1][17];		/*** taxon labels ***/
char tokenlabel[MAXTAXA+1][17];		/*** token labels ***/
char charlabel[16][17];				/*** area labels ***/
char intflag;				/*** interrupt flags **/
float *vicfreq[5][256];				/*** summary variable ***/
float limit[6];					/*** segment bounds variable ***/
unsigned short *d[MAXLINE/2];		/*** optimize variables ***/
unsigned char *ld[MAXLINE/2];		/*** optimize variables ***/
short nrep;							/*** variable for rarefy ***/
unsigned short aflag, sumaflag;		/*** areas in dist, areas summed ***/
char tokenset, taxlabelset;		/*** indicate whether taxon- or tokenlabels have been set ***/
char transpose;					/*** indicate transpose mode in DATA block ***/
char dimensionsset;				/*** indicate dimensionsset ***/
char inblock;					/*** indicate block ***/
short hitnodes;					/*** number of nodes for which hold limit hit ***/

