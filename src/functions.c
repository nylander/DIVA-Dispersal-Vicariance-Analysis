/*JN 03/14/2007 11:06:43 PM CET*/
/*Bayes-DIVA*/


#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include "diva.h"

short bits (unsigned short x);

/*** bits: count bits in an unsigned short, inline function to improve speed ***/
short bits (unsigned short x)
{

	short b;
	
	for (b=0;x!=0;b++)
		
		x&=(x-1);
	
	return b;

}

/***** pause: check for interrupt *****/
/*JN 03/14/2007 02:41:07 PM CET. Function version taken from divamac.c*/
char pause (FILE *fp)
{
	extern char intflag,allowint;
	extern FILE *infile;

	if (interrupt())
		OKStop();

	return intflag;	
}

/**** interrupt: return 1 if interrupt requested, else 0 ****/
/*03/14/2007 11:14:08 PM CET*/
/*Workaround for Linux: always return 0*/
/* To Do: Listen to 'Ctrl+C', and 'Ctrl+D' */
char interrupt (void)
{
	return 0;
}


/*** stringtotree: convert treespec s to tree, return 0 if failure else ntax in tree ***/
char stringtotree (char *s)
{
	extern short *root,*smal,*large,*l,*r;
	extern float *nodeage;
	extern short ntax;
	extern char tokenset, inblock, taxlabelset;
	extern char taxlabel[][17],tokenlabel[][17];
	
	short b3,b2,b1,ia,ib,c;

	char label[MAXTAXA+1][17], match[MAXTAXA+1];
	char copy[MAXLINE],*t,u[5],*v;
	char istoken=0, istaxlabel=0;
	
	/*** set labels, count number taxa ***/
	
	strcpy (copy,s);	
	
	t=strtok (copy,"(),");
	
	for (ia=1; t!=NULL;ia++) {
		strncpy(label[ia],t,16);
		label[ia][16]='\0';
		t=strtok(NULL,"(),");
	}
	
	/*** if new set of taxa, reset taxonlabels ***/
	if ((ia-1)!=ntax) {
		if (inblock==TREEBLOCK)
			return 0;
		taxlabelset=FALSE;
		for (ib=1; ib<=MAXTAXA; ib++)
			intoa (ib, taxlabel[ib]);
	}		
	
	ntax=ia-1;
	
	if (ntax==0)
		return 0;
	
	/*** check if labels match tokenlabels; default is taxonnumbers ***/
	istoken=1;
	for (ia=1; ia<=ntax; ia++)
		match[ia]=0;
	for (ia=1; ia<=ntax; ia++) {
		for (ib=1; ib<=ntax; ib++)
			if (strcmp(label[ia],tokenlabel[ib])==0)
				break;
		if (ib>ntax || match[ib]==1) {
			istoken=0;
			break;
		}
		match[ib]=1;
	}
	if (!istoken && tokenset)
		return 0;
			
	/*** check if labels match specified taxlabels if not tokenlabels ***/
	if (!istoken && taxlabelset==TRUE) {
		istaxlabel=1;
		for (ia=1; ia<=ntax; ia++)
			match[ia]=0;
		for (ia=1; ia<=ntax; ia++) {
			for (ib=1; ib<=ntax; ib++)
				if (strcmp(label[ia],taxlabel[ib])==0)
					break;
			if (ib>ntax || match[ib]==1) {
				istaxlabel=0;
				break;
			}
			match[ib]=1;
		}
	}
	
	/*** set taxlabels from tree, reset remaining taxlabels ***/
	if (!istoken && !istaxlabel) {
		if (inblock==TREEBLOCK)
			return 0;
		istaxlabel=1;
		for (ia=1; ia<=ntax; ia++) {
			strncpy (taxlabel[ia],label[ia],16);	/** set labels **/
			taxlabel[ia][16]='\0';		
		}
		taxlabelset=TRUE;
	}
	

	/*** convert labels to taxon numbers in copy if taxlabels or preset tokens ***/
	if (tokenset || istaxlabel) {
		for (ia=0; ia<=MAXLINE; ia++)		/** empty string needed **/
			copy[ia]='\0';
		
		ia=1;			/*** label number being replaced ***/
		ib=0;			/*** position in copy **/
		b1=0;			/*** not taxlabel part **/
		
		for (t=s; *t!='\0'; t++) {
			if (*t=='(' || *t==')' || *t==',') {
				if (b1==1) {
					for (c=1; c<=ntax; c++) {			/*** find number to replace with **/
						if (istoken && strcmp(tokenlabel[c],label[ia])==0)
							break;
						if (istaxlabel && strcmp(taxlabel[c],label[ia])==0)
							break;
					}
					intoa(c,u);
					strcat(copy,u);
					ib=ib+strlen(u);
					ia++;			/*** next label up ***/
					b1=0;			/*** out of taxlabel part **/
				}
				copy[ib++]=*t;
			}
			else
				b1=1;		/*** in taxlabel part **/
		}
		
		strcpy (s,copy);			/*** copy result to s ***/
	}
		
	for (ia=1;ia<=ntax;ia++)		/*** reset nodeages ***/
		nodeage[ia]=0;
	
	/*** convert string to tree ***/
	
	strcpy (copy,s);		/*** work with copy and not with s ***/
	
	for (b3=ntax+1; b3<2*ntax; b3++) {
	
		for (ia=1; (c=copy[ia])!=')'; ia++)
			if (c=='\0')
				return 0;			/*** can't find right parenthesis **/
		
		for (ib=ia-1; (c=copy[ib])!='('; ib--)
			if (ib==0)
				return 0;			/*** can't find left parenthesis **/
		
		s=&copy[ib+1];				/*** pointer to first position **/
				
		for (t=s; isdigit (*t);t++)
			u[t-s]=*t;
		u[t-s]='\0';

		if (u[0]=='\0')
			return 0;				/*** can't find first number ***/
		
		b1=atoi(u);

		if (*t!=',')
			return 0;				/*** bad format ***/
			
		for (v=++t; isdigit (*t);t++)
			u[t-v]=*t;
		u[t-v]='\0';

		if (u[0]=='\0')
			return 0;				/*** can't find second number ***/

		b2=atoi(u);

		if (*t!=')')
			return 0;				/*** bad format ***/		
		
		/*** pair b1 and b2 in b3 ***/	
		root[b1]=root[b2]=b3;	
		l[b3]=b1;
		r[b3]=b2;
		
		/*** set min. time for node***/
		nodeage[b3]=((nodeage[b1]>nodeage[b2]) ? nodeage[b1] : nodeage[b2])+1; 
		
		
		if (b1<=ntax)				/**** set range of terminals subtended by node **/
			smal[b3]=b1;
		else
			smal[b3]=smal[b1];
		if (b2<=ntax)
			large[b3]=b2;
		else
			large[b3]=large[b2];
		
		/*** update copy, s and t points to locations in copy***/
		
		*(--s)='\0';					/*** truncate s ***/
		
		intoa(b3,u);
		
		strcat(copy,u);
		strcat(copy,++t);				
		
	}

	if (strlen(copy)!=strlen(u))
		return 0;
	
	return ntax;
}

/*** distribution: check string, fill di[i], return no. filled, MAXTAXA+1 if string ***/
/*** too long, MAXTAXA+2 if duplicated or erroneous taxonlabels ***/
short distribution (char *s)

{
	
	extern unsigned short *di,aflag;
	extern long amax;
	extern char nareas;
	
	short n;	
	unsigned short j,x;
	long k;
				
	/*** reset di[] ***/
	for (j=1; j<=MAXTAXA; j++)
		di[j]=0;
	
	for (j=1; j<=MAXTAXA; j++) {
		
		k=j;
		
		if (s==NULL)
			break;
			
		if (isdigit(*s)) {		/*** numeric label **/
			k=atol(s);
			if (k<1 || k>ntax || di[k]!=0)
				return MAXTAXA+2;

			s=strtok(NULL," ;");
			
			if (s==NULL)
				return MAXTAXA+2;
		}
		
		/*** distribution spec **/
		di[k]=disttono(s);		/*** set di using s **/
					
		if (di[k]==0)
			return MAXTAXA+1;			/*** erroneous distribution ***/
				
		x=di[k];	
		
		for (n=0;x!=0;n++)
			x>>=1;
						
		if (n>nareas)
			nareas=n;
		
		s=strtok(NULL," ;");	
				
	}
	
	k=j-1;
	
	if (s!=NULL)
		return MAXTAXA+1;				/*** too many characters in string ***/
	
	amax=1<<nareas;				/*** success, set amax ***/
	
	aflag=0;						/*** set aflag **/
	for (j=1; j<=k; j++)
		aflag|=di[j];
	
	return k;				
}	 									

/*** disttono: convert area string to distribution number ***/
unsigned short disttono (char *s)
{
	unsigned short d,x;
	short c;
	
	d=0;
	
	while ((c=*s++)!='\0') {
		c=c-'a';
		if (c<0 || c>=15)
			return 0;			/** erroneous char. in spec. **/
		x=1<<c;
		if (d&x)
			return 0;			/** repeated areas **/
		d|=x;
	}
	
	return d;
}
		

/*** printanc: print results from search on file fp ***/
void printanc (FILE *fp, unsigned short *d[])
{
	extern short ntax,*smal,*large;
	extern short min;
	extern char maxareas, taxlabel[][17];
	extern unsigned short keep;
	extern short setbound;
	extern float age,weight;
	extern char printrecs, isheuristic;
	extern short hitnodes;
	
	unsigned short a,i,j;
	char s[20];
	short lines,position;
	
	fprintf (fp,"\tsettings: maxareas=%d, bound=%d, hold=%d, weight=%.3f, age=%.3f\n",maxareas,setbound,keep,weight,age);
	
	if (isheuristic) {
		fprintf (fp, "\thold limit hit for %d node%s\n",hitnodes, (hitnodes==1)? "":"s"); 

		fprintf (fp, "\tbest reconstruction found requires %d dispersal%s\n",min, min==1? "":"s"); 
	}
	else
		fprintf (fp, "\toptimal reconstruction requires %d dispersal%s\n",min, min==1? "":"s"); 
	
	fprintf (fp,"\noptimal distributions at each node:\n");
	
	lines=8;
	for (j=ntax+1; j<2*ntax; j++) {
	
		if (lines%20==0 && pause(fp))
			return;
		
		position=fprintf (fp,"node %d (anc. of terminals %s-%s):", j, taxlabel[smal[j]],taxlabel[large[j]]);

        /* JN 03/14/2007 11:05:33 PM CET: Prin all reconstructions on one line */
		for (a=0; (i=d[j][a])!=0; a++) {
			notodist (i,s);
            position+=fprintf (fp," %s",s);
		}
/*
		for (a=0; (i=d[j][a])!=0; a++) {
			notodist (i,s);
			if (position+strlen(s)+2>80) {
				fprintf (fp,"\n");
				position=fprintf (fp," %s",s);
				lines++;
			}
			else if (a==0)
				position+=fprintf (fp," %s",s);
			else
				position+=fprintf (fp," %s",s);
		}
*/
		fprintf(fp,"\n");
		
	}

	if (printrecs==1)
		pause(fp);
	else
		fprintf (fp,"\n");

	return;	
}

/*** printrec: print results of reconstruction on file fp ***/
void printrec (FILE *fp)
{
	extern char taxlabel[][17];
	extern short ntax,*smal,*large;
	extern unsigned short nposs;
	
	short j;
	char s[20];

	fprintf (fp, "\noptimal reconstruction number %i:\n",nposs);
	
	for (j=ntax+1; j<2*ntax; j++) {
	
		if ((j-ntax+6)%20==0 && pause(fp))
			return;
		
		notodist(di[j],s);
		fprintf (fp,"node %d (anc. of terminals %s-%s): %s\n", j, taxlabel[smal[j]],taxlabel[large[j]],s);
		
	}
	
	pause(fp);
	
	return;	
}

/*** command: return no. of command, 0 if illegal ***/
char command (char *s)

{
	const char *clist [] = {"illegal command","help","distribution","echo","noecho",
						"optimize","output","proc","return","tree","quit","reset",
						"sum","rarefy","nodeage","begin","endblock","end","dimensions",
						"matrix","format","charlabels","taxlabels","translate","showareas",
						"showtaxa","utree",""};
	
	char i;
	
	for (i=1; clist[i][0]!='\0';i++) {
		if (!strcmp(clist[i],s))
			return i;
	}
	return 0;
}

/*** option: return no. of option, 0 if illegal ***/
char option (char*s)
{
	const char *clist []= {"illegal option","areas","ambiguous","unambiguous","relative",
							"absolute","classes","interval","sumareas","maxareas",
							"bound","hold","age","weight","status","all","output",
							"seed","nrep","bounds","printrecs","none","keep",""};
	char i;
	
	for (i=1; clist[i][0]!='\0';i++) {
		if (!strcmp(clist[i],s))
			return i;
	}
	return 0;
}
	                  			 
/*** eprint: print string, add newline character, echo & ask cont. if status is set ***/
char eprint (char *s)
{
	extern FILE *outfile,*infile;
	extern char status,proc;
	
	if (outfile!=stdout)
		fprintf (outfile,"%s\n",s);
	
	if (status==1) {
		printf ("%s\n",s);
		
		if (infile!=stdin || proc==13) {
			printf ("\ado you want to continue (y/n)?\n");
			if (!yes()) {
				if (infile==stdin)
					printf ("rarefaction aborted\n");
				else {
					printf ("procedure aborted - control returned to console\n");
					fclose (infile);
					infile=stdin;
				}
				return 1;
			}
			else			
				printf ("procedure continues\n");	
		}
	}		
	return 0;
}

/*** statprint: print string, add newline character, echo if status is set ***/
void statprint (char *s)
{
	extern FILE *outfile,*infile;
	extern char status;

	if (outfile!=stdout)
		fprintf (outfile,"%s\n",s);
	
	if (status==1)
		printf ("%s\n",s);
	
	return;
}
		  						
/***** optimize: exact optimization, return stat:0 for success*****/
char optimize (void)
{
	/*** MAX=max poss length, keep=no. alternatives to keep at each node **/
	/*** setbound = set max length ***/
	
	extern char maxareas,echo,sumareas;
	extern short setbound;
	extern unsigned short keep;
	extern short summax;
	extern long amax;
	extern unsigned short *di;
	extern unsigned short nposs;
	extern short *l,*r;
	extern char ambiguous,isheuristic,classes;
	extern short *tempfreq[];
	extern float *tempvicfreq[][256];
	extern char input[];
	extern short hitnodes;			/*** number of nodes for which hold limit is hit ***/
	extern char intflag;		/** interrupt flag **/
	
	extern unsigned short *d[MAXLINE/2];
	extern unsigned char *ld[MAXLINE/2];
	
	unsigned char *lleft, *lright, *ltemp;
	unsigned short j,a,allocatedtaxa,ekeep;
	
	unsigned short ia,ir,il,nright,nleft,nanc,b;
	short c;
	long ncomb=0;
	unsigned char bound,upper,lower,maximal,improvedBound,minnode[2*MAXTAXA];
	
	char stat;

	/** initialize terminals ***/
	for (j=1; j<=ntax; j++) {
		minnode[j]=0;
		d[j]=(unsigned short *) calloc (2,sizeof(unsigned short));
		*d[j]=di[j];
		ld[j]=(unsigned char *) calloc (1,sizeof (unsigned char));
	}
	allocatedtaxa=ntax;
	
	if (ld[ntax]== NULL) {
		stat=3;
		goto freedld;
	}
	
	isheuristic=0;    /* mark status of search result ***/
	stat=0; 		 /* assume success */
	hitnodes=0;			/* hold limit not hit */
	
	ekeep=((keep>KEEPMAX-15)? KEEPMAX : keep+15);	/*** update keep to hold extra one-areas **/

	/****** downpass *****/
	
	improvedBound=setbound;
	
	ltemp=(unsigned char*) calloc (amax,sizeof(unsigned char));
	
	if (ltemp==NULL) {
		stat=3;
		free (ltemp);
		goto freedld;
	}

	if (status==1)
		printf ("optimizing...\n"); /*JN 03/14/2007 11:12:28 PM CET: Include here the 'Ctrl+C' when the 'interrupt()' function is correctly written*/
	if (echo==1)
		printf ("\tdown on node:         ");
	
	for (j=ntax+1; j<2*ntax;j++) {
		if (echo==1) {
			printf ("\b\b\b\b\b\b\b\b%-3d %3d%%",j,0);
			fflush (stdout);
		}
		for (ia=1; ia<amax; ia++)
			ltemp[ia]=MAX;
		
		maximal=nposs=0;
		min=MAX;
		improvedBound=improvedBound+minnode[l[j]]+minnode[r[j]];
		
/*		if (intflag || interrupt()) {					*/
/*			if (echo==1)*/
/*				printf ("\n");*/
/*			if (OKStop()) {*/
/*				stat=5;*/
/*				free (ltemp);*/
/*				goto freedld;*/
/*			}*/
/*			if (echo==1) {*/
/*				printf ("\tdown on node: %-3d %3d%%",j,100*ia/32767L);*/
/*				fflush(stdout);*/
/*			}*/
/*		}*/

	    for (ia=1; ia<amax; ia++) {
	    	a=bits(ia);
	    	if (a<=maxareas) {
	    		if (ncomb>50000) {			/*** suitable interrupt-asking frequency ***/
	    			ncomb=0;
/*	    			if (intflag || interrupt()) {		*/
/*		    			if (echo==1)*/
/*		    				printf ("\n");*/
/*		    			if (OKStop()) {*/
/*							stat=5;*/
/*							free (ltemp);*/
/*							goto freedld;*/
/*						}*/
/*						if (echo==1) {*/
/*							printf ("\tdown on node: %-3d %3d%%",j,100*ia/32767L);*/
/*							fflush(stdout);*/
/*						}*/
/*					}*/
	    			if (echo==1) {
	    				printf ("\b\b\b\b%3d%%",100*ia/32767L);
	    				fflush (stdout);
	    			}
				}
				for (nleft=0; (il=d[l[j]][nleft])!=0; nleft++) {
	    			if ((il&ia)!=0) {
						for (nright=0; (ir=d[r[j]][nright])!=0; nright++)
	    					if ((ir&ia)!=0 && (ia&(ir|il))==ia) {  	 /* possible comb */
	    						ncomb++;
	    						b=bits(ir)+bits(il)-a;		/*  b=cost  */
	    						if (a==1)
	    							b--;
	    						b = b + ld[l[j]][nleft] + ld[r[j]][nright];
	    						c=ltemp[ia];
	    						if (b<c && b<=improvedBound) {
	    							ltemp[ia]=b;
	    							if (b<min)
	    								min=b;
	    						}
	    					}
	    			}
	    		}
	    	}
	    }
		
		minnode[j]=min;
	    improvedBound-=min;
		
		nposs=getno(ltemp,MAX);
	    
	    if (nposs==0) {			/*** setbound or MAX too low ***/
			if (echo==1)
				printf ("\n");
	    	stat=1;
	    	free (ltemp);
	    	goto freedld;
	    }
	    	
	    if (j<2*ntax-1) {	/*** not for root node  ***/
	    	bound=MAX;
	    	
	    	if (nposs>ekeep) {         /*** does not fit ****/
	    		isheuristic=1;
	    		hitnodes++;
	    		for (ia=1; ia<amax; ia++) {      /** find max ***/
	    			c=ltemp[ia];
					if (c<MAX && c>maximal)
	    				maximal=c;
	    		}	
	    		bound=min+ekeep*(maximal-min)/nposs;
	    		upper=maximal;
	    		lower=min;
	    		
	    		while (bound>lower) {
	    			nposs=getno(ltemp,bound+1);
	    			if (nposs>ekeep) {
	    				upper=bound;
	    				bound=(bound+lower)/2;
	    			}
	 				else {
	 					lower=bound;
	 					bound=(bound+upper)/2;
	 				}
	 			}
	 			
	 			nposs=getno(ltemp,lower+1); 
	 			
	 			if (nposs>ekeep) {	        /*** not possible to fit ***/
	    			if (echo==1)
	    				printf ("\n");
	 				stat=2;
	 				free (ltemp);
	 				goto freedld;
	 			}

	 			bound=lower+1;		/*** will fit with this bound ***/
	 		}
	 		
	 		/***** assign distr and lengths to node ****/
	 		d[j]=(unsigned short *) calloc (nposs+1,sizeof(unsigned short));
	 		ld[j]=(unsigned char *) calloc (nposs,sizeof(unsigned char));
	 		allocatedtaxa=j;
	 		
			if (ld[j]== NULL) {
    			if (echo==1)
    				printf ("\n");
				stat=3;
				free (ltemp);
				goto freedld;
			}

	 		a=0;
	 		for (ia=1; ia<amax; ia++) {
	    		c=ltemp[ia];
	 			if (c<bound) {
	 				*(ld[j]+a)=c;
	 				*(d[j]+a)=ia;
	 				a++;
	 			}
	 		}
			/*** always keep one-areas ***/
			for (ia=1; ia<amax; ia=ia<<1) {
	    		c=ltemp[ia];
				if (c>=bound && c<MAX) {
	 				*(ld[j]+a)=c;
	 				*(d[j]+a)=ia;
	 				a++;
	 			}
	 		}				
	 		*(d[j]+a)=0;	/*** terminate d[j] ***/
	 	}	 							
	    		
	}   /****** next node ****/
	
	/**** rootsort ******/							
	
	a=0;
	j=2*ntax-1;
	
	for (ia=1; ia<amax; ia++) {
	    c=ltemp[ia];
		if (c==min)
			a++;
	}
	
	d[j]=(unsigned short *) calloc (a+1,sizeof(unsigned short));
	ld[j]=(unsigned char *) calloc (a,sizeof(unsigned char));
		    						
	allocatedtaxa=j;
	
	if (ld[j]== NULL) {
		stat=3;
		free (ltemp);
		goto freedld;
	}
	 		
	a=0;
	
	for (ia=1; ia<amax; ia++) {
	   	c=ltemp[ia];
		if (c==min) {
			*(d[j]+a)=ia;			/*** ld is already initialized by calloc */
			a++;   							
	    }
	}
	
	*(d[j]+a)=0;
	
	
	/**** up and final ****/
	
	free (ltemp);
	lright=(unsigned char *) calloc (ekeep,sizeof(unsigned char));
	lleft=(unsigned char *) calloc (ekeep,sizeof(unsigned char));

	if (lleft== NULL) {
		stat=3;
		free(lright);
		free(lleft);
		goto freedld;
	}
	 		
	if (echo==1)
		printf ("\n\tup & final on node:    ");
		
	ncomb=0;
	
	for (j=2*ntax-1; j>ntax; j--) {

		if (echo==1) {
			printf ("\b\b\b%-3d",j);
			fflush(stdout);
		}
			
		for (ia=0; ia<ekeep; ia++)
			*(lleft+ia)=*(lright+ia)=MAX;
			
		for (nanc=0; (ia=*(d[j]+nanc))!=0; 	nanc++) {
    		if (ncomb>100000) {			/*** suitable interrupt-asking frequency ***/
    			ncomb=0;
/*	    		if (intflag || interrupt()) {			*/
/*	    			if (echo==1)*/
/*	    				printf ("\n");*/
/*	    			if (OKStop()) {*/
/*						stat=5;*/
/*						free (lright);*/
/*						free (lleft);*/
/*						goto freedld;*/
/*					}*/
/*					if (echo==1) {*/
/*						printf ("\tup & final on node: %-3d",j);*/
/*						fflush(stdout);*/
/*					}*/
/*				}*/
			}
			for (nleft=0; (il=*(d[l[j]]+nleft))!=0; nleft++)
				if ((il&ia)!=0 && *(ld[l[j]]+nleft)<=min)
					for (nright=0; (ir=*(d[r[j]]+nright))!=0;nright++)
						if ((ir&ia)!=0 && (ia&(ir|il))==ia && *(ld[r[j]]+nright)<=min) {
							/* poss. comb. */
							ncomb++;
							a=bits(ia);
							b=bits(ir)+bits(il)-a;
							if (a==1)
								b--;				/* decrease one	if anc in 1 area */					
	    					b = b + *(ld[j]+nanc);
	    					b = b + *(ld[l[j]]+nleft);
	    					b =	b + *(ld[r[j]]+nright);
	    					
	    					if (b==min) {
	    						*(lright+nright)=min;
	    						*(lleft+nleft)=min;
	    					}
	    				}	
	    }
	    							
	    nleft=nright=0;
	    
	    /*** update ldowns to lups and update distr. array ****/ 			

	    for (ia=0; *(d[l[j]]+ia)!=0; ia++)
	    	if (*(lleft+ia)==min) {
	    		*(d[l[j]]+nleft)=*(d[l[j]]+ia);
	    		*(ld[l[j]]+nleft)=min - *(ld[l[j]]+ia);
	    		nleft++;
			}
		*(d[l[j]]+nleft)=0;

	    for (ia=0; *(d[r[j]]+ia)!=0; ia++)
	    	if (*(lright+ia)==min) {
	    		*(d[r[j]]+nright)=*(d[r[j]]+ia);
	    		*(ld[r[j]]+nright)=min - *(ld[r[j]]+ia);
	    		nright++;
	    	}
		*(d[r[j]]+nright)=0;
		
	}	/***** next node ****/
		
	if (echo==1)
		printf ("\n");
		
	free (lleft);
	free (lright);
	
	/*** write results ***/
	
	if (isheuristic)
		statprint ("optimization successful - heuristic solution");
	else
		statprint ("optimization successful - exact solution");
	
	if (outfile!=stdout)
		printanc (outfile,d);
	
	if (echo==1) {
		printanc (stdout,d);
	}
	
	if (intflag) {
		stat=5;
		goto freedld;
	}
		
	nposs=0;
		
	if (stat==0 && (sumareas>0 || printrecs)) {			/*** only if necessary  **/
	
		/**** add results to counters ***/
	
		if (status==1)
			printf ("retrieving alternative reconstructions...\n\tpress command-period (Mac) or 'B' (Win) to stop\n");
			
		if (ambiguous)
			equis();
		else
			unambiguous();

		if (intflag)
			stat=5;
		else {
			transfer();
			
			if (nposs==1)
				sprintf(input,"\nthere is only one optimal reconstruction");
			else
				sprintf(input,"\nthere are %i alternative, equally optimal reconstructions",nposs);

			if (outfile!=stdout)
				fprintf (outfile,"%s\n",input);

			if (echo)
				printf ("%s\n",input);

			statprint ("results added to counters");
		}
	}
		
freedld:

	for (j=1; j<=allocatedtaxa; j++) {
		free (d[j]);
		free (ld[j]);
	}

	return stat;

}				


/**** getno: get no. of lengths shorter than bound   *****/
unsigned short getno(unsigned char *p, unsigned char bound)
{
	extern long amax;
	unsigned short i,n;
	short a;
	
	n=0;
	
	for (i=1;i<amax;i++) {
		a=*(p+i);
		if (a<bound)
			n++;
	}
	
	for (i=1;i<amax;i<<=1) {
		a=*(p+i);
		if (a>=bound && a<MAX)
			n++;
	}
	return n;
}	    							


/***** equis: get alternative reconstructions ***/
void equis (void)
{
	extern unsigned short nposs;
	extern unsigned short *di;
	extern short min;
	extern unsigned short *d[];
	extern unsigned char *ld[];
	extern char intflag;
	
	short len,c,step[MAXLINE/2];
	unsigned short a,j,pos[MAXLINE/2];
	
	j=ntax+1;
	pos[j]=0;
	step[j-1]=len=0;
	
start:
	
	len=step[j-1];
	
	if ((a=d[j][pos[j]++])==0) {
		if (j==ntax+1)
			return;				/** no more possibilities **/
		j--;
		goto start;
	}
	
	if ((a&di[l[j]])!=0 && (a&di[r[j]])!=0 && (a&(di[l[j]]|di[r[j]]))==a) {
		c=bits(di[l[j]])+bits(di[r[j]])-bits(a);
		if (bits(a)==1)
			c--;
		if ((c+len)<=min) {
			di[j]=a;
			len+=c;
			step[j]=len;
			if (j==(2*ntax-1)) {
				nposs++;
				if (sum())
					return;
				if (printrecs) {
					if (outfile!=stdout)
						printrec (outfile);
					if (echo==1) {
						printrec (stdout);
						if (intflag)
							return;
					}	
				}
				j--;
				goto start;
			}
			pos[++j]=0;
			goto start;
		}
	
	}
	goto start;				
}

/*** area: return position of rightmost bit set in b, from 0 to 15; 0 if b=0 ***/
char area (unsigned short b)
{
	char c;
	
	for (c=0; b!=0; c++) {
		if (b&1)
			break;
		b>>=1;
	}
	
	return c;
}

/*** sum: sum results for one possible reconstruction, ambiguous or unambiguous ***/
char sum(void)
{
	extern unsigned short *di;
	extern short ntax,summax;
	extern short *r,*l;
	extern short *tempfreq[];
	extern char *nodeclass;
	
	unsigned short j;
	char a;
	
	for (j=2*ntax-1; j>ntax; j--) {
		
		a=nodeclass[j];
								
		if (di[j]<summax) {
			tempfreq[a][di[j]]++;		/*** count distribution frequency ***/
			if (bits(di[j])>2 && bits(di[l[j]])+bits(di[r[j]])>0)
				vicariance (a-1,di[j],di[l[j]],di[r[j]]);	/** vicariance freq. **/
		}
			
		if (bits(di[l[j]])==2)		/*** count dispersals on left branch **/
			moveadd(a-1,di[j],di[l[j]],di[r[j]]);
						
		if (bits(di[r[j]])==2)		/*** count dispersals on right branch **/
			moveadd(a-1,di[j],di[r[j]],di[l[j]]);
			
	}
/*	if (interrupt() && OKStop())			*/
/*		return 1;*/
/*	else*/
		return 0;
}

/**** moveadd: add dispersals if any between a (anc) and b (desc) ***/
/**** considering c (other descendant of a) segment=age segment of a *****/
void moveadd(char segment,unsigned short a, unsigned short b, unsigned short c)
{
	extern float *tempmove[16][16];
	char ub,uc,d,e;
	
	if (bits(a&b)==1)			/**** 1 common area with ancestor ***/
		(*(tempmove[area(a&b)][area(b&(a^b))]+segment))++;
	else {    					 /**** 2 common areas ***/
		ub=2-bits(a&b&c);		
		uc=bits(a&c)-bits(a&b&c);
		d=area(a&b&c);
		e=area(b&~(1<<d));
		switch (ub) {
		case 0:					/**** b's areas also in c ***/
			if (uc==0) {				/*** same areas in b and c ***/
				(*(tempmove[d][e]+segment))+=0.50;
				(*(tempmove[e][d]+segment))+=0.50;
			}
			else {							/*** additional areas shared by a&c ***/
				(*(tempmove[d][e]+segment))+=(1.0/3.0);	/*** both b areas inherited also possible **/				
				(*(tempmove[e][d]+segment))+=(1.0/3.0);
			}
			break;
		case 1:								/*** 1 area (=d) in b also in c ***/		
			if (uc==0)						/*** d only area shared by a&c ***/
				(*(tempmove[e][d]+segment))+=1.0;
			if (uc!=0)			/*** additional areas shared by a&c ***/
				(*(tempmove[e][d]+segment))+=0.5;
		
		}	/*** if case 2 (not in c or ambiguous in a, b or c) then no dispersal ***/
	}
}



/*** transfer: transfer temporary results to counters,initialize temps ***/
void transfer (void)
{
	extern short ntax,summax;
	extern unsigned short nposs;
	extern float *move[16][16],*freq[],*vicfreq[][256];
	extern float weight,terms,ancs,disps,ntrees;
	extern float *tempmove[16][16], *tempvicfreq[][256];
	extern short *tempfreq[];
	extern unsigned short aflag,sumaflag;
	extern long amax;
	extern char nareas;
	
	char a,b;	
	unsigned short c,d,i,j;
			
	for (a=0; a<nareas;a++)
		for (b=0;b<nareas;b++) 
			for (j=0;j<classes;j++) {
				move[a][b][j]+=(tempmove[a][b][j]*weight/(float)nposs);
				tempmove[a][b][j]=0;
			}
			
	c=(amax<summax)? amax : summax;
	
	for (j=1;j<=classes;j++)				/*** transfer **/
		for (d=0; d<c; d++) {
			freq[j][d]+=(tempfreq[j][d]*weight/(float)nposs);
			tempfreq[j][d]=0;
			if ((b=bits(d))>2)
				for (i=1; i<(1<<(b-1));i++) {
					vicfreq[j-1][d][i]+=(tempvicfreq[j-1][d][i]*weight/(float)nposs);
					tempvicfreq[j-1][d][i]=0;
				}
		}					
					
	for (j=1;j<=ntax;j++)
		if (di[j]<summax)
			freq[0][di[j]]+=weight;
	
	sumaflag|=aflag;
	
	ntrees+=weight;
	ancs+=(weight*(ntax-1));
	terms+=(weight*ntax);
	disps+=(weight*min);
}

/*** printsum: print summary statistics on file fp for x areas ***/
void printsum(FILE *fp,char x)
{
	extern float terms,ancs,disps,limit[];
	extern float *move[16][16],*freq[];
	extern FILE *outfile, *infile;
	extern long amaxhigh;
	extern char absolute,intflag,integer;
	
	unsigned short a,b,i,j,nomax,left,right;
	char c,lines, s[20],sright[20],sleft[20],temp[100];
	float f,g[7],h[7],m[7],n[7][8];
		
	nomax=1<<x;
	
	fprintf (fp,"Summary options:\n");
	fprintf (fp,"\tNumber of areas considered: %i\n",x);
	if (nomax<amaxhigh)
		fprintf (fp,"\tNB! additional areas present in terminals\n");
	if (ambiguous)
		fprintf (fp,"\tAmbiguous events included\n");
	else
		fprintf (fp,"\tOnly unambiguous events included");
	fprintf (fp,"\tNumber of time classes: %i\n",classes);
	if (classes>1) {
		if (absolute)
			fprintf(fp,"\tAbsolute time classes (speciations from terminals)\n");
		else
			fprintf(fp,"\tRelative time classes\n");
			
		fprintf(fp,"\tClass bounds: ");
		
		if (integer) {
			for (c=1;c<classes;c++) {
				if (limit[c]-limit[c-1]==1)
					fprintf(fp,"%i, ",(int)limit[c]);
				else
					fprintf(fp,"%i-%i, ",(int)limit[c-1]+1,(int)limit[c]);
			}
			fprintf (fp, "%i-\n",(int)limit[c-1]+1);
		}
		else {										/*** floats ***/
			for (c=1;c<classes;c++)
				fprintf (fp,"<%.2f, ",limit[c]);
			fprintf (fp, ">%.2f\n",limit[c-1]);
		}
	}
	if (pause(fp))
		goto endsum;
	
	lines=3;
	
	fprintf (fp,"\nFREQUENCY OF DISTRIBUTIONS\n");
	fprintf (fp,"Distribution\t");
	if (classes==1)
		fprintf (fp," Term.\t  Anc.\n");
	else {
		lines++;
		if (absolute) {
			fprintf (fp,"Frequency in time segment (absolute time)\n            \t");
			fprintf (fp," Term.\t");
		}
		else {					/*** relative ***/
			fprintf (fp,"Frequency in time segment (relative time)\n            \t");
			fprintf (fp," Term.\t");
		}
		if (integer) {
			for (c=1;c<classes;c++) {
				if (limit[c]-limit[c-1]==1)
					fprintf(fp,"%4i\t",(int)limit[c]);
				else
					fprintf(fp,"%3i-%i\t",(int)limit[c-1]+1,(int)limit[c]);
			}
			fprintf (fp, "%3i- \tAll anc.\n",(int)limit[c-1]+1);
		}
		else {					/*** floats ***/
			for (c=1;c<classes;c++)
				fprintf (fp," <%.2f\t",limit[c]);
			fprintf (fp, " >%.2f \tAll anc.\n",limit[c-1]);
		}
	}

	f=0;							/*** reset summary variables ***/

	for (c=0;c<=(classes+1);c++) {
		g[c]=h[c]=m[c]=0;
		for (i=0;i<x;i++)
			n[c][i]=0;
	}
	
	for (a=1;a<nomax;a++) {
		notodist (a,s);
		sprintf(temp,"%-12s\t%7.3f\t",s,freq[0][a]);
		g[0]+=freq[0][a];
		h[0]+=freq[0][a]*bits(a);
		if (bits(a)>1)
			m[0]+=freq[0][a];
		for (i=j=1;i<=x;i++) {
			if (a&j)
				n[0][i-1]+=freq[0][a];
			j<<=1;
		}
		for (c=1;c<=classes;c++) {
			sprintf (s,"%7.3f\t",freq[c][a]);
			strcat (temp,s);
			g[c]+=freq[c][a];
			h[c]+=freq[c][a]*bits(a);
			if (bits(a)>1)
				m[c]+=freq[c][a];
			for (i=j=1;i<=x;i++) {
				if (a&j)
					n[c][i-1]+=freq[c][a];
				j<<=1;
			}
			f+=freq[c][a];
		}
		if (classes>1) {
			sprintf (s,"%7.3f",f);
			strcat (temp,s);
		}
		g[c]+=f;
		h[c]+=f*bits(a);
		if (bits(a)>1)
			m[c]+=f;
		for (i=j=1;i<=x;i++) {
			if (a&j)
				n[c][i-1]+=freq[c][a];
			j<<=1;
		}
		if (f!=0 || freq[0][a]!=0) {
			fprintf (fp,"%s\n",temp);
			lines++;
		}
		
		f=0;
		if (lines>=20) {
			if (pause (fp))
				goto endsum;
			lines=0;
		}
	}
	
	fprintf (fp,"Total      \t");
	for (c=0; c<=classes; c++)
		fprintf (fp,"%7.3f\t",g[c]);
	if (classes>1)
		fprintf (fp,"%7.3f\t",g[c]);
	fprintf (fp,"\n");

	fprintf (fp,"Mean # areas\t");
	for (c=0; c<=classes; c++)
		if (g[c]==0)
			fprintf (fp,"%7.3f\t",0);
		else
			fprintf (fp,"%7.3f\t",h[c]/g[c]);
	if (classes>1)
		if (g[c]==0)
			fprintf (fp,"%7.3f",0);
		else
			fprintf (fp,"%7.3f",h[c]/g[c]);
	fprintf (fp,"\n");

	fprintf (fp,"Vicar. ratio\t   -   \t");
	for (c=1; c<=classes; c++)
		if (g[c]==0)
			fprintf (fp,"%7.3f\t",0);
		else
			fprintf (fp,"%7.3f\t",m[c]/g[c]);
	if (classes>1)
		if (g[c]==0)
			fprintf (fp,"%7.3f",0);
		else
			fprintf (fp,"%7.3f",m[c]/g[c]);
	fprintf (fp,"\n");

	
	if (pause(fp))
		goto endsum;
	
	f=i=0;
	lines=3;

	fprintf (fp,"\nSPECIES RICHNESS\n");
	fprintf (fp,"Area \t");
	if (classes==1)
		fprintf (fp," Term.\t  Anc.\n");
	else {
		lines++;
		if (absolute) {
			fprintf (fp,"Frequency in time segment (absolute time)\n     \t");
			fprintf (fp," Term.\t");
		}
		else {					/*** relative ***/
			fprintf (fp,"Frequency in time segment (relative time)\n     \t");
			fprintf (fp," Term.\t");
		}
		if (integer) {
			for (c=1;c<classes;c++) {
				if (limit[c]-limit[c-1]==1)
					fprintf(fp,"%4i\t",(int)limit[c]);
				else
					fprintf(fp,"%3i-%i\t",(int)limit[c-1]+1,(int)limit[c]);
			}
			fprintf (fp, "%3i- \tAll anc.\n",(int)limit[c-1]+1);
		}
		else {					/*** floats ***/
			for (c=1;c<classes;c++)
				fprintf (fp," <%.2f\t",limit[c]);
			fprintf (fp, " >%.2f \tAll anc.\n",limit[c-1]);
		}
	}

				
	f=0;
	for (c=0;c<=classes+1;c++)
		g[c]=0;
	
	for (a=0;a<x;a++) {
		sprintf (temp,"%c    \t%7.3f\t",a+'A',n[0][a]);
		g[0]+=n[0][a];
		for (c=1;c<=classes;c++) {
			sprintf (s,"%7.3f\t",n[c][a]);
			strcat (temp,s);
			g[c]+=n[c][a];
			f+=n[c][a];
		}
		if (classes>1) {
			sprintf (s,"%7.3f",f);
			strcat (temp,s);
		}
		g[c]+=f;
		
		if (f!=0 || n[0][a]!=0) {
			fprintf (fp,"%s\n",temp);
			lines++;
		}

		f=0;
		if (lines>=20) {
			if (pause (fp))
				goto endsum;
			lines=0;
		}
	}
	
	fprintf (fp,"Total\t");
	for (c=0; c<=classes; c++)
		fprintf (fp,"%7.3f\t",g[c]);
	if (classes>1)
		fprintf (fp,"%7.3f\t",g[c]);
	fprintf (fp,"\n");

	if (pause(fp))
		goto endsum;	
	
	fprintf (fp,"\nDISPERSALS BETWEEN SINGLE AREAS\n");
	fprintf (fp,"From\tTo\t");
	if (classes==1)
		fprintf (fp,"Frequency\n");
	else {
		lines++;
		if (absolute)
			fprintf (fp,"Frequency in time segment (absolute time)\n   \t \t");
		else					/*** relative ***/
			fprintf (fp,"Frequency in time segment (relative time)\n   \t \t");
		if (integer) {
			for (c=1;c<classes;c++) {
				if (limit[c]-limit[c-1]==1)
					fprintf(fp,"%4i\t",(int)limit[c]);
				else
					fprintf(fp,"%3i-%i\t",(int)limit[c-1]+1,(int)limit[c]);
			}
			fprintf (fp, "%3i- \tAll anc.\n",(int)limit[c-1]+1);		
		}
		else {					/*** floats ***/
			for (c=1;c<classes;c++)
				fprintf (fp," <%.2f\t",limit[c]);
			fprintf (fp, " >%.2f \tAll anc.\n",limit[c-1]);
		}
	}
				
	f=0;
	for (c=0;c<=classes;c++)
		g[c]=0;
	
	for (a=0;a<x;a++)
		for (b=0;b<x;b++)
			if (b!=a) {
				sprintf (temp,"%c   \t%c \t",a+'A',b+'A');
				for (c=0;c<classes;c++) {
					sprintf (s,"%7.3f\t",move[a][b][c]);
					strcat (temp,s);
					g[c]+=move[a][b][c];
					f+=move[a][b][c];
				}
				if (classes>1) {
					sprintf (s,"%7.3f",f);
					strcat (temp,s);
				}
				
				g[c]+=f;
				
				if (f!=0) {
					fprintf (fp,"%s\n",temp);
					lines++;
				}
				
				f=0;

				if (lines>=20) {
					if (pause (fp))
						goto endsum;
					lines=0;
				}
			}
	
	if (g[classes]==0)
		fprintf (fp,"No dispersals recorded\n");
	else {
		fprintf (fp,"Total   \t");
		for (c=0; c<classes; c++)
			fprintf (fp,"%7.3f\t",g[c]);
		if (classes>1)
			fprintf (fp,"%7.3f\t",g[c]);
		fprintf (fp,"\n");
	}
	
	if (pause(fp))
		goto endsum;	
	
	f=i=0;
	lines=3;

	if (x>2) {
	
	fprintf (fp,"\nFREQUENCY OF VICARIANCE EVENTS (INVOLVING MORE THAN 2 AREAS)\n");
	fprintf (fp,"Event    \t");
	if (classes==1)
		fprintf (fp,"Frequency\n");
	else {
		lines++;
		if (absolute)
			fprintf (fp,"Frequency in time segment (absolute time)\n         \t");
		else					/*** relative ***/
			fprintf (fp,"Frequency in time segment (relative time)\n         \t");
		if (integer) {
			for (c=1;c<classes;c++) {
				if (limit[c]-limit[c-1]==1)
					fprintf(fp,"%4i\t",(int)limit[c]);
				else
					fprintf(fp,"%3i-%i\t",(int)limit[c-1]+1,(int)limit[c]);
			}
			fprintf (fp, "%3i- \tAll anc.\n",(int)limit[c-1]+1);
		}
		else {					/*** floats ***/
			for (c=1;c<classes;c++)
				fprintf (fp," <%.2f\t",limit[c]);
			fprintf (fp, " >%.2f \tAll anc.\n",limit[c-1]);
		}
	}
				
	f=0;
	for (c=0;c<=classes;c++)
		g[c]=0;
	
	for (a=1;a<nomax;a++)
		if (bits(a)>2) {
			for (b=1;b<(1<<(bits(a)-1));b++) {
				left=vnotoa(b,a);
				right=a^left;
				notodist(left,sleft);
				notodist(right,sright);
				strcat(sleft,"-");
				strcat(sleft,sright);
				sprintf (temp,"%-9s\t",sleft);
				for (c=0;c<classes;c++) {
					sprintf (s,"%7.3f\t",vicfreq[c][a][b]);
					strcat (temp,s);
					g[c]+=vicfreq[c][a][b];
					f+=vicfreq[c][a][b];
				}
				if (classes>1) {
					sprintf (s,"%7.3f",f);
					strcat (temp,s);
				}
				
				g[c]+=f;
				
				if (f!=0) {
					fprintf (fp,"%s\n",temp);
					lines++;
				}
				
				f=0;

				if (lines>=20) {
					if (pause (fp))
						goto endsum;
					lines=0;
				}
			}
		}
		
	if (g[classes]==0)
		fprintf (fp,"No vicariance events recorded\n");
	else {
		fprintf (fp,"Total    \t");
		for (c=0; c<classes; c++)
			fprintf (fp,"%7.3f\t",g[c]);
		if (classes>1)
			fprintf (fp,"%7.3f\t",g[c]);
		fprintf (fp,"\n");
	}
	if (pause (fp))
		goto endsum;	
	}
	fprintf (fp,"\nSummary tree statistics\n");
	fprintf (fp,"\tTotal number of trees: %.3f\n",ntrees);
	fprintf (fp,"\tTotal number of terminals: %.3f\n",terms);
	fprintf (fp,"\tTotal number of speciations: %.3f\n",ancs);
	fprintf (fp,"\tTotal number of dispersals: %.3f\n",disps);

endsum:

	return;
}

/**** unambiguous: put unambiguous anc distr in di[j], zero if ambiguous ***/
void unambiguous (void)
{
	extern unsigned short *di;
	extern unsigned short nposs;
	
	short j;
	
	nposs=1;
	
	for (j=ntax+1; j<2*ntax; j++) {
		if (d[j][1]==0)
			di[j]=d[j][0];		/** unambiguous assignment ***/
		else
			di[j]=0;			/** ambiguous assignment, set di[j] to zero **/
	}
	
	sum();
			
	return;
}

/*** notodist: transform distribution number in a to distribution string in s ***/
void notodist (unsigned short a, char * s)
{
	char b;
	unsigned short c;
	
	c=1;

	for (b=0; b<15; b++) {
		
		if (a&c)
			
			*s++='A'+b;
		
		c<<=1;
	}
	
	*s='\0';
	
	return;
}

/*** vicariance: record frequency of vicariance events for >2areas ***/
void vicariance (char segment,unsigned short a, unsigned short b, unsigned short c)
{
	extern float *tempvicfreq[][256];
	
	short vicposs,common;
	unsigned short i,j,k[16],distmin, uniqueleft, uniqueright, left, right,shared;
	
	uniqueleft=a&b&(~c);
	uniqueright=a&c&(~b);
	shared=a&b&c;
	
	common=bits(shared);
	
	if (common==0) {
		distmin=(uniqueleft<uniqueright) ? uniqueleft : uniqueright;
		distmin=atovno(distmin,a);
		tempvicfreq[segment][a][distmin]+=1.0;
	}
	else {
	
		vicposs=(1<<common);
	
		if ((bits(b)-bits(b&shared))==0)
			vicposs--;
		if ((bits(c)-bits(c&shared))==0)
			vicposs--;
	
		i=0;
		for (j=shared; j!=0; j&=(j-1))
			k[i++]=(1<<area(j));
		
		for (i=0; i<(1<<common);i++) {
						
			j=0;
			
			left=uniqueleft;
			
			for (c=i; c!=0; c>>=1) {
				if (c&1)
					left=left|k[j];
				j++;
			}
			
			right=a^left;
			
			distmin=(left<right)? left:right;
			distmin=atovno(distmin,a);
			tempvicfreq[segment][a][distmin]+=(1.0/(float) vicposs);		

		}
		
	}
	
	return;
}

/*** intoa: convert n to characters in s, n must be positive ***/		
void intoa (int n, char *s)
{
	char i,j,c;
	
	i=0;
	
	do {
		s[i++]=n%10+'0';
	} while ((n/=10)>0);
	
	s[i--]='\0';
	
	for (j=0;j<(i+1)/2;j++) {
		c=s[j];
		s[j]=s[i-j];
		s[i-j]=c;
	}	
	
}
	
/*** OKStop: ask if stop, set and return intflag if OK ***/
int OKStop(void)
{
	extern char intflag;
	extern FILE *infile;
	
	if (infile==stdin)
		printf( "\ado you want to stop processing of the command? (y or n)\n");
	else
		printf( "\ado you want to abort the procedure? (y or n)\n");
		
	if (yes())
		intflag=1;
	else {
		intflag=0;
		if (infile==stdin)
			printf ("processing of command continues\n");
		else
			printf ("procedure continues\n");
	}
		
	return intflag;
}

/*** OKExit: ask if exit, return if not, otherwise exit ***/
void OKExit(int sig)
{
	signal (sig, ignore);
		
	printf( "\ado you want to exit the program? (y or n)\n");
	if (yes())
		exit(0);
	
	printf("execution resumed\n");
	
	signal (sig,OKExit);

	return;
}

/*** vnotoa: convert vic. event n of ancestral dist. a to smalest vic. dist **/
unsigned short vnotoa (unsigned short n, unsigned short a)
{
	unsigned short j,area[16];
	
	char i;
	
	i=0;
	for (j=1; j<amax;j<<=1) {
		if (a&j)
			area[i++]=j;
	}
	
	for (i=j=0;n!=0;n>>=1) {
		if (n&1)
			j+=area[i];
		i++;
	}
	
	return j;
}

/*** atovno: convert smalest vic. dist. of ancestral dist. a to vic. event no. **/		
unsigned short atovno (unsigned short b, unsigned short a)
{
	unsigned short x,n;
	
	n=0;
	for (x=1;b!=0;b>>=1) {
		if (a&1) {
			if (b&1)
				n|=x;
			x<<=1;
		}
		a>>=1;
	}
	return n;
}

/*** agetoclass: convert nodeage to nodeclass ***/
void agetoclass (void)
{
	extern float *nodeage;
	extern char *nodeclass,absolute;
	extern float age;
	extern float limit[];
	
	short j;
	char b;
	
	for (j=ntax+1; j<2*ntax; j++) {
		if (absolute) {
			for (b=classes-1; nodeage[j]<=limit[b];b--)
				;
			nodeclass[j]=b+1;
		}
		else {
			for (b=classes-1; (nodeage[j]*age/(nodeage[2*ntax-1]))<=limit[b];b--)
				;
			nodeclass[j]=b+1;
		}
	}
}		

/*** matrix: set distributions according to matrixstring, return 0 if success, 1 if failure ***/
char matrix (char *s)
{
	extern unsigned short *di;
	extern char taxlabel[][17], charlabel[][17];
	extern char nareas;
	extern short ntax;
	extern char transpose;
	extern long amax;
	
	unsigned short x;
	short j,k,n;
	
	/*** reset di[] ***/
	
	for (j=1; j<=ntax; j++)
		di[j]=0;
	
	if (transpose==TRUE) {
	
		/*** transpose ***/
		
		for (j=1; j<=nareas; j++) {
		
			s=strtok(s," ;");

			if (s==NULL)
				return 1;
			
			strncpy (charlabel[j],s,16);
			charlabel[j][16]='\0';
			
			s=strtok(NULL,";");
			
			if (s==NULL)
				return 1;
			
			x=1<<(j-1);
				
			for (k=1; k<=ntax; k++) {
			
				while (!isdigit(*s)) {
					if (isgraph(*s))
						return 1;
					s++;
				}
			
				if (*s!='0' && *s!='1')
					return 1;
				
				if (*s=='1')
					di[k]=di[k]+x;
					
				s++;
			}
		}
	}
	else {
		
		/** notranspose ***/
	
		for (j=1; j<=ntax; j++) {

			s=strtok(s," ;");

			if (s==NULL || !isalpha(*s))
				return 1;
			
			strncpy (taxlabel[j],s,16);
			taxlabel[j][16]='\0';
			
			s=strtok(NULL,";");
			
			if (s==NULL)
				return 1;
				
			for (k=1; k<=nareas; k++) {
			
				while (!isdigit(*s)) {
					if (isgraph(*s))
						return 1;
					s++;
				}
			
				if (*s!='0' && *s!='1')
					return 1;
				
				x=1<<(k-1);
				
				if (*s=='1')
					di[j]=di[j]+x;
					
				s++;
			}
		}
	}
	
	/*** find nareas; check if distribution empty ***/
	
	for (j=1;j<=ntax;j++) {
		x=di[j];	
		
		if (x==0)
			return 1;			/*** distribution is empty !!! ***/
			
		for (n=0;x!=0;n++)
			x>>=1;
						
		if (n>nareas)
			nareas=n;
		
	}
	
	if ((s=strtok(s," ;"))!=NULL)
		return 1;
	else {
		amax=1<<nareas;				/*** success, set amax ***/
	
		aflag=0;						/*** set aflag **/
		for (j=1; j<=k; j++)
			aflag|=di[j];
	}

	return 0;		/*** success ***/
}

/*** charlabels: set charlabels, return no. set ***/
char charlabels (char *s)
{
	extern char charlabel[][17];
	
	char j;
	
	for (j=1; j<=15; j++) {

		if (j==1)
			s=strtok(s," ;");
		else
			s=strtok(NULL," ;");

		if (s==NULL)
			return j-1;
		
		strncpy (charlabel[j],s,16);
		charlabel[j][16]='\0';
	}
	
	if ((s=strtok(NULL," ;"))!=NULL)
		return 16;
	else
		return 16;
}

/*** taxlabels: set taxlabels, return no. set ***/
short taxlabels (char *s)
{
	extern char taxlabel[][17];
	
	short j;

	for (j=1; j<=MAXTAXA; j++) {

		if (j==1)
			s=strtok(s," ;");
		else
			s=strtok(NULL," ;");

		if (s==NULL)
			return j-1;
		
		if (!isalpha(*s))
			return MAXTAXA+1;
					
		strncpy (taxlabel[j],s,16);
		taxlabel[j][16]='\0';
		
	}
	
	if ((s=strtok(NULL," ;"))!=NULL)
		return MAXTAXA+1;
	else
		return MAXTAXA+1;
}

/*** translate: set tokenlabels, return no. set, MAXTAXA+1 if failure ***/
short translate (char *s)
{
	extern char tokenlabel[][17];
	extern char taxlabel[][17];
	extern short ntax;
	
	short j,k;
	char *t;
	
	/*** clear tokenlabels ***/
	
	for (j=1; j<=ntax;j++)
		tokenlabel[j][0]='\0';
			
	for (j=1; j<=ntax; j++) {

		if (j==1)
			s=strtok(s," ;,");
		else
			s=strtok(NULL," ;,");

		if (s==NULL)
			break;
		
		t=strtok(NULL," ;,");
		
		if (t==NULL)
			return MAXTAXA+1;
			
		for (k=1;k<=ntax; k++) {
			if (strncmp(taxlabel[k],t,16)==0)
				break;
		}
		if (k>ntax || *tokenlabel[k]!='\0')
			return MAXTAXA+1;
		
		strncpy (tokenlabel[k],s,16);
		tokenlabel[k][16]='\0';
		
	}
	
	if (s!=NULL && (s=strtok(NULL," ;"))!=NULL)
		return MAXTAXA+1;
	
	for (k=1; k<=ntax; k++)				/*** set remaining, if any, to default ***/
		if (*tokenlabel[k]=='\0')
			intoa(k,tokenlabel[k]);
	return j-1;
}

/*** showareas ***/
void showareas (FILE *fp)
{
	char c;
	
	fprintf (fp,"Area\tLabel\n");
	for (c=0;c<15;c++)
		fprintf (fp,"%c    \t%s\n",'A'+c,charlabel[c+1]);
	fprintf (fp,"\n");
}

/*** showtaxa ***/
void showtaxa (FILE *fp)
{
	short i;
	
	fprintf (fp,"No.\tTaxon\n");
	for (i=1;i<=ntax;i++) {
		fprintf (fp,"%i\t%s\n",i,taxlabel[i]);
		if ((i+1)%20==0 && pause(fp))
			return;
	}		
	fprintf(fp,"\n");
	
}

/*** ignore ***/
void ignore(int sig)
{
	signal (sig,ignore);
	intflag=1;
	
	return;
}

/*** yes ***/
char yes (void)
{
    char c;

    while(1) {
        c=getchar();
        if (c=='y' || c=='Y')
            return 1;
        else
            if (c=='n' || c=='N')
                return 0;
    }
}
