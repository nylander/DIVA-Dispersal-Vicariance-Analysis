/*JN 12/11/2007 10:07:08 PM CET*/
/*Bayes-DIVA*/

#include <stdio.h>

#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <limits.h>
#include <float.h>
#include "diva.h"

short bitcount (unsigned short x);

short bitcount (unsigned short x)
{

	short b;
	
	for (b=0;x!=0;b++)
		
		x&=(x-1);
	
	return b;

}

/*void main (void)*/
main ()
{
	unsigned short a,b,c,i,j;
	char distrname[20], treename[20],*s, *word;
/*	char distrname[MAXDISTRNAME], treename[MAXTREENAME],*s, *word;*/ /* JN */
	char temp[MAXLINE+1]; 
	long g;
	float f;
	char m;
	short nlinesread=0;
	
	/* set up interrupt response */
	signal (SIGINT,OKExit);
	intflag=0;		/*** used to signal when getting keyboard events ***/

	/** set up console if necessary ***/
	/*setup();*/
	
	/*** initialize ***/

	echo=1;					/*** default input-output option values ***/
	status=1;
	printrecs=0;
	infile=stdin;
	outfile=stdout;
	
	absolute=1;				/*** default summary option values ***/
	ambiguous=1;
	classes=1;
	interval=5;
	sumareas=0;
	summax=1<<sumareas;
	amaxhigh=nareashigh=sumaflag=0;
	ancs=terms=disps=ntrees=0;
	integer=1;
	for (j=0;j<6;j++)
		limit[j]=0;
	
	inblock=OUTOFBLOCK;			/*** default input block mode ***/
	transpose=FALSE;
	
	tokenset=0;
	taxlabelset=FALSE;		/*** default for tree input ***/
	
	for (j=1; j<=MAXTAXA; j++) {
		intoa(j,tokenlabel[j]);	/*** default tokenlabels ***/
		intoa(j,taxlabel[j]);	/*** default taxonlabels ***/
	}
	
	for (i=1;i<=15;i++) {
		charlabel[i][0]='A'+i-1;		/*** default charlabels ***/
		charlabel[i][1]='\0';
	}
	
	/*** allocate summary variables ***/
	for (j=0; j<=classes+1;j++) {
		freq[j]=(float *) calloc (summax, sizeof (float));
		tempfreq[j]=(short *) calloc (summax, sizeof (short));
	}	
	for (j=0; j<classes; j++)
		for (a=1; a<summax; a++) {
			vicfreq[j][a]=(float *) calloc (1<<(bitcount(a)-1),sizeof(float));
			tempvicfreq[j][a]=(float *) calloc (1<<(bitcount(a)-1),sizeof(float));
		}
	
	for (i=0;i<16;i++)
		for(j=0;j<16;j++) {
			move[i][j]=(float *) calloc (classes, sizeof (float));
			tempmove[i][j]=(float *) calloc (classes, sizeof (float));
		}
	
	/*** allocate tree variables **/		
	root= (short *) calloc (2*MAXTAXA,sizeof(short));
	r= (short *) calloc (2*MAXTAXA,sizeof(short));
	l= (short *) calloc (2*MAXTAXA,sizeof(short));
	nodeage= (float*) calloc (2*MAXTAXA,sizeof(float));
	smal= (short *) calloc (2*MAXTAXA,sizeof(short));
	large= (short *) calloc (2*MAXTAXA,sizeof(short));
	nodeclass=(char *) calloc (2*MAXTAXA,sizeof(char));
	
	/*** allocate distribution array **/
	di= (unsigned short *) calloc (2*MAXTAXA,sizeof(unsigned short));
	
	ntax=tree=dist=0;		/*** no taxa, no tree, no distribution set  ***/

	printf ("\nWARNING: Experimental version!\n");
	printf ("\nDIVA version 1.2.JN -- 12/11/2007 10:07:08 PM CET\n");
	printf ("Johan Nylander\n");
    printf ("This version is based on DIVA WIN/MAC version:\n\n");
	printf ("DIVA version 1.2.    4 SEPT 2001\n");
	printf ("(c) Fredrik Ronquist.\n\n");
	printf ("Type 'help;' or 'help <command>;' for help.\n");
    printf ("Type 'quit;' or use 'Ctrl+C' to quit.\n\n");
	
new_input:

	if (infile==stdin) {
		printf(">");
		inblock=OUTOFBLOCK;
		signal (SIGINT,OKExit);		/*** signal handling out of execution ***/
	}
	if (infile!=stdin) {			/*** check for interrupt ***/
		nlinesread++;
		if (nlinesread>100) {					/*** suitable interval ***/
			nlinesread=0;
/*			if (interrupt() && OKStop()) {*/ /* JN 03/14/2007 02:37:20 PM CET*/
/*				printf("procedure aborted - control returned to console\n");*/
/*				fclose (infile);*/
/*				infile=stdin;*/
/*			}*/
			intflag=0;
		}
	}
					
	s=fgets (temp,MAXLINE,infile);
	
	if (s==NULL) {
		if (infile!=stdin) {
			if (feof(infile))
				statprint("end of file - control returned to console");
			else
				statprint ("file error - control returned to console");

			fclose (infile);
			infile=stdin;
		}
		goto new_input;
	}

	for (word=&temp[0]; *word!='\0';word++)
		if (isgraph(*word))
			break;
	
	if (*word=='*' || *word=='/') {		/*** echo comments before reading next ***/
		if (infile!=stdin)
			printf ("%s",temp);
		goto new_input;
	}	
	
	if (*word=='#' || *word=='\0')		/*** ignore NEXUS statements and null lines ***/
		goto new_input;
	
	/*** copy string from temp to input, remove comments, ignore case ***/
	m=1;
	for (i=j=0; temp[i]!=';' && temp[i]!='\n' && temp[i]!='\0'; i++) {
		if (temp[i]=='[')
			m=0;
		if (m==1) {
			if (isgraph(temp[i]))
				input[j++]=tolower(temp[i]);
			else
				input[j++]=' ';
		}	
		if (temp[i]==']')
			m=1;
	}	
	
	input[j]='\0';
		
	while (temp[i]!=';' && j<(MAXLINE-1)) {			/*** get rest of command if long **/
		input[j++]=' ';
		if (infile==stdin)
			printf(">");
		s=fgets (temp,MAXLINE,infile);
		if (s==NULL) {
			eprint("incomplete instruction at end of file - control returned to console");
			infile=stdin;
			goto new_input;
		}
		for (i=0; temp[i]!=';' && temp[i]!='\n' && j<(MAXLINE-1); i++) {
			if (temp[i]=='[')
				m=0;
			if (m==1) {
				if (isgraph(temp[i]))
					input[j++]=tolower(temp[i]);
				else
					input[j++]=' ';
			}	
			if (temp[i]==']')
				m=1;
		}	
	}
	input[j]='\0';
	 	
	if (temp[i]!=';' && j>=(MAXLINE-1)) {
		eprint("error - too many characters in input");
		goto new_input;
	}
	
	word=strtok (input," ");
	
	proc=command(word);
	
	if (proc==26) /*** equate utree with tree ***/
		proc=9;

	/** block sensitive commands ***/
	if (proc>=18 && proc <=22 && inblock!=DATABLOCK)
		proc=0;
	if (proc==23 && inblock!=TREEBLOCK)
		proc=0;
	if (proc<=14 && proc!=9 && inblock!=DIVABLOCK && inblock!=OUTOFBLOCK)
		proc=0;
	if (proc==9 && inblock!=TREEBLOCK && inblock!=DIVABLOCK && inblock!=OUTOFBLOCK)
		proc=0;
	if (proc>=24 && proc<=25 && inblock!=DIVABLOCK && inblock!=OUTOFBLOCK)
		proc=0;
			
	/** warn only if appropriate, otherwise ignore ***/
	if (proc==0 && (inblock==DIVABLOCK || inblock==OUTOFBLOCK)) {
		eprint ("error - bad command");
		goto new_input;
	}
	
	
	/*** correct input ***/
	
	
	signal (SIGINT,ignore);			/*** leave signal handling to functions ***/
	
	switch (proc) {
	
	case 1:					/***** help ***/
		
		j=g=0;
		
		while ((word=strtok(NULL," ;"))!=NULL) {
			i=command(word);
			if (i==0 || (i>14 && i<24)) {
				if (eprint ("warning - illegal help option"))
					break;
			}
			else
				g|=(1<<i);
		}
		
		if (g==0)
			g=(1<<26)-1;

		if (outfile!=stdout)
			j=help(g,outfile);
		if (j==0 && echo==1)
			j=help(g,stdout);
		if (j==0) {
			if (outfile!=stdout && status==1)
				printf ("help information printed\n");
		}
		else
			statprint ("output of help information interrupted");
			
		break;
		
	case 2: 				/**** distribution ****/
		 
		dist=nareas=0;
				
		word=strtok(NULL," ;");
		
		temp[0]='\0';
		temp[16]='\0';
				
		if (word!=NULL && *word=='+') {
			strncpy (temp,word,16);
			strcat (temp," ");
			word=strtok(NULL," ;");
		}
		
		if (word==NULL) {
			eprint ("error - lacking distribution specification");
			goto new_input;
		}
		
		i=distribution (word);
		
		if (i==MAXTAXA+2) {
			eprint ("error - taxon labels duplicated or erroneous");
			goto new_input;
		}
		
		if (i==0 || i==MAXTAXA+1) {
			eprint ("error - distribution specification erroneous");
			goto new_input;
		}
		
		if (i!=ntax) { /* Here is the bug! i is not equal to ntax and tree is set to 0. JN*/
			tree=0;
            printf ("\ni is not equal to ntax. set tree to 0\n"); /* JN */
			ntax=i;
		}
		
		strcpy(distrname,temp);
		statprint (strcat(temp,"distribution read successfully"));
		dist=1;
		if (amax>amaxhigh)
			amaxhigh=amax;
		if (nareas>nareashigh)
			nareashigh=nareas;
		break;
	
	case 3:				/***  echo ***/	

		echo=status=1;		/*** default ***/
		while ((word=strtok(NULL," ;"))!=NULL) {
			i=option(word);
			if (i==14)
				echo=0;
			else if (i==15)
				echo=1;
			else if (i==21)
				echo=status=0;
			else if (eprint ("warning - illegal echo option"))
					break;
		}

		if (echo)
			statprint ("echo on");
		else if (status)
			statprint ("echo status on");
		else
			statprint ("echo off");
		break;
		
	case 4:				/***  noecho ***/
		statprint ("echo off");
		echo=status=0;
		break;
		
	case 5:				/***  optimize  ***/
	
		if (tree==0) {
			eprint("error - no tree specified"); /* Here is the error message. JN */
			goto new_input;
		}
		
		if (dist==0) {
			eprint("error - no distribution specified");
			goto new_input;
		}

		weight=1;				/*** set options to default values ***/
		setbound=MAX;
		keep=1000;
		maxareas=nareas;
		age=1;
		printrecs=0;
		
		while ((word=strtok(NULL,"= ;"))!=NULL) {
			i=option(word);

			switch (i) {
			case 9:				/*** maxareas ***/
		
				word=strtok(NULL," ;");
		
				if (word!=NULL)
					g=atol(word);
				else
					g=0;
		
				if (g<2 || g>nareas) {
					if (eprint ("warning - bad maxareas argument"))
						goto new_input;
				}
				else
					maxareas=g;
				break;
				
			case 10:			/*** bound ***/
		
				word=strtok(NULL," ;");
			
				if (word!=NULL)
					g=atol(word);
				else
					g=-1;
		
				if (g<0 || g>MAX) {
					if (eprint ("warning - bad bound argument"))
						goto new_input;
				}
				else
					setbound=g;
		
				break;
			
			case 11:			/*** hold ***/
			case 22:			/*** keep ***/
				
				word=strtok(NULL," ;");
		
				if (word!=NULL)
					g=atol(word);
				else
					g=0;
		
				if (g<1 || g>KEEPMAX) {
					if (eprint ("warning - bad hold argument"))
						goto new_input;
				}
				else
					keep=g;
		
				break;
				
			case 12:			/*** age ***/
				
				word=strtok(NULL," ;");
				
				if (word!=NULL)
					f=atof(word);
				else
					f=0;
			
				if (f<=0 || f>FLT_MAX) {
					if (eprint ("warning - bad age argument"))
						goto new_input;
				}
				else
					age=f;
					
				break;
			
			case 13:			/*** weight ***/
	
				word=strtok(NULL," ;");
			
				if (word==NULL)
					f=0;
				else
					f=atof(word);
	
				if (f>0 && f<=1)
					weight=f;
				else {
					if (eprint ("warning - bad weight argument"))
						goto new_input;		 
				}
				break;

			case 20:		/*** printrecs ***/
				
				printrecs=1;
				
				break;
			
			default:
				if (eprint ("warning - illegal option"))
					goto new_input;
				break;	
			}
	
		}	
		
		if (!absolute && age!=1.0)
			agetoclass ();
		
		intflag=0;
		
		i=optimize ();
				
		if (i) {
			if (i==1 && setbound<MAX)
				eprint("error - no reconstruction shorter than bound");
			if (i==1 && setbound==MAX)
				eprint("error - optimal reconstruction requires too many dispersals");
			if (i==2 && keep==KEEPMAX)
				eprint("error - not possible to fit alternatives within boundaries");
			if (i==2 && keep<KEEPMAX)
				eprint("error - overflow of alternatives, try to increase keep");
			if (i==3)
				eprint("error - out of memory");
			if (i==4)
				eprint("error - out of memory, results not added to counters");
			if (i==5 && infile==stdin) {
				if (sumareas==0)
					eprint("optimization aborted");
				else
					eprint("optimization aborted - counters unaffected");
			}
			if (i==5 && infile!=stdin) {
				fclose (infile);
				infile=stdin;
			}			
		}
		else
			if (outfile==stdout && infile!=stdin && echo==1) {
/*JN 03/14/2007 02:37:50 PM CET*/
/*				if (halt()) {*/ 
/*					printf("procedure aborted - control returned to console\n");*/
/*					fclose (infile);*/
/*					infile=stdin;*/
/*				}*/
			}
			
		intflag=0;
		
		break;
	
	case 6:				/*** output ***/
	
		word=strtok(NULL," ;");
		
		if (word==NULL) {
			eprint ("error - lacking filename");
			goto new_input;
		}	
		
		if (outfile!=stdout) {
			fclose (outfile);
			outfile=stdout;
			statprint("previous output file closed");
		}
			
		if (strcmp(word,"--")!=0) {
			outfile=fopen(word,"r");
			if (outfile!=NULL)	{
				fclose (outfile);
				eprint ("error - this file already exists");
				break;
			}
			outfile=fopen(word,"w");
		}
		
		if (outfile==NULL) {
			outfile=stdout;
			eprint ("error - could not open file");
		}		
		else
			if (strcmp(word,"--")!=0)
				statprint (strcat(word," opened as output file"));
		break;
			
	case 7:				/*** proc ***/
		
		word=strtok(NULL," ;");
		
		if (word==NULL) {
			eprint ("error - lacking filename");
			goto new_input;
		}
		
		if (strcmp(word,"--")==0 && infile!=stdin) {
			infile=stdin;
			statprint ("control returned to console");
			break;
		}
		
		if (infile!=stdin)
			fclose (infile);
		
		infile=fopen(word,"r");
		
		if (infile==NULL) {
			infile=stdin;
			eprint ("error - could not open file");
		}
		else
			statprint (strcat(word," opened as proc file"));
				
		break;
	
	case 8:				/*** return ***/
	
		if (infile!=stdin)
			fclose (infile);
		infile=stdin;
		eprint ("control returned to console");
		
		break;
	
	case 9:				/*** tree ***/
		
		if (inblock==TREEBLOCK && dimensionsset==FALSE) {
			eprint ("error - lacking dimensions statement");
			break;
		}
		
		word=strtok(NULL," =;");
		
		if (*word=='*')
			word=strtok (NULL," =;");
		
		temp[0]=temp[16]='\0'; /* Check this '16'. JN */
		
		if (word!=NULL && *word!='(') {
			strncpy (temp,word,16); /* Check this '16'. JN */
			strcat (temp, " ");
			word=strtok(NULL," =;");
		}

		if (word==NULL) {
			eprint ("error - lacking tree specification");
			goto new_input;
		}
		
		for (;(s=strtok(NULL," ;"))!=NULL;strcat(word,s))	/*** skip spaces in spec **/
			;
		b=ntax;

		ntax=stringtotree (word);
		
		if (ntax!=b)			/*** new tree - no distribution ***/
			dist=0;
		
		if (ntax==0) {
			eprint ("error - bad tree specification");
			goto new_input;
		}
		
		age=1.0;
		agetoclass();
		
		strcpy (treename,temp);
		statprint (strcat(temp,"tree read successfully"));
		
		tree=1;					/**** tree specified ***/
		
		break;
	
	case 10:				/*** quit ***/
		
		exit (0);
		break;

	case 11:				/*** reset ***/
		
		for (j=0; j<=classes;j++) {			/*** free summary variables ***/
			free (freq[j]);
			free (tempfreq[j]);
		}
		for (j=0; j<classes; j++)
			for (a=1; a<summax; a++) {
				free (vicfreq[j][a]);
				free (tempvicfreq[j][a]);
		}
		for (i=0;i<16;i++)
			for(j=0;j<16;j++) {
				free (move[i][j]);
				free (tempmove[i][j]);
			}
	
		absolute=1;				/*** set options to default values ***/
		ambiguous=1;
		classes=1;
		interval=5;
		sumareas=0;
		integer=1;
		b=0;			/*** no interval set ***/

		terms=ancs=disps=ntrees=0;
		amaxhigh=nareashigh=sumaflag=0;
		
		for (j=0;j<6;j++)
			limit[j]=0;
		
		while ((word=strtok(NULL,"= ;"))!=NULL) {

			switch (option(word)) {

			case 2:				/*** ambiguous ***/
		
				ambiguous=1;
				break;
			
			case 3:				/*** unambiguous ***/
				
				ambiguous=0;
				break;
			
			case 4:				/*** relative ***/
				
				absolute=0;
				break;
			
			case 5:				/*** absolute ***/
			
				absolute=1;
				break;
			
			case 6:				/*** classes ***/
			
				word=strtok(NULL," ;");
		
				if (word!=NULL)
					g=atol(word);
				else
					g=0;
		
				if (g<1 || g>5) {
					if (eprint ("warning - bad classes argument"))
						goto reallocate;
				}
				else
					classes=g;
				break;
				
			case 7:				/*** interval ***/
		
				word=strtok(NULL," ;");
			
				if (word!=NULL)
					f=atof(word);
				else
					f=0;
		
				if (f<=0 || f>FLT_MAX) {
					if (eprint ("warning - bad interval argument"))
						goto reallocate;
				}
				else {
					interval=f;
					b=1;
				}		
				break;
			
			case 8:				/*** sumareas ***/
				
				word=strtok(NULL," ;");
		
				if (word!=NULL)
					g=atol(word);
				else
					g=0;
		
				if (g<0 || g>8) {
					if (eprint ("warning - bad sumareas argument"))
						goto reallocate;
				}
				else
					sumareas=g;
		
				break;
			case 19:			/*** bounds ***/
			
				word=strtok(NULL," ;");
		
				if (word!=NULL)
					g=atol(word);
				else
					g=0;
		
				if (g<1 || g>4) {
					if (eprint ("warning - bad bounds argument"))
						goto reallocate;
				}
				else
					a=g;
				
				for (j=1; j<=a; j++) {
					word=strtok(NULL," ;");
		
					if (word!=NULL)
						f=atof(word);
					else
						f=0;
					
					if (f<=limit[j-1] || f>FLT_MAX) {
						if (eprint ("warning - bad bounds argument"))
							goto reallocate;
						break;
					}
					
					else
						limit[j]=f;
						
				}
				
				for (j=a+1;j<6;j++)
					limit[j]=0;
					
				break;

			default: 			/*** illegal option ***/
				
				if (eprint ("warning - illegal option"))
					goto reallocate;
				
				break;

			}	
		
		}
		
		if (!absolute) {				/*** set interval & integer ***/
			integer=0;
			if (b==0)
				interval=0.5;
		}
			
		for (j=1;j<classes;j++) {				/*** fix bounds && check integer ***/
			if (limit[j]<limit[j-1]+1)
				limit[j]=limit[j-1]+interval;
			if (limit[j]-(int)limit[j]!=0)
				integer=0;
		}
		
	reallocate:
		summax=1<<sumareas;		/*** reallocate summary variables ***/

		for (j=0; j<=classes;j++) {	
			freq[j]=(float *) calloc (summax, sizeof (float));
			tempfreq[j]=(short *) calloc (summax, sizeof (short));
		}
		for (j=0; j<classes; j++)
			for (a=1; a<summax; a++) {
				vicfreq[j][a]=(float *) calloc (1<<(bitcount(a)-1),sizeof(float));
				tempvicfreq[j][a]=(float *) calloc (1<<(bitcount(a)-1),sizeof(float));
		}
		for (i=0;i<16;i++)
			for(j=0;j<16;j++) {
				move[i][j]=(float *) calloc (classes, sizeof (float));
				tempmove[i][j]=(float *) calloc (classes, sizeof (float));
			}
	
		statprint ("counters for summary statistics reset");
		break;
	
	case 12:				/**** sum ****/
		
		a=bitcount(sumaflag&(summax-1));
		
		if (a<2) {
			if (a==0)
				eprint ("error - no reconstructions summarized");
			else
				eprint ("error - only one area in reconstruction(s) summarized");
			goto new_input;
		}
		
		if (echo!=1 && outfile==stdout) {
			eprint ("error - no output stream selected; set echo on or specify output file");
			break;
		}
			
		while ((word=strtok(NULL,"= ;"))!=NULL) {
			switch(option(word)){
			
			case 1:							/*** areas ***/	

				word=strtok(NULL," ;");
				if (word!=NULL)
					g=atol(word);
				else
					g=0;
				
				if (g<2 || g>sumareas) {
					if(eprint ("warning - bad areas argument"))
						goto new_input;
				}
				else
					a=g;
				break;
				
			default:						/*** illegal option ***/
				
				if (eprint ("warning - illegal option"))
					goto new_input;
				break;	
			
			}
		}
		
		intflag=0;
		
		if (outfile!=stdout)
			printsum (outfile,a);

		if (echo==1 && intflag==0)
			printsum(stdout,a);

		if (intflag)
			if (infile==stdin)
				statprint ("output of summary statistics aborted");
			else {
				fclose (infile);
				infile=stdin;
			}
		else {
			if (outfile!=stdout && status)
				printf ("summary statistics calculated");
/*			if (outfile==stdout && infile!=stdin && echo==1)*/
/*				if (halt()) {*/
/*					printf("procedure aborted - control returned to console\n");*/
/*					fclose (infile);*/
/*					infile=stdin;*/
/*				}*/
		}

		intflag=0;
		
		break;
	
	case 13:								/**** rarefy ***/
		
		word=strtok(NULL," ;");
		
		if (word==NULL) {
			eprint ("error - lacking input filename");
			goto new_input;
		}
		
		rarein=fopen(word,"r");
		
		if (rarein==NULL) {
			eprint ("error - could not open file");
			goto new_input;
		}
		
		b=0;				/*** set defaults ***/
		srand(1);
		nrep=1;
		temp[0]='\0';

		while ((word=strtok(NULL,"= ;"))!=NULL) {
			switch(option(word)){

			case 1:							/*** areas ***/	

				word=strtok(NULL," ;");
				if (word!=NULL)
					i=disttono(word);
				else
					i=0;
				
				if (bitcount(i)<2) {
					if(eprint ("warning - bad areas argument")) {
						fclose (rarein);
						goto new_input;
					}
				}
				else
					b=i;
				break;
				
			case 16:						/*** output ***/
				
				word=strtok(NULL," ;");
				
				strcpy(temp,word);
				
				break;
			
			case 17:						/*** seed ***/
			
				word=strtok(NULL," ;");
				if (word!=NULL)
					g=atol(word);
				else
					g=0;
			
				if (g==0) {
					if(eprint ("warning - bad seed argument")) {
						fclose (rarein);
						goto new_input;
					}
				}
				else
					srand(g);
				
				break;
				
			case 18:						/*** nrep ***/
				
				word=strtok(NULL," ;");
				if (word!=NULL)
					g=atol(word);
				else
					g=0;
			
				if (g<=0 || g>SHRT_MAX) {
					if(eprint ("warning - bad nrep argument")) {
						fclose (rarein);
						goto new_input;
					}
				}
				else
					nrep=g;
				
				break;
				
				
			default:						/*** illegal option ***/
				
				if(eprint ("warning - illegal option")) {
					fclose (rarein);
					goto new_input;
				}
				break;	
			}
			
		}
		
		proc=0; 			/*** for correct eprint handling ***/
		
		if (b==0) {
			fclose (rarein);
			eprint ("error - no specification of areas");
			goto new_input;
		}
		
		if (temp[0]=='\0') {
			fclose (rarein);
			eprint ("error - lacking output specification");
			goto new_input;
		}
		
		rareout=fopen(temp,"r");
		if (rareout!=NULL)	{
			fclose (rareout);
			fclose (rarein);
			eprint ("error - this file already exists");
			goto new_input;
		}
		
		rareout=fopen(temp,"w");
		
		if (rareout==NULL) {
			fclose (rarein);
			eprint ("error - could not open outfile");
			goto new_input;
		}
				
		intflag=0;
		i=rarefy(b);
		
		if (i==0)
			statprint ("rarefaction completed successfully");
		
		if (i==2)
			eprint ("error - no occurrences in one of the specified areas");
		
		if (i==3)
			eprint ("error - no relevant content in file to be rarefied");
			 	
		if (intflag) {
			statprint ("rarefaction aborted");
			intflag=0;
		}

		fclose (rarein);
		fclose (rareout);
		
		intflag=0;
		
		break;
	
	case 14:			/*** nodeage ***/
		
		if (tree==0) {
			eprint ("error - no tree specified");
			goto new_input;
		}
		
		f=0;
		
		for (j=ntax+1;j<2*ntax;j++) {
		
			word=strtok(NULL," ;");
		
			f=atof(word);
			
			if (f<nodeage[l[j]] || f<nodeage[r[j]] ||
							f>FLT_MAX || f==0)
				break;
		
			nodeage[j]=f;
		}
		
		if (j!=2*ntax || strtok(NULL," ;")!=NULL) {
			eprint ("error - nodeage specification erroneous");
			goto new_input;
		}
		
		agetoclass();
				
		statprint ("nodeages set successfully");
		break;		
		
	case 15:		/*** BEGIN ***/
			
		word = strtok(NULL," ;");
		
		if (strcmp(word,"diva")==0)
			inblock=DIVABLOCK;
		else if (strcmp(word,"data")==0) {
			inblock=DATABLOCK;
			dimensionsset=FALSE;
			transpose=FALSE;
			taxlabelset=FALSE;
		}
		else if (strcmp(word,"trees")==0) {
			inblock=TREEBLOCK;
			if (dimensionsset==FALSE) {
				eprint ("error - tree block encountered without preceding dimensions statement");
				goto new_input;
			}
			tokenset=0;
		}
		else
			inblock=FOREIGNBLOCK;
		
		
		break;
				
	case 16:		/*** ENDBLOCK ***/
			
		if (inblock==TREEBLOCK) {
			for (i=1;i<=ntax; i++)
				intoa(i,tokenlabel[i]);
			tokenset=0;
		}
		inblock=OUTOFBLOCK;
			
		break;
			
	case 17:		/*** END ***/
			
		if (inblock==TREEBLOCK) {
			for (i=1;i<=ntax; i++)
				intoa(i,tokenlabel[i]);
			tokenset=0;
		}
		inblock=OUTOFBLOCK;
			
		break;
			
	case 18:		/*** DIMENSIONS ***/
	
		dist=0;
		a=b=0;

		/*** get ntax and nchar ***/
		
		while ((word=strtok(NULL," =;"))!=NULL) {
			if (strcmp(word,"ntax")==0)
				c=1;
			else if (strcmp(word,"nchar")==0)
				c=2;
			else {
				eprint ("error - incorrect dimensions statement");
				goto new_input;
			}
			
			if ((word = strtok (NULL," =;"))==NULL) {
				eprint ("error - incorrect dimensions statement");
				goto new_input;
			}

			g=atol(word);	
				
			if (c==1) {					/*** set NTAX ***/
				if (g>MAXTAXA) {
				 	eprint("error - too many taxa");
					goto new_input;
				}
				if (g<2) {
				 	eprint("error - too few taxa");
					goto new_input;
				}
				a=g;
			}	
			else {						/*** set NCHAR ***/
			
				if (g<1) {
					eprint("error - too few areas");
					goto new_input;
				}
				if (g>15) {
					eprint ("error - too many areas");
					goto new_input;
				}
				b=g;
			}
		}

		/*** ntax and nchar OK? ***/
		
		if (a==0 || b==0) {
			eprint ("error - incorrect dimensions statement");
			goto new_input;
		}
		
		if (a!=ntax) {
			ntax=a;
			tree=0;
		}
		
		nareas=b;
		
		dimensionsset=TRUE;

		break;
				
	case 19:			/***MATRIX ***/
	
		if (dimensionsset==FALSE) {
			eprint ("error - no dimensions set");
			break;
		}
		
		word=strtok(NULL,";");
		
		if (matrix(word)) {
			eprint ("error - incorrect matrix specification");
			taxlabelset=FALSE;
			for (i=0;i<15;i++) {
				charlabel[i][0]='A'+i;		/*** reset charlabels ***/
				charlabel[i][1]='\0';
			}
			for (i=1; i<=ntax; i++)
				intoa(i,taxlabel[i]);		/*** reset taxlabels ***/
			goto new_input;
		}
		else {
			statprint ("distribution matrix read successfully");
			dist=1;
			if (amax>amaxhigh)
				amaxhigh=amax;
			if (nareas>nareashigh)
				nareashigh=nareas;
			if (transpose==FALSE)
				taxlabelset=TRUE;
		}
		break;

	case 20:			/*** FORMAT ***/
	
		while ((word=strtok (NULL, " =;"))!=NULL) {
			if (strcmp(word,"transpose")==0)
				transpose=TRUE;
		}
		break;

	case 21:			/*** CHARLABELS ***/
	
		if (dimensionsset==FALSE) {
			eprint ("error - no dimensions set");
			goto new_input;
		}

		word=strtok(NULL,";");
		
		a=charlabels(word);
		
		if (a>nareas || a<=0) {
			eprint ("error - incorrect charlabels statement");
			for (i=0;i<15;i++) {
				charlabel[i][0]='A'+i;		/*** reset charlabels ***/
				charlabel[i][1]='\0';
			}
			goto new_input;
		}
			
		break;
		
	case 22:			/*** TAXLABELS ***/

		if (dimensionsset==FALSE) {
			eprint ("error - no dimensions set");
			goto new_input;
		}

		word=strtok(NULL,";");
		
		if (taxlabels(word)>ntax) {
			eprint ("error - too many taxlabels");
			for (i=1; i<=ntax; i++)
				intoa(i,taxlabel[i]);		/*** reset taxlabels ***/
			goto new_input;
		}
		
		taxlabelset=TRUE;
		
		break;
			
	case 23:			/*** TRANSLATE ***/
		
		if (dimensionsset==FALSE)
			break;
		
		tokenset=0;
		
		word=strtok(NULL,";");
		
		if (translate(word)==MAXTAXA+1) {
			eprint ("error - duplicated or incorrect number of tokens");
			for (i=1;i<=MAXTAXA; i++)
				intoa(i,tokenlabel[i]); 	/*** set tokenlabels to default ***/
		}
		else
			tokenset=1;
		break;
	
	case 24:			/*** showareas ***/
		
		if (outfile!=stdout)
			showareas (outfile);
		
		if (echo==1)
			showareas (stdout);
		break;

	case 25:			/*** showtaxa ***/
		
		if (taxlabelset==FALSE) {
			eprint ("error - no taxlabels set");
			break;
		}
		
		if (outfile!=stdout)
			showtaxa (outfile);
		
		if (echo==1)
			showtaxa (stdout);
		break;
		
	default:
		break;
		
	}
	
	goto new_input;
	
	return;
	
}

/*** rarefy: rarefaction of b's areas in rarein nrep times to rareout return 1 if ***
****  interrupted, 2 if no occurrences in one of the specified areas 3 if rarein **
****  without relevant content ******/
char rarefy (unsigned short b)
{
	extern short nrep;
	extern FILE *rarein,*rareout;
	extern char tree,dist,nareas;
	extern short ntax;
	
	unsigned short i,j,a,c,d,n,neu,ditemp[2*MAXTAXA],area[16],max;
	unsigned long threshold[16];
	char *s,*word, optimizecopy[MAXLINE],*t;
	char temp[MAXLINE],nodeagecopy[MAXLINE],nodeageset; 
	double f;
	long lg;
	float occ[16],rarest,g[2*MAXTAXA];
	
	/** initialize **/
	max=0;
	for (i=1; i<=b; i<<=1)
		if (b&i)
			area[max++]=i;
	
	/*** calculate frequency ***/
	
	a=0;		/*** no relevant command encountered ***/

	dist=tree=ntax=0;		/*** reset variables ***/
	for (c=0;c<max;c++)
		occ[c]=0;
	
	/**** read until EOF, return, quit or proc  encountered ***/
	while ((i=instring (rarein,input))!=100 && i!=8 && i!=10 && i!=7) {
		
		if (intflag)
			return 1;
		
		if (a && i==11)
			break;			/*** break if reset encountered after an optimize ***/

		word=strtok(input," ;");		/** get command **/
		
		switch (i) {
		
		case 2:			/**** distribution ***/
									
			j=dist=0;
	
			if (tree==0)
				if (eprint ("warning - distribution command encountered but no tree specified"))
					return 1;
			
			if (tree!=0) {
				word=strtok(NULL," ;");
						
				if (*word=='+')
					word=strtok(NULL," ;");
				if (word!=NULL)
					j=distribution (word);
			}
	
			if (j>MAXTAXA) {
				if (eprint("warning - erroneous distribution specification"))
					return 1;
			}
			else {
				if (j!=ntax) {
					ntax=j;
					tree=0;
				}
				dist=1;
			}
			
			break;
								
		case 5:		/*** optimize ***/
						
			a=1;
			
			if (tree==0) {
				if (eprint("warning - optimize command encountered but no tree specified"))
					return 1;
				else
					break;
			}
			if (dist==0) {
				if (eprint("warning - optimize command encountered but no distribution specified"))
					return 1;
				else
					break;
			}
				
			weight=1;
			while ((word=strtok(NULL,"= ;"))!=NULL) {
				switch (option(word)) {
						
				case 9:				/*** maxareas ***/
	
					word=strtok(NULL," ;");
	
					if (word!=NULL)
						lg=atol(word);
					else
						lg=0;
		
					if (lg<2 || lg>nareas)
						if (eprint ("warning - bad maxareas argument"))
							return 1;
					break;
				
				case 10:			/*** bound ***/
		
					word=strtok(NULL," ;");
			
					if (word!=NULL)
						lg=atol(word);
					else
						lg=MAX+1;
		
					if (lg<0 || lg>MAX)
						if (eprint ("warning - bad bound argument"))
							return 1;
					break;
			
				case 11:			/*** hold ***/
				
					word=strtok(NULL," ;");
		
					if (word!=NULL)
						lg=atol(word);
					else
						lg=0;
		
					if (lg<1 || lg>KEEPMAX)
						if (eprint ("warning - bad hold argument"))
							return 1;
					break;
				
				case 12:			/*** age ***/
				
					word=strtok(NULL," ;");
				
					if (word!=NULL)
						f=atof(word);
					else
						f=0;
			
					if (f<=0 || f>FLT_MAX)
						if (eprint ("warning - bad age argument"))
							return 1;
					break;
			
				case 13:			/*** weight ***/
	
					word=strtok(NULL," ;");
				
					if (word==NULL)
						f=0;
					else
						f=atof(word);
	
					if (f>0 && f<=1)
						weight=f;
					else
						if (eprint ("warning - bad weight argument"))
							return 1;
					break;

				case 20:		/*** printrecs ***/
					
					break;
			
				default:
					if (eprint ("warning - illegal optimize option"))
						return 1;
					break;	
				}
			}
			
			/*** summarize frequency of areas ***/
				
			for (j=1;j<=ntax;j++) {
				for(c=0;c<max;c++) {
					if(area[c]&di[j])
						occ[c]+=weight;
				}
			}
			break;
			
		case 9:				/*** tree ***/
			
			d=ntax;
			tree=ntax=0;
			nodeageset=0;
		
			word=strtok(NULL," ;");
	
			if (*word!='(')
				word=strtok(NULL," ;");
	
			if (word!=NULL) {
				/*** skip spaces in spec **/
				for (;(s=strtok(NULL," ;"))!=NULL;strcat(word,s))
					;
		
				ntax=stringtotree (word);
				
			}
			
			if (ntax==0) {
				if (eprint("warning - erroneous tree specification"))
					return 1;
			}
			else
				tree=1;
			
			if (ntax!=d)
				dist=0;
					
			break;
				
		case 14:			/*** nodeage ***/
			if (tree==0 && eprint ("warning - no tree specified"))
				return 1;
		
			for (j=1; j<=ntax; j++)
				g[j]=0;

			for (j=ntax+1;j<2*ntax;j++) {
		
			word=strtok(NULL," ;");
		
			f=atof(word);
			
			if (f<g[l[j]] || f<g[r[j]] ||
							f>FLT_MAX || f==0)
				break;
			}
			
			g[j]=f;
			
			if (j!=2*ntax || strtok(NULL," ;")!=NULL)
				if (eprint ("warning - erroneous nodeage specification"))
					return 1;
			
			nodeageset=1;
			
			break;		
		
		default:					/*** do nothing if not necessary ***/	
			break;
		}
	}
	
	if (a==0)			/*** nothing relevant in rarein ***/
		return 3;
	
	rarest=occ[0];
	
	for (c=1;c<max;c++) {
		if (occ[c]<rarest)
			rarest=occ[c];
	}
	if (rarest==0)
		return 2;
	
	for (c=0; c<max; c++)
		threshold[c]=RAND_MAX*rarest/occ[c];
	
	for (n=0; n<nrep; n++) {
		rewind(rarein);
		dist=tree=nodeageset=0;
		a=0; 	/*** no optimize command issued ***/
		
		while ((i=instring (rarein,input))!=100 && i!=8 && i!=10 && i!=7) {
			if (intflag)
				return 1;
			
			if (a && i==11)
				break;			/*** break if reset after an optimize **/
				
			word=strtok(input," ;");		/** get command **/
		
			switch (i) {
			
			case 2:					/*** distribution ***/
									
				d=dist=nareas=0;
	
				if (tree!=0) {
					word=strtok(NULL," ;");
						
					if (*word=='+')
						word=strtok(NULL," ;");
					if (word!=NULL)
						d=distribution (word);
	
					if (d==ntax)
						dist=1;
				}
				break;
					
			case 9:					/*** tree ***/
		
				d=ntax;
				tree=ntax=nodeageset=0;
		
				word=strtok(NULL," ;");
	
				if (*word!='(')
					word=strtok(NULL," ;");
	
				if (word!=NULL) {
					/*** skip spaces in spec **/
					while ((s=strtok(NULL," ;"))!=NULL)
						strcat(word,s);
		
					ntax=stringtotree (word);
				
				}
			
				if (ntax!=0)
					tree=1;
				if (ntax!=d)
					dist=0;
					
				break;
				
			case 14:				/*** nodeage ***/
			
				strcpy (nodeagecopy,"nodeage ");
				strcat(nodeagecopy,strtok(NULL,";"));
				nodeageset=1;
				break;
			
			case 5:					/*** optimize ***/
						
				a=1;
				
				if (tree==0 || dist==0)
					break;

				strcpy (optimizecopy,"optimize");
				if ((t=strtok(NULL,";"))!=NULL)
					strcat(strcat(optimizecopy," "),t);
				
				/*** exterminate ***/
				for (j=1;j<=ntax;j++)	/** copy in case dist used again **/
					ditemp[j]=di[j];
					
				for (j=1;j<=ntax;j++)
					for(c=0;c<max;c++)
						if ((area[c]&ditemp[j]) && rand()>threshold[c])
							ditemp[j]^=area[c];
				
				i=0;
				for (j=1;j<=ntax;j++)
					if (ditemp[j]!=0)
						i++;
				
				if (i<2)
					break;			/*** ignore if less than two taxa left **/
					
				for (j=ntax+1; j<2*ntax; j++)		/*** set extinct ancestors to 0 ***/
					ditemp[j]=(ditemp[l[j]]|ditemp[r[j]]);
						
				/*** write new tree ***/
				fprintf(rareout,"tree ");
				
				i=0;		/*** largest terminal accounted for**/
				j=2*ntax-1;		/*** the node where we are **/
				neu=0;		/*** no. taxa written to new tree **/
				
				while (j<2*ntax) {
					while (j>ntax) {		/*** pass up tree **/
						if (ditemp[r[j]]!=0 && ditemp[l[j]]!=0)
							fprintf (rareout,"(");
						if (smal[j]>i && ditemp[l[j]]!=0)
							j=l[j];
						else
							j=r[j];
					}
					i=j;
					intoa(++neu,temp);
					fprintf(rareout,"%s",temp);
					j=root[j];
					while (large[j]==i || ditemp[r[j]]==0) {	/** pass down **/
						if (ditemp[r[j]]!=0 && ditemp[l[j]]!=0)
							fprintf (rareout,")");
						if (j<2*ntax-1)
							j=root[j];
						else {
							j++;		/** passed out of tree=finished **/
							break;
						}
					}
					if (j<2*ntax) {
						fprintf (rareout,",");
						j=r[j];
					}
				}
					
				fprintf (rareout,";\n");
				
				/*** write distribution ***/
				fprintf (rareout,"distribution");
				for (j=1; j<=ntax; j++) {
					if (ditemp[j]!=0) {
						notodist(ditemp[j],temp);
						fprintf (rareout," %s",temp);
					}
				}
				fprintf (rareout,";\n");
				
				if (nodeageset)							/** write nodeage if set **/
					fprintf (rareout,"%s;\n",nodeagecopy);
				
				fprintf (rareout,"%s;\n",optimizecopy);	/** write optimize **/
				break;
			
			default:
				break;						/** do not write other lines **/
			
			}								
		}
	}		
		
	
	fprintf (rareout,"return;\n");
	
	ntax=tree=dist=0;
	
	return 0;
		
}	


/*** instring: input string from infile, return command (100 if EOF), string in s ***/
char instring (FILE *fp, char *s)
{
	char temp[MAXLINE+1],t[35],*word;
	short i,j,index;


/*JN 03/14/2007 02:38:42 PM CET*/
/*	if (interrupt() && OKStop())*/
/*		return 101;*/
	
	if (fgets(temp,MAXLINE,fp)==NULL)
		return 100;
	
	for (i=0; temp[i]!=';' && temp[i]!='\n'; i++)
		s[i]=tolower(temp[i]);
	
	s[i]='\0';
	
	t[30]='\0';
	strncpy (t,s,30);
	
	word=strtok(t," ;");
	
	j=0;
	
	if (word==NULL || (j=command(word))==0)
		return 0;
		
	index=i;
			
	/*** get rest of command if long **/
			
	while (temp[i]!=';' && index<MAXLINE) {
						
		input[index++]=' ';
		
		if (fgets(temp,MAXLINE,fp)==NULL)
			return 100;
					
		for (i=0; temp[i]!=';' && temp[i]!='\n' && index<MAXLINE;i++) {
			s[index]=tolower(temp[i]);
			index++;
		}
	}
	
	s[index]='\0';
			
	if (temp[i]!=';' || index==MAXLINE)
		return 0;
			
	return j;
}

