/*JN 03/15/2007 10:10:15 AM CET*/
/*Bayes-DIVA*/


#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


char pause (FILE *);
char help(unsigned long b,FILE *fp);
		
static char *line [17][37]={

	{"General information",
	"\tDIVA is entirely controlled by commands typed from the keyboard or",
	"\tread from a batch file. This gives the program flexibility, while",
	"\tkeeping it small and portable. When entering commands, it is necessary",
	"\tto enter the full name of the command and finish with semicolon (;).",
	"\tIf it is not possible to fit a command with its arguments and options",
	"\tin a single line, just finish the line with carriage return and con-",
	"\ttinue typing on the next line. DIVA will not process the command until",
	"\tsemicolon ';' is encountered in a line fed to the program by hitting",
	"\tcarriage return. Case is ignored in all input. A good way of running",
	"\tDIVA is to prepare a batch file (text only) using your favourite",
	"\tword processor and then run it in DIVA with the 'proc' command.",
	"\tThe 'proc' command can also be used to read distribution matrices",
	"\tand trees from NEXUS files, e.g., files created with MacClade.","",
	"\r"},
	
	{"help [commands];",
	"\tPrints help information about the specified command(s). If no command",
	"\tis specified, the entire helpfile will be output.","",
	"\r" },

	{"tree [treename] treespec;",
	"\tThis command sets a tree structure. 'Treename' is an optional label,",
	"\tmaximally 16 characters long. 'Treespec' is a treespecification in",
	"\tNEXUS (parenthetical) format. The tree must be fully bifurcate and",
	"\tcontain less than 180 taxa. Taxa may be labelled using any combination",
	"\tof printing characters except parentheses and comma, but a consecutive",
	"\tseries of integers between 1 and the number of taxa is recommended.","",
	"\tExamples:","",
	"\t>tree example ((1,2),(3,(4,5));",
	"\t>tree (5,(2,(4,(1,3))));",		
	"\t>tree (apple,(pear,(1,(orange,banana))));","",		
	"\tIf the labels in the last example do not correspond to previously",
	"\tused taxon labels, new labels are set assuming that 'apple' is",
	"\ttaxon 1, 'pear' is taxon 2, '1' is taxon 3, etc.","",		
	"\r"},

	{"distribution [+distributionname] distributionspec;",
	"\tThis command sets the distributions of a set of terminal taxa.",
	"\t'Distributionname' is an optional label, maximally 16 characters long.",
	"\tThe label must begin with '+', otherwise it will be interpreted as a",
	"\tdistribution. 'Distributionspec' is a list of the distributions of the",
	"\tterminal taxa in terms of unit areas (max 15). Call the areas 'A', 'B',",
	"\t'C', etc. and specify multiple-area distributions like 'BD' or 'ACE'. ",
	"\tLetters from A to O must be used. The distributions may be preceded by",
	"\tnumeric labels corresponding to taxon numbers. If such labels are",
	"\tmissing, the first distribution will be assumed to belong to taxon 1,",
	"\tthe second to taxon 2, etc.","",
	"\tExamples (two equivalent distribution specifications):","",
	"\t>distribution +test  A AB C;",
	"\t>distribution",
	"\t>3 c",
	"\t>1 a",
	"\t>2 ab;","",
	"\r"},
	
	{"optimize [options];",
	"\tReconstructs the distribution(s) of the ancestral nodes in the last",
	"\ttree specified. Depending on the setting of parameters and the",
	"\tdifficulty of the problem, the optimization is either exact or",
	"\theuristic (signalled in output). The following options are available:",
	"    bound=x",
	"\tSets an upper bound x to the length of the optimal reconstruction.",
	"\tThis will in most cases speed up the optimization and increase the",
	"\tchances of finding an exact solution. The value of x must be smaller",
	"\tthan 250, which is the default value.",
	"    maxareas=x",
	"\tConstrains ancestral distributions to contain maximally x unit areas.",
	"\tThe value of x must be in the range from 2 to 15. The default value is",
	"\tthe total number of unit areas inhabited by the terminals. The speed",
	"\tof the optimization is strongly dependent on the value of maxareas.",
	"\tThe smaller the value, the faster the optimization.",
	"    hold=x",
	"\tSets the maximum number of alternative reconstructions that will be",
	"\tkept at a node. The value of x must be smaller than or equal to 32767.",
	"\tThe default setting is 1000. If x is set to 32767, the optimization",
	"\tis guaranteed to be exact.",
	"    keep=x",
	"\tEquivalent to hold=x.",
	"    age=x",
	"\tSets the age of the deepest node in the tree. This value is used in",
	"\tthe calculation of summary statistics if relative age classes are",
	"\tchosen. The default value is 1.0.",
	"    weight=x",
	"\tSets the weight for a particular optimization to x, which must be",
	"\tbetween 0 and 1. The weight is used in the calculation of summary",
	"\tstatistics. The default value is 1.0.",
	"    printrecs",
	"\tPrints all alternative, equally optimal reconstructions. If printrecs",
	"\tis not requested, output is restricted to a summary of the optimal",
	"\t(most parsimonious) distributions at each node","",
	"\r"},
		
	{"proc filename;",
	"\tChanges control from screen to the batch file named 'filename'. The",
	"\tbatch file must specify a series of commands as they would have been",
	"\ttyped in from the keyboard. It must be a text file, and it must be in",
	"\tthe same folder as DIVA unless a correct file path is specified.",
	"\tControl is returned to the screen console when 'proc --;' or 'return;',",
	"\tis encountered or the end of the file is reached. The 'proc' command",
	"\t can also be used to process NEXUS files containing a presence/absence",
	"\tdistribution matrix and a tree specification. If the NEXUS file",
	"\tcontains several tree descriptions, only the last will be retained.","",
	"\r"},

	{"return;",
	"\tIf encountered in a batch file, returns control to the screen console.",
	"\tIf no return statement is encountered, the batch file is read until ",
	"\tthe end of the file is reached.","",
	"\r"},

	{"output filename;",
	"\tRedirects output to the specified file. The file will be created in", 
	"\tthe same folder as DIVA. DIVA can only create new files, it refuses",
	"\tto overwrite existing files. The setting of echo determines whether",
	"\tthe output will be echoed to screen. If filename is '--', the output",
	"\tfile is closed and output is redirected to the screen console.","",
	"\r"},
		
	{"echo [option];",
	"\tDetermines whether output will be echoed to screen. Three options",
	"\tare available:",
	"    all",
	"\tCauses all output to be echoed to screen (default setting).",
	"    status",
	"\tOnly status reports and error messages will be echoed to screen.",
	"    none",
	"\tNothing will be echoed to screen.","",
	"\r"},

	{
	"\r"},

	{"showareas;",
	"\tDisplays the area labels (character labels) set in a NEXUS file.", "",
	"\r"},
		
	{"showtaxa;",
	"\tDisplays the taxon labels set in a NEXUS file or in a DIVA tree", 
	"\tdescription.","",
	"\r"},
		
	{"reset [options];",
	"\tResets the counters for summary statistics. Six options are available:",
	"    ambiguous/unambiguous",
	"\tDetermines whether or not ancestral nodes for which the optimal dis-",
	"\ttribution is ambiguous will be included in summary statistics. The",
	"\tdefault setting is to include ambiguous ancestral distributions.",
	"    relative/absolute",
	"\tDetermines the type of time classes used. If absolute is specified",
	"\tthe node ages are used as is. If relative is specified, the node ages",
	"\tare scaled to lie between 0 (for terminal nodes) and the age of the",
	"\tgroup. The default setting is absolute.",
	"    classes=x",
	"\tDefines the number of time classes used. The number must be smaller",
	"\tthan or equal to 5. The default value is 1.",
	"    interval=x",
	"\tSets the width of all time classes except the oldest one. The default",
	"\tsetting is 5 for absolute time classes and 0.5 for relative time",
	"\tclasses.",
	"    bounds n x1 [x2 .. xn]",
	"\tAllows the user to specify time classes of unequal size. The number of",
	"\tbound values is given by an integer n, which is followed by a list of",
	"\tinteger or floating point numbers (x1 etc.) defining the upper bound",
	"\tof time classes from the youngest to the oldest. If bounds are not",
	"\tgiven for all classes, the interval value is used to obtain the missing",
	"\tbound values.",
	"    sumareas=x",
	"\tConstrains the summation to only consider the first x areas. The value",
	"\tof x must be in the range from 0 to 8. The default value is 0.","",
	"\r"},

	{"nodeage nodeagespec;",
	"\tSets the ages of the ancestral nodes in a previously specified tree",
	"\taccording to the list of values in 'nodeagespec'. The order of the",
	"\tancestral nodes must follow the standard used by the program (see",
	"\tmanual for further information). If no nodeage command is issued, ",
	"\tthe age of a node is calculated as the maximum number of nodes,",
	"\tincluding the node itself, separating the node from any descendant",
	"\tterminal.","",
	"\tExample (for a six-taxon tree with five ancestral nodes):","",
	"\tnodeage 0.02 0.03 0.10 1.0 2.3;","",
	"\r"},
	
	{"sum [option];",
	"\tPrints summary statistics for the optimizations performed since the",
	"\tlast 'reset' statement. One option:",
	"    areas=x",
	"\tRestricts output of summary statistics to the first x areas. The value",
	"\tof x must be smaller than or equal to the number of areas summed. The",
	"\tdefault value is either the number of areas summed or the number of",
	"\tareas encountered in the terminals, depending on which is smaller.","",
	"\r"},
	
	{"rarefy filename1 output=filename2 areas=distributionspec [options];",
	"\tThis command is used to examine the effects of random extinction in",
	"\tcertain areas. First, the frequency of occurrence in the areas in",
	"\tdistributionspec is calculated for the optimize commands in the",
	"\tspecified batch file ('filename1'). Second, occurrences in the areas",
	"\tin distributionspec are randomly deleted such that the frequencies",
	"\tin the areas become equal. The result is written to a new batch file",
	"\t('filename2'). Options:",
	"    nrep=x;",
	"\tSets the number of replications to x. The default setting is 1.",
	"    seed=x;",
	"\tFeeds the pseudorandom number generator the seed x. The default is 1.","",
	"\r"},
	
	{"quit;",
	"\tTerminates the program.","",
	"\r"}
};
	
/*** help: print help for commands flagged in b, return 1 if interrupt, 0 if success ***/
char help(unsigned long b,FILE *fp)
{
	static char ord [17] ={ 0,1,9,2,5,7,8,6,3,4,24,25,11,14,12,13,10};
	
	unsigned long a;
	short j,k,lines;
	
	/*** help print generator **/
	
	lines=0;
	
	for (j=0;j<=16;j++) {
		a=1<<ord[j];
		if (a&b)
			for (k=0;line[j][k][0]!='\r';k++) {
				if ((lines>=16 && line[j][k+1][0]!='\t' && line[j][k+1][0]!='\0')
									 || lines>=21) {
					lines=0;
					if (pause(fp))
						return 1;
				}
				fprintf (fp,"%s\n",line[j][k]);
				lines++;
			}
	}
	
	return 0;
}
