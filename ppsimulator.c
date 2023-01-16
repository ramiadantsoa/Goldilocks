#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <math.h>

#define MAXSTRINGLEN 1024

#define OPTIONAL 0
#define REQUIRED 1

/* process types */
#define RATE_PER_AREA 1
#define RATE_PER_SPECIES_COUNT 2
#define SPECIES_CONNECTIVITY 3
char *PROCESS_TYPE[4] = {"NONE", "RATE_PER_AREA", "RATE_PER_SPECIES_COUNT", "SPECIES_CONNECTIVITY"};

/* process point types */
#define NONE 0
#define PRODUCT 1
#define REACTANT 2
#define CATALYST 3
char *POINT_TYPE[4] = {"NONE", "PRODUCT", "REACTANT", "CATALYST"};

/* kernel families */
#define CONSTANT_KERNEL 1
#define TOPHAT_KERNEL 2
#define TRUNCATED_GAUSSIAN_KERNEL 3
char *KERNEL_FAMILY[4] = {"NONE", "CONSTANT", "TOPHAT", "TRUNCATED_GAUSSIAN"};

/* kernel types */
#define CONNECTIVITY_KERNEL 1
#define OTHER_KERNEL 2
char *KERNEL_TYPE[3] = {"NONE", "CONNECTIVITY_KERNEL", "OTHER_KERNEL"};

typedef struct {
  int species;
  double *position;
  double *cactivity; /* for each unique connectivity kernel defined in a model */
  double *pactivity; /* for each process defined in a model */
} Point;

typedef struct {
  char *species;
  char *coord;
} SymbolPoint;

typedef struct {  
  char *fullname; /* with parameter values, but excluding point distance argument */
  int type, family;
  int species1, species2;
  double integral, maxdensity, radius, sigma, halfivar, trunceffect;
  int *pind; /* array of process indices where this Kernel is involved */
  int num_pind; 
} Kernel;

typedef struct {  
  char *name;
  char *coord1, *coord2;  
  int ptype1, ptype2; /* CATALYST, REACTANT, PRODUCT */
  int argindex; /* index of ProcessDefinition argument list */
  int type;
} SymbolKernel;

typedef struct {
  char *name;
  char **args;
  int num_args;
  SymbolPoint *catalysts, *reactants, *products;
  int num_catalysts, num_reactants, num_products;

  /* number of individuals per species/argument index, length of array num_args */
  int *argind_catalyst_count; 
  int *argind_reactant_count; 
  int *argind_product_count; 

  int type;
  int source_argindex; /* index of process source in process call argument list */

  /* for changeSpeciesInPlace type processes, [i][0]: reactant species argument index, [i][1]: product species argument index */
  int **same_coordinate_reactant_product_pair; 
  int num_same_coordinate_reactant_product_pairs;
  
  SymbolKernel *kernels;
  int num_kernels;
} ProcessModule;

typedef struct {
  char *name;
  char **args;
  int num_args;
} ModelComponent;

typedef struct {
  ProcessModule *pm;
  int *cind, num_cind; /* indexes to unique connectivity kernels */
  int *oind, num_oind; /* indexes to other kernels */
  /* ptype1: kernel->species1, ptype2: kernel->species2, ptype PRODUCT, REACTANT, CATALYST */
  int *ctype1, *ctype2; 
  int *otype1, *otype2; 
  int source_species; /* -1 if process has no input */
  int *arg_species, num_args; /* arguments converted to species ids (integers), -1 if kernel */
} Process;

typedef struct LinkedList {
  Point *p;
  struct LinkedList *prev,*next;
} LinkedList;

typedef struct {
  /* bottom left corner */
  double *position;
  double *pactivity; /* for each process defined in a model, sum of point values within cell */
  LinkedList *first, *last;
  int num_points;
  int *species_count;
} Cell;

typedef struct {
  int dimension; 
  Cell *cell;
  int num_cells, L; /* 1D: U, 2D: U*U */
  double cell_width;
  double U; /* length of space in 1D: U=L*cell_width, 2D: U-by-U area */
  double *pactivity; /* for each process defined in a model, sum of point values within entire space */
  int *species_count; /* species_count[0] contains number of all points, species-specific start from 1 */
  int num_species_count; /* length of species_count vector in each cell */
} SimulationSpace;

char **extractFunctionArguments(char *s, int *num_args) {
  char **args;
  int i,j,slen,num,bracketcount;
  char *thisfunction = "extractFunctionArguments";

  /* some arguments can be functions themselves, i.e. 'a,b,f[c,d]' has only 3 arguments */

  /* remove trailing '_' in symbol names */

  bracketcount=0;
  num=1;
  for (i=0; s[i] != '\0'; i++) {
    if (s[i] == '[')
      bracketcount++;
    else if (s[i] == ']')
      bracketcount--;
    else if ((s[i] == ',') && (bracketcount == 0))
      num++;
  }

  if (bracketcount != 0) {
    fprintf(stderr,"ERROR (%s): mismatching brackets (bracketcount %d, num_args %d).\n",thisfunction,bracketcount,num);
    return (NULL);
  }

  *num_args = num;

  if ((args = (char **) malloc(num * sizeof(char *))) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot malloc %d char ptrs.\n",thisfunction,num);
    perror(""); return(NULL);
  }

  bracketcount=0;
  num=0;
  j=0;
  for (i=0; s[i] != '\0'; i++) {
    if (s[i] == '[')
      bracketcount++;
    else if (s[i] == ']')
      bracketcount--;
    else if ((s[i] == ',') && (bracketcount == 0)) {
      slen = i-j;
      /* remove trailing '_' in symbol names */
      if ((i>0) && (s[i-1] == '_')) slen--;
      if ((args[num] = strndup(s+j,slen)) == NULL) {
	fprintf(stderr,"ERROR (%s): cannot malloc %d chars for symbol %d/%d.\n",thisfunction, slen,num+1,*num_args);
	perror(""); return(NULL);
      }
      j=i+1;
      num++;
    }
  }
  /* last item */
  slen = i-j;
  /* remove trailing '_' in symbol names */
  if ((i>0) && (s[i-1] == '_')) slen--;
  if ((args[num] = strndup(s+j,slen)) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot malloc %d chars for symbol %d/%d.\n",thisfunction, slen,num+1,*num_args);
    perror(""); return(NULL);
  }
  
  return (args);
}

int readProcessDefinition(FILE *fp, char **names, char *filename) {
  int c,prevc=0,found=0;
  int comment=0, nameind=0, i;
  int bracketcount;
  char *thisfunction = "readProcessDefinition";

  /* find pattern: 'name1[args1] := name2[args2]' and store name1,args1,name2,args2 to names */

  for (i=0; i<4; i++)
    names[i][0] = '\0';

  i=0;
  while (!found && ((c = fgetc(fp)) != EOF)) {
    if (!comment && !isspace(c)) {

      names[nameind][i] = c;
      i++;

      if (i>=MAXSTRINGLEN) {
	fprintf(stderr,"ERROR (%s): Function definition exceeds maximum allowed characters (%d). Please check the input file '%s'.\n",thisfunction, MAXSTRINGLEN, filename);
	fprintf(stderr,"Error occurred in Module: '");
	int j;
	for (j=0; (j<i) && (names[0][j] != '\0'); j++) {
	  fprintf(stderr,"%c",names[0][j]);
	}
	fprintf(stderr,"...'\n");
	fprintf(stderr,"If the input file is ok, increase the value of MAXSTRINGLEN in line '#define MAXSTRINGLEN %d' in c-code and re-compile. Sorry for the inconvenience.\n",MAXSTRINGLEN);
	fclose(fp); exit(-1);
      }

      if ((nameind == 0) && (c == '[')) {
	i--;
	names[nameind][i] = '\0';
	nameind++;
	i=0;
      }
      else if ((nameind == 1) && (c == ']')) {
	i--;
	names[nameind][i] = '\0';
	nameind++;
	i=0;
      }
      else if ((nameind == 2) && (c == '[')) {
	i--;
	names[nameind][i] = '\0';
	nameind++;
	i=0;	
	bracketcount=1;
      }
      else if (nameind == 3) {
	if (c == '[') bracketcount++;
	else if (c == ']') bracketcount--;
	if (bracketcount == 0) {
	  i--;
	  names[nameind][i] = '\0';
	  found=1;
	}
      }
    }

    if ((c == '*') && (prevc == '(')) {
      comment = 1;
    }
    if ((c == ')') && (prevc == '*')) {
      comment = 0;
      i=0;
    }
    if (!isspace(c)) {
      prevc = c;
    }
  }

  /* read away space and potential ';' character */
  i=0;
  while (!i && ((c = fgetc(fp)) != EOF)) {
    if (!isspace(c) && (c != ';')) {
      i=1;
      ungetc(c,fp);
    }
  }
  
  return(found);
}

SymbolPoint *extractPointList(char *s, int *num_points) {
  int i,num, bracketcount, j,k, commacount, slen;
  SymbolPoint *p;
  char *thisfunction = "extractPointList";

  if (s[0] != '{') {
    fprintf(stderr,"ERROR (%s): first character '%c' instead of expected '{'.\n",thisfunction,s[0]);
    return(NULL);
  }
  if ((s[1] != '{') && (s[1] != '}')) {
    fprintf(stderr,"ERROR (%s): second character '%c' instead of expected '{' or '}'.\n",thisfunction,s[1]);
    return(NULL);
  }
  i=0;
  bracketcount=1;
  num=0;
  while (bracketcount != 0) {
    i++;
    if (s[i] == '\0') {
      fprintf(stderr,"ERROR (%s): string terminated before reading entire list (bracketcount %d).\n",thisfunction, bracketcount);
      return (NULL);
    }
    if (s[i] == '{') {
      bracketcount++;
      num++;
    }
    else if (s[i] == '}') 
      bracketcount--;
  }

  *num_points = num;

  if ((p = (SymbolPoint *) malloc (num * sizeof(SymbolPoint))) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot malloc %d SymbolPoints.\n",thisfunction,num);
    return (NULL);
  }
  
  bracketcount=1;
  commacount = 0;
  num=0;
  i=0;
  while (bracketcount != 0) {
    i++;
    if (s[i] == '\0') {
      fprintf(stderr,"ERROR (%s): string terminated before reading entire list (bracketcount %d).\n",thisfunction, bracketcount);
      return (NULL);
    }
    if (s[i] == '{') {
      bracketcount++;
      if (bracketcount > 2) {
	fprintf(stderr,"ERROR (%s): check for missing/extra brackets '{}' and commas (bracketcount %d).\n",thisfunction,bracketcount);
	return(NULL);
      }
      j=i+1;
      commacount=0;
    }
    else if (s[i] == '}') {
      bracketcount--;
      if (bracketcount != 0) {
	if (commacount != 1) {
	  fprintf(stderr,"ERROR (%s): check for missing/extra brackets '{}' and commas (commacount %d).\n",thisfunction,commacount);
	  return (NULL);
	}
	slen = i-j;
	if ((p[num].coord = strndup(s+j,slen)) == NULL) {
	  fprintf(stderr,"ERROR (%s): cannot malloc %d chars for symbol.\n",thisfunction, slen);
	  perror(""); return(NULL);
	}
      }
      commacount=0;
      num++;
    }
    else if (s[i] == ',') {
      commacount++;
      if (commacount > 1) {
	fprintf(stderr,"ERROR (%s): check for missing/extra brackets '{}' and commas (commacount %d).\n",thisfunction,commacount);
	return (NULL);
      }
      if (bracketcount > 1) {
	slen = i-j;
	if ((p[num].species = strndup(s+j,slen)) == NULL) {
	  fprintf(stderr,"ERROR (%s): cannot malloc %d chars for symbol.\n",thisfunction, slen);
	  perror(""); return(NULL);
	}
      }
      j=i+1;
    }
  }

  return(p);
}

SymbolKernel *extractKernels(char *s, int *num_kernels) {
  SymbolKernel *k;
  char *c;
  int i,j, bracketcount, num, slen;
  char *thisfunction = "extractKernels";

  if ((c = strstr(s, ":=")) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot find ':=' in rate/kernel definition.\n",thisfunction);
    return(NULL);
  }
  i=2;
  bracketcount=0;
  num=0;
  while ((c[i] != ';') && (c[i] != '\0')) {
    i++;
    if (c[i] == '[') 
      bracketcount++;
    else if (c[i] == ']') {
      bracketcount--;
      num++;
    }
  }
  
  if (c[i] != ';') {
    fprintf(stderr,"ERROR (%s): check rate/kernel definition, couldn't find ';' character in the end.\n",thisfunction);
    return(NULL);
  }
  
  if (bracketcount != 0) {
    fprintf(stderr,"ERROR (%s): check rate/kernel definition, brackets '[]' don't match (bracketcount %d).\n",thisfunction, bracketcount);
    return(NULL);
  }

  if (num == 0) {
    /* constant rate, set coord2 */
    num = 1;
    *num_kernels = num;

    if ((k = (SymbolKernel *) malloc (num * sizeof(SymbolKernel))) == NULL) {
      fprintf(stderr,"ERROR (%s): cannot malloc %d SymbolKernels.\n",thisfunction,num);
      perror(""); return(NULL);
    }
    slen = i-2;
    if ((k[0].name = strndup(c+2,slen)) == NULL) {
      fprintf(stderr,"ERROR (%s): cannot malloc %d chars.\n",thisfunction,slen);
      perror(""); return(NULL);
    }

    /* get coordinate argument from the beginning */
    i=0;
    while ((s[i] != ',') && (s[i] != ']')) {
      i++;
    }
    /* remove trailing '_' in symbol names */
    if ((i>0) && (s[i-1] == '_')) i--;
    if ((k[0].coord2 = strndup(s,i)) == NULL) {
      fprintf(stderr,"ERROR (%s): cannot malloc %d chars.\n",thisfunction,i);
      perror(""); return(NULL);
    }
    k[0].coord1 = NULL;
  }
  else {
    /* one or more kernel functions a[..] */
    *num_kernels = num;

    if ((k = (SymbolKernel *) malloc (num * sizeof(SymbolKernel))) == NULL) {
      fprintf(stderr,"ERROR (%s): cannot malloc %d SymbolKernels.\n",thisfunction,num);
      perror(""); return(NULL);
    }

    i=2;
    j=i;
    num=0;
    while (c[i] != ';') {
      if (c[i] == '[') {
	slen = i-j;
	if ((k[num].name = strndup(c+j,slen)) == NULL) {
	  fprintf(stderr,"ERROR (%s): cannot malloc %d chars.\n",thisfunction,slen);
	  perror(""); return(NULL);
	}
	j=i+1;
      }
      else if (c[i] == '-') {
	slen = i-j;
	if ((k[num].coord1 = strndup(c+j,slen)) == NULL) {
	  fprintf(stderr,"ERROR (%s): cannot malloc %d chars.\n",thisfunction,slen);
	  perror(""); return(NULL);
	}
	j=i+1;
      }
      else if (c[i] == ']') {
	slen = i-j;
	if ((k[num].coord2 = strndup(c+j,slen)) == NULL) {
	  fprintf(stderr,"ERROR (%s): cannot malloc %d chars.\n",thisfunction,slen);
	  perror(""); return(NULL);
	}
	j=i+1;
	num++;
      }
      i++;
    } 
  }
  
  return(k);
}

int getKernelType(SymbolKernel *k) {  

  if ( ((k->ptype1 == CATALYST) || (k->ptype1 == REACTANT)) && ((k->ptype2 == CATALYST) || (k->ptype2 == REACTANT)) )
    return (CONNECTIVITY_KERNEL);

  if ( ((k->ptype1 != NONE) && (k->ptype2 == NONE)) || ((k->ptype1 == NONE) && (k->ptype2 != NONE)) )
    return (OTHER_KERNEL);
  
  if ((k->ptype1 == PRODUCT) || (k->ptype2 == PRODUCT)) 
    return (OTHER_KERNEL);

  
  return(0);
}

int getProcessSourceArgIndex(ProcessModule *p) {
  /* define process source: one species from reactants or catalysts which has connections to all other inputs */
  /* source_argindex is the index of the source species in process call argument list */

  int i,j,ind,found,ok, num_connectivity_kernels;
  char *source_species, *coord;
  char *thisfunction = "getProcessSourceArgIndex";

  if (p->num_catalysts + p->num_reactants + p->num_products == 0) {
    fprintf(stderr,"ERROR (%s): process has no catalysts, reactants, or products.\n",thisfunction);
    return(-2);
  }
  if (p->num_catalysts + p->num_reactants == 0) {  
    return(-1);    
  }

  /* prefer reactants, so try them first */

  found=0;
  if (p->num_reactants > 0) {
    for (j=0; (j<p->num_reactants) && !found; j++) {
      source_species = p->reactants[j].species;
      coord = p->reactants[j].coord;
      ok=1;
      for (i=0; (i<p->num_kernels) && ok; i++) {
	if ((p->kernels[i].coord1 != NULL) && (p->kernels[i].coord2 != NULL)) {
	  if (strcmp(coord, p->kernels[i].coord1) && strcmp(coord, p->kernels[i].coord2))
	    ok=0;
	}
      }
      if (ok)
	found=1;
    }
  }

  if (!found) {
    /* look in catalysts list */
    if (p->num_catalysts > 0) {
      for (j=0; (j<p->num_catalysts) && !found; j++) {
	source_species = p->catalysts[j].species;
	coord = p->catalysts[j].coord;
	ok=1;
	for (i=0; (i<p->num_kernels) && ok; i++) {
	  if ((p->kernels[i].coord1 != NULL) && (p->kernels[i].coord2 != NULL)) {
	    if (strcmp(coord, p->kernels[i].coord1) && strcmp(coord, p->kernels[i].coord2))
	      ok=0;
	  }
	}
	if (ok)
	  found=1;
      }
    }
  }

  if (!found) {
    /* process is not "central-point" process */
    fprintf(stderr,"ERROR (%s): couldn't find a central point in process input (if process has input (Reactants or Catalysts), it must have a point which has connectivity to all other points of the process).\n",thisfunction);
    return(-2);    
  }

  ind=-1;
  for (i=0; (i<p->num_args) && (ind<0); i++) {
    if (!strcmp(source_species,p->args[i]))
      ind=i;
  }

  return (ind);
}

int **getSwapSpeciesArray(ProcessModule *p, int *alen) {
  /* get species pairs where Reactant and Product coordinates are the same */
  /* return indices of the source species in process call argument list */

  int i,j,k,num,ok,found,**a;
  char *thisfunction = "getSwapSpeciesArray";


  num=0;
  for (i=0; i<p->num_reactants; i++) {
    for (found=0, j=0; (j<p->num_products) && !found; j++) {
      if (!strcmp(p->reactants[i].coord, p->products[j].coord))
	found=1;
    }
    if (found)
      num++;
  }
  
  if (num == 0) {
    *alen=0;
    return(NULL);
  }

  if ((a = (int **) malloc (num * sizeof(int *))) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot malloc %d*2 int array.\n",thisfunction,num);
    perror(""); *alen=-1; return(NULL);
  }
  for (i=0; i<num; i++) {
    if ((a[i] = (int *) malloc (2 * sizeof(int))) == NULL) {
      fprintf(stderr,"ERROR (%s): cannot malloc %d*2 int array.\n",thisfunction,num);
      perror(""); *alen=-1; return(NULL);
    }
  }

  num=0;
  for (i=0; i<p->num_reactants; i++) {
    for (found=0, j=0; (j<p->num_products) && !found; j++) {
      if (!strcmp(p->reactants[i].coord, p->products[j].coord)) {
	/* find argindex for reactants[i] and products[j] */
	for (ok=0, k=0; (k<p->num_args) && !ok; k++) {
	  if (!strcmp(p->reactants[i].species, p->args[k])) {
	    a[num][0] = k;
	    ok=1;
	  }
	}
	if (!ok) {
	  fprintf(stderr,"ERROR (%s): cannot find Reactant species '%s' in argument list.\n",thisfunction,p->reactants[i].species);
	  *alen=-1; return(NULL);
	}

	for (ok=0, k=0; (k<p->num_args) && !ok; k++) {
	  if (!strcmp(p->products[j].species, p->args[k])) {
	    a[num][1] = k;
	    ok=1;
	  }
	}
	if (!ok) {
	  fprintf(stderr,"ERROR (%s): cannot find Product species '%s' in argument list.\n",thisfunction,p->reactants[i].species);
	  *alen=-1; return(NULL);
	}
	found=1;
	num++;
      }
    }
  }

  *alen = num;
  return (a);
}

int *getArgIndexSpeciesPointCount(SymbolPoint *sp, int num_sp, char **args, int num_args) {
  int i,j,k,*count;
  char *thisfunction = "getArgIndexSpeciesPointCount";

  if ((count = (int *) malloc (num_args * sizeof(int))) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot malloc %d ints.\n",thisfunction,num_args);
    perror(""); return(NULL);
  }
  for (j=0; j<num_args; j++)
    count[j]=0;
  
  /* find argindex for each species and update count */

  for (i=0; i<num_sp; i++) {
    k=-1;
    for (k=-1,j=0; (j<num_args) && (k<0); j++) {
      if (!strcmp(sp[i].species, args[j]))
	k=j;
    }
    if (k<0) {
      fprintf(stderr,"ERROR (%s): cannot find species '%s' (item %d/%d) in process argument list (length %d).\n",thisfunction,sp[i].species, i+1, num_sp, num_args);
      return(NULL);
    }
    count[k]++;
  }
  
  return (count);
}

int getKernelArgumentIndices(ProcessModule *pm, char *coord1, char *coord2, int *i1, int *i2, int *o1, int *o2) {
  int i;
  char *s1=NULL,*s2=NULL;

  /* pm contains {species,coord} as character strings, return species argument index based on coord */
  /* i1,i2 if found in input (Reactants,Catalysts), o1,o2 if found in output (Products) */

  *i1 = -1;
  *i2 = -1;
  *o1 = -1;
  *o2 = -1;
  if (coord2 == NULL) 
    return (-1);

  if (coord1 == NULL) {
    /* constant kernel, use only coord2 */
    for (i=0; i<pm->num_catalysts; i++) 
      if (!strcmp(pm->catalysts[i].coord, coord2))
	s2 = pm->catalysts[i].species;

    for (i=0; i<pm->num_reactants; i++) 
      if (!strcmp(pm->reactants[i].coord, coord2)) 
	s2 = pm->reactants[i].species;

    if (s2 != NULL) {
      for (i=0; i<pm->num_args; i++)
	if (!strcmp(pm->args[i], s2))
	  *i2 = i;
    }

    for (i=0; i<pm->num_products; i++) 
      if (!strcmp(pm->products[i].coord, coord2))
	s2 = pm->products[i].species;

    if (s2 != NULL) {
      for (i=0; i<pm->num_args; i++)
	if (!strcmp(pm->args[i], s2))
	  *o2 = i;
    }
    
    if ((*i2 < 0) && (*o2 < 0))
      return(-1);
    
    return (0);
  }

  /* kernel with two arguments */

  for (i=0; i<pm->num_catalysts; i++) {
    if (!strcmp(pm->catalysts[i].coord, coord1))
      s1 = pm->catalysts[i].species;
    if (!strcmp(pm->catalysts[i].coord, coord2))
      s2 = pm->catalysts[i].species;
  }

  for (i=0; i<pm->num_reactants; i++) {
    if (!strcmp(pm->reactants[i].coord, coord1)) 
      s1 = pm->reactants[i].species;
    if (!strcmp(pm->reactants[i].coord, coord2)) 
      s2 = pm->reactants[i].species;
  }

  if (s1 != NULL) {
    for (i=0; i<pm->num_args; i++) 
      if (!strcmp(pm->args[i], s1))
	*i1 = i;
  }
  if (s2 != NULL) {
    for (i=0; i<pm->num_args; i++) 
      if (!strcmp(pm->args[i], s2))
	*i2 = i;
  }
  
  for (i=0; i<pm->num_products; i++) {
    if (!strcmp(pm->products[i].coord, coord1))
      s1 = pm->products[i].species;
    if (!strcmp(pm->products[i].coord, coord2))
      s2 = pm->products[i].species;
  }

  if (s1 != NULL) {
    for (i=0; i<pm->num_args; i++) 
      if (!strcmp(pm->args[i], s1))
	*o1 = i;
  }
  if (s2 != NULL) {
    for (i=0; i<pm->num_args; i++) 
      if (!strcmp(pm->args[i], s2))
	*o2 = i;
  }

  if (((*i1 < 0) && (*o1 < 0)) || ((*i2 < 0) && (*o2 < 0)))
    return(-1);

  return (0);
}

int getProcessType(ProcessModule *p) {
  int num_wsckernels,num_ckernels, i, j;
  int i1,i2,o1,o2;
  char *thisfunction = "getProcessType";

  if ((p->num_reactants == 0) && (p->num_catalysts == 0))
    return (RATE_PER_AREA);

  if (p->num_reactants + p->num_catalysts <= 1)
    return (RATE_PER_SPECIES_COUNT);

  num_wsckernels=0;
  num_ckernels=0;
  for (i=0; i<p->num_kernels; i++) {
    if (p->kernels[i].type == CONNECTIVITY_KERNEL) {
      num_ckernels++;
      if (getKernelArgumentIndices(p, p->kernels[i].coord1, p->kernels[i].coord2, &i1, &i2, &o1, &o2) < 0) {
	fprintf(stderr,"ERROR (%s): cannot get kernel argument indices for kernel %d/%d.\n",thisfunction,i+1,p->num_kernels);
	return(-1);
      }
      if ((i1 < 0) && (i2 < 0)) {
	fprintf(stderr,"ERROR (%s): cannot get kernel argument indices from input (Catalysts,Reactants) for kernel %d/%d.\n",thisfunction,i+1,p->num_kernels);
	return(-1);
      }
      if (!strcmp(p->args[i1],p->args[i2]))
	num_wsckernels++;
    }
  }
  if (num_wsckernels > 1) {
    fprintf(stderr,"ERROR (%s): max 1 within-species connectivity kernel allowed (found %d).\n",thisfunction,num_wsckernels);
    return(-1);
  }

  return (SPECIES_CONNECTIVITY);
}

int extractProcessDefinition(char **names, ProcessModule *p) {
  int i,j,patlen,found;
  char *c, *pat;
  char *thisfunction = "extractProcessDefinition";

  if ((p->name = strdup(names[0])) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot copy name of ProcessModule.\n",thisfunction);
    perror(""); return(-1);
  }
  
  /* get list of arguments */

  if ((p->args = extractFunctionArguments(names[1], &(p->num_args))) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot extract arguments of ProcessModule.\n",thisfunction);
    return (-1);
  }

  /* get Products={...} */

  pat = "Products=";
  patlen = strlen(pat);
  if ((c = strstr(names[3], pat)) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot find pattern '%s'.\n",thisfunction,pat);
    return(-1);
  }
  if ((p->products = extractPointList(c+patlen, &(p->num_products))) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot extract arguments {...} for '%s'.\n",thisfunction,pat);
    return(-1);
  }

  /* get Reactants={...} */

  pat = "Reactants=";
  patlen = strlen(pat);
  if ((c = strstr(names[3], pat)) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot find pattern '%s'.\n",thisfunction,pat);
    return(-1);
  }
  if ((p->reactants = extractPointList(c+patlen, &(p->num_reactants))) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot extract arguments {...} for '%s'.\n",thisfunction,pat);
    return(-1);
  }

  /* get Catalysts={...} */

  pat = "Catalysts=";
  patlen = strlen(pat);
  if ((c = strstr(names[3], pat)) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot find pattern '%s'.\n",thisfunction,pat);
    return(-1);
  }
  if ((p->catalysts = extractPointList(c+patlen, &(p->num_catalysts))) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot extract arguments {...} for '%s'.\n",thisfunction,pat);
    return(-1);
  }

  /* calculate number of points per species for each process point type */

  if ((p->argind_catalyst_count = getArgIndexSpeciesPointCount(p->catalysts, p->num_catalysts, p->args, p->num_args)) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot calculate number of points per species for Catalysts.\n",thisfunction);
    return(-1);
  }
  if ((p->argind_reactant_count = getArgIndexSpeciesPointCount(p->reactants, p->num_reactants, p->args, p->num_args)) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot calculate number of points per species for Reactants.\n",thisfunction);
    return(-1);
  }
  if ((p->argind_product_count = getArgIndexSpeciesPointCount(p->products, p->num_products, p->args, p->num_args)) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot calculate number of points per species for Products.\n",thisfunction);
    return(-1);
  }

  /* get function[args] := xyz; where xyz is constant or kernel1[args1]..kernelN[argsN] */

  pat = "function[";
  patlen = strlen(pat);
  if ((c = strstr(names[3], pat)) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot find pattern '%s'.\n",thisfunction,pat);
    return(-1);
  }
  if ((p->kernels = extractKernels(c+patlen, &(p->num_kernels))) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot extract arguments for '%s'.\n",thisfunction,pat);
    return(-1);
  }

  /* get position of kernel in Process function call */

  for (i=0; i<p->num_kernels; i++) {
    p->kernels[i].argindex = -1;
    for (j=0; j<p->num_args; j++) {
      if (!strcmp(p->kernels[i].name, p->args[j])) {
	p->kernels[i].argindex = j;
      }
    }
    if (p->kernels[i].argindex < 0) {
      fprintf(stderr,"ERROR (%s): cannot find kernel '%s' in ProcessDefinition arguments.\n",thisfunction,p->kernels[i].name);
      return(-1);
    }
  }
  
  /* set ptypes for kernel arguments */
  for (i=0; i<p->num_kernels; i++) {
    p->kernels[i].ptype1 = 0;
    p->kernels[i].ptype2 = 0;
    if (p->kernels[i].coord1 != NULL) {
      found=0;
      for (j=0; (j<p->num_catalysts) && !found; j++) {
	if (!strcmp(p->kernels[i].coord1, p->catalysts[j].coord)) {
	  p->kernels[i].ptype1 = CATALYST;
	  found=1;
	}
      }
      for (j=0; (j<p->num_reactants) && !found; j++) {
	if (!strcmp(p->kernels[i].coord1, p->reactants[j].coord)) {
	  p->kernels[i].ptype1 = REACTANT;
	  found=1;
	}
      }
      for (j=0; (j<p->num_products) && !found; j++) {
	if (!strcmp(p->kernels[i].coord1, p->products[j].coord)) {
	  p->kernels[i].ptype1 = PRODUCT;
	  found=1;
	}	  
      }
    }
    if (p->kernels[i].coord2 != NULL) {
      found=0;
      for (j=0; (j<p->num_catalysts) && !found; j++) {
	if (!strcmp(p->kernels[i].coord2, p->catalysts[j].coord)) {
	  p->kernels[i].ptype2 = CATALYST;
	  found=1;
	}
      }
      for (j=0; (j<p->num_reactants) && !found; j++) {
	if (!strcmp(p->kernels[i].coord2, p->reactants[j].coord)) {
	  p->kernels[i].ptype2 = REACTANT;
	  found=1;
	}
      }
      for (j=0; (j<p->num_products) && !found; j++) {
	if (!strcmp(p->kernels[i].coord2, p->products[j].coord)) {
	  p->kernels[i].ptype2 = PRODUCT;
	  found=1;
	}	  
      }
    }    
  }

  /* get kernel type */
  
  for (i=0; i<p->num_kernels; i++) {
    p->kernels[i].type = getKernelType(&(p->kernels[i]));
  }

  /* define process source: one species from catalysts or reactants */
  /* source_argindex is the index of the source species in process call argument list */
  /* kernel type must be defined before this */

  if ((p->source_argindex = getProcessSourceArgIndex(p)) < -1) {
    fprintf(stderr,"ERROR (%s): cannot assign process source.\n",thisfunction);
    return(-1);
  }

  p->same_coordinate_reactant_product_pair = getSwapSpeciesArray(p, &(p->num_same_coordinate_reactant_product_pairs));
  if (p->num_same_coordinate_reactant_product_pairs < 0) {
    fprintf(stderr,"ERROR (%s): cannot assign same_coordinate_reactant_product_pair array.\n",thisfunction);
    return(-1);
  }

  /* determine process type */

  if ((p->type = getProcessType(p)) < 0) {
    fprintf(stderr,"ERROR (%s): cannot getProcessType.\n",thisfunction);
    return(-1);
  }

  return (0);
}

ProcessModule *readProcessFile(char *filename, int *num_processes) {
  FILE *fp;
  char **names;
  ProcessModule *p;
  int count, i;
  char *thisfunction = "readProcessFile";
  
  if ((names = (char **) malloc (4*sizeof(char *))) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot malloc 4*%d char.\n",thisfunction, MAXSTRINGLEN);
    perror(""); exit(-1);
  }
  for (i=0; i<4; i++) {
    if ((names[i] = (char *) malloc (MAXSTRINGLEN*sizeof(char))) == NULL) {
      fprintf(stderr,"ERROR (%s): cannot malloc 4*%d char.\n",thisfunction, MAXSTRINGLEN);
      perror(""); exit(-1);
    }
  }

  if ((fp = fopen(filename,"r")) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot read file '%s'.\n",thisfunction,filename);
    perror(""); exit(-1);
  }

  count=0;
  while (!feof(fp)) {
    if ((i = readProcessDefinition(fp, names, filename)) >0) {
      count++;
      if (strcmp(names[2], ":=Module")) {
	fprintf(stderr,"ERROR (%s): check processDefinitionFile '%s', Module %d.\n", thisfunction,filename,count);
	fprintf(stderr,"expected ':=Module' but found '%s'.\n",names[2]);
	fclose(fp); exit(-1);
      }
    }
  }
  if ((i == 0) && (names[0][0] != '\0')) {
    fprintf(stderr,"ERROR (%s): check processDefinitionFile '%s', read %d Modules after which there were extra characters.\n", thisfunction,filename,count);
    fclose(fp); exit(-1);
  }

  if ((p = (ProcessModule *) malloc (count * sizeof(ProcessModule))) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot malloc %d ProcessModules.\n",thisfunction,count);
    perror(""); exit(-1);
  }

  rewind(fp);

  i=0;
  while (!feof(fp)) {
    if (readProcessDefinition(fp, names, filename)) {
      if (extractProcessDefinition(names,&(p[i])) < 0) {
	fprintf(stderr,"ProcessModule %d/%d '%s'.\n",i+1,count,names[0]);
	perror(""); exit(-1);	
      }
      i++;
    }
  }
  fclose(fp);
  
  *num_processes = count;
  return(p);
}

int printProcessModules (ProcessModule *p, int num_processes) {
  int k,i;

  for (k=0; k<num_processes; k++) {
    printf("ProcessDefinition %d/%d '%s', type %d '%s':\n",k+1,num_processes,p[k].name, p[k].type, PROCESS_TYPE[p[k].type]);
    printf("  num_arguments %d:",p[k].num_args);
    for (i=0; i<p[k].num_args; i++) 
      printf(" '%s'",p[k].args[i]);
    printf(" (source argindex %d)\n", p[k].source_argindex);
    printf("  num_products %d:",p[k].num_products);
    for (i=0; i<p[k].num_products; i++) 
      printf(" (%s,%s)",p[k].products[i].species,p[k].products[i].coord);
    printf("\n");
    printf("  num_reactants %d:",p[k].num_reactants);
    for (i=0; i<p[k].num_reactants; i++) 
      printf(" (%s,%s)",p[k].reactants[i].species,p[k].reactants[i].coord);
    printf("\n");
    printf("  num_catalysts %d:",p[k].num_catalysts);
    for (i=0; i<p[k].num_catalysts; i++) 
      printf(" (%s,%s)",p[k].catalysts[i].species,p[k].catalysts[i].coord);
    printf("\n");
    printf("  num_kernels %d:",p[k].num_kernels);
    for (i=0; i<p[k].num_kernels; i++) 
      printf(" %s[%s %s,%s %s](%s)",p[k].kernels[i].name,p[k].kernels[i].coord1,POINT_TYPE[p[k].kernels[i].ptype1],p[k].kernels[i].coord2,POINT_TYPE[p[k].kernels[i].ptype2], KERNEL_TYPE[p[k].kernels[i].type]);
    printf("\n");

    printf("  argind_product_count :");
    for (i=0; i<p[k].num_args; i++) 
      printf(" %d",p[k].argind_product_count[i]);
    printf("\n");
    printf("  argind_reactant_count:");
    for (i=0; i<p[k].num_args; i++) 
      printf(" %d",p[k].argind_reactant_count[i]);
    printf("\n");
    printf("  argind_catalyst_count:");
    for (i=0; i<p[k].num_args; i++) 
      printf(" %d",p[k].argind_catalyst_count[i]);
    printf("\n");
    printf("  %d same coordinate reactant,product pairs:",p[k].num_same_coordinate_reactant_product_pairs);
    for (i=0; i<p[k].num_same_coordinate_reactant_product_pairs; i++)
      printf(" (%d,%d)",p[k].same_coordinate_reactant_product_pair[i][0],p[k].same_coordinate_reactant_product_pair[i][1]);
    printf("\n");
  }
  return(0);
}

ModelComponent *readModelFile(char *filename, int *num_components) {
  FILE *fp;
  char line[MAXSTRINGLEN];
  ModelComponent *m;
  int i,c,num, bracketcount;
  char *thisfunction = "readModelFile";

  if ((fp = fopen(filename,"r")) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot open file '%s'.\n",thisfunction,filename);
    perror(""); exit(-1);
  }

  /* count number of process functions */

  num=0;
  bracketcount=0;
  while ((c = fgetc(fp)) != EOF) {
    if (c == '[') {
      bracketcount++;
    }
    else if (c == ']') {
      bracketcount--;
      if (bracketcount == 0) {
	num++;
      }
    }
  }

  if ((num < 1) || (bracketcount != 0)) {
    fprintf(stderr,"ERROR (%s): check model file '%s' (%d process components, bracketcount %d).\n", thisfunction, filename, num, bracketcount);
    fclose(fp); exit(-1);
  }

  *num_components = num;
  if ((m = (ModelComponent *) malloc (num * sizeof(ModelComponent))) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot malloc %d ModelComponents, file '%s'.\n",thisfunction,num,filename);
    perror(""); exit(-1);
  }

  /* read process functions: name[args] */

  rewind(fp);
  num=0;
  bracketcount=0;
  i=0;
  while ((c = fgetc(fp)) != EOF) {    
    if (!isspace(c)) {
      if (c == '[') { 
	bracketcount++;
	if (bracketcount == 1) {
	  if ((m[num].name = strndup(line,i)) == NULL) {
	    fprintf(stderr,"ERROR (%s): cannot malloc %d chars for function name %d/%d, file '%s'.\n",thisfunction,i,num+1, *num_components,filename);
	    perror(""); exit(-1);
	  }
	  i=0;
	}
	else {
	  line[i] = c;
	  i++;
	}
      }
      else if (c == ']') {
	bracketcount--;
	if (bracketcount == 0) {
	  line[i] = '\0';
	  if ((m[num].args = extractFunctionArguments(line, &(m[num].num_args))) == NULL) {
	    fprintf(stderr,"ERROR (%s): cannot extract arguments for function %d/%d, file '%s'.\n",thisfunction,num+1,*num_components,filename);
	  }
	  num++;
	  i=0;
	}
	else {
	  line[i] = c;
	  i++;
	}
      }
      else {
	line[i] = c;
	i++;
      }

      if (i >= MAXSTRINGLEN) {
	fprintf(stderr,"ERROR (%s): Model component string exceeds maximum allowed characters (%d). Please check the input file '%s'.\n",thisfunction, MAXSTRINGLEN, filename);
	fprintf(stderr,"Error occurred in function %d/%d: '",num+1,*num_components);
	int j;
	for (j=0; j<i; j++) {
	  fprintf(stderr,"%c",line[j]);
	}
	fprintf(stderr,"...'\n");
	fprintf(stderr,"If the input file is ok, increase the value of MAXSTRINGLEN in line '#define MAXSTRINGLEN %d' in c-code and re-compile. Sorry for the inconvenience.\n",MAXSTRINGLEN);
	fclose(fp); exit(-1);
      }
    }
  }

  fclose(fp);
  return (m);
}

int printModelFunctions(ModelComponent *m, int num_components) {
  int i,k;

  for (k=0; k<num_components; k++) {
    printf("ModelFunction %d/%d: '%s'\n",k+1,num_components,m[k].name);
    printf("  num_arguments %d:",m[k].num_args);
    for (i=0; i<m[k].num_args; i++) {
      printf(" '%s'",m[k].args[i]);
    }
    printf("\n");
  }
  return (0);
}

ProcessModule *findProcessModule(ProcessModule *pm, int num_pm, char *name) {
  int i;
  ProcessModule *ptr=NULL;
  for (i=0; (i<num_pm) && (ptr==NULL); i++) {
    if (!strcmp(pm[i].name, name))
      ptr = &(pm[i]);
  }
  return(ptr);
}

Process *createProcessArray(ProcessModule *pm, int num_pm, ModelComponent *m, int num_m, Kernel **c_ptr, int *num_ckernels, Kernel **o_ptr, int *num_okernels) {
  Kernel *ckernels, *rkernels, *okernels;
  ProcessModule *pk;
  Process *pa;
  int i,j,kc,ko,i1,i2,o1,o2,s1,s2,t,v,found,num,swap;
  char *kname;
  char *thisfunction = "createProcessArray";

  /* number of different processes is the same as number of defined model components */

  if ((pa = (Process *) malloc (num_m * sizeof(Process))) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot malloc %d Processes.\n",thisfunction,num_m);
    perror(""); return(NULL);
  }

  /* some processes may use several kernels and some kernels may be used by more than one process */

  /* count all kernels */

  num=0;
  for (i=0; i<num_m; i++) {
    if ((pk = findProcessModule(pm, num_pm, m[i].name)) == NULL) {
      fprintf(stderr,"ERROR (%s): cannot find ProcessDefinition for '%s'.\n",thisfunction,m[i].name);
      return (NULL);
    }
    num += pk->num_kernels;
  }

  /* number of unique c and o kernels is <= all kernels, but this is ok for malloc */

  if ((ckernels = (Kernel *) malloc (num * sizeof(Kernel))) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot malloc %d Kernels.\n",thisfunction,num);
    perror(""); return(NULL);
  }
  if ((okernels = (Kernel *) malloc (num * sizeof(Kernel))) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot malloc %d Kernels.\n",thisfunction,num);
    perror(""); return(NULL);
  }

  /* get unique kernels in model definition */
  /* species arguments (coord->symbol) from pm and parameters of kernel function from m */

  kc=0;
  ko=0;
  for (i=0; i<num_m; i++) {
    if ((pk = findProcessModule(pm, num_pm, m[i].name)) == NULL) {
      fprintf(stderr,"ERROR (%s): cannot find ProcessDefinition for '%s'.\n",thisfunction,m[i].name);
      return (NULL);
    }
    pa[i].pm = pk;
    if ((pa[i].cind = (int *) malloc (pk->num_kernels * sizeof(int))) == NULL) {
      fprintf(stderr,"ERROR (%s): cannot malloc %d ints for cind of Process %d/%d.\n",thisfunction,pk->num_kernels,i+1,num_m);
      perror(""); return(NULL);
    }
    if ((pa[i].oind = (int *) malloc (pk->num_kernels * sizeof(int))) == NULL) {
      fprintf(stderr,"ERROR (%s): cannot malloc %d ints for oind of Process %d/%d.\n",thisfunction,pk->num_kernels,i+1,num_m);
      perror(""); return(NULL);
    }
    pa[i].num_cind = 0; 
    pa[i].num_oind = 0;

    if ((pa[i].ctype1 = (int *) malloc (pk->num_kernels * sizeof(int))) == NULL) {
      fprintf(stderr,"ERROR (%s): cannot malloc %d ints for ptype1 of Process %d/%d.\n",thisfunction,pk->num_kernels,i+1,num_m);
      perror(""); return(NULL);
    }
    if ((pa[i].ctype2 = (int *) malloc (pk->num_kernels * sizeof(int))) == NULL) {
      fprintf(stderr,"ERROR (%s): cannot malloc %d ints for ptype2 of Process %d/%d.\n",thisfunction,pk->num_kernels,i+1,num_m);
      perror(""); return(NULL);
    }
    if ((pa[i].otype1 = (int *) malloc (pk->num_kernels * sizeof(int))) == NULL) {
      fprintf(stderr,"ERROR (%s): cannot malloc %d ints for ptype1 of Process %d/%d.\n",thisfunction,pk->num_kernels,i+1,num_m);
      perror(""); return(NULL);
    }
    if ((pa[i].otype2 = (int *) malloc (pk->num_kernels * sizeof(int))) == NULL) {
      fprintf(stderr,"ERROR (%s): cannot malloc %d ints for ptype2 of Process %d/%d.\n",thisfunction,pk->num_kernels,i+1,num_m);
      perror(""); return(NULL);
    }

    if ((pa[i].arg_species = (int *) malloc (m[i].num_args * sizeof(int))) == NULL) {
      fprintf(stderr,"ERROR (%s): cannot malloc %d ints for arg_species of Process %d/%d.\n",thisfunction,m[i].num_args,i+1,num_m);
      perror(""); return(NULL);
    }
    pa[i].num_args = m[i].num_args;
    for (j=0; j<pa[i].num_args; j++) {
      for (found=0, t=0; (t<pk->num_kernels) && !found; t++) {
	if (j == pk->kernels[t].argindex)
	  found=1;
      }
      if (found)
	pa[i].arg_species[j] = -1;
      else 
	pa[i].arg_species[j] = atoi(m[i].args[j]);
    }

    /* set process source species */
    
    if (pk->source_argindex < 0)
      pa[i].source_species = -1;
    else
      pa[i].source_species = atoi(m[i].args[pk->source_argindex]);

    for (j=0; j<pk->num_kernels; j++) {
      if (getKernelArgumentIndices(pk, pk->kernels[j].coord1, pk->kernels[j].coord2, &i1, &i2, &o1, &o2) < 0) {
	fprintf(stderr,"ERROR (%s): cannot find Kernel %d/%d arguments (species symbol based on coord symbol) in ProcessDefinition '%s'.\n",thisfunction,j+1,pk->num_kernels,m[i].name);
	return (NULL);
      }
      /* model arguments are in the same order as in ProcessModule */

      /* choose c,o based on pk->kernels[j].type */

      if (i1 >= 0) 
	s1 = atoi(m[i].args[i1]);
      else if (o1 >= 0)
	s1 = atoi(m[i].args[o1]);
      else
	s1 = -1;

      if (i2 >= 0) 
	s2 = atoi(m[i].args[i2]);
      else if (o2 >= 0) 
	s2 = atoi(m[i].args[o2]);
      else
	s2 = -1;

      /* use sorted indices (so that ukernel doesn't depend on argument order) */
      swap=0;
      if (s2 < s1) {
	t=s1; s1=s2; s2=t;
	swap=1;
      }
      kname = m[i].args[pk->kernels[j].argindex];

      /* check if already exists, k is num of currently existing unique kernels */
      /* note: kernel function can be identical but used for type 1 or type 2 purpose */

      if (pk->kernels[j].type == CONNECTIVITY_KERNEL) {
	found=0;
	for (t=0; (t<kc) && !found; t++) {
	  if (!strcmp(ckernels[t].fullname, kname) && (ckernels[t].species1 == s1) && (ckernels[t].species2 == s2)) {
	    found=1;
	    v = pa[i].num_cind;
	    ckernels[t].num_pind++;
	    pa[i].cind[v] = t;
	    if (swap) {
	      pa[i].ctype1[v] = pk->kernels[j].ptype2;
	      pa[i].ctype2[v] = pk->kernels[j].ptype1;
	    }
	    else {
	      pa[i].ctype1[v] = pk->kernels[j].ptype1;
	      pa[i].ctype2[v] = pk->kernels[j].ptype2;
	    }
	    pa[i].num_cind++;
	  }
	}
	if (!found) {
	  ckernels[kc].fullname = strdup(kname);
	  ckernels[kc].species1 = s1;
	  ckernels[kc].species2 = s2;
	  ckernels[kc].type = pk->kernels[j].type;	
	  ckernels[kc].num_pind = 1;
	  v = pa[i].num_cind;
	  pa[i].cind[v] = kc;
	  if (swap) {
	    pa[i].ctype1[v] = pk->kernels[j].ptype2;
	    pa[i].ctype2[v] = pk->kernels[j].ptype1;
	  }
	  else {
	    pa[i].ctype1[v] = pk->kernels[j].ptype1;
	    pa[i].ctype2[v] = pk->kernels[j].ptype2;
	  }
	  pa[i].num_cind++;
	  kc++;
	}
      }
      else if (pk->kernels[j].type == OTHER_KERNEL) {
	found=0;
	for (t=0; (t<ko) && !found; t++) {
	  if (!strcmp(okernels[t].fullname, kname) && (okernels[t].species1 == s1) && (okernels[t].species2 == s2)) {
	    found=1;
	    v = pa[i].num_oind;
	    okernels[t].num_pind++;
	    pa[i].oind[v] = t;
	    if (swap) {
	      pa[i].otype1[v] = pk->kernels[j].ptype2;
	      pa[i].otype2[v] = pk->kernels[j].ptype1;
	    }
	    else {
	      pa[i].otype1[v] = pk->kernels[j].ptype1;
	      pa[i].otype2[v] = pk->kernels[j].ptype2;
	    }
	    pa[i].num_oind++;
	  }
	}
	if (!found) {
	  okernels[ko].fullname = strdup(kname);
	  okernels[ko].species1 = s1;
	  okernels[ko].species2 = s2;
	  okernels[ko].type = pk->kernels[j].type;	
	  okernels[ko].num_pind = 1;
	  v = pa[i].num_oind;
	  pa[i].oind[v] = ko;
	  if (swap) {
	    pa[i].otype1[v] = pk->kernels[j].ptype2;
	    pa[i].otype2[v] = pk->kernels[j].ptype1;
	  }
	  else {
	    pa[i].otype1[v] = pk->kernels[j].ptype1;
	    pa[i].otype2[v] = pk->kernels[j].ptype2;
	  }
	  pa[i].num_oind++;
	  ko++;
	}
      }
      else {
	fprintf(stderr,"ERROR (%s): unrecognized kernel type %d, kernel %d/%d '%s', process %d/%d.\n",thisfunction,pk->kernels[j].type,j+1,pk->num_kernels,kname,i+1,num_m);
	return(NULL);
      }

      /* printf("model %d/%d kname '%s' %d %d, already exists: %d, %d unique\n",i+1,num_m,kname,s1,s2,found,k); */      
    }
    
  }

  /* create links (index array) from ukernels to processes */
  /* kc/ko is number of unique connectivity/other kernels in all models */

  for (t=0; t<kc; t++) {
    if ((ckernels[t].pind = (int *) malloc (ckernels[t].num_pind * sizeof(int))) == NULL) {
      fprintf(stderr,"ERROR (%s): cannot malloc %d ints for ckernel %d/%d.\n",thisfunction,ckernels[t].num_pind,t+1,kc);
      perror(""); return(NULL);
    }

    /* scan ckernel indices of processes */
    s1=0;
    for (i=0; i<num_m; i++) {
      for (j=0; j<pa[i].num_cind; j++) {
	if (pa[i].cind[j] == t) {
	  ckernels[t].pind[s1] = i;
	  s1++;
	}
      }
    }    
  }
  for (t=0; t<ko; t++) {
    if ((okernels[t].pind = (int *) malloc (okernels[t].num_pind * sizeof(int))) == NULL) {
      fprintf(stderr,"ERROR (%s): cannot malloc %d ints for okernel %d/%d.\n",thisfunction,okernels[t].num_pind,t+1,ko);
      perror(""); return(NULL);
    }

    /* scan ckernel indices of processes */
    s1=0;
    for (i=0; i<num_m; i++) {
      for (j=0; j<pa[i].num_oind; j++) {
	if (pa[i].oind[j] == t) {
	  okernels[t].pind[s1] = i;
	  s1++;
	}
      }
    }    
  }
  
  *num_ckernels = kc; 
  *c_ptr = ckernels;
  *num_okernels = ko; 
  *o_ptr = okernels;
  return(pa);
}

int printProcesses(Process *p, int num_processes, Kernel *ckernel, Kernel *okernel) {
  int i,j;
  Kernel *u;

  for (i=0; i<num_processes; i++) {
    printf("process %d/%d '%s' source species %d (%d ckernels, %d okernels) type: %s\n",i+1,num_processes,p[i].pm->name,p[i].source_species,p[i].num_cind,p[i].num_oind, PROCESS_TYPE[p[i].pm->type]);
    printf("  arg_species:");
    for (j=0; j<p[i].num_args; j++)
      printf(" %d",p[i].arg_species[j]);
    printf("\n");
    for (j=0; j<p[i].num_cind; j++) {
      u = &(ckernel[p[i].cind[j]]);
      printf("  c '%s' species (%d,%d) ptypes (%s,%s)\n",u->fullname, u->species1, u->species2, POINT_TYPE[p[i].ctype1[j]], POINT_TYPE[p[i].ctype2[j]]);
    }
    for (j=0; j<p[i].num_oind; j++) {
      u = &(okernel[p[i].oind[j]]);
      printf("  o '%s' species (%d,%d) ptypes (%s,%s)\n",u->fullname, u->species1, u->species2, POINT_TYPE[p[i].otype1[j]], POINT_TYPE[p[i].otype2[j]]);
    }
  }
  return(0);
}

int assignKernelFamily(Kernel *kernels, int num_kernels, int dimension) {
  int i,j,k,found,num_args;
  char str[MAXSTRINGLEN], **args;
  char *thisfunction = "assignKernelFamily";

  for (k=0; k<num_kernels; k++) {
    if ((kernels[k].species1 < 0) || (kernels[k].species2 < 0)) {
      kernels[k].family = CONSTANT_KERNEL;
      kernels[k].integral = atof(kernels[k].fullname);
      kernels[k].maxdensity = kernels[k].integral;
      kernels[k].radius = -1.0;
      kernels[k].sigma = 0.0;
      kernels[k].halfivar = 0.0;
      kernels[k].trunceffect = 0.0;
    }
    else {
      /* find pattern kernelname[arg1,arg2,...,argN] */
      j=0;
      found=0;
      for (i=0; kernels[k].fullname[i] != '\0'; i++) {
	str[i] = kernels[k].fullname[i];
	if (str[i] == '[') {
	  str[i] = '\0';
	  if (!strcmp(str, "tophat")) {
	    kernels[k].family = TOPHAT_KERNEL;	    
	    found=1;
	  }
	  else if (!strcmp(str, "truncatedGaussian")) {
	    kernels[k].family = TRUNCATED_GAUSSIAN_KERNEL;
	    found=1;
	  }
	  else {
	    fprintf(stderr,"ERROR (%s): unrecognized kernel family '%s' (fullname '%s'), ukernel %d/%d.\n",thisfunction,str, kernels[k].fullname, k+1, num_kernels);
	    return (-1);
	  }
	  j=i+1;
	}
      }
      if (!found) {
	fprintf(stderr,"ERROR (%s): cannot find '[' in kernel name, fullname '%s', ukernel %d/%d.\n",thisfunction,kernels[k].fullname, k+1, num_kernels);
	return (-1);
      }

      num_args=0;
      for (i=0; kernels[k].fullname[j+i] != '\0'; i++) {
	str[i] = kernels[k].fullname[j+i];
	if (str[i] == ']') {
	  str[i] = '\0';
	  if ((args = extractFunctionArguments(str, &num_args)) == NULL) {
	    fprintf(stderr,"ERROR (%s): cannot extract kernel function arguments, fullname '%s', ukernel %d/%d.\n",thisfunction,kernels[k].fullname, k+1, num_kernels);
	  }
	}
      }
      if ((kernels[k].fullname[j+i] != '\0') || (num_args == 0)) {
	fprintf(stderr,"ERROR (%s): extra characters in kernel function (num_args %d), fullname '%s', ukernel %d/%d.\n",thisfunction,num_args,kernels[k].fullname, k+1, num_kernels);
      }
      
      if (kernels[k].family == TOPHAT_KERNEL) {
	if (num_args != 2) {
	  fprintf(stderr,"ERROR (%s): number of arguments (integral, radius) should be 2 for tophat kernel (was %d), fullname '%s', ukernel %d/%d.\n",thisfunction,num_args,kernels[k].fullname, k+1, num_kernels);
	  return(-1);
	}
	kernels[k].integral = atof(args[0]);
	kernels[k].radius = atof(args[1]);
	if (dimension == 1) {
	  kernels[k].maxdensity = kernels[k].integral/(2.0 * kernels[k].radius);
	}
	else if (dimension == 2) {
	  kernels[k].maxdensity = kernels[k].integral / (3.141593 * kernels[k].radius * kernels[k].radius);
	}
	else {
	  fprintf(stderr,"ERROR (%s): only 1D or 2D space supported (dimension was defined %d).\n",thisfunction,dimension);
	  return(-1);
	}

	kernels[k].sigma = 0.0;
	kernels[k].halfivar = 0.0;
	kernels[k].trunceffect = 0.0;
      }
      else if (kernels[k].family == TRUNCATED_GAUSSIAN_KERNEL) {
	if (num_args != 2) {
	  fprintf(stderr,"ERROR (%s): number of arguments (integral, sigma) should be 2 for truncatedGaussian kernel (was %d), fullname '%s', ukernel %d/%d.\n",thisfunction,num_args,kernels[k].fullname, k+1, num_kernels);
	  return(-1);
	}
	kernels[k].integral = atof(args[0]);
	kernels[k].sigma = atof(args[1]);
	kernels[k].halfivar = 0.5/(kernels[k].sigma*kernels[k].sigma);
	kernels[k].radius = 3.0 * kernels[k].sigma;
	/* maxdensity value is stored without the effect of truncated area of Gaussian */
	if (dimension == 1) {
	  /* normalizing factor: (2*pi*sigma^2)^-1/2 = 0.3989423/sigma */
	  /* area outside truncation -> uniform -3sigma...3sigma : 2*Integral_r=3sigma^Inf(Gaussian) / (2*3sigma) = 0.000449966/sigma */
	  kernels[k].maxdensity = kernels[k].integral*0.3989423/kernels[k].sigma;
	  kernels[k].trunceffect = kernels[k].integral*0.000449966/kernels[k].sigma;
	}
	else if (dimension == 2) {
	  /* normalizing factor: (2*pi*sigma^2)^-1 = 0.1591549/(sigma*sigma) */
	  /* area outside truncation -> uniform area within radius 3sigma : Integral_r=3sigma^Inf(2D-Gaussian) / (pi*(3sigma)^2) = 0.03536777 * exp(-1.5/sigma)/(sigma*sigma) */
	  kernels[k].maxdensity = kernels[k].integral*0.1591549/(kernels[k].sigma*kernels[k].sigma);
	  kernels[k].trunceffect = kernels[k].integral*0.03536777 * exp(-1.5/kernels[k].sigma)/(kernels[k].sigma*kernels[k].sigma);
	}
	else {
	  fprintf(stderr,"ERROR (%s): only 1D or 2D space supported (dimension was defined %d).\n",thisfunction,dimension);
	  return(-1);
	}
      }
      else {
	fprintf(stderr,"ERROR (%s): unrecognized kernel family %d, fullname '%s', ukernel %d/%d.\n",thisfunction,kernels[k].family,kernels[k].fullname, k+1, num_kernels);
	return(-1);
      }
    }
  }
  return (0);
}

int printUniqueKernels(Kernel *k, int num, Process *p) {
  int i,j;
  for (i=0; i<num; i++) {
    printf("Kernel %d/%d: '%s' with species (%d,%d) type: %s\n",i+1,num,k[i].fullname,k[i].species1,k[i].species2, KERNEL_TYPE[k[i].type]);
    printf("  integral %f, maxdensity %f, trunceffect %f, sigma %f, radius %f\n",k[i].integral, k[i].maxdensity, k[i].trunceffect, k[i].sigma, k[i].radius);
    printf("  num_processes %d:",k[i].num_pind);
    for (j=0; j<k[i].num_pind; j++) 
      printf(" '%s' (source species %d)",p[k[i].pind[j]].pm->name, p[k[i].pind[j]].source_species);
    printf("\n");
  }
 
  return (0);
}

Point *readPointFile(char *filename, int *num_points, int *dimension) {
  FILE *fp;
  Point *p;
  int i,num,linecount, count,prevcount,dim;
  char line[MAXSTRINGLEN], *tok;
  char *thisfunction = "readPointFile";

  if ((fp = fopen(filename,"r")) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot open file '%s'.\n",thisfunction, filename);
    perror(""); exit(-1);
  }
  count=0;
  linecount=0;
  num=0;
  prevcount=0;
  while (fgets(line, MAXSTRINGLEN, fp) != NULL) {
    linecount++;
    count=0;
    if (strtok(line," \t\r\n")) {
      count++;
      while (strtok(NULL," \t\r\n")) {
	count++;
      }
    }
    if (count >0) {
      num++;
      if (prevcount > 0) {
	if (count != prevcount) {
	  fprintf(stderr,"ERROR (%s): number of columns in two consequtive non-empty lines don't match (%d vs %d), line: %d (non-empty line %d), file: '%s'.\n",thisfunction,count,prevcount,linecount,num,filename);
	  fclose(fp); exit(-1);
	}
      }
      prevcount = count;
    }
  }
  
  if (count < 2) {
    fprintf(stderr,"ERROR (%s): minimum point coordinate dimension is 1 (line format: species coord), now only %d value per line in file '%s'.\n",thisfunction,count, filename);
    fclose(fp); exit(-1);
  }

  dim = count-1;

  if ((p = (Point *) malloc (num * sizeof(Point))) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot malloc %d Points.\n",thisfunction,num);
    perror(""); fclose(fp); exit(-1);
  }
  for (i=0; i<num; i++) {
    if ((p[i].position = (double *) malloc (dim * sizeof(double))) == NULL) {
      fprintf(stderr,"ERROR (%s): cannot malloc %d doubles for position of point %d/%d.\n",thisfunction,dim,i+1,num);
      perror(""); fclose(fp); exit(-1);
    }
  }

  *num_points = num;
  *dimension = dim;

  rewind(fp);

  num=0;
  while (fgets(line, MAXSTRINGLEN, fp) != NULL) {
    if ((tok = strtok(line," \t\r\n")) != NULL) {
      p[num].species = atoi(tok);
      for (i=0; i<dim; i++) {
	if ((tok = strtok(NULL," \t\r\n")) == NULL) {
	  fprintf(stderr,"ERROR (%s): cannot read coordinate %d/%d for point %d/%d, file '%s'.\n",thisfunction,i+1,dim,num,*num_points, filename);
	  perror(""); fclose(fp); exit(-1);
	}
	p[num].position[i] = atof(tok);
      }
      num++;
    }
  }
  
  fclose(fp);

  return (p);
}

int printPoints(Point *p, int num_points, int dimension) {
  int i,j;

  for (i=0; i<num_points; i++) {
    printf("Point %d/%d, species %d, position (%f",i+1,num_points,p[i].species,p[i].position[0]);
    for (j=1; j<dimension; j++) 
      printf(",%f",p[i].position[j]);
    printf(")\n");
  }
      
  return (0);
}

int addPointToList(Point *p, Cell *c) {
  LinkedList *item;
  char *thisfunction = "addPointToList";

  if ((item = (LinkedList *) malloc (sizeof(LinkedList))) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot malloc LinkedList item.\n",thisfunction);
    perror(""); exit(-1);
  }
  item->p = p;
  item->next = NULL;

  if (c->first == NULL) {
    c->first = item;
    c->last = item;
    item->prev = NULL;
  }
  else {
    item->prev = c->last;
    c->last->next = item;
    c->last = item;
  }

  c->num_points++;
  return (0);
}

SimulationSpace *createSimulationSpace(double U, double w, int dimension, Point *point, int num_points) {
  SimulationSpace *s;
  int i,j,k,t,d,ci,cj;
  char *thisfunction = "createSimulationSpace";

  if ((s = (SimulationSpace *) malloc (sizeof(SimulationSpace))) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot malloc SimulationSpace.\n",thisfunction);
    perror(""); exit(-1);
  }
  s->dimension = dimension;
  /* ceil U to integer multiple of w so that area of each cell is constant */
  s->L = (int) (U/w + 0.5);
  s->cell_width = w;
  s->U = s->L * w;

  if (s->L < 1) {
    fprintf(stderr,"ERROR (%s): L is less than 1, happens when cell width exceeds space size U (here cell_width %f, U %f).\n",thisfunction,w,U);
    exit(-1);
  }

  U=s->U;

  switch (dimension) {
  case 1:
    s->num_cells = s->L;
    break;
  case 2:
    s->num_cells = s->L * s->L;
    break;
  default:
    fprintf(stderr,"ERROR (%s): only 1D and 2D spaces supported (input point coordinate dimension was %d, num_points %d).\n",thisfunction,dimension,num_points);
    return(NULL);
    break;
  }

  if ((s->cell = (Cell *) malloc (s->num_cells * sizeof(Cell))) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot malloc %d cells (dimension %d, L %d).\n",thisfunction,s->num_cells,dimension,s->L);
    perror(""); return(NULL);
  }

  for (k=0; k<s->num_cells; k++) {
    if ((s->cell[k].position = (double *) malloc (dimension * sizeof(double))) == NULL) {
      fprintf(stderr,"ERROR (%s): cannot malloc position vector for cell (dimension %d, %d cells).\n",thisfunction, dimension, s->num_cells);
      perror(""); exit(-1);
    }    
  }

  switch (dimension) {
  case 1:
    for (k=0; k<s->num_cells; k++) {
      s->cell[k].position[0] = k * s->cell_width;
      s->cell[k].num_points = 0;
      s->cell[k].first = NULL;
      s->cell[k].last = NULL;
    }

    for (i=0; i<num_points; i++) {
      while (point[i].position[0] < 0.0) 
	point[i].position[0] += s->U;
      while (point[i].position[0] >= s->U) 
	point[i].position[0] -= s->U;

      k = (int) (point[i].position[0]/s->cell_width);
      addPointToList(&(point[i]), &(s->cell[k]));
    }

    break;
  case 2:
    k=0;
    for (i=0; i<s->L; i++) {
      for (j=0; j<s->L; j++) {
	s->cell[k].position[0] = i * s->cell_width;
	s->cell[k].position[1] = j * s->cell_width;
	s->cell[k].num_points = 0;
	s->cell[k].first = NULL;
	s->cell[k].last = NULL;
	k++;
      }
    }

    for (i=0; i<num_points; i++) {
      for (d=0; d<dimension; d++) {
	while (point[i].position[d] < 0.0) 
	  point[i].position[d] += s->U;
	while (point[i].position[d] >= s->U) 
	  point[i].position[d] -= s->U;
      }
      ci = (int) (point[i].position[0]/s->cell_width);
      cj = (int) (point[i].position[1]/s->cell_width);
      k = ci * s->L + cj;
      addPointToList(&(point[i]), &(s->cell[k]));
    }
    break;
  default:
    fprintf(stderr,"ERROR (%s): dimension %d.\n",thisfunction,dimension);
    exit(-1);
    break;
  }

  return(s);
}

int printSimulationSpace(SimulationSpace *s, int num_ckernels, int num_okernels, int num_processes) {
  int i,k,d,j;
  LinkedList *li;

  printf("SimulationSpace: dimension %d, U %f, %d cells, cell width %f.\n",s->dimension,s->U,s->num_cells, s->cell_width);

  if (num_processes >0) {
    printf("pactivity: (%f",s->pactivity[0]);
    for (j=1; j<num_processes; j++)
      printf(",%f",s->pactivity[j]);
    printf(")\n");
  }

  for (k=0; k<s->num_cells; k++) {
    printf("Cell %d/%d:",k+1,s->num_cells);
    printf(" position (%f",s->cell[k].position[0]);
    for (d=1; d<s->dimension; d++) 
      printf(",%f",s->cell[k].position[d]);
    printf(") num_points: %d species count (%d",s->cell[k].num_points,s->cell[k].species_count[0]);
    for (i=1; i<s->num_species_count; i++)
      printf(",%d",s->cell[k].species_count[i]);
    printf(")\n");

    printf("pactivity: (%f",s->cell[k].pactivity[0]);
    for (j=1; j<num_processes; j++)
      printf(",%f",s->cell[k].pactivity[j]);
    printf(")\n");

    for (li=s->cell[k].first; li != NULL; li=li->next) {
      printf("  point (%f",li->p->position[0]);
      for (d=1; d<s->dimension; d++) 
	printf(",%f",li->p->position[d]);
      printf(") species %d\n",li->p->species);

      if (num_ckernels > 0) {
	printf("   cactivity: (%f",li->p->cactivity[0]);
	for (j=1; j<num_ckernels; j++)
	  printf(",%f",li->p->cactivity[j]);
	printf(") pactivity: (%f",li->p->pactivity[0]);
	for (j=1; j<num_processes; j++)
	  printf(",%f",li->p->pactivity[j]);
	printf(")\n");
      }
      
    }
  }

  return (0);
}

double distance2(double *p1, double *p2, int dim, double U) {

  /* return squared Euclidean distance in a torus topology */

  int i;
  double d1,d2=0.0,halfU;

  halfU = U/2.0;
  for (i=0; i<dim; i++) {
    if ((d1 = p2[i]-p1[i]) < 0.0)
      d1 = -d1;
    if (d1 > halfU)
      d1 = U - d1;
    d2 += d1*d1;
  }

  return(d2);
}

int *getCellNeighbors(double *p, int dim, int L, double cell_width, double r, int *num_ind) {
  /* note: argument r is radius instead of squared radius */

  static int maxnum=0, *index=NULL;
  int i,j,k,ci,cj,ki,kj,nc,num,t,found;
  double U;
  char *thisfunction = "getCellNeighbors";

  nc = (int) ceil(r/cell_width);
  U = L * cell_width;

  /* calculate number of cells within box */
  /* no need to cover more than entire space once */
  if (2*nc+1 > L)
    nc = L/2;

  num = 2*nc +1;
  for (i=1; i<dim; i++)
    num *= 2*nc +1;

  if (num > maxnum) {
    if ((index = (int *) realloc(index, num * sizeof(int))) == NULL) {
      fprintf(stderr,"ERROR (%s): cannot realloc integer vector from size %d to %d.\n",thisfunction,maxnum,num);
      perror(""); exit(-1);
    }
    maxnum=num;
  }

  num=0;
  switch (dim) {
  case 1:
    ci = (int) (p[0]/cell_width);
    for (i=ci-nc; i<=ci+nc; i++) {
      k = i;
      while (k<0)
	k += L;
      while (k>=L)
	k -= L;
      /* one cell index only once in index array */
      for (found=0, t=0; (t<num) && (!found); t++) {
	if (index[t] == k)
	  found =1;
      }
      if (!found) {
	index[num] = k;
	num++;
      }
    }
    break;
  case 2:
    ci = (int) (p[0]/cell_width);
    cj = (int) (p[1]/cell_width);
    for (i=ci-nc; i<=ci+nc; i++) {
      ki=i;
      while (ki<0)
	ki += L;
      while (ki>=L)
	ki -= L;      
      for (j=cj-nc; j<=cj+nc; j++) {
	kj=j;
	while (kj<0)
	  kj += L;
	while (kj>=L)
	  kj -= L;      

	k = ki*L + kj;
	/* one cell index only once in index array */
	for (found=0, t=0; (t<num) && (!found); t++) {
	  if (index[t] == k)
	    found =1;
	}
	if (!found) {
	  index[num] = k;
	  num++;
	}
      }
    }
    break;
  default:
    fprintf(stderr,"ERROR (%s): only 1D and 2D lattices supported (dim %d).\n",thisfunction,dim);
    exit(-1);
    break;
  }

  *num_ind = num;
  return(index);
}

int initializeActivities(SimulationSpace *s, Process *p, int num_p, Kernel *ckernels, int num_ckernels, Kernel *okernels, int num_okernels) {
  int i,k,t,sk,st,*ind,num_ind;
  LinkedList *li, *ti;
  double r2, dist2;
  char *thisfunction = "initializeActivities";

  if ((s->pactivity = (double *) malloc (num_p * sizeof(double))) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot malloc %d doubles for pactivity.\n",thisfunction,num_p);
    perror(""); exit(-1);
  }

  for (k=0; k<s->num_cells; k++) {
    if ((s->cell[k].pactivity = (double *) malloc (num_p * sizeof(double))) == NULL) {
      fprintf(stderr,"ERROR (%s): cannot malloc %d doubles for pactivity, cell %d/%d.\n",thisfunction,num_p,k+1,s->num_cells);
      perror(""); exit(-1);
    }
  }

  for (k=0; k<s->num_cells; k++) {
    for (li=s->cell[k].first; li != NULL; li=li->next) {
      if ((li->p->cactivity = (double *) malloc (num_ckernels * sizeof(double))) == NULL) {
	fprintf(stderr,"ERROR (%s): cannot malloc %d doubles for cactivity, point of cell %d/%d",thisfunction,num_ckernels,k+1,s->num_cells);
	perror(""); exit(-1);
      }
      for (i=0; i<num_ckernels; i++)
	li->p->cactivity[i] = 0.0;

      if ((li->p->pactivity = (double *) malloc (num_p * sizeof(double))) == NULL) {
	fprintf(stderr,"ERROR (%s): cannot malloc %d doubles for pactivity, point of cell %d/%d.\n",thisfunction,num_p,k+1,s->num_cells);
	perror(""); exit(-1);
      }
      for (i=0; i<num_p; i++)
	li->p->pactivity[i] = 0.0;
    }
  }

  /* this allows indexing species_count directly by point->species */

  s->num_species_count = 0;
  for (i=0; i<num_p; i++) {
    for (k=0; k<p[i].num_args; k++) {
      if (p[i].arg_species[k] +1 > s->num_species_count)
	s->num_species_count = p[i].arg_species[k] +1;
    }
  }

  if ((s->species_count = (int *) malloc (s->num_species_count * sizeof(int))) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot malloc %d ints for species_count.\n",thisfunction,s->num_species_count);
    perror(""); exit(-1);
  }
  for (i=0; i<s->num_species_count; i++)
    s->species_count[i] = 0;

  for (k=0; k<s->num_cells; k++) {
    if ((s->cell[k].species_count = (int *) malloc (s->num_species_count * sizeof(int))) == NULL) {
      fprintf(stderr,"ERROR (%s): cannot malloc %d ints for species_count, cell %d/%d.\n",thisfunction,s->num_species_count,k+1,s->num_cells);
      perror(""); exit(-1);
    }
    for (i=0; i<s->num_species_count; i++)
      s->cell[k].species_count[i] = 0;
  }

  for (k=0; k<s->num_cells; k++) {
    for (li=s->cell[k].first; li != NULL; li=li->next) {
      s->cell[k].species_count[li->p->species]++;
    }
  }

  for (i=0; i<s->num_species_count; i++) {
    for (k=0; k<s->num_cells; k++)
      s->species_count[i] += s->cell[k].species_count[i];
  }
  for (i=1; i<s->num_species_count; i++) 
    s->species_count[0] += s->species_count[i];

  /* calculate cactivities */

  for (k=0; k<s->num_cells; k++) {
    for(li=s->cell[k].first; li != NULL; li=li->next) {
      for (i=0; i<num_ckernels; i++) {

	r2 = ckernels[i].radius * ckernels[i].radius;
	sk=-1;
	if (li->p->species == ckernels[i].species1) {
	  sk = ckernels[i].species1;
	  st = ckernels[i].species2;
	}
	else if (li->p->species == ckernels[i].species2) {
	  sk = ckernels[i].species2;
	  st = ckernels[i].species1;
	}
	if (sk >= 0) {
	  ind = getCellNeighbors(li->p->position, s->dimension, s->L, s->cell_width, ckernels[i].radius, &num_ind);
	  switch (ckernels[i].family) {
	  case TOPHAT_KERNEL:
	    for (t=0; t<num_ind; t++) {
	      if (ind[t] >= k) { 
		/* not yet processed */
		/* if points within the same cell (t==k) update only one direction */
		if (s->cell[ind[t]].species_count[st] > 0) {
		  /* calculate distances */
		  for (ti=s->cell[ind[t]].first; ti != NULL; ti=ti->next) {
		    if ((ti != li) && (ti->p->species == st) && (distance2(li->p->position, ti->p->position, s->dimension, s->U) <= r2)) {
		      li->p->cactivity[i] += ckernels[i].maxdensity;
		      if (ind[t] != k) {
			ti->p->cactivity[i] += ckernels[i].maxdensity;
		      }
		    }
		  }
		}
	      }
	    }
	    break;
	  case TRUNCATED_GAUSSIAN_KERNEL:
	    for (t=0; t<num_ind; t++) {
	      if (ind[t] >= k) {
		/* not yet processed */
		/* if points within the same cell (t==k) update only one direction */
		if (s->cell[ind[t]].species_count[st] > 0) {
		  for (ti=s->cell[ind[t]].first; ti != NULL; ti=ti->next) {
		    if ((ti != li) && (ti->p->species == st) && ((dist2=distance2(li->p->position, ti->p->position, s->dimension, s->U)) <= r2)) {
		      li->p->cactivity[i] += ckernels[i].maxdensity * exp(-ckernels[i].halfivar*dist2) + ckernels[i].trunceffect;
		      if (ind[t] != k)
			ti->p->cactivity[i] += ckernels[i].maxdensity * exp(-ckernels[i].halfivar*dist2) + ckernels[i].trunceffect;
		    }
		  }
		}
	      }
	    }
	    break;
	  default:
	    fprintf(stderr,"ERROR (%s): unrecognized/forbidden kernel family %d used as connectivity kernel, kernel %d/%d.\n",thisfunction,ckernels[i].family,i+1,num_ckernels);
	    exit(-1);
	    break;
	  }	    
	}
      }
    }
  }

  /* calculate pactivities, first for processes with source species (i.e. process has input) */

  for (i=0; i<num_p; i++) 
    s->pactivity[i] = 0.0; 

  for (k=0; k<s->num_cells; k++) {
    for (i=0; i<num_p; i++) 
      s->cell[k].pactivity[i] = 0.0; 
    
    for (li=s->cell[k].first; li != NULL; li=li->next) {
      for (i=0; i<num_p; i++) {
	li->p->pactivity[i] = 0.0;
	/* store values only in process sources (which has connections to all other process inputs) */
	if (li->p->species == p[i].source_species) {
	  li->p->pactivity[i] = 1.0;
	  for (t=0; t<p[i].num_cind; t++) {
	    li->p->pactivity[i] *= li->p->cactivity[ p[i].cind[t] ];
	  }
	  for (t=0; t<p[i].num_oind; t++) {
	    li->p->pactivity[i] *= okernels[ p[i].oind[t] ].integral;
	  }
	}
      }
      for (i=0; i<num_p; i++) {
	s->cell[k].pactivity[i] += li->p->pactivity[i]; 
      }
    }
    for (i=0; i<num_p; i++) {
      s->pactivity[i] += s->cell[k].pactivity[i]; 
    }
  }
  
  /* calculate o and p activities for processes without input (rate_per_area events) */
  /* there can be only one okernel (storing rate_per_area coefficient) related to process */

  r2 = s->cell_width;
  for (i=1; i<s->dimension; i++)
    r2 *= s->cell_width;

  for (k=0; k<s->num_cells; k++) {
    for (i=0; i<num_p; i++) {
      if (p[i].pm->type == RATE_PER_AREA)
	s->cell[k].pactivity[i] = r2 * okernels[ p[i].oind[0] ].integral;
    }
  }

  r2 = s->U;
  for (i=1; i<s->dimension; i++)
    r2 *= s->U;

  for (i=0; i<num_p; i++) {
    if (p[i].pm->type == RATE_PER_AREA)
      s->pactivity[i] = r2 * okernels[ p[i].oind[0] ].integral;
  }

  return(0);
}

int getPactivityArrayIndex(double th, double *a, int n) {
  int i=0;
  double sum;

  sum=a[i];
  while ((sum < th) && (i<n-1)) {
    i++;
    sum+=a[i];
  }  
  return(i);
}

int getPactivityCellIndex(double th, int pi, Cell *c, int n) {
  int i=0;
  double sum;

  sum=c[i].pactivity[pi];
  while ((sum < th) && (i<n-1)) {
    i++;
    sum+=c[i].pactivity[pi];
  }  
  return(i);
}

LinkedList *getPactivityListItem(double th, int pi, LinkedList *first) {
  LinkedList *li;
  double sum;

  if (first == NULL)
    return (NULL);
  
  li = first;
  sum = li->p->pactivity[pi];
  while ((sum < th) && (li->next != NULL)) {
    li = li->next;
    sum += li->p->pactivity[pi];
  }
  return (li);
}

LinkedList *getOnePointWithinKernel(int target_species, double *center, Kernel *k, LinkedList *exclude_item, SimulationSpace *s, int *cell_ind) {
  int *ind, num_ind, i, t, count,count_th,found;
  double r2, dist2, dsum, dsum_th;
  LinkedList *ti, *ri=NULL;
  char *thisfunction = "getPointWithinKernel";

  if (k->family == TOPHAT_KERNEL) {
    r2 = k->radius * k->radius;  
    count=0;
    ind = getCellNeighbors(center, s->dimension, s->L, s->cell_width, k->radius, &num_ind);
    for (t=0; t<num_ind; t++) {
      /* have to calculate distances */
      for (ti=s->cell[ind[t]].first; ti != NULL; ti=ti->next) {
	if ((ti->p->species == target_species) && (distance2(center, ti->p->position, s->dimension, s->U) <= r2) && (ti != exclude_item))
	  count++;
      }
    }
    
    count_th = (int) (count*drand48());
    count=0;
    found=0;
    for (t=0; (t<num_ind) && !found; t++) {
      /* have to calculate distances */
      for (ti=s->cell[ind[t]].first; ti != NULL; ti=ti->next) {
	if ((ti->p->species == target_species) && (distance2(center, ti->p->position, s->dimension, s->U) <= r2) && (ti != exclude_item)) {
	  count++;
	  if (count >= count_th) {
	    ri = ti;
	    *cell_ind = ind[t];
	    found = 1;
	  }
	}
      }
    }
  }
  else if (k->family == TRUNCATED_GAUSSIAN_KERNEL) {
    r2 = k->radius * k->radius;  
    dsum=0;
    ind = getCellNeighbors(center, s->dimension, s->L, s->cell_width, k->radius, &num_ind);
    for (t=0; t<num_ind; t++) {
      /* have to calculate distances */
      for (ti=s->cell[ind[t]].first; ti != NULL; ti=ti->next) {
	if ((ti->p->species == target_species) && ((dist2 = distance2(center, ti->p->position, s->dimension, s->U)) <= r2) && (ti != exclude_item))
	  dsum += k->maxdensity * exp(-k->halfivar*dist2) + k->trunceffect;
      }
    }

    dsum_th = dsum*drand48();
    dsum=0;
    found=0;
    for (t=0; (t<num_ind) && !found; t++) {
      /* have to calculate distances */
      for (ti=s->cell[ind[t]].first; ti != NULL; ti=ti->next) {
	if ((ti->p->species == target_species) && ((dist2 = distance2(center, ti->p->position, s->dimension, s->U)) <= r2) && (ti != exclude_item)) {
	  dsum += k->maxdensity * exp(-k->halfivar*dist2) + k->trunceffect;
	  if (dsum >= dsum_th) {
	    ri = ti;
	    *cell_ind = ind[t];
	    found = 1;
	  }
	}
      }
    }
  }
  else {
    fprintf(stderr,"ERROR (%s): only TOPHAT and TRUNCATED_GAUSSIAN possible (asked for %s).\n",thisfunction,KERNEL_TYPE[k->family]);
    return(NULL);
  }

  return(ri);
}

LinkedList *newPoint(LinkedList **garbage, int dim, int num_ckernels, int num_processes) {
  int i;
  LinkedList *item;
  char *thisfunction = "addPointToList";

  /* check first if there are any points in garbage list, if not, malloc new */

  /* take item from the end of the list */
  if (*garbage != NULL) {
    item = *garbage;
    *garbage = (*garbage)->prev;
    if (*garbage != NULL)
      (*garbage)->next = NULL;
  }
  else {
    if ((item = (LinkedList *) malloc (sizeof(LinkedList))) == NULL) {
      fprintf(stderr,"ERROR (%s): cannot malloc LinkedList item.\n",thisfunction);
      perror(""); exit(-1);
    }
    if ((item->p = (Point *) malloc (sizeof(Point))) == NULL) {
      fprintf(stderr,"ERROR (%s): cannot malloc Point.\n",thisfunction);
      perror(""); exit(-1);
    }
    if ((item->p->position = (double *) malloc (dim * sizeof(double))) == NULL) {
      fprintf(stderr,"ERROR (%s): cannot malloc %d doubles for position of point.\n",thisfunction,dim);
      perror(""); exit(-1);
    }
    if ((item->p->cactivity = (double *) malloc (num_ckernels * sizeof(double))) == NULL) {
      fprintf(stderr,"ERROR (%s): cannot malloc %d doubles for cactivity of point.\n",thisfunction,num_ckernels);
      perror(""); exit(-1);
    }
    if ((item->p->pactivity = (double *) malloc (num_processes * sizeof(double))) == NULL) {
      fprintf(stderr,"ERROR (%s): cannot malloc %d doubles for pactivity of point.\n",thisfunction,num_processes);
      perror(""); exit(-1);
    }    
  }
  /* reset cactivity and pactivity, but don't touch position */
  for (i=0; i<num_ckernels; i++)
    item->p->cactivity[i] = 0.0;
  for (i=0; i<num_processes; i++)
    item->p->pactivity[i] = 0.0;

  return (item);
}

int removePoint(LinkedList *item, Cell *c, LinkedList **garbage) {

  /* don't free memory but add point to garbage list */

  if (item->prev == NULL) {
    /* first item in list */
    c->first = item->next;
    if (item->next == NULL) {
      /* the only item in the list */
      c->last = NULL;
    }
    else {
      c->first->prev = NULL;
    }
  }
  else {
    if (item->next == NULL) {
      /* last item in the list */
      item->prev->next = NULL;
      c->last = item->prev;
    }
    else {
      /* in the middle of the list */
      item->prev->next = item->next;
      item->next->prev = item->prev;
    }
  }

  /* garbage points to the end of unused LinkedList items */
  if (*garbage != NULL) {
    /* add to the end of garbage list */
    (*garbage)->next = item;
  }
  item->prev = *garbage;
  item->next = NULL;
  *garbage = item;

  c->num_points--;
  c->species_count[item->p->species]--;

  return(0);
}

int printList(LinkedList *item) {
  LinkedList *li;

  for (li=item->prev; li != NULL; li=li->prev)
    printf(" Before-Listitem %d %f\n",li->p->species,li->p->position[0]);
  printf(" ListItem %d %f\n",item->p->species,item->p->position[0]);
  for (li=item->next; li != NULL; li=li->next)
    printf(" After-ListItem %d %f\n",li->p->species,li->p->position[0]);
  return(0);
}

double randGaussian(void) {
  double u,s;
  static double v;
  static int ready = 0;

  /* Marsaglia's polar method */
  
  if (ready) {
    ready=0;
    return (v);
  }
  ready=1;

  do {
    u = 2.0*drand48() - 1.0;
    v = 2.0*drand48() - 1.0;
    s = u*u + v*v;
  }
  while ((s>=1.0) || (s==0.0));

  s=sqrt(-2.0 * log(s)/s);
  v *= s;
  u *= s;

  return(u);
}

int randomPosition(double *center, Kernel *k, int dim, double U, double *p) {
  int i;
  double a,r;

  if (k->family == CONSTANT_KERNEL) {
    for (i=0; i<dim; i++) {
      p[i] = U * drand48();
    }
  }
  else if (k->family == TOPHAT_KERNEL) {
    if (dim == 1) {
      p[0] = center[0] + k->radius * (2.0 * drand48() - 1.0);
    }
    else if (dim == 2) {
      a = 6.283185 * drand48();
      r = k->radius * sqrt(drand48());
      p[0] = center[0] + r * cos(a);
      p[1] = center[1] + r * sin(a);      
    }
    else {
      fprintf(stderr,"ERROR (randomPosition): dimension must be 1 or 2 (was %d).\n",dim);
      return(-1);    
    }
  }
  else if (k->family == TRUNCATED_GAUSSIAN_KERNEL) {
    if (dim == 1) {
      r = randGaussian();
      if ((r < -3.0) || (r > 3.0))
	r = 6.0*drand48() - 3.0;
      p[0] = center[0] + k->sigma * r;    
    }
    else if (dim == 2) {
      /* first (r,a)=(x,y) from randGaussian, but if beyond truncation border then (radius,angle) of new location */
      r = randGaussian();
      a = randGaussian();
      if (r*r + a*a > 9.0) {
	a = 6.283185 * drand48();
	r = 3.0 * k->sigma * sqrt(drand48());
	p[0] = center[0] + r * cos(a);
	p[1] = center[1] + r * sin(a);
      }
      else {
	p[0] = center[0] + k->sigma * r;
	p[1] = center[1] + k->sigma * a;
      }
    }
    else {
      fprintf(stderr,"ERROR (randomPosition): dimension must be 1 or 2 (was %d).\n",dim);
      return(-1);    
    }
  }
  else {
    fprintf(stderr,"ERROR (randomPosition): unsupported kernel %s.\n",KERNEL_FAMILY[k->family]);
    return(-1);
  }
  
  /* ensure point is within space limits */
  for (i=0; i<dim; i++) {
    while (p[i] < 0.0) 
      p[i] += U;
    while (p[i] >= U) 
      p[i] -= U;
  }

  return(0);
}

int getKernelIndex(int target_species, int point_type, Kernel *k, int *kind, int num_kind, int *ptype1, int *ptype2) {
  int i;

  /* note: kind is an index array, this function returns the index of original kernel array, not the index of the index array */

  for (i=0; i<num_kind; i++) {
    if ( ((k[kind[i]].species1 == target_species) && (ptype1[i] == point_type)) || ((k[kind[i]].species2 == target_species) && (ptype2[i] == point_type)) ) {
      return(kind[i]);
    }
  }
  return (-1);
}

int getCellIndex(double *p, SimulationSpace *s) {
  int k;
  char *thisfunction = "getCellIndex";

  switch (s->dimension) {
  case 1:
    while (p[0] < 0.0) 
      p[0] += s->U;
    while (p[0] >= s->U) 
      p[0] -= s->U;
    k = (int) (p[0]/s->cell_width);
    break;
  case 2:
    for (k=0; k<s->dimension; k++) {
      while (p[k] < 0.0) 
	p[k] += s->U;
      while (p[k] >= s->U) 
	p[k] -= s->U;
    }
    /* i*L + j */
    k = ((int) (p[0]/s->cell_width)) * s->L + ((int) (p[1]/s->cell_width));
    break;
  default:
    fprintf(stderr,"ERROR (%s): dimension %d.\n",thisfunction,s->dimension);
    exit(-1);
    break;
  }
  return(k);
}

int appendList(LinkedList *item, Cell *c) {

  item->next = NULL;

  if (c->first == NULL) {
    c->first = item;
    c->last = item;
    item->prev = NULL;
  }
  else {
    item->prev = c->last;
    c->last->next = item;
    c->last = item;
  }

  c->num_points++;
  c->species_count[item->p->species]++;

  return (0);
}

int addEffects(Point *point, int ci,  SimulationSpace *s, Kernel *ckernel, int num_ckernels, Kernel *okernel, int num_okernels, Process *process, int num_processes) {
  int i,j,t,maxi, *ind, num_ind, update;
  double dist2;
  LinkedList *li;
  char *thisfunction = "addEffects";

  /* find largest ckernel radius applicable to point->species */
  maxi=-1;
  for (i=0; i<num_ckernels; i++) {
    if (maxi < 0) {
      if ((ckernel[i].species1 == point->species) || (ckernel[i].species2 == point->species))
	maxi=i;
    }
    else if (((ckernel[i].species1 == point->species) || (ckernel[i].species2 == point->species)) && (ckernel[i].radius > ckernel[maxi].radius))
      maxi = i;
  }

  if (maxi >= 0) {
    /* this species has connectivity effect, update them before calculating pactivity */

    ind = getCellNeighbors(point->position, s->dimension, s->L, s->cell_width, ckernel[maxi].radius, &num_ind);
    for (j=0; j<num_ind; j++) {
      for (li=s->cell[ind[j]].first; li != NULL; li=li->next) {
	/* check which kernel(s) apply to point->species, li->p->species pair within r2 */
	update=0;
	for (i=0; i<num_ckernels; i++) {
	  if (((ckernel[i].species1 == point->species) && (ckernel[i].species2 == li->p->species)) || ((ckernel[i].species2 == point->species) && (ckernel[i].species1 == li->p->species))) {
	    if ((dist2=distance2(point->position, li->p->position, s->dimension, s->U)) <= ckernel[i].radius * ckernel[i].radius) {
	      if (ckernel[i].family == TOPHAT_KERNEL) {
		li->p->cactivity[i] += ckernel[i].maxdensity;
		point->cactivity[i] += ckernel[i].maxdensity;
	      }
	      else if (ckernel[i].family == TRUNCATED_GAUSSIAN_KERNEL) {
		li->p->cactivity[i] += ckernel[i].maxdensity * exp(-ckernel[i].halfivar*dist2) + ckernel[i].trunceffect;
		point->cactivity[i] += ckernel[i].maxdensity * exp(-ckernel[i].halfivar*dist2) + ckernel[i].trunceffect;
	      }
	      else {
		fprintf(stderr,"ERROR (%s): unsupported/forbidden kernel family %d (ckernel %d/%d).\n",thisfunction,ckernel[i].family,i+1,num_ckernels);
		exit(-1);
	      }
	      update=1;
	    }
	  }
	}

	if (update) {
	  /* update pactivity of li->p (in principle would be enough just to multiply by new_cactivity/old_cactivity...but now calculating all) */
	  for (i=0; i<num_processes; i++) {
	    /* store values only in process sources (which has connections to all other process inputs) */
	    if (li->p->species == process[i].source_species) {
	      /* subtract effect of old value */
	      s->cell[ind[j]].pactivity[i] -= li->p->pactivity[i];
	      s->pactivity[i] -= li->p->pactivity[i];
	      
	      li->p->pactivity[i] = 1.0;
	      for (t=0; t<process[i].num_cind; t++) 
		li->p->pactivity[i] *= li->p->cactivity[ process[i].cind[t] ];
	      for (t=0; t<process[i].num_oind; t++)
		li->p->pactivity[i] *= okernel[ process[i].oind[t] ].integral;
	      
	      /* add effect of new value */
	      s->cell[ind[j]].pactivity[i] += li->p->pactivity[i];
	      s->pactivity[i] += li->p->pactivity[i];
	    }
	  }
	}
      }
    }
  }

  /* calculate pactivity of new point */
  for (i=0; i<num_processes; i++) {
    /* store values only in process sources (which has connections to all other process inputs) */
    if (point->species == process[i].source_species) {
      point->pactivity[i] = 1.0;
      for (t=0; t<process[i].num_cind; t++) 
	point->pactivity[i] *= point->cactivity[ process[i].cind[t] ];
      for (t=0; t<process[i].num_oind; t++)
	point->pactivity[i] *= okernel[ process[i].oind[t] ].integral;

      /* add effect of new value */
      s->cell[ci].pactivity[i] += point->pactivity[i];
      s->pactivity[i] += point->pactivity[i];
    }
  }

  return(0);
}

int removeEffects(Point *point, int ci,  SimulationSpace *s, Kernel *ckernel, int num_ckernels, Kernel *okernel, int num_okernels, Process *process, int num_processes) {
  int i,j,t,maxi, *ind, num_ind, update;
  double dist2;
  LinkedList *li;
  char *thisfunction = "removeEffects";

  /* find largest ckernel radius applicable to point->species */
  maxi=-1;
  for (i=0; i<num_ckernels; i++) {
    if (maxi < 0) {
      if ((ckernel[i].species1 == point->species) || (ckernel[i].species2 == point->species))
	maxi=i;
    }
    else if (((ckernel[i].species1 == point->species) || (ckernel[i].species2 == point->species)) && (ckernel[i].radius > ckernel[maxi].radius))
      maxi = i;
  }
  if (maxi >= 0) {
    /* this species has connectivity effect */

    ind = getCellNeighbors(point->position, s->dimension, s->L, s->cell_width, ckernel[maxi].radius, &num_ind);
    for (j=0; j<num_ind; j++) {
      for (li=s->cell[ind[j]].first; li != NULL; li=li->next) {
	/* check which kernel(s) apply to point->species, li->p->species pair within r2 */
	update=0;
	for (i=0; i<num_ckernels; i++) {
	  if (((ckernel[i].species1 == point->species) && (ckernel[i].species2 == li->p->species)) || ((ckernel[i].species2 == point->species) && (ckernel[i].species1 == li->p->species))) {
	    if ((dist2=distance2(point->position, li->p->position, s->dimension, s->U)) <= ckernel[i].radius * ckernel[i].radius) {
	      if (ckernel[i].family == TOPHAT_KERNEL) {
		li->p->cactivity[i] -= ckernel[i].maxdensity;
	      }
	      else if (ckernel[i].family == TRUNCATED_GAUSSIAN_KERNEL) {
		li->p->cactivity[i] -= ckernel[i].maxdensity * exp(-ckernel[i].halfivar*dist2) + ckernel[i].trunceffect;
	      }
	      else {
		fprintf(stderr,"ERROR (%s): unsupported/forbidden kernel family %d (ckernel %d/%d).\n",thisfunction,ckernel[i].family,i+1,num_ckernels);
		exit(-1);
	      }
	      update=1;
	    }
	  }
	}
	  
	if (update) {
	  /* update pactivity of li->p (in principle would be enough just to multiply by new_cactivity/old_cactivity...but now calculating all) */
	  for (i=0; i<num_processes; i++) {
	    /* store values only in process sources (which has connections to all other process inputs) */
	    if (li->p->species == process[i].source_species) {
	      /* subtract effect of old value */
	      s->cell[ind[j]].pactivity[i] -= li->p->pactivity[i];
	      s->pactivity[i] -= li->p->pactivity[i];
	      
	      li->p->pactivity[i] = 1.0;
	      for (t=0; t<process[i].num_cind; t++)
		li->p->pactivity[i] *= li->p->cactivity[ process[i].cind[t] ];
	      for (t=0; t<process[i].num_oind; t++)
		li->p->pactivity[i] *= okernel[ process[i].oind[t] ].integral;
	      
	      /* add effect of new value */
	      s->cell[ind[j]].pactivity[i] += li->p->pactivity[i];
	      s->pactivity[i] += li->p->pactivity[i];
	    }
	  }
	}
      }
    }
  }

  /* remove pactivity of old point */
  for (i=0; i<num_processes; i++) {
    /* store values only in process sources (which has connections to all other process inputs) */
    if (point->species == process[i].source_species) {
      s->cell[ci].pactivity[i] -= point->pactivity[i];
      s->pactivity[i] -= point->pactivity[i];
    }
  }
  
  return(0);
}

int updateConfiguration(SimulationSpace *s, int ci, LinkedList *li, int pi, Process *processes, int num_processes, Kernel *ckernel, int num_ckernels, Kernel *okernel, int num_okernels) {
  int i,j, cj, ki, hi, target_species, processed_one;
  ProcessModule *pm;
  LinkedList *item;
  Process *p;
  static LinkedList *garbage=NULL;
  static Point **pptr=NULL;
  static int num_pptrs=0;
  static LinkedList **hold_item_ptr = NULL;
  static int num_hold_item_ptrs;
  char *thisfunction = "updateConfiguration";

  p = &(processes[pi]);
  pm = p->pm;

  if (pm->num_products > num_hold_item_ptrs) {
    num_hold_item_ptrs = pm->num_products;
    if ((hold_item_ptr = (LinkedList **) realloc(hold_item_ptr, num_hold_item_ptrs * sizeof(LinkedList *))) == NULL) {
      fprintf(stderr,"ERROR (%s): cannot realloc %d LinkedList ptrs.\n",thisfunction,num_hold_item_ptrs);
      perror(""); exit(-1);
    }
  }

  /* first check if any products need exact coordinates */

  if (pm->num_same_coordinate_reactant_product_pairs > 0) {
    if (pm->num_same_coordinate_reactant_product_pairs > num_pptrs) {
      num_pptrs = pm->num_same_coordinate_reactant_product_pairs;
      if ((pptr = (Point **) realloc(pptr, num_pptrs * sizeof(Point *))) == NULL) {
	fprintf(stderr,"ERROR (%s): cannot realloc %d Point ptrs.\n",thisfunction,num_pptrs);
	perror(""); exit(-1);
      }
    }
    /* get reactants of point pairs, store reactant coordinate */
    for (i=0; i<pm->num_same_coordinate_reactant_product_pairs; i++) {
      /* target_species is Reactant species, first check source point, if not ok, then look within kernel radius */
      target_species = p->arg_species[pm->same_coordinate_reactant_product_pair[i][0]];
      if ((pm->source_argindex >= 0) && (pm->argind_reactant_count[pm->source_argindex] >0) && (p->arg_species[pm->source_argindex] == target_species)) {
	pptr[i] = li->p;
      }
      else {
	/* find kernel and store the position of Reactant (assumption there is only one kernel where target_species is reactant) */
	if ((ki = getKernelIndex(target_species, REACTANT, ckernel, p->cind, p->num_cind, p->ctype1, p->ctype2)) < 0) {
	  fprintf(stderr,"ERROR (%s): cannot find ckernel for Reactant-Product pair %d/%d, species %d.\n",thisfunction,i+1,pm->num_same_coordinate_reactant_product_pairs,target_species);
	  return (-1);
	}
	item = getOnePointWithinKernel(target_species, li->p->position, &(ckernel[ki]), li, s, &cj);
	pptr[i] = item->p;
      }
    }
  }

  /* Products: add point to Configuration */
  /* process first products, in case there are same_coordinate_reactant_product_pairs, reactant coordinate has not yet been changed */
  /* all products are first put into waiting list and moved to configuration only after reactants have been removed */

  hi=0;
  for (i=0; i<pm->num_args; i++) {
    if (pm->argind_product_count[i] > 0) {
      target_species = p->arg_species[i];

      item = newPoint(&garbage, s->dimension, num_ckernels, num_processes);
      item->p->species = target_species;

      /* get position of new point */

      /* check if it is position of reactant */
      cj = -1;
      for (j=0; (j<pm->num_same_coordinate_reactant_product_pairs) && (cj<0); j++) {
	if (pm->same_coordinate_reactant_product_pair[j][1] == i) {
	  cj=j;	    
	}
      }
      if (cj >= 0) {
	/* exact position based on reactant */ 
	for (j=0; j<s->dimension; j++)
	  item->p->position[j] = pptr[cj]->position[j];
      }
      else {
	if ((ki = getKernelIndex(target_species, PRODUCT, okernel, p->oind, p->num_oind, p->otype1, p->otype2)) < 0) {
	  fprintf(stderr,"ERROR (%s): cannot find okernel for Product %d/%d, species %d.\n",thisfunction,i+1,pm->num_products,target_species);
	  return (-1);
	}

	if (li == NULL) {
	  /* position anywhere in space */
	  randomPosition(NULL, &(okernel[ki]), s->dimension, s->U, item->p->position);	
	}
	else {
	  /* random position within kernel */
	  randomPosition(li->p->position, &(okernel[ki]), s->dimension, s->U, item->p->position);
	}
      }

      /* put new item to waiting list, add to configuration after Reactants have been removed */

      hold_item_ptr[hi] = item;
      hi++;
    }
  }
  
  /* Reactants: remove point from Configuration */  
  /* first check source species (no need to calculate distances to find the point) */

  processed_one = -1;
  if ((pm->source_argindex >= 0) && (pm->argind_reactant_count[pm->source_argindex] > 0)) {
    removePoint(li, &(s->cell[ci]), &garbage);
    s->species_count[li->p->species]--;
    s->species_count[0]--;
    removeEffects(li->p, ci, s, ckernel, num_ckernels, okernel, num_okernels, processes, num_processes);
    processed_one = pm->source_argindex;
  }

  for (i=0; i<pm->num_args; i++) {
    
    if (pm->argind_reactant_count[i] > 0) {
      /* get species point within kernel, first determine which kernel to use */
      /* simplifying assumption for Reactants: max 1 point per species */
      target_species = p->arg_species[i];
      if (i != processed_one) {
	for (ki=-1, j=0; (j<p->num_cind) && (ki<0); j++) {
	  if ( ((ckernel[p->cind[j]].species1 == target_species) && (p->ctype1[j] == REACTANT)) || ((ckernel[p->cind[j]].species2 == target_species) && (p->ctype2[j] == REACTANT)) ) {
	    ki = p->cind[j];
	  }
	}
	if (ki < 0) {
	  fprintf(stderr,"ERROR (%s): cannot find ckernel for Reactant %d/%d, species %d.\n",thisfunction,i+1,pm->num_reactants,target_species);
	  return (-1);
	}
	item = getOnePointWithinKernel(target_species, li->p->position, &(ckernel[ki]), li, s, &cj);
	removePoint(item, &(s->cell[cj]), &garbage);
	s->species_count[item->p->species]--;
	s->species_count[0]--;
	removeEffects(item->p, cj, s, ckernel, num_ckernels, okernel, num_okernels, processes, num_processes);
      }
    }
  }
  
  /* add Products to Configuration */

  for (i=0; i<hi; i++) {
    /* 1) find cell, 2) add effects, and 3) add item to the end of point list, NOTE: do this in this order!!! */

    item = hold_item_ptr[i];
    cj = getCellIndex(item->p->position, s);
    addEffects(item->p, cj, s, ckernel, num_ckernels, okernel, num_okernels, processes, num_processes);
    appendList(item, &(s->cell[cj]));
    s->species_count[item->p->species]++;
    s->species_count[0]++;
  }

  return (0);
}

int writePointFile(char *ofile, SimulationSpace *s) {
  FILE *fp;
  int i,d;
  LinkedList *li;
  char *thisfunction = "writePointFile";

  if ((fp = fopen(ofile,"w")) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot open file '%s' for writing.\n",thisfunction,ofile);
    perror(""); exit(-1);
  }

  for (i=0; i<s->num_cells; i++) {
    for (li=s->cell[i].first; li != NULL; li=li->next) {
      fprintf(fp,"%d",li->p->species);
      for (d=0; d<s->dimension; d++) 
        fprintf(fp," %f",li->p->position[d]);
      fprintf(fp,"\n");
    }
  }

  if (ferror(fp)) {
    fprintf(stderr,"ERROR (%s): cannot write to '%s'.\n",thisfunction,ofile);
    perror(""); exit(-1);
  }
  fclose(fp);

  return(0);
}

int writeConfigurationToFile(FILE *fp, double t, int num_events, SimulationSpace *s) {
  int i,d;
  LinkedList *li;

  fprintf(fp,"%f %d",t,num_events);
  for (i=0; i<s->num_cells; i++) {
    for (li=s->cell[i].first; li != NULL; li=li->next) {
      fprintf(fp," %d",li->p->species);
      for (d=0; d<s->dimension; d++) 
        fprintf(fp," %f",li->p->position[d]);
    }
  }
  fprintf(fp,"\n");

  if (ferror(fp))
    return(-1);
  
  return(0);
}

int writeSpeciesCountToFile(FILE *fp, double t, int num_events, SimulationSpace *s) {
  int i,sum;

  fprintf(fp,"%f\t%d",t,num_events);
  for (i=0; i<s->num_species_count; i++) 
    fprintf(fp,"\t%d",s->species_count[i]);
  fprintf(fp,"\n");

  if (ferror(fp))
    return(-1);
  
  return(0);
}

int writeNumberOfEventsToFile(FILE *fp, double t, int num_events, int *process_count, int num_processes) {
  int i;

  fprintf(fp,"%f\t%d",t,num_events);
  for (i=0; i<num_processes; i++) 
    fprintf(fp,"\t%d",process_count[i]);
  fprintf(fp,"\n");

  if (ferror(fp))
    return(-1);  
  return(0);
}

int gillespie(SimulationSpace *s, Process *process, int num_processes, ModelComponent *m, Kernel *ckernel, int num_ckernels, Kernel *okernel, int num_okernels, double max_time, int max_num_events, long int randomseed, char *ofile, double dT, int stop_empty, char *ofile_separate, int outc, int outn, int oute, int info) {
  int num_events, pi, ci, ok, *process_count,i, ofile_s_count=0;
  double t,tau, r1,r2, a0, tout;
  LinkedList *li;
  FILE *fpc, *fpn, *fpe;
  char ofile_c[MAXSTRINGLEN], ofile_n[MAXSTRINGLEN], ofile_e[MAXSTRINGLEN], ofile_s[MAXSTRINGLEN];
  char *thisfunction = "gillespie";

  printf("gillespie!!! domain size U=%f (dimension %d), cell_width=%f\n",s->U, s->dimension, s->cell_width);
  if ((process_count = (int *) malloc (num_processes * sizeof(int))) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot malloc %d ints for process count.\n",thisfunction,num_processes);
    perror(""); exit(-1);
  }
  for (pi=0; pi<num_processes; pi++)
    process_count[pi] = 0;

  if (oute) {
    /* write event explanation file */
    sprintf(ofile_e,"%s.eventnames.txt",ofile);
    if ((fpe = fopen(ofile_e,"w")) == NULL) {
      fprintf(stderr,"ERROR (%s): cannot open file '%s' for writing.\n",thisfunction,ofile_e);
      perror(""); exit(-1);
    }
    for (pi=0; pi<num_processes; pi++) {
      if (fprintf(fpe,"e%d\t%s\t[%s",pi+1,m[pi].name,m[pi].args[0]) < 0) {
	fprintf(stderr,"ERROR (%s): cannot write to file '%s'.\n",thisfunction,ofile_e);
	perror(""); exit(-1);
      }
      for (i=1; i<m[pi].num_args; i++) {
	if (fprintf(fpe,", %s",m[pi].args[i]) < 0) {
	  fprintf(stderr,"ERROR (%s): cannot write to file '%s'.\n",thisfunction,ofile_e);
	  perror(""); exit(-1);
	}
      }
      if (fprintf(fpe,"]\n") < 0) {
	fprintf(stderr,"ERROR (%s): cannot write to file '%s'.\n",thisfunction,ofile_e);
	perror(""); exit(-1);
      }
    }
    fclose(fpe);
  }

  if (oute) {
    sprintf(ofile_e,"%s.events",ofile);
    if ((fpe = fopen(ofile_e,"w")) == NULL) {
      fprintf(stderr,"ERROR (%s): cannot open file '%s' for writing.\n",thisfunction,ofile_e);
      perror(""); exit(-1);
    }
    if (fprintf(fpe,"time\tnum.events") < 0) {
      fprintf(stderr,"ERROR (%s): cannot write to file '%s'.\n",thisfunction,ofile_e);
      perror(""); exit(-1);
    }
    for (i=0; i<num_processes; i++) {
      if (fprintf(fpe,"\te%d",i+1) < 0) {
	fprintf(stderr,"ERROR (%s): cannot write to file '%s'.\n",thisfunction,ofile_e);
	perror(""); exit(-1);
      }
    }
    if (fprintf(fpe,"\n") < 0) {
      fprintf(stderr,"ERROR (%s): cannot write to file '%s'.\n",thisfunction,ofile_e);
      perror(""); exit(-1);
    }
  }

  if (outc) {
    sprintf(ofile_c,"%s.points",ofile);
    if ((fpc = fopen(ofile_c,"w")) == NULL) {
      fprintf(stderr,"ERROR (%s): cannot open file '%s' for writing.\n",thisfunction,ofile_c);
      perror(""); exit(-1);
    }
  }

  if (outn) {
    sprintf(ofile_n,"%s.counts",ofile);
    if ((fpn = fopen(ofile_n,"w")) == NULL) {
      fprintf(stderr,"ERROR (%s): cannot open file '%s' for writing.\n",thisfunction,ofile_n);
      perror(""); exit(-1);
    }
    if (fprintf(fpn,"time\tnum.events\tnum.points") < 0) {
      fprintf(stderr,"ERROR (%s): cannot write to file '%s'.\n",thisfunction,ofile_n);
      perror(""); exit(-1);
    }
    for (i=1; i<s->num_species_count; i++) {
      if (fprintf(fpn,"\ts%d",i) < 0) {
	fprintf(stderr,"ERROR (%s): cannot write to file '%s'.\n",thisfunction,ofile_n);
	perror(""); exit(-1);
      }
    }
    if (fprintf(fpn,"\n") < 0) {
      fprintf(stderr,"ERROR (%s): cannot write to file '%s'.\n",thisfunction,ofile_n);
      perror(""); exit(-1);
    }
  }

  srand48(randomseed);
  t=0.0;
  num_events=0;
  tout=t;
  if (dT >= 0.0) {
    if (outc) {
      if (writeConfigurationToFile(fpc, tout, num_events, s) < 0) {
	fprintf(stderr,"ERROR (%s): cannot write initial state to file '%s' (event %d).\n",thisfunction,ofile_c, num_events);
	perror(""); exit(-1);
      }
    }
    if (outn) {
      if (writeSpeciesCountToFile(fpn, tout, num_events, s) < 0) {
	fprintf(stderr,"ERROR (%s): cannot write species count to file '%s' (event %d).\n",thisfunction,ofile_n, num_events);
	perror(""); exit(-1);
      }
    }
    if (oute) {
      if (writeNumberOfEventsToFile(fpe, tout, num_events, process_count, num_processes) < 0) {
	fprintf(stderr,"ERROR (%s): cannot write to file '%s' (event %d).\n",thisfunction,ofile_e, num_events);
	perror(""); exit(-1);
      }
    }
    if (ofile_separate) {
      ofile_s_count++;
      sprintf(ofile_s,"%s%06d",ofile_separate,ofile_s_count);
      writePointFile(ofile_s, s);
    }
  }
  
  ok=1;
  if (stop_empty && (s->species_count[0] <= 0))
    ok=0;

  while ((t < max_time) && (num_events < max_num_events) && ok) {

    /* when next event occurs */

    r1=drand48();
    for (a0=0.0, pi=0; pi<num_processes; pi++) 
      a0 += s->pactivity[pi];

    if (a0 < 1e-10) {
      ok = 0;
      printf("%s: time %f event %d, a0 below threshold, %g\n",thisfunction, t,num_events,a0);
    }

    if (ok) {
      tau = -log(r1)/a0;
      t=t+tau;

      /* tout moves by fixed time steps dT */
      if (dT > 0.0) {
	while (tout+dT < t) {
	  tout += dT;
	  if (outc) {
	    if (writeConfigurationToFile(fpc, tout, num_events, s) < 0) {
	      fprintf(stderr,"ERROR (%s): cannot write to file '%s' (event %d).\n",thisfunction,ofile_c, num_events);
	      perror(""); exit(-1);
	    }
	  }
	  if (outn) {
	    if (writeSpeciesCountToFile(fpn, tout, num_events, s) < 0) {
	      fprintf(stderr,"ERROR (%s): cannot write to file '%s' (event %d).\n",thisfunction,ofile_n, num_events);
	      perror(""); exit(-1);
	    }
	  }
	  if (oute) {
	    if (writeNumberOfEventsToFile(fpe, tout, num_events, process_count, num_processes) < 0) {
	      fprintf(stderr,"ERROR (%s): cannot write to file '%s' (event %d).\n",thisfunction,ofile_e, num_events);
	      perror(""); exit(-1);
	    }
	  }
	  if (ofile_separate) {
	    ofile_s_count++;
	    sprintf(ofile_s,"%s%06d",ofile_separate,ofile_s_count);
	    writePointFile(ofile_s, s);
	  }
	}
      }

      num_events++;
      
      /* which process */
      
      r2=a0*drand48();
      pi=getPactivityArrayIndex(r2,s->pactivity,num_processes);
      
      if (info) {
	printf("time %f event %d process %d '%s' (tau %g, a0 %g, r1 %g) a:",t, num_events,pi+1,process[pi].pm->name,tau,a0,r1);
	for (i=0; i<num_processes; i++)
	  printf(" %g",s->pactivity[i]);
	printf("\n");
      }      

      /* where process occurs: cell and individual */
      
      li = NULL;
      if (process[pi].source_species >= 0) {
	r2=s->pactivity[pi]*drand48();
	ci=getPactivityCellIndex(r2,pi,s->cell,s->num_cells);
	
	r2=s->cell[ci].pactivity[pi]*drand48();
	li=getPactivityListItem(r2,pi,s->cell[ci].first);
      }
      
      /* update configuration according to process[pi] */
      process_count[pi]++;
      updateConfiguration(s,ci,li, pi, process, num_processes, ckernel, num_ckernels, okernel, num_okernels);
      
      if (stop_empty && (s->species_count[0] <= 0))
	ok=0;
      
      if (info>1)
	printSimulationSpace(s, num_ckernels, num_okernels, num_processes);
    }
  }

  if (outc) {
    /* write final state of configuration */
    if (writeConfigurationToFile(fpc, t, num_events, s) < 0) {
      fprintf(stderr,"ERROR (%s): cannot write final state to file '%s' (event %d).\n",thisfunction,ofile_c, num_events);
      perror(""); exit(-1);
    }
    fclose(fpc);
  }
  if (outn) {
    if (writeSpeciesCountToFile(fpn, t, num_events, s) < 0) {
      fprintf(stderr,"ERROR (%s): cannot write to file '%s' (event %d).\n",thisfunction,ofile_n, num_events);
      perror(""); exit(-1);
    }
    fclose(fpn);
  }  
  if (oute) {
    if (writeNumberOfEventsToFile(fpe, t, num_events, process_count, num_processes) < 0) {
      fprintf(stderr,"ERROR (%s): cannot write to file '%s' (event %d).\n",thisfunction,ofile_e, num_events);
      perror(""); exit(-1);
    }
    fclose(fpe);
  }
  if (ofile_separate) {
    ofile_s_count++;
    sprintf(ofile_s,"%s%06d",ofile_separate,ofile_s_count);
    writePointFile(ofile_s, s);
  }

  /* show process event histogram */
  printf("Number of events per process:\n");
  for (pi=0; pi<num_processes; pi++) {
    printf("%d)\t%s[%s",pi+1,m[pi].name,m[pi].args[0]);
    for (i=1; i<m[pi].num_args; i++) 
      printf(",%s",m[pi].args[i]);
    printf("]\t%d\n",process_count[pi]);
  }

  return (0);
}

char *get_argument(char *key, int mandatory, char *default_value, char **argv, int argc, char *usage) {
  int i;
  char *ret = NULL;

  /* argv[0] is name of executable */
  for (i=1; i<argc; i++) {
    if (!strcmp(argv[i], key)) {
      if (i+1 < argc) {
        ret = argv[i+1];
      }
      break;
    }
  }
  if (ret == NULL) {
    if (mandatory) {
      fprintf(stderr,"ERROR: cannot find required argument '%s'.\nusage: %s\n",key, usage);
      exit(-1);
    }
    else {
      ret = default_value;
    }
  }

  return(ret);
}

int main (int argc, char **argv) {
  char *file_p, *file_m, *file_i, *file_o, *ofile_separate;
  ProcessModule *p;
  ModelComponent *m;
  Kernel *ckernels, *okernels;
  int num_ckernels, num_okernels;
  Process *processes;
  int num_process_modules, num_model_functions, num_processes, i;
  Point *ipoints;
  int num_ipoints, dimension, max_num_events, stop_empty, outc, outn, oute;
  double max_time, w, U, dT;
  SimulationSpace *s;
  long int random_seed;
  int info;
  char *usage = "PointProcessSimulator version 1.0\nrequired arguments:\n -p filename process definitions\n -m filename model components\n -i filename initial configuration\n -o filename basename of output\n -U integer/float space length\n -T integer max time\noptional arguments:\n -dT integer/float interval for saving the state of configuration, negative value means saving only final output (default -1)\n -osep filename basename for separate point coordinate files (default none, if defined, this will set -outc 0)\n -outc 1/0 output merged point coordinate file '.points', 1: yes, 0: no (default 1)\n -outn 1/0 output point count file '.counts', 1: yes, 0: no (default 1)\n -oute 1/0 output event count file '.events', 1: yes, 0: no (default 1)\n -w integer/float cell width (default min kernel radius), note: w will be rounded up so that U is its multiple\n -r integer seed value for random number generator (default 1)\n -E integer max number of events (default 2147483647)\n -s 1/0 stop if extinction of points, 1: yes, 0: no (default 1)\n -info integer (default 0)\n";
  
  if (argc < 4) {
    printf("%s",usage);
    exit(0);
  }

  file_p = get_argument("-p",REQUIRED,NULL,argv,argc,usage);
  file_m = get_argument("-m",REQUIRED,NULL,argv,argc,usage);
  file_i = get_argument("-i",REQUIRED,NULL,argv,argc,usage);
  file_o = get_argument("-o",REQUIRED,NULL,argv,argc,usage);
  U = atof(get_argument("-U",REQUIRED,NULL,argv,argc,usage));
  max_time = atof(get_argument("-T",REQUIRED,NULL,argv,argc,usage));
  dT = atof(get_argument("-dT",OPTIONAL,"-1",argv,argc,usage));
  ofile_separate = get_argument("-osep",OPTIONAL,NULL,argv,argc,usage);
  outc = atoi(get_argument("-outc",OPTIONAL,"1",argv,argc,usage));
  outn = atoi(get_argument("-outn",OPTIONAL,"1",argv,argc,usage));
  oute = atoi(get_argument("-oute",OPTIONAL,"1",argv,argc,usage));
  w = atof(get_argument("-w",OPTIONAL,"-1",argv,argc,usage));
  max_num_events = atoi(get_argument("-E",OPTIONAL,"2147483647",argv,argc,usage));
  stop_empty = atoi(get_argument("-s",OPTIONAL,"1",argv,argc,usage));
  random_seed = atol(get_argument("-r",OPTIONAL,"1",argv,argc,usage));
  info = atoi(get_argument("-info",OPTIONAL,"0",argv,argc,usage));

  if (ofile_separate != NULL)
    outc = 0;

  p = readProcessFile(file_p, &num_process_modules);
  m = readModelFile(file_m, &num_model_functions);
  if ((processes = createProcessArray(p,num_process_modules,m,num_model_functions, &ckernels, &num_ckernels, &okernels, &num_okernels)) == NULL) {
    fprintf(stderr,"processDefinitionFile: '%s'\nmodelFile: '%s'\n",file_p,file_m);
    exit(-1);
  }

  ipoints = readPointFile(file_i, &num_ipoints, &dimension);

  if (assignKernelFamily(ckernels, num_ckernels, dimension) < 0) {
    fprintf(stderr,"ckernels, processDefinitionFile: '%s'\nmodelFile: '%s'\n",file_p,file_m);
    exit(-1);
  }
  if (assignKernelFamily(okernels, num_okernels, dimension) < 0) {
    fprintf(stderr,"okernels, processDefinitionFile: '%s'\nmodelFile: '%s'\n",file_p,file_m);
    exit(-1);
  }
  
  if (w < 0.0) {
    if (num_ckernels == 0) {
      w = U;
    }
    else {
      /* find smallest connectivity kernel radius */
      w = U;
      for (i=0; i<num_ckernels; i++) {
	if (ckernels[i].radius < w)
	  w = ckernels[i].radius;
      }
    }
  }
  /* round up w so that U is its multiple */
  if (w < U) {
    w = U/floor(U/w);
  }
  else {
    w = U;
  }
  
  if ((s = createSimulationSpace(U,w, dimension, ipoints, num_ipoints)) == NULL) {
    fprintf(stderr,"inputPointFile: '%s'.\n",file_i);
    exit(-1);
  }

  initializeActivities(s, processes, num_model_functions, ckernels, num_ckernels, okernels, num_okernels);

  if (info > 1) {
    printProcessModules(p,num_process_modules);
    printModelFunctions(m, num_model_functions);
    printf("Connectivity kernels:\n");
    printUniqueKernels(ckernels, num_ckernels, processes);
    printf("Other kernels:\n");
    printUniqueKernels(okernels, num_okernels, processes);
    printPoints(ipoints, num_ipoints, dimension);
    printSimulationSpace(s, num_ckernels, num_okernels, num_model_functions);
    printProcesses(processes, num_model_functions, ckernels, okernels);
  }

  gillespie(s, processes, num_model_functions, m, ckernels, num_ckernels, okernels, num_okernels, max_time, max_num_events, random_seed, file_o, dT, stop_empty, ofile_separate, outc, outn, oute, info);

  return(0);
}
