/*
-------------------------------------------------------------------  
 A program for fast parsing any valid mmcif files. (HY,2013-01-10)
-------------------------------------------------------------------
      How to use the cif parser
      
Here is an example of how to use the cifparse:
A). Include the module in your c program  (like #include "cifparse.c")

B). To parse non-looped item:
1. declare the cif table by adding a line cifparse(ciffile, "_cell.");
2. parse the value by adding a line, c=parse_value("_cell.length_c"));
Here, c is a string pointer (must be declared as  char *c)

C). To parse the looped item:
1. declare the cif table by adding a line cifparse(ciffile, "_atom_site.");
2. parse the value by adding a line, x=parse_values("_atom_site.Cartn_x",&n);
Here, x is a double string pointer (must be declared as char **x).

To display the values 
for (i=0; i<n; i++){ printf("row %d %s\n",  i,x[i]);}

If you don't know it is a loop or non-loop item, you can print parse.loop
after declaring a table. parse.loop=0 for non-looped, parse.loop=1 looped.

For a successful parsing, cifparse(ciffile, "_cell.")=1, else =0.

-------------------------------------------------------------------
              Defination of valid mmcif format:
1. for non-looped items, below are the valid syntaxes
_citation.title "The Structure and Evolution of the Major Capsid Protein"
or 
_citation.title 'The Structure and Evolution of the Major Capsid Protein'
or 
_citation.title
;The Structure and Evolution of
the Major Capsid Protein
;

2. for looped items, below are the valid syntaxes
loop_
_citation_author.citation_id 
_citation_author.name 
_citation_author.ordinal 
primary '???, N.' 1 
primary '???, A.'     2 

or
primary
"???, N." 1 
primary
'???, A.'     2 

or
primary
;???, N.
;
1 
primary
'???, A.'     2 

-------------------------------------------------------------------
*/


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

struct
{
    int nitem, nss, nrow, loop;
    char **items, **ss; 
}parse;

int  cifparse(char *file,  const char *cate);
char *parse_value(const char *item);
char **parse_values(const char *item, int *nrow);
char **get_ciftables(char *file, int *ntable);

/* for testing 
int main(int argc, char **argv)
{
    char **tables,  **tmp=NULL, inpfile[512], outfile[512],command[512];
    int i, ntmp;
    
    strcpy(outfile, "CIF2PDB.pdb");
    if (argc==3)strcpy(outfile, argv[2]);
    sprintf(command,"rm -f %s", outfile);
    system(command);

 
    tables=get_ciftables(argv[1],&ntmp);
    for (i=0; i<ntmp; i++){
        cifparse(argv[1], tables[i]);
    }
   
    cifparse(argv[1], "_citation.");
    tmp=parse_values("_citation.journal_abbrev", &ntmp);
    for (i=0; i<ntmp; i++){
        printf("row %d %d %s\n", parse.nrow, i,tmp[i]);
    }
    cifparse(argv[1], "_cell.");
    printf("cell-id=%s\n",parse_value("_cell.entry_id"));
    printf("cell-c=%s\n",parse_value("_cell.length_c"));
    printf("cell-a=%s\n",parse_value("_cell.angle_Alpha"));

    cifparse(argv[1], "_entity.");
    tmp=parse_values("_entity.type", &ntmp);
    for (i=0; i<ntmp; i++){
        printf("row %d %d %s\n", parse.nrow, i,tmp[i]);
    }

    cifparse(argv[1], "_pdbx_database_related.");
    tmp=parse_values("_pdbx_database_related.details", &ntmp);
    for (i=0; i<ntmp; i++){ printf("row %d %d %s\n", parse.nrow, i,tmp[i]);}

    cifparse(argv[1], "_entity_name_com.");
    tmp=parse_values("_entity_name_com.name", &ntmp);
    for (i=0; i<ntmp; i++){ printf("row %d %d %s\n", parse.nrow, i,tmp[i]);}

    cifparse(argv[1], "_entity_poly.");
    tmp=parse_values("_entity_poly.pdbx_seq_one_letter_code", &ntmp);
    for (i=0; i<ntmp; i++){ printf("row %d %d %s\n", parse.nrow, i,tmp[i]);}

    cifparse(argv[1], "_atom_type.");
    tmp=parse_values("_atom_type.symbol", &ntmp);
    for (i=0; i<ntmp; i++){ printf("row %d %d %s\n", parse.nrow, i,tmp[i]);}

    cifparse(argv[1], "_atom_site.");
    tmp=parse_values("_atom_site.Cartn_x", &ntmp);
    for (i=0; i<ntmp; i++){ printf("row %d %d %s\n", parse.nrow, i,tmp[i]);}

    cifparse(argv[1], "_refln.");
    tmp=parse_values("_refln.F_meas_sigma_au", &ntmp);
    for (i=0; i<ntmp; i++){ printf("row %d %d %s\n", parse.nrow, i,tmp[i]);}

}
*/  

void rid_of_char(char *str, char ch)
/* get rid of space*/
{
    int i=0, j=0, len=0;

    len=strlen(str);
    for(i=0; i<len; i++){
        if(str[i]==ch) continue;
        str[j++]= str[i];
    }
    str[j]='\0';
    if(iscntrl(str[j-1]))str[j-1]='\0';
}

char ** string_token(char *str, const char *token, int *nstr)
/* use characters in 'token' to parse. return a double pointer
*/
{
    char *tokenPtr, **line=NULL, *tmp=NULL;
    int i=0;

    rid_of_char(str, '\n');

    tmp = (char *) malloc( (strlen(str)+1) * sizeof (char));
    strcpy(tmp,str);

    tokenPtr = strtok(tmp, token); //get number of token
    i = 0;
    while(tokenPtr != NULL){
        tokenPtr = strtok(NULL, token);
        i++;
    }
    free(tmp);

    line = (char **) malloc( (i+1) * sizeof (char *));

    tokenPtr = strtok(str, token);
    i = 0;
    while(tokenPtr != NULL){
        line[i] = (char *) malloc( (strlen(tokenPtr)+1) * sizeof (char));
        strcpy(line[i], tokenPtr);
        tokenPtr = strtok(NULL, token);
        i++;
    }
    *nstr = i;
    return (line);
}

int is_space_line(const char *str)
// n=0, is a space line; n!=0 NOT.
{
    int i, n=0, len;

    len = strlen(str);
    for(i=0; i<len; i++){
        if(isspace(str[i])) continue;
        n++;
    }
    return n;
}

void rid_of_front_space(char *str)
/* get rid of first space*/
{
    int i=0, j=0, len=0, n=0;

    len=strlen(str);
    for(i=0; i<len; i++){
        if(isspace(str[i]))
            continue;
        else{
            j=i;
            break;
        }
    }
    
    for(i=j; i<len; i++){
        str[n++]=str[i];
    }
    str[n]='\0';
}

char **get_lines_from_file(char *inpfile, int *nline)
/*Input a file and return a double pointer, no space between lines
*/
{
    char str[2000], **all_line=NULL;
    int m;
    FILE *fp=NULL;

    if((fp=fopen(inpfile, "r"))==NULL){
        printf("Can not open file (%s) in get_lines_from_file\n", inpfile);
        *nline=0;
        return NULL;
    }

    m=0;
    while (fgets(str, sizeof str, fp) != NULL) {
        if(is_space_line(str)<=0)continue;
        m++;
    }
    rewind(fp);
    all_line= (char **) malloc( (m+2) * sizeof (char *)); 

    m=0;
    while (fgets(str, sizeof str, fp) != NULL) {
        if(is_space_line(str)<=0)continue;
        rid_of_front_space(str);
        all_line[m] = (char *) malloc( (strlen(str)+1) * sizeof (char )); 
        strcpy(all_line[m],str);
        m++;
    }
    *nline =m;
    fclose(fp);
    return (all_line);
}

void rid_of_front_end_space(char *str)
{
    int i=0, j=0, len=0, n=0;

    len=strlen(str);
    for(i=0; i<len; i++){  //front space
        if(isspace(str[i]))
            continue;
        else{
            j=i;
            break;
        }
    }
    n=0;
    for(i=j; i<len; i++){
        str[n++]=str[i];
    }
    str[n]='\0';

    if(str[0]=='\0') return;

    n=0;
    len=strlen(str);
    for(i=len-1; i>=0; i--){ //end space
        if(isspace(str[i]))
            n++;
        else
            break;
    }
    str[len-n]='\0';
    if(iscntrl(str[len-n-1]))str[len-n-1]='\0';
}


int is_cif(const char *inpfile)
/* all the possible categories are separated by a space */
{
    char str[512], tmp[512];
    int ncat=0,  n=0;
    FILE *fp;

    fp = fopen(inpfile, "r");

    while(fgets(str, sizeof str, fp)!= NULL){
        rid_of_front_end_space(str);

        sscanf(str,"%s",tmp);
        if ((tmp[0]=='_' && strchr(tmp,'.'))|| !strncmp(str,"data_",5) ||
            !strncmp(str,"_loop",5)){
            ncat++;
            if(ncat>2){
                n=1;
                break;
            }
        }
    }
    fclose(fp);
    return n;

}

void string_between_char(char *tmp1, char cc, char *strnew)
{
        /*parse string between char cc*/
    
    int m, nlen;

    strcpy(strnew, strchr(tmp1, cc)+1);
    
    nlen=strlen(strnew);
    for (m=0;  m < nlen; m++) {
        if (m>0 && strnew[m]==cc){
            strnew[m]='\0';
            break;
        }
    }
    
}

void string_after_char(int *k, char *tmp1, char cc, char *strnew)
{
        /*parse string between char cc*/
    
    int m,n, nlen;

    nlen=strlen(tmp1);
    strcpy(strnew, "");
    
    if(isspace(cc)){
        n=0;
        for (m=*k;  m < nlen; m++) {
            if(n>0 && isspace(tmp1[m])){
                *k=m;
                strnew[n]='\0';
                break;
            }else if (!isspace(tmp1[m])){
                strnew[n]=tmp1[m];
                n++;
            }
        }
        
    }else{
        n=0;
        for (m=*k+1;  m < nlen; m++) {
            if (n>0 && tmp1[m]==cc ){
                *k=m+1;
                strnew[n]='\0';
                break;
            }else if (tmp1[m]!=cc){
                strnew[n]=tmp1[m];
                n++;
            }
        }
    }
}

void strings_between_char(int nline, char **line, int *k, char cc, char *ssnew)
{
    int m,n=0;

    strcpy(ssnew, "");
    for (m=*k;  m < nline; m++) {
        if (n>0 && line[m][0]==cc ){
            *k=m;
            break;
        }
        strcat(ssnew, line[m]);
        n++;

        if (!strncmp(line[m], "loop_",5)){
            printf("Error: can not find the second ;. Stop here!\n");
            *k=m;
            break;
        }
        
    }
    if (ssnew[0]==cc) ssnew[0]=' ';
    
}


void clean_str(char *ss)
{
    int i,n=0, slen;

    slen = strlen(ss);
    
    for (i=0;  i < slen; i++) {
        if (isspace(ss[i])) ss[i]=' ';
    }
    
    for (i=0;  i < slen; i++) {
        if (i>0 && ss[i-1]==' ' && ss[i]== ' ') continue;
        ss[n] = ss[i];
        n++;
    }
    ss[n]='\0';
    
}
    
char **get_ciftables(char *file, int *ntable)
{
    char **line, **tables, tmp[8000],stab[500];
    int nline, ntab=0, i, n=0;

    tables = (char **) malloc(5000 * sizeof (char *));

    line=get_lines_from_file(file, &nline);

    for (i=0;  i < nline; i++) {
        sscanf(line[i], "%s", tmp);
        rid_of_front_end_space(tmp);
        if (!(tmp[0] == '_' && strchr(tmp,'.')))continue;
        n=strchr(tmp, '.')-tmp+1;
        strncpy(stab, tmp, n);
        stab[n]='\0';
        if (ntab==0){
            tables[ntab]=(char *) malloc(strlen(stab)+5 * sizeof (char ));
            strcpy(tables[ntab], stab);
            printf("table=%3d  %s\n",ntab, tables[ntab]);
            ntab++;
        }else{
            if (strcmp(stab,tables[ntab-1])){
                tables[ntab]=(char *) malloc(strlen(stab)+5 * sizeof (char ));
                strcpy(tables[ntab], stab);
                printf("table=%3d  %s\n",ntab, tables[ntab]);
                ntab++;
            }
        }
    }
    *ntable=ntab;
    return tables;
}
int  cifparse(char *file,  const char *cate) 
{/* parse cate */
    char **line, tmp[800000],tmp1[800], tmp2[800], **ss, **items,cc;
    int nline, i,j, k, start=0,n=0, m, nt=0, nlen=0, nitem=0,loop=0,clen;


    line=get_lines_from_file(file, &nline);
    clen=strlen(cate);

    ss = (char **) malloc(500*nline * sizeof (char *));
    items = (char **) malloc(200 * sizeof (char *));
   
    m=0; 
    for (i=0;  i < nline; i++) {
        sscanf(line[i], "%s", tmp);
        rid_of_front_end_space(tmp);
        if (!strncmp(tmp, cate, clen)){
            if (i>0 && strstr(line[i-1], "loop_")){
                loop=1;
            }
            start=i;
            m++;
            break; /*find start position and loop or nonloop*/
        }
    }
    if(m==0){
        printf("cif table (%s) not found.\n", cate);
        return 0;
    }
    
    if (loop==0){/*for nonlooped items*/
        for (i=start;  i < nline; i++) {
            m=sscanf(line[i], "%s%s",tmp1, tmp2);
            rid_of_front_end_space(tmp1);
            if (m==2 && strstr(tmp1, cate)){
                if (tmp2[0]=='\''){
                    string_between_char(line[i], '\'', tmp);
                }else if(tmp2[0]=='"'){
                    string_between_char(line[i], '"', tmp);
                }else{
                    strcpy(tmp,tmp2);
                }
                ss[nitem]=(char *) malloc((strlen(tmp)+2) * sizeof (char ));
                items[nitem]=(char *) malloc((strlen(tmp1)+2) * sizeof (char ));
                strcpy(items[nitem], tmp1);
                strcpy(ss[nitem], tmp);
                nitem++;
               
            }else if (m==1 && i<nline-1 && strstr(tmp1, cate)){
                n=0;
                strcpy(tmp,"");
                if (line[i+1][0]==';'){
                    for (j=i+1;  j < nline; j++) {
                        if(n>0 && line[j][0]==';') {
                            i=j;
                            break;
                        }
                        strcat(tmp,line[j]);
                        n++;
                        if ((strlen(line[j])<6 && !strncmp(line[j], "loop_",5)) ||
                            !strncmp(line[j], cate,clen)){
                            printf("Error: can not find the second ;. Stop here!\n");
                            i=j;
                            break;
                        }
                        
                    }
                    if (tmp[0]==';') tmp[0]=' ';
                    
                }else if (line[i+1][0]=='\'' || line[i+1][0]=='"'){
                    cc=line[i+1][0];
                    strcpy(tmp,line[i+1]);
                    rid_of_front_end_space(tmp);
                    k=strlen(tmp);
                    if (tmp[0]==cc && tmp[k-1]==cc) {
                        tmp[0]=' ';
                        tmp[k-1]=' ';
                    }
                }else if (strncmp(line[i+1], cate, clen)){
                    strcpy(tmp,line[i+1]);
                }
//                printf("tmp=%s|%s|%s\n", tmp, tmp1,tmp2);
                
                
                ss[nitem]=(char *) malloc((strlen(tmp)+2) * sizeof (char ));
                items[nitem]=(char *) malloc((strlen(tmp1)+2) * sizeof (char ));
                strcpy(items[nitem], tmp1);
                strcpy(ss[nitem], tmp);
                nitem++;
                
            }else if (line[i][0]=='#' || !strncmp(line[i], "loop_",5) ||
                      (line[i][0]=='_' && !strstr(line[i], cate))){
                break;
            }
        }
        nt=nitem;
        
    } else {  /*for looped items*/
        nitem=0;
        for (i=start;  i < nline; i++) {
            m=sscanf(line[i], "%s%s",tmp1, tmp2);
            rid_of_front_end_space(tmp1);
            if (!strncmp(tmp1, cate, clen)){
                items[nitem]=(char *) malloc((strlen(tmp1)+2) * sizeof (char ));
                strcpy(items[nitem], tmp1);
                nitem++;
                
            }else if (nitem>0 && i<nline-1 && !strstr(line[i+1], cate)){
                start=i;
                break;
            }
        }
        
        for (j=start;  j < nline; j++) {
            if (line[j][0]==';'){
                strings_between_char(nline, line, &j, ';', tmp);
                ss[nt]=(char *) malloc((strlen(tmp)+2) * sizeof (char ));
                strcpy(ss[nt], tmp);
                nt++;
                
            }else if(line[j][0]=='#' || !strncmp(line[j], "loop_",5) ||
                     (line[j][0]=='_' && strchr(line[j],'.') )){
                break;
                
            }else{
                strcpy(tmp1,line[j]);
                clean_str(tmp1);
                nlen=strlen(tmp1);
                
                for (k=0;  k < nlen; k++) {
                    if ((k>0 && (tmp1[k]=='\'' ) && tmp1[k-1] ==' ') ||
                        (k==0 && tmp1[0]=='\'') ){
                        string_after_char(&k, tmp1, '\'', tmp);

                    }else if ((k>0 && (tmp1[k]=='"') && tmp1[k-1] ==' ') ||
                              (k==0 && tmp1[0]=='"') ){
                        string_after_char(&k, tmp1, '"', tmp);
                        
                    }else{
                        string_after_char(&k, tmp1, ' ', tmp);
                        rid_of_front_end_space(tmp);
                        if (strlen(tmp)==0) continue; 
                    }
                    ss[nt]=(char *) malloc((strlen(tmp)+2) * sizeof (char ));
                    strcpy(ss[nt], tmp);
                    nt++;
                }
            }
        }
    }

    parse.nss=nt;
    parse.nitem=nitem;
    parse.nrow=nt/nitem;
    parse.loop=loop;
    parse.ss=ss;
    parse.items=items;
    n=nt%nitem;
    m=(nt/nitem);
    if (n==0) {
        printf("\nparsing (%s)  ntot=%d : Row=%d : data_item=%d : loop=%d\n",cate,nt, m,nitem,loop); 
        return 1;
    }else{
        printf("\nProblem parsing (%s)  ntot=%d : Row=%d : data_item=%d : loop=%d\n",cate,nt, m,nitem,loop);
        return 0;
    } 

/*
// below for testing/printing

    for (k=0;  k < parse.nitem; k++) {
        printf("%s\n",parse.items[k]);
    }
    
    n=0;
    for (k=0;  k < parse.nss; k++) {
        printf("%s ", parse.ss[k]);
        if(++n==parse.nitem){
            printf("\n");
            n=0;
        }
    }
    n=nt/nitem;

    for (k=0;  k < parse.nrow; k++) {
        printf("%s %s \n", ss[0 + k*parse.nitem], ss[2 + k*parse.nitem]);
    }
*/
    //return loop;
    
}
   
char **parse_values(const char *item, int *nrow)
{
// parse looped values of the input item
    char **ss=NULL;
    int k, m=-1;

    printf("\n%s", item);
    for (k=0;  k < parse.nitem; k++) {
        if (!strcmp(item, parse.items[k])){
            m=k;
            break;
        }
    }
    if (m==-1){
        printf("the item (%s) is not found\n",item);
        *nrow=0;
        return NULL;
    }
    ss=(char **) malloc ((parse.nrow + 10) * sizeof (char *));
    for (k=0;  k < parse.nrow; k++) {
        ss[k]=(char *) malloc ((strlen(parse.ss[m + k*parse.nitem])+2) * sizeof (char ));
        strcpy(ss[k], parse.ss[m + k*parse.nitem]);
    }
    *nrow = parse.nrow;
    printf("\n%d", *nrow);
    return ss;
}

char *parse_value(const char *item)
{
// parse a non-loop value of the input item
    int k, m=-1;
    for (k=0;  k < parse.nitem; k++) {
        if (!strcmp(item, parse.items[k])){
            m=k;
            break;
        }
    }
    if (m==-1){
        printf("the item (%s) is not found\n",item);
        return "";
    }
    return parse.ss[m];
}

/////// below are some applications  /////
 


void cif2pdb(char *file, char *fpdb)
{ /*parse cif and convert to SIMPLE pdb*/
    char **line,tmp[200], tmp1[2000], tmp2[2000], pdbid[10], *sym, **mtype ;
    char **m11,**m12,**m13,**t1,**m21,**m22,**m23,**t2, **m31,**m32,**m33,**t3;
    int nline, i,j, k=0, nn, n=0, m=0;
    float cell[6];
    FILE *fw=NULL;


    fw=fopen(fpdb, "w");
// use _pdbx_struct_oper_list, not _pdbx_struct_legacy_oper_list

    if (cifparse(file, "_pdbx_struct_oper_list.")){
        mtype=parse_values("_pdbx_struct_oper_list.type", &nn);
        m11=parse_values("_pdbx_struct_oper_list.matrix[1][1]", &nn);
        m12=parse_values("_pdbx_struct_oper_list.matrix[1][2]", &nn);
        m13=parse_values("_pdbx_struct_oper_list.matrix[1][3]", &nn);
        t1 =parse_values("_pdbx_struct_oper_list.vector[1]", &nn);
        m21=parse_values("_pdbx_struct_oper_list.matrix[2][1]", &nn);
        m22=parse_values("_pdbx_struct_oper_list.matrix[2][2]", &nn);
        m23=parse_values("_pdbx_struct_oper_list.matrix[2][3]", &nn);
        t2 =parse_values("_pdbx_struct_oper_list.vector[2]", &nn);
        m31=parse_values("_pdbx_struct_oper_list.matrix[3][1]", &nn);
        m32=parse_values("_pdbx_struct_oper_list.matrix[3][2]", &nn);
        m33=parse_values("_pdbx_struct_oper_list.matrix[3][3]", &nn);
        t3 =parse_values("_pdbx_struct_oper_list.vector[3]", &nn);
        strcpy(tmp,"REMARK 350   BIOMT");
        j=0;
        for (i=0;  i < nn; i++) {
            
            if (!strstr(mtype[i], "point symmetry operation") &&
                !strstr(mtype[i], "helical symmetry operation") &&
                !strstr(mtype[i], "general operation") ) continue;
        
            j++;
            
            
            fprintf(fw,"%s1%4d%10.6f%10.6f%10.6f%15.5f \n",
                    tmp,j,atof(m11[i]),atof(m12[i]),atof(m13[i]),atof(t1[i]) );
            fprintf(fw,"%s2%4d%10.6f%10.6f%10.6f%15.5f \n",
                    tmp,j,atof(m21[i]),atof(m22[i]),atof(m23[i]),atof(t2[i]) );
            fprintf(fw,"%s3%4d%10.6f%10.6f%10.6f%15.5f \n",
                    tmp,j,atof(m31[i]),atof(m32[i]),atof(m33[i]),atof(t3[i]) );
        }
    }else{
        printf("Note:The infor. from_pdbx_struct_oper_list is not parsed to pdb\n");
    }

    cifparse(file, "_cell.");
    cell[0] = atof(parse_value("_cell.length_a"));
    cell[1] = atof(parse_value("_cell.length_b"));
    cell[2] = atof(parse_value("_cell.length_c"));
    cell[3] = atof(parse_value("_cell.angle_alpha"));
    cell[4] = atof(parse_value("_cell.angle_beta"));
    cell[5] = atof(parse_value("_cell.angle_gamma"));

    cifparse(file, "_symmetry.");
    sym=parse_value("_symmetry.space_group_name_H-M");

    fprintf(fw,"CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %10s \n",
            cell[0],cell[1], cell[2],cell[3],cell[4], cell[5], sym);

    char **group, **natm, **symbol, **atom, **asym, **comp, **seq, **biso, **occ;
    char **x, **y, **z, **ins,**alt, alter[10], inst[10], atomname[10];
    cifparse(file, "_atom_site.");
    group=parse_values("_atom_site.group_PDB", &nn);
    natm=parse_values("_atom_site.id", &nn);
    symbol=parse_values("_atom_site.type_symbol", &nn);
    atom=parse_values("_atom_site.auth_atom_id", &nn);
    asym=parse_values("_atom_site.auth_asym_id", &nn);
    comp=parse_values("_atom_site.auth_comp_id", &nn);
    seq=parse_values("_atom_site.auth_seq_id", &nn);
    biso=parse_values("_atom_site.B_iso_or_equiv", &nn);
    occ=parse_values("_atom_site.occupancy", &nn);
    x=parse_values("_atom_site.Cartn_x", &nn);
    y=parse_values("_atom_site.Cartn_y", &nn);
    z=parse_values("_atom_site.Cartn_z", &nn);
    ins=parse_values("_atom_site.pdbx_PDB_ins_code", &nn);
    alt=parse_values("_atom_site.label_alt_id", &nn);

    for (i=0;  i < nn; i++) {

        if(strlen(symbol[i])>1){ //
             sprintf(tmp1,"%s    ", atom[i]);
        }else {
             sprintf(tmp1," %s    ", atom[i]);
             if (strlen(atom[i]) == 4) sprintf(tmp1,"%s    ", atom[i]);
        }
        strncpy(atomname, tmp1, 4);
        atomname[4] = '\0';

        (strchr(ins[i],'.') || strchr(ins[i],'?'))? strcpy(inst," "): strcpy(inst,ins[i]);
        (strchr(alt[i],'.') || strchr(alt[i],'?'))? strcpy(alter," "): strcpy(alter,alt[i]);

        fprintf(fw,"%-6s",group[i]);
        fprintf(fw, "%5s %4s%1s%3s%1s%1s%4s    %8s%8s%8s%6s%6s",natm[i], atomname,
                alter,comp[i], inst,asym[i],seq[i], x[i], y[i],z[i],occ[i], biso[i]);
        fprintf(fw,"%12s  \n",symbol[i]);


    }
    fclose(fw);

}



void cif2pdb_old(char *file, char *fpdb)
{ /*parse cif and convert to SIMPLE pdb*/
    char **line, tmp1[2000], tmp2[2000], pdbid[10], sym[20] ;
    int nline, i,j, k=0, n=0, m=0;
    float cell[6];
    FILE *fw=NULL;

/*    sprintf(fpdb, "%s_tmp.pdb", file); */
    fw=fopen(fpdb, "w");

    line=get_lines_from_file(file, &nline);

    for (i=0;  i < nline; i++) {

        n=sscanf(line[i], "%s %s", tmp1, tmp2);

        if (!strcmp(tmp1, "_cell.entry_id")){
            strcpy(pdbid, tmp2);
        }else if (!strcmp(tmp1, "_cell.length_a")){
            cell[0]=atof(tmp2);
        }else if (!strcmp(tmp1, "_cell.length_b")){
            cell[1]=atof(tmp2);
        }else if (!strcmp(tmp1, "_cell.length_c")){
            cell[2]=atof(tmp2);
        }else if (!strcmp(tmp1, "_cell.angle_alpha")){
            cell[3]=atof(tmp2);
        }else if (!strcmp(tmp1, "_cell.angle_beta")){
            cell[4]=atof(tmp2);
        }else if (!strcmp(tmp1, "_cell.angle_gamma")){
            cell[5]=atof(tmp2);

        }else if (!strcmp(tmp1, "_symmetry.space_group_name_H-M")){
            strcpy(sym, strstr(line[i], "'"));
            rid_of_char(sym, '\'');

        }else if (strstr(tmp1,"loop_") && strstr(line[i+1], "_atom_site.")){
            k=i;
            break;
        }
    }
    fprintf(fw,"CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %10s \n",
            cell[0],cell[1], cell[2],cell[3],cell[4], cell[5], sym);


    int group=0,natm=0, symbol=0, atom=0,asym=0,comp=0,seq=0,biso=0,occ=0,x=0,y=0,z=0,ins=0,alt=0;



    n=0;
    m=0;
    for (j=k+1;  j < nline; j++) {
 /*       printf("%s", line[j]);*/
        rid_of_front_end_space(line[j]);

        if (!strcmp(line[j],"_atom_site.group_PDB")){
            group = n;
        }else if (!strcmp(line[j], "_atom_site.id")){
            natm = n;
        }else if (!strcmp(line[j], "_atom_site.type_symbol")){
            symbol = n;
        }else if (!strcmp(line[j],"_atom_site.auth_atom_id")){
            atom= n;
        }else if (!strcmp(line[j],"_atom_site.auth_asym_id")){
            asym= n;
        }else if (!strcmp(line[j],"_atom_site.auth_comp_id")){
            comp= n;
         }else if (!strcmp(line[j],"_atom_site.auth_seq_id")){
            seq= n;
        }else if (!strcmp(line[j],"_atom_site.B_iso_or_equiv")){
            biso= n;
        }else if (!strcmp(line[j],"_atom_site.occupancy")){
            occ= n;
        }else if (!strcmp(line[j],"_atom_site.Cartn_x")){
            x= n;
        }else if (!strcmp(line[j],"_atom_site.Cartn_y")){
            y= n;
        }else if (!strcmp(line[j],"_atom_site.Cartn_z")){
            z= n;
        }else if (!strcmp(line[j],"_atom_site.pdbx_PDB_ins_code")){
            ins= n;

        }else if (!strcmp(line[j],"_atom_site.label_alt_id")){
            alt= n;


        }else if (n>5 && !strstr( line[j],"_atom_site.")){
            m=j;
            break;
        }
        n++;
    }
    char **ss,atomname[10], inst[2], alter[2];
    int nn;

    for (j=m;  j < nline; j++) {
        ss=string_token(line[j], " " , &nn);
        if (nn>n-2){


            if (ss[atom][0]=='"'){
                rid_of_char(ss[atom], '"');
            }else if (ss[atom][0]=='\''){
                rid_of_char(ss[atom], '\'');
            }

            if(strlen(ss[symbol])>1){ //
                sprintf(tmp1,"%s    ", ss[atom]);
            }else {
                sprintf(tmp1," %s    ", ss[atom]);

                if (strlen(ss[atom]) == 4) sprintf(tmp1,"%s    ", ss[atom]);
            }
            strncpy(atomname, tmp1, 4);
            atomname[4] = '\0';


            (strchr(ss[ins],'.') || strchr(ss[ins],'?'))? strcpy(inst," "): strcpy(inst,ss[ins]);
            (strchr(ss[alt],'.') || strchr(ss[alt],'?'))? strcpy(alter," "): strcpy(alter,ss[alt]);

            fprintf(fw,"%-6s",ss[group]);
            fprintf(fw, "%5s %4s%1s%3s%1s%1s%4s    %8s%8s%8s%6s%6s",ss[natm], atomname,
                    alter,ss[comp], inst,ss[asym],ss[seq], ss[x], ss[y],ss[z], ss[occ], ss[biso]);

            fprintf(fw,"%12s  \n",ss[symbol]);

        }

    }

}

