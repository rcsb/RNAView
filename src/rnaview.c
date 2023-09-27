/* the input coordinate file either is PDF format or CIF */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <math.h>
#include <sys/types.h>
#include "nrutil.h"
#include "rna.h"

long PS=1, VRML=0, ANAL=0, XML=0, HETA=1, VERB=1; /*globle variables */
char FILEOUT[BUF512];



void base_edge_stat(char *pdbfile, long *A, long *U,long *G,long *C,long *T,
                    long *P, long *I);
void print_edge_stat(FILE *fs,  long *A, long *U, long *G,long *C, long *T,
                      long *P,  long *I);

void print_statistic(FILE *fstat, long *type_stat_tot, long **pair_stat);
void sixteen_pair_statistics(long num_pair_tot,long **bs_pairs_tot, char *bseq,
                             char **pair_type,long **type_stat);
void write_single_Hbond_stat(char *pdbfile, char *bseq, long **pair_stat);
void delete_file(char *pdbfile, char *extension);

long clean_inpfile(char *inpfile, char *chain_in, float reso_in, char *file_new);


int main(int argc, char *argv[])
{
    clock_t start, finish;
    char inpfile[BUF512], pdbfile_new[BUF512], chain[80], str[BUF512];
    
    long i, j,n,  key, base_all, npdb, nxml;
    long type_stat[20]; /* maxmum 20 different pairs */
    long **pair_stat; /* maxmum 20 different pairs */
    static long A[4],U[4],G[4],C[4],T[4],P[4],I[4];
    float reso;
    
    FILE  *fstat, *fp;

    start = clock();

    fstat=fopen("base_pair_statistics.out", "w");
    
    if(argc<=1 || (argc==2 && strstr(argv[1], "-h"))) usage();

    strcpy(chain,"");
    for(i=1; i<argc; i++){
        upperstr(argv[i]);        
        if(!strcmp(argv[i],"-F")|| !strcmp(argv[i],"-P")
           ||!strcmp(argv[i],"-PX")){
            strcpy(inpfile,argv[i+1]);
            strcpy(FILEOUT,argv[i+1]);
            
        }else if(!strcmp(argv[i],"-C")){
            strcpy(chain,argv[i+1]);
            upperstr(chain);
            
        }else if(!strcmp(argv[i],"-V")){
            VRML=1;
        }else if(!strcmp(argv[i],"-A")){
            ANAL=1;
        }else if(!strcmp(argv[i],"-B")){
            VERB=0;
        }else if(!strcmp(argv[i],"-R")){
            reso=atof(argv[i+1]);
        }
    }

    if((fp=fopen(inpfile, "r"))==NULL){
        printf("Error! can not open file (%s)\n", inpfile);
        return 0;
    }

    pair_stat = lmatrix(0, 20, 0, 40);
    for (i = 0; i < 20; i++){
        type_stat[i] =0;        
    }
    
 /* i for A-A ... pairs (16);  j for Leontis-Westhof base-pairs */
    for (i = 0; i <20; i++) 
        for (j = 0; j <40; j++)
            pair_stat[i][j]=0;

    npdb=0;
    nxml=0;
    while (fgets(str, sizeof str, fp) != NULL) {
        if(!strncmp(str, "ATOM", 4) || !strncmp(str, "HETATM", 4) ){
            if(npdb++>10){
                break;
            }
        }else if(strstr(str, "<!DOCTYPE rnaml") || strstr(str, "<rnaml") ||
                 strstr(str, "<molecule") || strstr(str, "<base-pair")){
            if(nxml++>10){
                XML=1;
                break;
            }
        }
    }
    
    if(nxml>10){ /* read information from RNAML and write into PS */
        xml2ps(inpfile, 0, 1); 
        fclose(fp);
        return 0;
        
    }else if(npdb>=10 && nxml<10){
        printf("Processing a single file ...\n");

        sprintf(pdbfile_new, "%s_new", inpfile);
        n=clean_inpfile(inpfile,chain, reso, pdbfile_new);
        rna(pdbfile_new, type_stat, pair_stat, &base_all);
        base_edge_stat(pdbfile_new, A, U, G, C, T, P, I);
    }else{
        printf("Processing multiple files ...\n");
        rewind(fp);
        
/* file contains a list of PDB file */        
        while (fgets(str, sizeof str, fp) != NULL) {
            sprintf(pdbfile_new, "%s_new", str);
            n=clean_inpfile(str,chain, reso, pdbfile_new);
            rna(pdbfile_new, type_stat, pair_stat, &base_all);
            base_edge_stat(pdbfile_new, A, U, G, C, T, P, I);
        }
    }
    fclose(fp);


    fprintf(fstat,"\nNumber of the total bases = %d\n", base_all);
    print_statistic(fstat, type_stat, pair_stat);
    print_edge_stat(fstat, A, U, G, C, T, P, I);
    fclose(fstat);
    
    free_lmatrix(pair_stat,0, 20, 0, 40);

    
    delete_file("", "pattern_tmp.out");
    delete_file("", "best_pair.out");
    delete_file(inpfile, "_new"); 
    delete_file(inpfile, "_new_sort.out");
    delete_file(inpfile, "_patt_tmp.out");

    
    finish = clock();
    printf( "\nTime used: %.2f seconds\n",
            ((double) (finish - start)) / CLOCKS_PER_SEC);
    fprintf(stderr, "\nJOB finished! Time used: %.2f seconds\n",
            ((double) (finish - start)) / CLOCKS_PER_SEC);

    
    return 0;    
}

long clean_inpfile(char *inpfile, char *chain_in, float reso_in, char *inpfile_new)
/* Get single pdb file with one alt conformer.
   It is selected by chain, resolution, only one model.
*/
{
    char str[256];
    int npdb=1, n;
    float reso;
    FILE *fr, *fw;

    fr=fopen(inpfile, "r");
    fw=fopen(inpfile_new, "w");

    while (fgets(str, sizeof str, fr) != NULL) {
        if(!strncmp(str, "ENDMDL", 6))break;
        if(reso_in>0.1){
            if(strstr(str, "REMARK   3   RESOLUTION RANGE HIGH (ANGSTROMS) :")){
                n=sscanf(strstr(str, ":")+1, "f", &reso);
                if(n==1 && reso>reso_in){
                    npdb=0;/*not satisfy resolution*/
                    break;
                }
            }
        }
        if(!strncmp(str, "ANISOU", 6)) continue;
        
        if(!strncmp(str, "ATOM  ", 6) || !strncmp(str, "HETATM", 6) ){
            if(strlen(chain_in)>0){ /*rid of chain out of range*/
                if(!strchr(chain_in, str[21])){
                    continue;
                }
            }

            if(str[17]==' ' && str[18]=='D'){
                if(str[19]=='A' || str[19]=='T' || str[19]=='G' || str[19]=='C'){
                    str[17]=' ';
                    str[18]=' ';
                }
            }
        }
        fprintf(fw,"%s", str);
        
    }
    fclose(fr);
    fclose(fw);
    
    return npdb;
}

void delete_file(char *pdbfile, char *extension)
{
    char command[512];
    
    strcpy(command,"rm -f ");
    strcat(command,pdbfile);
    strcat(command,extension);
    system(command);
    return;
}


void usage(void)
{
    printf("\nUsage: rnaview -f filename \n");

    printf(  "--------------------------------------------------------------\n"
             "        Options of the rnaview program\n");
    printf(
            "+-------------------------------------------------------------+\n"
            "| (1) -f followed by the file name.                           |\n"
            "|     The filename is either a single PDB file or a single XML|\n"
            "|     file or a file containing a list of many PDB files      |\n"
            "|     Example:    rnaview  -f filename                        |\n"
            "|                                                             |\n"
            "| (2) -c followed by a string of chain IDs                    |\n"
            "|     For a complicated structure, you may only select several|\n"
            "|     chains for clear display of the 2D graphic.             |\n"
            "|                                                             |\n"
            "|     Example:    rnaview  -f pdbfile -c AB                   |\n"
            "|     Here only chain A & B are selected for calculations.    |\n"
            "|     There should be no space between chain A & B.           |\n"
            "|                                                             |\n"
            "| (3) -r followed by a high resolution limit                  |\n"
            "|     This option is useful to calculate a list of pdbfiles.  |\n"
            "|     The PDB file with resolution lower than the given value |\n"
            "|     will be filtered out.                                   |\n"
            "|                                                             |\n"
            "|     Example:    rnaview  -f file.list  -r 3.0               |\n"
            "|     It means that only the pdbfiles with resolution value   |\n"
            "|     < 3.0A are selected for calculation.                    |\n"
            "|                                                             |\n"
            "| (4) -v  to generate a 3D structure in VRML format.          |\n"
            "|     It can be displayed on internet (with VRML plug in).    |\n"
            "|     Example:    rnaview  -f pdbfile  -v                     |\n"
            "|                                                             |\n"
            "| (5) -a  to calculate  Shear  Stretch  Stagger Buckle ...    |\n"
            "|     The program 3DNA is used to get all the parameters.     |\n"
            "|     Example:    rnaview  -f  pdbfile -a                     |\n"
            "|                                                             |\n"
            "| For further information please contact:                     |\n"
            "| ndbadmin@ndbserver.rutgers.edu                              |\n"
            "+-------------------------------------------------------------+\n");
    exit (1);
}

double  get_torsion(long o3p_i0, long p, long o5p, long c5p, double **xyz)
{
    long  n;
    double  **xyz4;
    
    xyz4 = dmatrix(1, 4, 1, 3);
    
    for (n = 1; n <= 3; n++){
        xyz4[1][n] = xyz[o3p_i0][n];
        xyz4[2][n] = xyz[p][n];
        xyz4[3][n] = xyz[o5p][n];
        xyz4[4][n] = xyz[c5p][n];
    }
    return torsion(xyz4);
}



void rna_dna_bb_angle(char *pdbfile, long nchain, long **chain_idx, long **seidx, long *RY,
                       char **AtomName,char *ChainID,long *ResSeq, char **ResName, double **xyz)
{
    char c2c4[5],  n1n9[5], c6c8[5], *idmsg="", torout[512];
    long i,j, k,m, n=0, ib, ie, nb, ne, dna_rna=0;
    long o3p_i0=0, p=0, o5p=0, c5p=0, c4p=0, c3p=0, o3p=0, p_i1=0;
    long o5p_i1=0, n9=0, c2=0, c6=0, o4p=0, c1p=0;
    double alpha, beta, gamma, delta,epsilon, zeta , chi, ini=-9999.99;
    FILE *fw=NULL;
    

/*     
Note: alpha:   O3'(i-1)-P-O5'-C5'
      beta:    P-O5'-C5'-C4'
      gamma:   O5'-C5'-C4'-C3'
      delta:   C5'-C4'-C3'-O3'
      epsilon: C4'-C3'-O3'-P(i+1)
      zeta:    C3'-O3'-P(i+1)-O5'(i+1)

      chi for pyrimidines(Y): O4'-C1'-N1-C2
          chi for purines(R): O4'-C1'-N9-C4

*/
    sprintf(torout, "%s_torsion.out", pdbfile);
    fw=fopen(torout,"w");
    fprintf(fw, "The backbone torsional angle for nucleic acid \n\n");
    
    fprintf(fw, " \n\
Note: alpha:   O3'(i-1)-P-O5'-C5'\n\
      beta:    P-O5'-C5'-C4'\n\
      gamma:   O5'-C5'-C4'-C3'\n\
      delta:   C5'-C4'-C3'-O3'\n\
      epsilon: C4'-C3'-O3'-P(i+1)\n\
      zeta:    C3'-O3'-P(i+1)-O5'(i+1)\n\
\n\
      chi for pyrimidines(Y): O4'-C1'-N1-C2\n\
          chi for purines(R): O4'-C1'-N9-C4\n\
\n\
    (The number -9999.99 means that the angle can not be calculated) \n\
\n");
    fprintf(fw, "Ch Res  Num    alpha    beta     gamma    delta  epsilon    zeta     chi\n");
    
    for (i=1; i<=nchain; i++){ /* rid of ligand */
        if((chain_idx[i][2] - chain_idx[i][1]) <= 0)continue;
    
        for (k=chain_idx[i][1]; k<=chain_idx[i][2]; k++){ /*residue #*/
            alpha= ini;
            beta= ini;
            gamma= ini;
            delta= ini;
            epsilon= ini;
            zeta= ini;
            
            ib = seidx[k][1];
            ie = seidx[k][2];
            dna_rna = residue_ident(AtomName, xyz, ib, ie);
            if (dna_rna <0) continue; /*only DNA/RNA */

            
            
            p = find_1st_atom(" P  ", AtomName, ib, ie, idmsg);
            o5p = find_1st_atom(" O5'", AtomName, ib, ie, idmsg);
            c5p = find_1st_atom(" C5'", AtomName, ib, ie, idmsg);
            c4p = find_1st_atom(" C4'", AtomName, ib, ie, idmsg);
            c3p = find_1st_atom(" C3'", AtomName, ib, ie, idmsg);
            o3p = find_1st_atom(" O3'", AtomName, ib, ie, idmsg);

        /* chi(R): O4'-C1'-N9-C4; chi(Y): O4'-C1'-N1-C2 */
            if (RY[k] == 1) {
                strcpy(n1n9, " N9 ");
                strcpy(c2c4, " C4 ");
                strcpy(c6c8, " C8 ");
            } else if (RY[k] == 0) {
                strcpy(n1n9, " N1 ");
                strcpy(c2c4, " C2 ");
                strcpy(c6c8, " C6 ");
            }
            o4p = find_1st_atom(" O4'", AtomName, ib, ie, idmsg);
            c1p = find_1st_atom(" C1'", AtomName, ib, ie, idmsg);
            n9 = find_1st_atom(n1n9, AtomName, ib, ie, "");
            c2 = find_1st_atom(c2c4, AtomName, ib, ie, "");
            c6 = find_1st_atom(c6c8, AtomName, ib, ie, "");
            


            if(k>1){
                nb = seidx[k-1][1];
                ne = seidx[k-1][2];
                o3p_i0 =  find_1st_atom(" O3'", AtomName, nb, ne, idmsg);
            }
            if(k<chain_idx[i][2]){
                nb = seidx[k+1][1];
                ne = seidx[k+1][2];
                p_i1 = find_1st_atom(" P  ", AtomName, nb, ne, idmsg);
                o5p_i1 = find_1st_atom(" O5'", AtomName, nb, ne, idmsg);
            }
            
            if(o3p_i0 >0 &&  p>0 && o5p>0  && c5p>0 )alpha = get_torsion(o3p_i0, p, o5p, c5p,xyz); 
            if(p>0 && o5p>0 && c5p>0 && c4p>0 ) beta = get_torsion(p, o5p, c5p,c4p, xyz); 
            if(o5p>0 && c5p>0 && c4p>0 && c3p>0) gamma = get_torsion(o5p, c5p,c4p,c3p, xyz); 
            if(c5p>0 && c4p>0 && c3p>0 && o3p>0) delta = get_torsion(c5p,c4p,c3p,o3p, xyz); 
            if(c4p>0 && c3p>0 && o3p>0 && p_i1>0) epsilon = get_torsion(c4p,c3p,o3p,p_i1, xyz); 
            if(c3p>0 && o3p>0 && p_i1>0 && o5p_i1>0) zeta = get_torsion(c3p,o3p,p_i1,o5p_i1, xyz); 

            if(o4p>0 && c1p>0 && n9>0 && c2>0) chi= get_torsion(o4p,c1p,n9,c2, xyz);
            
            fprintf(fw,"%c %3s %4ld  %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n",
                   ChainID[ib],ResName[ib],ResSeq[ib],
                   alpha, beta, gamma , delta, epsilon, zeta, chi);
        
/*           
              
            printf("i,k, %ld, %4ld, %4ld, %4ld, %4ld, %4ld, %4ld, %4ld, %4ld, \
 %4ld, %4ld, %4ld, %4ld, %4ld, %4ld, %4ld, %ld\n",
                   i, k, o3p_i0, p,o5p,c5p,c4p,c3p,o3p,
                   p_i1,o5p_i1, o4p, c1p, n9,c2,c6,RY[k] );
       
            for(n = ib; n <= ie; n++){
                printf("%5ld %4s %3s %c %4ld  %8.3lf%8.3lf%8.3lf\n", n, AtomName[n], 
                       ResName[n], ChainID[n], ResSeq[n],  xyz[n][1],  xyz[n][2], xyz[n][3]);
            }

*/
            
        }
    }
    fclose(fw);
    printf("\nThe backbone torsion angles are in %s\n\n", torout);
    
}



    
void rna(char *pdbfile, long *type_stat, long **pair_stat, long *bs_all)
/* do all sorts of calculations */
{
    char outfile[BUF512], str[BUF512];
    char HB_ATOM[BUF512], ALT_LIST[BUF512];
    char *ChainID, *bseq, **AtomName, **ResName, **Miscs;
    long i, j, k,m,n, ie, ib, dna_rna, num, num_residue, nres, bs_atoms;
    long *ResSeq, *RY, **seidx, num_modify, *modify_idx, nprot_atom=0;
    long **chain_idx,nchain;
    double HB_UPPER[2], **xyz;
    static long base_all;    
    FILE *fout;
    
/*    sprintf(outfile, "%s.out", pdbfile);*/
    sprintf(outfile, "%s.out", FILEOUT);
    fout=fopen(outfile, "w");
    
/* read in H-bond length upper limit etc from <misc_rna.par> */    
    hb_crt_alt(HB_UPPER, HB_ATOM, ALT_LIST);

/* read in the PDB file */
    num = number_of_atoms(pdbfile);
    AtomName = cmatrix(1, num, 0, 4);
    ResName = cmatrix(1, num, 0, 3);
    ChainID = cvector(1, num);
    ResSeq = lvector(1, num);
    xyz = dmatrix(1, num, 1, 3);
    Miscs = cmatrix(1, num, 0, NMISC);
    
    printf("\nPDB data file name: %s\n",  pdbfile);
    fprintf(fout,"PDB data file name: %s\n",  pdbfile);
    num = read_pdb(pdbfile,AtomName, ResName, ChainID, ResSeq, xyz, Miscs,
                   ALT_LIST);

/* get the numbering information of each residue.
   seidx[i][j]; i = 1-num_residue  j=1,2
*/
    seidx = residue_idx(num, ResSeq, Miscs, ChainID, ResName, &num_residue);
        
/* Below is only for nucleic acids ie RY >= 0*/  
    bs_atoms = 0;
    for(i = 1; i <= num_residue; i++){
        ib = seidx[i][1];
        ie = seidx[i][2];
        dna_rna = residue_ident(AtomName, xyz, ib, ie);
       
        if (dna_rna >= 0){
            for(j = seidx[i][1]; j <= seidx[i][2]; j++){
                bs_atoms++;
                strcpy(AtomName[bs_atoms], AtomName[j]);
                strcpy(ResName[bs_atoms], ResName[j]);
                ChainID[bs_atoms] = ChainID[j];
                ResSeq[bs_atoms] = ResSeq[j];
/*                printf("%c %4s %4ld\n", ChainID[j], ResName[bs_atoms], ResSeq[bs_atoms]);*/
                

                       
                for(k = 0 ; k <=NMISC; k++)
                    Miscs[bs_atoms][k] = Miscs[j][k];
                for(k = 1 ; k <=3; k++)
                    xyz[bs_atoms][k] = xyz[j][k];
                
            }
        }
    }
    
/* get the new  numbering information of each residue */
/* get base sequence, RY identification */
/* identifying a residue as follows:  RY[j]
 *  R-base  Y-base  amino-acid, others [default] 
 *   +1        0        -1        -2 [default]

 bseq[j] --> the sigle letter name for each residue. j = 1 - num_residue.
 */

    bseq = cvector(1, num_residue);
    nres=0;    
    seidx=residue_idx(bs_atoms, ResSeq, Miscs, ChainID, ResName, &nres);

    chain_idx = lmatrix(1,500 , 1, 2);  /* # of chains max = 200 */    
    get_chain_idx(nres, seidx, ChainID, &nchain, chain_idx);

/*    
    for (i=1; i<=nchain; i++){ 
        for (k=chain_idx[i][1]; k<=chain_idx[i][2]; k++){
            printf("!nchain, chain_idx %4d %4d %4d %4d\n",
                   i, k ,chain_idx[i][1], chain_idx[i][2] );
        }
    }
*/

    
    bs_atoms = 0;
    for (i=1; i<=nchain; i++){ /* rid of ligand */
        if((chain_idx[i][2] - chain_idx[i][1]) <= 0)continue;
        ib=chain_idx[i][1];
        ie=chain_idx[i][2];
        
        printf("RNA/DNA chain_ID:  %c  from residue %4d to %4d\n",
               ChainID[ seidx[ib][1]], ResSeq[ seidx[ib][1] ], ResSeq[ seidx[ie][1] ]);
        
        for (k=chain_idx[i][1]; k<=chain_idx[i][2]; k++){
            ib = seidx[k][1];
            ie = seidx[k][2];
            dna_rna = residue_ident(AtomName, xyz, ib, ie);
       
            if (dna_rna >= 0){          
                for(j = ib; j <= ie; j++){
                    bs_atoms++;
                    strcpy(AtomName[bs_atoms], AtomName[j]);
                    strcpy(ResName[bs_atoms], ResName[j]);
                    ChainID[bs_atoms] = ChainID[j];
                    ResSeq[bs_atoms] = ResSeq[j];
                    for(m = 0 ; m <=NMISC; m++)
                        Miscs[bs_atoms][m] = Miscs[j][m];
                    for(m = 1 ; m <=3; m++)
                        xyz[bs_atoms][m] = xyz[j][m];
                    n=bs_atoms;
                       
/*                    
         printf("%s%5ld %4s%c%3s %c%4ld%c   %8.3lf%8.3lf%8.3lf\n", 
                 "ATOM  ", n, AtomName[n], Miscs[n][1],
                ResName[n], ChainID[n], ResSeq[n], Miscs[n][2], xyz[n][1],
                xyz[n][2], xyz[n][3]);
                        
*/                    
                }
            }
            
        }
        
    }
    
    nres=0;    
    seidx=residue_idx(bs_atoms, ResSeq, Miscs, ChainID, ResName, &nres);

    RY = lvector(1, num_residue);
    modify_idx = lvector(1, num_residue);    
    get_seq(fout,nres, seidx, AtomName, ResName, ChainID, ResSeq, Miscs,
            xyz, bseq, RY, &num_modify,modify_idx);  /* get the new RY */
 
    rna_dna_bb_angle(pdbfile,nchain, chain_idx, seidx, RY, AtomName, ChainID,  ResSeq, ResName,xyz);

    
    work_horse(pdbfile, fout, nres, bs_atoms, bseq, seidx, RY, AtomName,
               ResName, ChainID, ResSeq, Miscs, xyz,num_modify, modify_idx,  
               type_stat, pair_stat);

    base_all=base_all+nres; /* acculate all the bases */
    *bs_all=base_all;
    
//    if(!(PS>0 || VRML>0 || XML>0 || ALL>0) )
    write_tmp_pdb(pdbfile,nres,seidx,AtomName,ResName,ChainID,ResSeq,xyz);
 
    
    free_cmatrix(AtomName, 1, num, 0, 4);
    free_cmatrix(ResName, 1, num, 0, 3);
    free_cvector(ChainID, 1, num);
    free_lvector(ResSeq, 1, num);
    free_dmatrix(xyz, 1, num, 1, 3);
    free_cmatrix(Miscs, 1, num, 0, NMISC);
    free_lmatrix(seidx, 1, num_residue, 1, 2);
    free_cvector(bseq, 1, num_residue);
    free_lmatrix(chain_idx, 1,500 , 1, 2);   
    free_lvector(RY, 1, num_residue);
    free_lvector(modify_idx, 1, num_residue);
    
}

void work_horse(char *pdbfile, FILE *fout, long num_residue, long num,
                char *bseq, long **seidx, long *RY, char **AtomName,
                char **ResName, char *ChainID, long *ResSeq,char **Miscs, 
                double **xyz,long num_modify, long *modify_idx, 
                long *type_stat,long **pair_stat)
/* perform all the calculations */

{    
    
    long **bs_pairs_tot, num_pair_tot=0, num_single_base=0,*single_base, ntot;
    long i, j, num_loop, **loop;
    long num_bp = 0, num_helix = 1, nout = 16, nout_p1 = 17;
    long pair_istat[17], pair_jstat[17];
    long *bp_idx, *helix_marker, *matched_idx;
    long **base_pairs, **helix_idx;
    long **bp_order, **end_list;
    long num_multi, *multi_idx, **multi_pair;
    long num_bp_best=0, **pair_num_best,*sugar_syn;
    long xml_nh, *xml_helix_len, **xml_helix, xml_ns, *xml_bases;
    
    double BPRS[7];
    double **orien, **org, **Nxyz, **o3_p, **bp_xyz, **base_xy;
    char **pair_type;

    multi_pair = lmatrix(1, num_residue, 1, 20); /* max 20-poles */
    multi_idx = lvector(1, num_residue);  /*max multipoles = num_residue */
    pair_type = cmatrix(1, num_residue*2, 0, 3); /* max base pairs */    
    bs_pairs_tot = lmatrix(1, 2*num_residue, 1, 2);
    single_base = lvector(1, num_residue);  /* max single base */   
    sugar_syn = lvector(1, num_residue); 
 
    orien = dmatrix(1, num_residue, 1, 9);
    org = dmatrix(1, num_residue, 1, 3);
    Nxyz = dmatrix(1, num_residue, 1, 3);     /* RN9/YN1 atomic coordinates */
    o3_p = dmatrix(1, num_residue, 1, 8);     /* O3'/P atomic coordinates */

/* get the  base information for locating possible pairs later
  orien --> rotation matrix for each library base to match each residue.
  org --> the fitted sxyz origin for each library base.
*/
    base_info(num_residue, bseq, seidx, RY, AtomName, ResName, ChainID,
              ResSeq, Miscs, xyz,  orien, org, Nxyz, o3_p, BPRS);    

/* find all the base-pairs */
    printf("Finding all the base pairs...\n");    
    all_pairs(pdbfile, fout, num_residue, RY, Nxyz, orien, org, BPRS,
              seidx, xyz, AtomName, ResName, ChainID, ResSeq, Miscs, bseq,
              &num_pair_tot, pair_type, bs_pairs_tot, &num_single_base,
              single_base, &num_multi, multi_idx, multi_pair, sugar_syn);
/*
	for(i=1; i<=num_pair_tot; i++){
		printf("pair-type %4d %4d %4d %s\n", i, bs_pairs_tot[i][1], bs_pairs_tot[i][2],
		pair_type[i]);
	}
*/


    
    fprintf(fout, "  The total base pairs =%4d (from %4d bases)\n",
           num_pair_tot,num_residue);
    printf("  The total base pairs =%4d (from %4d bases);\n",
           num_pair_tot,num_residue);
    
    if (!num_pair_tot) {        
        printf( "No base-pairs found for (%s) "
                "(May be a single strand!!)\n\n",pdbfile);
        fclose(fout);
        if(PS>0) printf( "No 2D structure plotted!\n");
        return; 
    }
    pair_type_statistics(fout, num_pair_tot, pair_type, type_stat); /*12 type*/
    sixteen_pair_statistics(num_pair_tot, bs_pairs_tot, bseq, pair_type, pair_stat);
     
        
    fprintf(fout, "------------------------------------------------\n");    
    fclose(fout);
  /* do statistics for the single H bonded pairs (base to base)  
	write_single_Hbond_stat(pdbfile, bseq, pair_stat); */

    ntot = 2*num_pair_tot;
/*    if(ARGC <=2){ */
    motif(pdbfile);  /*  write the RNA motif pattern */   
    	print_sorted_pair(ntot, pdbfile);/* sort pairs according to W-H-S*/
/*    	write_multiplets(pdbfile);  */
            /*
        return;
    }
            */   

/* find best base-pairs */    
    matched_idx = lvector(1, num_residue);
    base_pairs = lmatrix(1, num_residue, 1, nout_p1);

/* base_pairs[][17]: i, j, bpid, d, dv, angle, dNN, dsum, bp-org, normal1,normal2
 *             col#  1  2    3   4  5     6     7     8    9-11    12-14   15-17
 * i.e., add one more column for i from pair_stat [best_pair]
*/
    for (i = 1; i <= num_residue; i++) {
        best_pair(i, num_residue, RY, seidx, xyz, Nxyz, matched_idx, orien,
                  org, AtomName, bseq, BPRS, pair_istat);
        if (pair_istat[1]) {        /* with paired base */
            best_pair(pair_istat[1], num_residue, RY, seidx, xyz, Nxyz,
                      matched_idx, orien, org, AtomName, bseq, BPRS,
                      pair_jstat);
            if (i == pair_jstat[1]) { /* best match between i && pair_istat[1] */
                matched_idx[i] = 1;
                matched_idx[pair_istat[1]] = 1;
                base_pairs[++num_bp][1] = i;
                for (j = 1; j <= nout; j++)
                    base_pairs[num_bp][j + 1] = pair_istat[j];
            }
        }
    }
    
    bp_idx = lvector(1, num_bp);
    helix_marker = lvector(1, num_bp);
    helix_idx = lmatrix(1, num_bp, 1, 7);
   
    re_ordering(num_bp, base_pairs, bp_idx, helix_marker, helix_idx, BPRS,
                &num_helix, o3_p, bseq, seidx, ResName, ChainID, ResSeq,
                Miscs);

    bp_order = lmatrix(1, num_bp, 1, 3);
    end_list = lmatrix(1, num_bp, 1, 3);
    bp_xyz = dmatrix(1, num_bp, 1, 9); /*   bp origin + base I/II normals: 9 - 17 */
    
    for (i = 1; i <= num_bp; i++)
        for (j = 1; j <= 9; j++)
            bp_xyz[i][j] = base_pairs[i][j + 8] / MFACTOR;
    
    pair_num_best = lmatrix(1, 3, 1, num_bp);
       
    write_best_pairs(num_helix, helix_idx, bp_idx, helix_marker,
                     base_pairs, seidx, ResName, ChainID, ResSeq,
                     Miscs, bseq, BPRS, &num_bp_best, pair_num_best);    

  
    if(PS>0){
        if(num_pair_tot<2){
            printf("Too few base pairs to form a helix(NO 2D structure plotted)!\n");
            return;
        }
        if(num_helix==0){
            printf("No anti-parallel helix (NO 2D structure plotted)!\n");
            return;
        }

        base_xy = dmatrix(0, num_residue, 1, 2);   
        loop = lmatrix(1, num_helix*2, 1, 2);
        xml_helix= lmatrix(0, num_residue, 1, 2); 
        xml_helix_len= lvector(0, num_residue); 
        xml_bases= lvector(0, num_residue);

        process_2d_fig(num_residue, bseq, seidx, RY, AtomName, ResName,ChainID,
                       ResSeq, Miscs, xyz, num_pair_tot, pair_type, bs_pairs_tot,
                       num_helix, helix_idx, bp_idx, base_pairs, base_xy,
                       &num_loop, loop, &xml_nh, xml_helix_len, xml_helix,
                       &xml_ns, xml_bases);
        
        if(xml_nh==0){
            printf("No anti-parallel helix (NO 2D structure plotted, No XML file output)!\n");
            return;
        }
        
        write_xml(FILEOUT, num_residue, bseq, seidx, AtomName, ResName,
                  ChainID, ResSeq, Miscs, xyz, xml_nh, xml_helix,
                  xml_helix_len, xml_ns, xml_bases, num_pair_tot,
                  bs_pairs_tot,pair_type, base_xy, num_modify,
                  modify_idx,num_loop, loop, num_multi,multi_idx,
                  multi_pair,sugar_syn);

        free_dmatrix(base_xy , 0, num_residue, 1, 2);         
        free_lmatrix(loop , 1, num_helix*2, 1, 2);
        free_lmatrix(xml_helix, 0, num_residue, 1, 2); 
        free_lvector(xml_helix_len, 0, num_residue); 
        free_lvector(xml_bases, 0, num_residue); 
      
    	printf("Ploting 2D structure...\n");        
        xml2ps(pdbfile, num_residue, XML); /* read information from RNAML */

    }

    if(ANAL >0){
        bp_analyze(pdbfile, num, AtomName, ResName, ChainID, ResSeq, xyz, 
                   num_residue,Miscs, seidx, num_bp_best, pair_num_best);
    }
    
    if(VRML>0){
        process_3d_fig(pdbfile,num_residue, bseq, seidx, AtomName, ResName,
                       ChainID,xyz, num_pair_tot, pair_type, bs_pairs_tot);
    }
    
    
    free_lmatrix(multi_pair, 1, num_residue, 1, 20);  
    free_lvector(multi_idx , 1, num_residue);       
    free_cmatrix(pair_type , 1, num_residue*2, 0, 3);          
    free_lmatrix(bs_pairs_tot , 1, 2*num_residue, 1, 2);     
    free_lvector(single_base , 1, num_residue);         
    free_lvector(sugar_syn , 1, num_residue); 
    free_dmatrix(bp_xyz, 1, num_bp, 1, 9);
    free_lmatrix(pair_num_best, 1, 3, 1, num_bp);

    free_lmatrix(bp_order, 1, num_bp, 1, 3);
    free_lmatrix(end_list, 1, num_bp, 1, 3);

    free_lvector(bp_idx, 1, num_bp);
    free_lvector(helix_marker, 1, num_bp);
    free_lmatrix(helix_idx, 1, num_bp, 1, 7);

    free_lvector(matched_idx, 1, num_residue);
    free_lmatrix(base_pairs, 1, num_residue, 1, nout_p1);


    free_dmatrix(orien, 1, num_residue, 1, 9);
    free_dmatrix(org, 1, num_residue, 1, 3);
    free_dmatrix(Nxyz, 1, num_residue, 1, 3);
    free_dmatrix(o3_p, 1, num_residue, 1, 8);

}


void write_tmp_pdb(char *pdbfile,long nres, long **seidx, char **AtomName,
                   char **ResName, char *ChainID, long *ResSeq, double **xyz)
/* write a tmp pdb file for the web */
{
    char parfile[100];    
    long i, j, ib, ie;
    FILE *fp;


/*    sprintf(parfile, "%s_tmp.pdb", pdbfile);*/
    sprintf(parfile, "%s_tmp.pdb", FILEOUT);
    fp = fopen(parfile, "w");

    for (i = 1; i <= nres; i++) {
        ib=seidx[i][1];
        ie=seidx[i][2];
        for (j = ib; j <= ie; j++) {
            fprintf(fp, "%s%5d %4s %3s %c%4d    %8.3f%8.3f%8.3f\n", "ATOM  ",
                    i, AtomName[j], ResName[j], ChainID[j], ResSeq[j],xyz[j][1],
                    xyz[j][2], xyz[j][3]);
        }
    }
    fclose(fp);   
}
        
void print_sorted_pair(long ntot, char *pdbfile)
{
    long i, j, np=0, n, n1[25], **index;
    
    char str[200], **str_pair,**str_tmp ;
    char inpfile[80], outfile[80];
    
    FILE  *finp, *fout;
    str_pair = cmatrix(0, ntot, 0, 120);
    str_tmp = cmatrix(0, ntot, 0, 120);   
    index=lmatrix(0,25, 0, ntot); 
/*
    sprintf(inpfile, "%s.out", pdbfile);
    sprintf(outfile, "%s_sort.out", pdbfile);
*/
    sprintf(inpfile, "%s.out", FILEOUT);
    sprintf(outfile, "%s_sort.out", FILEOUT);
    
    fout = fopen(outfile, "w");
    finp = fopen(inpfile, "r");
    while(fgets(str, sizeof str, finp) !=NULL){        
        if (strstr(str, "BEGIN_base-pair")){
            np=0;
            if(np>=ntot){
                printf("Increase memory for str_pair(in print_sorted_pair)\n"); 
                return;
            }
            while(fgets(str, sizeof str, finp) !=NULL){
                if(strstr(str, "!") || strstr(str, "stack")) continue;                
                strcpy(str_pair[np], str);
                np++;					
                if (strstr(str, "END_base-pair"))               
                    break;
            }            
        }
    }
    fclose(finp);

    
    for (i=0; i<25; i++){
        n1[i]=0;
        for (j=0; j<ntot; j++)
            index[i][j]=0;
    }
  

    for (i=0; i<np; i++){
        if((strstr(str_pair[i], "+/+")||strstr(str_pair[i], "-/- ") )
           && strstr(str_pair[i], "cis ")){
            
            n1[0]++;
            n = n1[0];
            index[0][n] = i;
              
            
        } else if (strstr(str_pair[i], "W/W") && strstr(str_pair[i], "cis")){
            n1[1]++;
            n = n1[1];
            index[1][n] = i;
        } else if (strstr(str_pair[i], "W/W") && strstr(str_pair[i], "tran")){
            n1[2]++;
            n = n1[2];
            index[2][n] = i;
        } else if (strstr(str_pair[i], "W/H") && strstr(str_pair[i], "cis")){
            n1[3]++;
            n = n1[3];
            index[3][n] = i;
        } else if (strstr(str_pair[i], "W/H") && strstr(str_pair[i], "tran")){
            n1[4]++;
            n = n1[4];
            index[4][n] = i;
        } else if (strstr(str_pair[i], "W/S") && strstr(str_pair[i], "cis")){
            n1[5]++;
            n = n1[5];
            index[5][n] = i;
        } else if (strstr(str_pair[i], "W/S") && strstr(str_pair[i], "tran")){
            n1[6]++;
            n = n1[6];
            index[6][n] = i;
        } else if (strstr(str_pair[i], "H/W") && strstr(str_pair[i], "cis")){
            n1[7]++;
            n = n1[7];
            index[7][n] = i;
        } else if (strstr(str_pair[i], "H/W") && strstr(str_pair[i], "tran")){
            n1[8]++;
            n = n1[8];
            index[8][n] = i;
        } else if (strstr(str_pair[i], "H/H") && strstr(str_pair[i], "cis")){
            n1[9]++;
            n = n1[9];
            index[9][n] = i;
        } else if (strstr(str_pair[i], "H/H") && strstr(str_pair[i], "tran")){
            n1[10]++;
            n = n1[10];
            index[10][n] = i;
        } else if (strstr(str_pair[i], "H/S") && strstr(str_pair[i], "cis")){
            n1[11]++;
            n = n1[11];
            index[11][n] = i;
        } else if (strstr(str_pair[i], "H/S") && strstr(str_pair[i], "tran")){
            n1[12]++;
            n = n1[12];
            index[12][n] = i;
        } else if (strstr(str_pair[i], "S/W") && strstr(str_pair[i], "cis")){
            n1[13]++;
            n = n1[13];
            index[13][n] = i;
        } else if (strstr(str_pair[i], "S/W") && strstr(str_pair[i], "tran")){
            n1[14]++;
            n = n1[14];
            index[14][n] = i;
        } else if (strstr(str_pair[i], "S/H") && strstr(str_pair[i], "cis")){
            n1[15]++;
            n = n1[15];
            index[15][n] = i;
        } else if (strstr(str_pair[i], "S/H") && strstr(str_pair[i], "tran")){
            n1[16]++;
            n = n1[16];
            index[16][n] = i;
        } else if (strstr(str_pair[i], "S/S") && strstr(str_pair[i], "cis")){
            n1[17]++;
            n = n1[17];
            index[17][n] = i;
        } else if (strstr(str_pair[i], "S/S") && strstr(str_pair[i], "tran")){
            n1[18]++;
            n = n1[18];
            index[18][n] = i;
            
        } else if ((strstr(str_pair[i], "W/.")||strstr(str_pair[i], "./W"))&&
                   strstr(str_pair[i], "cis")|| strstr(str_pair[i], "tran") ){
            n1[19]++;
            n = n1[19];
            index[19][n] = i;

        } else if ((strstr(str_pair[i], "H/.")||strstr(str_pair[i], "./H"))&&
                   strstr(str_pair[i], "cis")|| strstr(str_pair[i], "tran") ){
            n1[20]++;
            n = n1[20];
            index[20][n] = i;

        } else if ((strstr(str_pair[i], "S/.")||strstr(str_pair[i], "./S"))&&
                   strstr(str_pair[i], "cis")|| strstr(str_pair[i], "tran") ){
            n1[21]++;
            n = n1[21];
            index[21][n] = i;
        } else if ((strstr(str_pair[i], "./.")||strstr(str_pair[i], "./."))&&
                   strstr(str_pair[i], "cis")|| strstr(str_pair[i], "tran") ){
            n1[22]++;
            n = n1[22];
            index[22][n] = i;
        } else if(strstr(str_pair[i], "__ ")){
            n1[23]++;
            n = n1[23];
            index[23][n] = i;
        }

    }
        
    finp = fopen(inpfile, "r");

    while(fgets(str, sizeof str, finp) !=NULL){
        fprintf(fout, "%s",str);
        if (strstr(str, "BEGIN_base-pair")){
            for(i=0; i<=23; i++){
                for(j=1; j<=n1[i]; j++){   
                    n = index[i][j];
                        
                    strcpy(str_tmp[j], str_pair[n]); 
                }
                n = n1[i];
                sort_by_pair(fout, n, str_tmp);
            }
            while(fgets(str, sizeof str, finp) !=NULL){
                if (strstr(str, "END_base-pair")){
                    fprintf(fout, "%s",str);
                    break;
                }
            }
        }
    }
    
    free_cmatrix(str_pair , 0, ntot, 0, 120);
    free_cmatrix(str_tmp , 0, ntot, 0, 120);
    free_lmatrix(index,0,25, 0, ntot); 
       
}

void sort_by_pair(FILE *fout, long nt, char **str)
/* sort the pairs by A>C>G>U */
{
    
    long i, j, n, n1[25], **index;
 
    index=lmatrix(0,25,0,nt);
    for (i=0; i<25; i++){
        n1[i]=0;
                
        for (j=0; j<=nt; j++)
            index[i][j]=0;
            
    }
    
    
    for (i=1; i<=nt; i++){
        if        (strstr(str[i], "A-A")){
            n1[1]++;
            n = n1[1];
            index[1][n] = i;
        } else if (strstr(str[i], "A-C") ){
            n1[2]++;
            n = n1[2];
            index[2][n] = i;
        } else if (strstr(str[i], "A-G") ){
            n1[3]++;
            n = n1[3];
            index[3][n] = i;
        } else if (strstr(str[i], "A-U") ){
            n1[4]++;
            n = n1[4];
            index[4][n] = i;
        } else if (strstr(str[i], "C-A") ){
            n1[5]++;
            n = n1[5];
            index[5][n] = i;
        } else if (strstr(str[i], "C-C") ){
            n1[6]++;
            n = n1[6];
            index[6][n] = i;
        } else if (strstr(str[i], "C-G") ){
            n1[7]++;
            n = n1[7];
            index[7][n] = i;
        } else if (strstr(str[i], "C-U") ){
            n1[8]++;
            n = n1[8];
            index[8][n] = i;
        } else if (strstr(str[i], "G-A") ){
            n1[9]++;
            n = n1[9];
            index[9][n] = i;
        } else if (strstr(str[i], "G-C ") ){
            n1[10]++;
            n = n1[10];
            index[10][n] = i;
        } else if (strstr(str[i], "G-G") ){
            n1[11]++;
            n = n1[11];
            index[11][n] = i;
        } else if (strstr(str[i], "G-U") ){
            n1[12]++;
            n = n1[12];
            index[12][n] = i;
        } else if (strstr(str[i], "U-A") ){
            n1[13]++;
            n = n1[13];
            index[13][n] = i;
        } else if (strstr(str[i], "U-C") ){
            n1[14]++;
            n = n1[14];
            index[14][n] = i;
        } else if (strstr(str[i], "U-G") ){
            n1[15]++;
            n = n1[15];
            index[15][n] = i;
        } else if (strstr(str[i], "U-U") ){
            n1[16]++;
            n = n1[16];
            index[16][n] = i;
        } else if (strstr(str[i], "__") ){
            n1[17]++;
            n = n1[17];
            index[17][n] = i;
        }
    }
        
    for(i=1; i<=17; i++){
        for(j=1; j<=n1[i]; j++){
            n = index[i][j];
            fprintf(fout, "%s",str[n]);
                
        }
    }
    
}







        
            
        
                
