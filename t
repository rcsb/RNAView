src/analyze.c:        backbone_torsion(ds, num_bp, bp_seq, backbone, sugar, chi, xyz,
src/analyze.c:void backbone_torsion(long ds, long num_bp, char **bp_seq, long **backbone,
src/analyze.c:                    alpha2zeta[i][ioffset + k] = torsion(xyz4);
src/analyze.c:            /* chi torsion angle */
src/analyze.c:                chi_angle[i][j] = torsion(xyz4);
src/analyze.c:            /* sugar ring torsion angles v0 to v4 */
src/analyze.c:                    sugar_angle[i][o7 + k] = torsion(xyz4);
src/analyze.c:    fprintf(fp, "Main chain and chi torsion angles: \n\n"
src/analyze.c:double torsion(double **d)
src/analyze.c:/* get torsion angle a-b-c-d in degrees */
src/fpair.c:            chi_angle = torsion(xyz4);
src/rnaview.c:double  get_torsion(long o3p_i0, long p, long o5p, long c5p, double **xyz)
src/rnaview.c:    return torsion(xyz4);
src/rnaview.c:  alpha = get_torsion(o3p_i0, p, o5p, c5p,xyz); 
src/rnaview.c:            if(o3p_i0 >0 &&  p>0 && o5p>0  && c5p>0 )alpha = get_torsion(o3p_i0, p, o5p, c5p,xyz); 
src/rnaview.c:            if(p>0 && o5p>0 && c5p>0 && c4p>0 ) beta = get_torsion(p, o5p, c5p,c4p, xyz); 
src/rnaview.c:            if(o5p>0 && c5p>0 && c4p>0 && c3p>0) gamma = get_torsion(o5p, c5p,c4p,c3p, xyz); 
src/rnaview.c:            if(c5p>0 && c4p>0 && c3p>0 && o3p>0) delta = get_torsion(c5p,c4p,c3p,o3p, xyz); 
src/rnaview.c:            if(c4p>0 && c3p>0 && o3p>0 && p_i1>0) epsilon = get_torsion(c4p,c3p,o3p,p_i1, xyz); 
src/rnaview.c:            if(c3p>0 && o3p>0 && p_i1>0 && o5p_i1>0) zeta = get_torsion(c3p,o3p,p_i1,o5p_i1, xyz); 
src/rnaview.c:            if(o4p>0 && c1p>0 && n9>0 && c2>0) chi= get_torsion(o4p,c1p,n9,c2, xyz);
