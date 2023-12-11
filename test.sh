rnaview -p --pdb test/pdb/pdb1nvy/pdb1nvy.pdb 1> /dev/null
rnaview -p --pdb test/pdb/test1/test1.pdb 1> /dev/null
rnaview -p --pdb test/pdb/tr0001/tr0001.pdb 1> /dev/null
rnaview -p --pdb test/pdb/url064/url064.pdb 1> /dev/null
rnaview -p --pdb test/pdb/urx053/urx053.pdb 1> /dev/null
rnaview -p --cif test/mmcif/insertion_code/1EFW/1EFW.cif 1> /dev/null
rnaview -p --cif test/mmcif/insertion_code/1VVJ/1VVJ.cif 1> /dev/null
rnaview -p --cif test/mmcif/insertion_code/4ARC/4ARC.cif 1> /dev/null
rnaview -p --cif test/mmcif/nmr_structure/8if5/8if5.cif 1> /dev/null
rnaview -p --cif test/mmcif/other/6pom/6pom.cif 1> /dev/null
rnaview -p --cif test/mmcif/x-ray/3P4J/assembly-1/3p4j-assembly1.cif 1> /dev/null
rnaview -p --cif test/mmcif/x-ray/434D/assembly-1/434d-assembly1.cif 1> /dev/null
rnaview -p --cif test/mmcif/x-ray/434D/assembly-2/434d-assembly2.cif 1> /dev/null
rnaview -p --cif test/mmcif/x-ray/4NMG/assembly-1/4nmg-assembly1.cif 1> /dev/null

diff test/pdb/pdb1nvy/pdb1nvy.pdb.ps test/pdb/pdb1nvy/pdb1nvy_gt.pdb.ps
diff test/pdb/test1/test1.pdb.ps test/pdb/test1/test1_gt.pdb.ps
diff test/pdb/tr0001/tr0001.pdb.ps test/pdb/tr0001/tr0001_gt.pdb.ps
diff test/pdb/url064/url064.pdb.ps test/pdb/url064/url064_gt.pdb.ps
diff test/pdb/urx053/urx053.pdb.ps test/pdb/urx053/urx053_gt.pdb.ps
diff test/mmcif/insertion_code/1EFW/1EFW.cif.ps test/mmcif/insertion_code/1EFW/1EFW_gt.cif.ps
diff test/mmcif/insertion_code/1VVJ/1VVJ.cif.ps test/mmcif/insertion_code/1VVJ/1VVJ_gt.cif.ps
diff test/mmcif/insertion_code/4ARC/4ARC.cif.ps test/mmcif/insertion_code/4ARC/4ARC_gt.cif.ps
diff test/mmcif/nmr_structure/8if5/8if5.cif.ps test/mmcif/nmr_structure/8if5/8if5_gt.cif.ps
diff test/mmcif/other/6pom/6pom.cif.ps test/mmcif/other/6pom/6pom_gt.cif.ps
diff test/mmcif/x-ray/3P4J/assembly-1/3p4j-assembly1.cif.ps test/mmcif/x-ray/3P4J/assembly-1/3p4j-assembly1_gt.cif.ps
diff test/mmcif/x-ray/434D/assembly-1/434d-assembly1.cif.ps test/mmcif/x-ray/434D/assembly-1/434d-assembly1_gt.cif.ps
diff test/mmcif/x-ray/434D/assembly-2/434d-assembly2.cif.ps test/mmcif/x-ray/434D/assembly-2/434d-assembly2_gt.cif.ps
diff test/mmcif/x-ray/4NMG/assembly-1/4nmg-assembly1.cif.ps test/mmcif/x-ray/4NMG/assembly-1/4nmg-assembly1_gt.cif.ps
