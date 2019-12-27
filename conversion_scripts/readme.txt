##############################################################################
##############################################################################

If running from crystal structure:

First, download pdb file.
Run mapping_atocg.py
	use .pdb file

Run tinker's pdbxyz.x
	use .pdb.pdb file

Run xyz2cg_xyz.py
	use .pdb.xyz file

output file .xyz.xyz is ready to run in TINKER. (can rename to just .xyz file)


#############################################################################
#############################################################################

If running from sequence:

Run tinker's nucleic.x program using the amber parameters (amber99.prm found
in params file). Note: can only create single stranded or duplexes; if more
than two strands or multihelical regions, will need to connect manually.
	use .xyz file

Run tinker's xyzpdb.x (amber parameters again)
	use .pdb file

Then, same flow as crystal structure:
Run mapping_atocg.py
        use .pdb file

Run tinker's pdbxyz.x
        use .pdb.pdb file

Run xyz2cg_xyz.py
        use .pdb.xyz file

output file .xyz.xyz is ready to run in TINKER. (can rename to just .xyz file)


#############################################################################
#############################################################################

Converting back to all-atom structure from coarse-grained structure:

./cg_xyz2aa_pdb.py cg_file.xyz            # where cg_file.xyz is the name of the coarse-grained structure

use the amber99.prm parameters in Tinker. (key file = Parameters amber99.prm)

Run tinker's pdbxyz.x 
	use cg_file.xyz.aa.pdb

Run tinker's minimize.x
	use cg_file.xyz.aa.xyz
	use gradient cutoff of 1.0

output cg_file.xyz.aa.xyz_2 is the minimized all atom xyz structure

If want to convert to pdb,

Run tinker's xyzpdb.x
	use cg_file.xyz.aa.xyz_2

	output cg_file.xyz.aa.pdb_2 is the minimized all atom pdb structure.




