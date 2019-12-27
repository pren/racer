07/25/17 Update of RACER conversion scripts
Author: Sara Cheng
Email: sara.y.cheng@utexas.edu

The previous version of xyz2cg_xyz.py, which outputs an xyz file with
connections between the pseduoatoms, the following order of pseudoatoms was
assumed.

For purines
1.phosphate
2.sugar
3.sugar connect
4.base

For pyrimidines
1. phosphate
2. sugar
3. base
4. sugar connect

For atomic PDBs with functional groups differing from the order of phsophate
coordinates, sugar coordinates, and base coordinates, the xyz2cg_xyz.py would
produce erroneous coarse-grained xyz files, with connections written
incorrectly.

The latest version of xyz2cg_xyz.py (xyz2cg_xyz_v2.py) is corrected. Now, the
newer version of the code can recognize purines with the base pseudoatoms
before the sugar connect.

Overall, the xyz2cg_xyz_v2.py can be used in place of the xyz2cg_xyz.py, with
no other changes to the other coversion scripts.

For additional questions, please contact Sara at sara.y.cheng@utexas.edu
