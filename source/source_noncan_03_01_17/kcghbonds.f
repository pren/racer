c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  module kcghbonds  --  Hydrogen bond forcefield parameters  ##
c     ##                                                            ##
c     ################################################################
c
c
c     maxncghb   maximum number of cg hydrogen bond parameter entries
c     cghb_types    atom type indices specifiying hydrogen bond
c     bcon    force constant parameters for harmonic bond stretch
c     blen    bond length parameters for harmonic bond stretch
c     kb      string of atom classes for harmonic bond stretch
c
c
      module kcghbonds
      implicit none
      integer maxncghb
      parameter (maxncghb=2000)
      integer, dimension(6,maxncghb) :: cghb_types
      real*8 cghb_ehbmax(maxncghb) 
      real*8 cghb_hbondmin(maxncghb) 
      real*8 cghb_cghbondcutoff(maxncghb) 
      save
      end
