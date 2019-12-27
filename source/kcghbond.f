c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine kcghbond  --  Hbond parameter assignment      ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "kcghbond" assigns a force constant, ideal bond length,
c     and cutoff to each hydrogen bond type
c
c
      subroutine kcghbond
      use sizes
      use atomid
      use atoms
      use vdwpot
      use couple
      use fields
      use inform
      use iounit
      use kcghbonds
      use keys
      use potent
      use usage
      implicit none
      integer i,j,next
      integer type_ia,type_ib
      integer type_ic,type_ja
      integer type_jb,type_jc
      real*8  t_ehbmax,t_hbondmin
      real*8  t_cghbondcutoff
      logical header
      
      character*20 keyword
      character*120 record
      character*120 string
c
c
c     process keywords containing bond stretch parameters
c
c initialize index type array to zeros
      cghb_types = 0 
      ncghb = 0

      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:8) .eq. 'CGHBOND ') then
            ncghb = ncghb + 1
            type_ia = 0
            type_ib = 0
            type_ic = 0 
            type_ja = 0 
            type_jb = 0
            type_jc = 0
            t_ehbmax = 0.0d0
            t_hbondmin = 0.0d0
            t_cghbondcutoff = 0.0d0
            string = record(next:120)
            read (string,*,err=10,end=10)  type_ia,type_ib,type_ic,
     &                                     type_ja,type_jb,type_jc,
     &                                     t_ehbmax,t_hbondmin,
     &                                     t_cghbondcutoff
   10       continue
            if (.not. silent) then
               if (header) then
                  header = .false.
                  write (iout,20)
   20             format (/,' Hydrogen Bond  Parameters :',
     &                    //,5x,'AtomClasses',9x,'maxpot',6x,
     &                    'equil dist','cutoff',/)
               end if
               write (iout,30)  type_ia,type_ib,type_ic,type_ja,
     &                          type_jb,type_jc,t_ehbmax,t_hbondmin,
     &                          t_cghbondcutoff
   30          format (6x,6i4,4x,f12.3,f12.4,f12.4)
            end if
            if (cghb_types(1,ncghb).eq. 0 ) then
               cghb_types(1,ncghb) = type_ia
               cghb_types(2,ncghb) = type_ib
               cghb_types(3,ncghb) = type_ic
               cghb_types(4,ncghb) = type_ja
               cghb_types(5,ncghb) = type_jb
               cghb_types(6,ncghb) = type_jc
               cghb_ehbmax(ncghb) = t_ehbmax
               cghb_hbondmin(ncghb) = t_hbondmin
               cghb_cghbondcutoff(ncghb) = t_cghbondcutoff
               goto 60
            end if
   60       continue
         end if
      end do

      return
      end
