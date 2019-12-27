c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine ecghbond1  --  CG H-bond potential energy    ##
c     ##                                                         ##
c     #############################################################
c
c
c     "ecghbond1" calculates the CG H-bond potential energy and
c     firsts derivatives
c
      subroutine ecghbond1
      use sizes
      use angbnd
      use angpot
      use atoms
      use bound
      use deriv
      use energi
      use group
      use math
      use usage
      use vdwpot
      use virial
      use kcghbonds
      implicit none
      integer i,ia,ib,ic,id
      integer j,ja,jb,jc,jd
      integer k
      integer NTHREADS,OMP_GET_NUM_THREADS
      real*8 rhbond,rhbond2
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xip,yip,zip
      real*8 xja,yja,zja
      real*8 xjb,yjb,zjb
      real*8 xjc,yjc,zjc
      real*8 xjp,yjp,zjp    
      real*8 rip,rjp            
      real*8 e,ecghbo
      real*8 doti,cosinei
      real*8 dotj,cosinej
      real*8 dangledcosinei,dangledcosinej
      real*8 danglekdangle
      real*8 anglei,anglej,anglek,fgrp
      real*8 rx,ry,rz
      real*8 xiab,yiab,ziab
      real*8 xicb,yicb,zicb
      real*8 xjab,yjab,zjab
      real*8 xjcb,yjcb,zjcb
      real*8 cghbondcutoff2
      real*8 r03
      real*8 ehbmax, hbondmin
      real*8 cghbondcutoff
      logical proceed
      real *8 dedr,dedt
      real *8 dedtxia,dedtyia,dedtzia
      real *8 dedtxib,dedtyib,dedtzib
      real *8 dedtxic,dedtyic,dedtzic
      real *8 dedtxja,dedtyja,dedtzja
      real *8 dedtxjb,dedtyjb,dedtzjb
      real *8 dedtxjc,dedtyjc,dedtzjc
      real *8 dedrxib,dedryib,dedrzib
      real *8 dedrxjb,dedryjb,dedrzjb
      real *8 dedxia,dedyia,dedzia 
      real *8 dedxib,dedyib,dedzib
      real *8 dedxic,dedyic,dedzic
      real *8 dedxja,dedyja,dedzja
      real *8 dedxjb,dedyjb,dedzjb
      real *8 dedxjc,dedyjc,dedzjc
      real*8, allocatable :: decghbo(:,:)

c
c     ## Constants   ehbmax,cghbondcutoff,hbondmin prmkey.f and vdwpot
 
c      r03 = hbondmin*hbondmin*hbondmin

c      cghbondcutoff2 = cghbondcutoff*cghbondcutoff
c      write (*,*) "Use CG H-bond, magnitude = ",  ehbmax, hbondmin,
c     &                cghbondcutoff
c      write (*,*) "Use CG H-bond"
c
c
c     zero out the angle bending energy component
c
      ecghb = 0.0d0
      do i = 1, n
         decghb(1,i) = 0.0d0
         decghb(2,i) = 0.0d0
         decghb(3,i) = 0.0d0
      end do
c
c
c     decide whether to compute the current interaction
c
         proceed = .true.
c         if (use_group)  call groups (proceed,fgrp,ia,ib,ic,0,0,0)
c         if (proceed)  proceed = (use(ia) .or. use(ib) .or. use(ic))         
c      

c
c     initialize local variables for OpenMP calculation
c
      ecghbo = ecghb
      allocate (decghbo(3,n))
      do i = 1, n
         decghbo(1,i) = 0.0d0
         decghbo(2,i) = 0.0d0
         decghbo(3,i) = 0.0d0
      end do

c
c     set OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(nangle,iang,x,y,z,
!$OMP& type,cghb_ehbmax,cghb_cghbondcutoff,cghb_types,
!$OMP& cghb_hbondmin,ncghb)
!$OMP& shared(ecghbo,decghbo)
!$OMP DO reduction(+:ecghbo,decghbo) schedule(guided)

c
c     calculate the bond angle bending energy term
c
      do i = 1, nangle
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
c     determine atom type is N6 N2 O6 O2 4 5 6 7  CG 3
         if ((type(ib) .GT. 2) .AND. (type(ib) .LT. 8)) then
         xib = x(ib)
         yib = y(ib)
         zib = z(ib)
c         write (*,*) "ia ib ic",ia,ib,ic
        do j = i+1, nangle
            jb = iang(2,j)
            if ((type(jb) .GT. 2) .AND. (type(jb) .LT. 8)) then
            xjb = x(jb)
            yjb = y(jb)
            zjb = z(jb)    
            rx = xib - xjb
            ry = yib - yjb
            rz = zib - zjb               
            rhbond2 = rx*rx + ry*ry + rz*rz
            if (rhbond2 .gt. 60) goto 30
            rhbond = sqrt(rhbond2)
            xia = x(ia)
            yia = y(ia)
            zia = z(ia)
            xic = x(ic)
            yic = y(ic)
            zic = z(ic)
            ja = iang(1,j)
            jc = iang(3,j)             
            xja = x(ja)
            yja = y(ja)
            zja = z(ja)
            xjc = x(jc)
            yjc = y(jc)
            zjc = z(jc) 
c     determine the atoms involved is H-bond atoms
c            print *," NCGHB is ", ncghb
            do k = 1, ncghb
                if (cghb_types(1,k) .NE. 0) then
                   if (((cghb_types(2,k) .eq. type(ib))
     &               .and. (cghb_types(5,k) .eq. type(jb))
     &               .and. (((cghb_types(1,k) .eq. type(ia))
     &               .and. (cghb_types(3,k) .eq. type(ic)))
     &               .or.  ((cghb_types(1,k) .eq. type(ic))
     &               .and.  (cghb_types(3,k) .eq. type(ia))))
     &               .and. (((cghb_types(4,k) .eq. type(ja))
     &               .and. (cghb_types(6,k) .eq. type(jc)))
     &               .or.  ((cghb_types(4,k) .eq. type(jc))
     &               .and.  (cghb_types(6,k) .eq. type(ja)))) 
     &               .or.
     &                     (cghb_types(2,k) .eq. type(jb))
     &               .and. (cghb_types(5,k) .eq. type(ib))
     &               .and. (((cghb_types(1,k) .eq. type(ja))
     &               .and. (cghb_types(3,k) .eq. type(jc)))
     &               .or.  ((cghb_types(1,k) .eq. type(jc))
     &               .and.  (cghb_types(3,k) .eq. type(ja))))
     &               .and. (((cghb_types(4,k) .eq. type(ia))
     &               .and. (cghb_types(6,k) .eq. type(ic)))
     &               .or.  ((cghb_types(4,k) .eq. type(ic))
     &               .and.  (cghb_types(6,k) .eq. type(ia)))))
     &     .and. (((ib - jb) .gt. 10) .or. ((jb-ib) .gt. 10))) then

                      ehbmax = cghb_ehbmax(k)
                      hbondmin = cghb_hbondmin(k)
                      cghbondcutoff = cghb_cghbondcutoff(k)
                      r03 = hbondmin*hbondmin*hbondmin
                      cghbondcutoff2 = cghbondcutoff*cghbondcutoff
                      if (rhbond2 .lt. cghbondcutoff2) then
                         goto 20
                      else
                         goto 30
                      endif
                   endif

                else
                      goto 30
                endif
            end do
C Get here if no hydrogen bond match. Go to next iteration
            goto 30

   20       continue

c     compute the value of the 1st angle i (alpha) 
              xiab = xia - xib
              yiab = yia - yib
              ziab = zia - zib
              xicb = xic - xib
              yicb = yic - yib
              zicb = zic - zib    
              xip = yiab*zicb - yicb*ziab
              yip = (-1.0d0)*xiab*zicb + xicb*ziab
              zip = xiab*yicb - xicb*yiab
              rip = sqrt(xip*xip + yip*yip + zip*zip)
              rip = max(rip,0.001d0)
              doti = xip*rx + yip*ry + zip*rz
              cosinei = doti / (rip*rhbond) 
              cosinei = min(1.0d0,max(-1.0d0,cosinei))
              anglei = acos(cosinei)
c     compute the second angle j (beta)         
              xjab = xja - xjb
              yjab = yja - yjb
              zjab = zja - zjb
              xjcb = xjc - xjb
              yjcb = yjc - yjb
              zjcb = zjc - zjb    
              xjp = yjab*zjcb - yjcb*zjab
              yjp = (-1.0d0)*xjab*zjcb + xjcb*zjab
              zjp = xjab*yjcb - xjcb*yjab
              rjp = sqrt(xjp*xjp + yjp*yjp + zjp*zjp)
              rjp = max(rjp,0.001d0)
              dotj = xjp*rx + yjp*ry + zjp*rz
              cosinej = dotj / (rjp*rhbond) 
              cosinej = min(1.0d0,max(-1.0d0,cosinej))
              anglej = acos(cosinej)
              anglek = anglei + anglej
              if (anglek .gt. (pi/2) .and. anglek .lt. (3*pi/2)) then
c               NTHREADS = OMP_GET_NUM_THREADS()
c               print *, '1Number of threads = ', NTHREADS
               if (anglek .gt. pi) anglek = 2*pi - anglek
               anglek = anglek*2 - pi
c      potential energy for CG Hbond               
               e = ((-1.0d0)*ehbmax/2.0d0)*(1-cos(anglek))*
     &                               (r03/(rhbond2*rhbond))
ccccccc         Added by David  12/6/16
cc               if (rhbond .lt. hbondmin) then
cc                e = (-1.0d0)*(ehbmax/2)*(1-cos(anglek))*
cc     &          (hbondmin/(rhbond))
cc               endif
cccccc
c      first derivatives for angle 
               dedt = ((-1.0d0)*ehbmax/2.0d0)*sin(anglek)*     
     &                          (r03/(rhbond2*rhbond)) 

ccccccc         Added by David  12/6/16
cc               if (rhbond .lt. hbondmin) then
cc                dedt = ((-1.0d0)*ehbmax/2.0d0)*sin(anglek)*
cc     &          (hbondmin/(rhbond))
cc               endif
ccccccc


               if ((anglei+anglej) .gt. pi) then
                  danglekdangle = -2
               else
                  danglekdangle = 2
               endif

               dangledcosinei = (-1 / sqrt(1-cosinei*cosinei)) 
               dangledcosinej = (-1 / sqrt(1-cosinej*cosinej))

c               write(*,*) "atoms anglei ", ia,ib,ic
c               write(*,*) "atoms anglej ", ja,jb,jc 
c               write(*,*) "dedt,anglek,rhbond " , dedt,anglek,rhbond
c               write(*,*) "cosinei,dangledcosinei " , cosinei,
c     &                                   dangledcosinei


c   #########  residue i ####### de = dedt*dnaglekdangle*dangledcosine*dcosinedx
 
               dedtxia = dedt * danglekdangle * dangledcosinei * 
c &                          (((-1.0d0)*zicb*ry + yicb*rz)/(rip*rhbond)
     &                          ((yicb*rz - zicb*ry)/(rip*rhbond)
c     &                          - doti*(yip*(-1.0d0)*zicb + zip*yicb)/
     &                          - doti*(zip*yicb -  yip*zicb)/
     &                          (rhbond*rip**3))

               dedtyia = dedt * danglekdangle * dangledcosinei * 
     &                          ((zicb*rx - xicb*rz)/(rip*rhbond)
     &                          - doti*(xip*zicb - zip*xicb)/
     &                          (rhbond*rip**3)) 

               dedtzia = dedt * danglekdangle * dangledcosinei *
c     &                          (((-1.0d0)*yicb*rx + xicb*ry)/(rip*rhbond)
     &                          ((xicb*ry - yicb*rx)/(rip*rhbond))


c For atom b, have dependence on anglei and anglej

               dedtxib = dedt * danglekdangle * ( dangledcosinei *
     &                          (((xip + (zicb-ziab)*ry + 
     &                          (yiab-yicb)*rz)/(rip*rhbond))
     &                          - (doti*(yip*(zicb-ziab) + 
     &                          zip*(yiab-yicb))/
     &                          (rhbond*rip**3)) - (doti * rx/
     &                          (rip * rhbond**3))) + dangledcosinej *
     &                          ((xjp/(rjp*rhbond)) - (dotj * rx/
     &                          (rjp*rhbond**3))))

c               write(*,*) "dedtxib,dedt,danglekdangle,dangledcosinei ",
c     &          dedtxib,dedt,danglekdangle,dangledcosinei
c               write(*,*) "xip,zicb,ziab,ry,yiab,yicb,rz,rip,rhbond",
c     &          xip,zicb,ziab,ry,yiab,yicb,rz,rip,rhbond
c               write(*,*) "doti,yip,zip,rx",
c     &          doti,yip,zip,rx


               dedtyib = dedt * danglekdangle * ( dangledcosinei *
     &                          ((((ziab-zicb)*rx + yip + 
     &                          (xicb-xiab)*rz)/(rip*rhbond))
     &                          - (doti*(xip*(ziab-zicb) +
     &                          zip*(xicb-xiab))/
     &                          (rhbond*rip**3)) - (doti * ry/
     &                          (rip * rhbond**3))) + dangledcosinej *
     &                          ((yjp/(rjp*rhbond)) - (dotj * ry/ 
     &                          (rjp*rhbond**3))))


               dedtzib = dedt * danglekdangle * ( dangledcosinei *
     &                          ((((yicb-yiab)*rx + 
     &                          (xiab-xicb)*ry + zip)/(rip*rhbond))
     &                          - (doti*(xip*(yicb-yiab) +
     &                          yip*(xiab-xicb))/
     &                          (rhbond*rip**3)) - (doti * rz/
     &                          (rip * rhbond**3))) + dangledcosinej *
     &                          ((zjp/(rjp*rhbond)) - (dotj * rz/ 
     &                          (rjp*rhbond**3))))


               dedtxic = dedt * danglekdangle * dangledcosinei *
     &                          ((ziab*ry - yiab*rz)/(rip*rhbond)
     &                          - doti*(yip*ziab - zip*yiab)/
     &                          (rhbond*rip**3))

               dedtyic = dedt * danglekdangle * dangledcosinei *
c     &                          (((-1.0d0)*ziab*rx + xiab*rz)/(rip*rhbond)
     &                          ((xiab*rz - ziab*rx)/(rip*rhbond)
     &                          - doti*((-1.0d0)*xip*ziab + zip*xiab)/
     &                          (rhbond*rip**3))

               dedtzic = dedt * danglekdangle * dangledcosinei *
     &                          ((yiab*rx - xiab*ry)/(rip*rhbond)
     &                          - doti*(xip*yiab - yip*xiab)/
     &                          (rhbond*rip**3))

c  ##########  residue j ########

               dedtxja = dedt * danglekdangle * dangledcosinej *
c     &                          (((-1.0d0)*zjcb*ry + yjcb*rz)/(rjp*rhbond)
c     &                          - dotj*((-1.0d0)*yjp*zjcb + zjp*yjcb)/
     &                          ((yjcb*rz - zjcb*ry)/(rjp*rhbond)
     &                          - dotj*(zjp*yjcb - yjp*zjcb)/
     &                          (rhbond*rjp**3))

               dedtyja = dedt * danglekdangle * dangledcosinej *
     &                          ((zjcb*rx - xjcb*rz)/(rjp*rhbond)
     &                          - dotj*(xjp*zjcb - zjp*xjcb)/
     &                          (rhbond*rjp**3))

               dedtzja = dedt * danglekdangle * dangledcosinej *
c     &                          (((-1.0d0)*yjcb*rx + xjcb*ry)/(rjp*rhbond)
c     &                          - dotj*((-1.0d0)*xjp*yjcb + yjp*xjcb)/
     &                          ((xjcb*ry - yjcb*rx)/(rjp*rhbond)
     &                          - dotj*(yjp*xjcb - xjp*yjcb)/
     &                          (rhbond*rjp**3))

               dedtxjb = dedt * danglekdangle * ( dangledcosinej *
c     &                          ((((-1.0d0)*xjp + (zjcb-zjab)*ry + 
     &                          ((((-1.0d0)*xjp + (zjcb-zjab)*ry + 
     &                          (yjab-yjcb)*rz)/(rjp*rhbond))
     &                          - (dotj*(yjp*(zjcb-zjab) +
     &                          zjp*(yjab-yjcb))/
     &                          (rhbond*rjp**3)) + (dotj * rx/
c     &                          (rjp * rhbond**3))) + dangledcosinei *
     &                          (rjp * rhbond**3))) - dangledcosinei *
c     &                          (((-1.0d0)*xip/(rip*rhbond)) + (doti * rx/ 
     &                          ((xip/(rip*rhbond)) + (doti * rx/ 
     &                          (rip*rhbond**3))))

               dedtyjb = dedt * danglekdangle * ( dangledcosinej *
     &                          ((((zjab-zjcb)*rx + (-1.0d0)*yjp + 
     &                          (xjcb-xjab)*rz)/(rjp*rhbond))
     &                          - (dotj*(xjp*(zjab-zjcb) +
     &                          zjp*(xjcb-xjab))/
     &                          (rhbond*rjp**3)) + (dotj * ry/
c     &                          (rjp * rhbond**3))) + dangledcosinei *
c     &                          (((-1.0d0)*yip/(rip*rhbond)) + (doti * ry/
     &                          (rjp * rhbond**3))) - dangledcosinei *
     &                          ((yip/(rip*rhbond)) + (doti * ry/
     &                          (rip*rhbond**3))))

               dedtzjb = dedt * danglekdangle * ( dangledcosinej *
     &                          ((((yjcb-yjab)*rx + (xjab-xjcb)*ry + 
     &                          (-1.0d0)*zjp)/(rjp*rhbond))
     &                          - (dotj*(xjp*(yjcb-yjab) +
     &                          yjp*(xjab-xjcb))/
     &                          (rhbond*rjp**3)) + (dotj * rz/
c     &                          (rjp * rhbond**3))) + dangledcosinei *
c     &                          (((-1.0d0)*zip/(rip*rhbond)) + (doti * rz/
     &                          (rjp * rhbond**3))) - dangledcosinei *
     &                          ((zip/(rip*rhbond)) + (doti * rz/
     &                          (rip*rhbond**3))))


               dedtxjc = dedt * danglekdangle * dangledcosinej *
     &                          ((zjab*ry - yjab*rz)/(rjp*rhbond)
     &                          - dotj*(yjp*zjab - zjp*yjab)/
     &                          (rhbond*rjp**3))

               dedtyjc = dedt * danglekdangle * dangledcosinej *
c     &                          (((-1.0d0)*zjab*rx + xjab*rz)/(rjp*rhbond)
c     &                          - dotj*((-1.0d0)*xjp*zjab + zjp*xjab)/
     &                          ((xjab*rz - zjab*rx)/(rjp*rhbond)
     &                          - dotj*(zjp*xjab - xjp*zjab)/
     &                          (rhbond*rjp**3))

               dedtzjc = dedt * danglekdangle * dangledcosinej *
     &                          ((yjab*rx - xjab*ry)/(rjp*rhbond)
     &                          - dotj*(xjp*yjab - yjp*xjab)/
     &                          (rhbond*rjp**3))


c                write(*,*) "dedtzjc,dedt,dangledcosinej,yjab,rx",
c     &          dedtzjc,dedt,dangledcosinej,yjab,rx
c                write(*,*) "xjab,ry,rjp,rhbond,dotj,xjp,yjp",
c     &          xjab,ry,rjp,rhbond,dotj,xjp,yjp
c                write(*,*) "danglekdangle,anglek", danglekdangle,anglek
c                write(*,*) "anglei,anglej", anglei,anglej

c               check

c               if (anglei .gt. (pi/2)) then
c               dedtbix = dedtbi*xip 
c               dedtbiy = dedtbi*yip
c               dedtbiz = dedtbi*zip
c               else
c               dedtbix = (-1.0d0)*dedtbi*xip 
c               dedtbiy = (-1.0d0)*dedtbi*yip
c               dedtbiz = (-1.0d0)*dedtbi*zip 
c               end if 
c               dedtbj = (-1.0d0)*ehbmax/2.0d0*sin(anglek)*
c     &                          (r03/rhbond2/rhbond)/rjp                
c               if (anglej .gt. (pi/2)) then
c               dedtbjx = (-1.0d0)*dedtbj*xjp 
c               dedtbjy = (-1.0d0)*dedtbj*yjp
c               dedtbjz = (-1.0d0)*dedtbj*zjp
c               else
c               dedtbjx = dedtbj*xjp 
c               dedtbjy = dedtbj*yjp
c               dedtbjz = dedtbj*zjp 
c               end if                  

     
c      first derivatives for distance 
               dedr = 1.5d0*ehbmax*(1-cos(anglek))*
     &             (r03/(rhbond2*rhbond2))

ccccccc         Added by David  12/6/16
cc               if (rhbond .lt. hbondmin) then
cc                dedr = (ehbmax/2)*(1-cos(anglek))*
cc     &             (hbondmin/rhbond2) 
cc               endif
ccccccc

               dedrxib = dedr * rx/rhbond  
               dedryib = dedr * ry/rhbond
               dedrzib = dedr * rz/rhbond

               dedrxjb = dedr * (-1.0d0)*rx/rhbond
               dedryjb = dedr * (-1.0d0)*ry/rhbond
               dedrzjb = dedr * (-1.0d0)*rz/rhbond
        
c               write(*,*) "dedr,anglek,rhbond",dedr,anglek,rhbond
c               write(*,*) "dedrxib,rx",dedrxib,rx

c               dedrbix = dedr * rx
c               dedrbiy = dedr * ry
c               dedrbiz = dedr * rz
c               dedrbjx = (-1.0d0)*dedr * rx
c               dedrbjy = (-1.0d0)*dedr * ry
c               dedrbjz = (-1.0d0)*dedr * rz 
c     total derivatives
c               dedbix = dedtbix + dedrbix 
c               dedbiy = dedtbiy + dedrbiy
c               dedbiz = dedtbiz + dedrbiz  
c               dedbjx = dedtbjx + dedrbjx 
c               dedbjy = dedtbjy + dedrbjy
c               dedbjz = dedtbjz + dedrbjz                                                     


c      total derivatives

                dedxia = dedtxia 
                dedyia = dedtyia
                dedzia = dedtzia

                dedxib = dedtxib + dedrxib
                dedyib = dedtyib + dedryib
                dedzib = dedtzib + dedrzib

                dedxic = dedtxic 
                dedyic = dedtyic
                dedzic = dedtzic


                dedxja = dedtxja 
                dedyja = dedtyja
                dedzja = dedtzja

                dedxjb = dedtxjb + dedrxjb
                dedyjb = dedtyjb + dedryjb
                dedzjb = dedtzjb + dedrzjb

                dedxjc = dedtxjc
                dedyjc = dedtyjc
                dedzjc = dedtzjc


              else 
               e = 0.0d0
               dedxia = 0.0d0
               dedyia = 0.0d0
               dedzia = 0.0d0
               dedxib = 0.0d0
               dedyib = 0.0d0
               dedzib = 0.0d0
               dedxic = 0.0d0
               dedyic = 0.0d0
               dedzic = 0.0d0
               dedxja = 0.0d0
               dedyja = 0.0d0
               dedzja = 0.0d0
               dedxjb = 0.0d0
               dedyjb = 0.0d0
               dedzjb = 0.0d0
               dedxjc = 0.0d0
               dedyjc = 0.0d0
               dedzjc = 0.0d0
              end if
              ecghbo = ecghbo + e
              decghbo(1,ia) = decghbo(1,ia) + dedxia
              decghbo(2,ia) = decghbo(2,ia) + dedyia
              decghbo(3,ia) = decghbo(3,ia) + dedzia
              decghbo(1,ib) = decghbo(1,ib) + dedxib
              decghbo(2,ib) = decghbo(2,ib) + dedyib
              decghbo(3,ib) = decghbo(3,ib) + dedzib
              decghbo(1,ic) = decghbo(1,ic) + dedxic
              decghbo(2,ic) = decghbo(2,ic) + dedyic
              decghbo(3,ic) = decghbo(3,ic) + dedzic

              decghbo(1,ja) = decghbo(1,ja) + dedxja
              decghbo(2,ja) = decghbo(2,ja) + dedyja
              decghbo(3,ja) = decghbo(3,ja) + dedzja
              decghbo(1,jb) = decghbo(1,jb) + dedxjb
              decghbo(2,jb) = decghbo(2,jb) + dedyjb
              decghbo(3,jb) = decghbo(3,jb) + dedzjb
              decghbo(1,jc) = decghbo(1,jc) + dedxjc
              decghbo(2,jc) = decghbo(2,jc) + dedyjc
              decghbo(3,jc) = decghbo(3,jc) + dedzjc

             
c              write(*,*) "final dedxib", dedxib 
c              write(*,*) "decghb", decghb(1,ib),decghb(1,jb), ib, jb

c              write (*,*) "ecghb", e, ib, jb, type(ib), type(jb), 
c     &             radian *anglei, radian *anglej, radian *anglek
           end if 
   30      continue
          end do
        end if
      end do
c
c     end OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
c
c     add local copies to global variables for OpenMP calculation
c
      ecghb = ecghbo
      do i = 1, n 
         decghb(1,i) = decghb(1,i) + decghbo(1,i) 
         decghb(2,i) = decghb(2,i) + decghbo(2,i) 
         decghb(3,i) = decghb(3,i) + decghbo(3,i) 
      end do


      deallocate(decghbo)

      return
      end
