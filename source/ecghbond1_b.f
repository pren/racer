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
      implicit none
      integer i,ia,ib,ic,id
      integer j,ja,jb,jc,jd
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
 
      r03 = hbondmin*hbondmin*hbondmin

      cghbondcutoff2 = cghbondcutoff*cghbondcutoff
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
      decghbo = decghb
c
c     set OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(nangle,iang,x,y,z,
!$OMP& ehbmax,cghbondcutoff2,r03)
!$OMP& shared(ecghbo,decghbo)
!$OMP DO reduction(+:ecghbo,decghbo) schedule(guided)


c
c     calculate the bond angle bending energy term
c
      do i = 1, nangle
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
c     determine atom type is N6 N2 O6 O2 4 5 6 7    
         if ((type(ib) .GT. 3) .AND. (type(ib) .LT. 8)) then
         xib = x(ib)
         yib = y(ib)
         zib = z(ib)
c         write (*,*) "ia ib ic",ia,ib,ic
        do j = i+1, nangle
            jb = iang(2,j)
            if ((type(jb) .GT. 3) .AND. (type(jb) .LT. 8)) then
            xjb = x(jb)
            yjb = y(jb)
            zjb = z(jb)    
            rx = xib - xjb
            ry = yib - yjb
            rz = zib - zjb               
            rhbond2 = rx*rx + ry*ry + rz*rz
            rhbond = sqrt(rhbond2)
            if (rhbond2 .lt. cghbondcutoff2) then
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
          if (((((type(ib) .eq. 6) .and. ((type(ia) .eq. 5) 
     &                      .or. (type(ic) .eq. 5))) .and.    
     &         ((type(jb) .eq. 4) .and. ((type(ja) .eq. 7) 
     &                      .or. (type(jc) .eq. 7)))) .or.      
     &        (((type(ib) .eq. 4) .and. ((type(ia) .eq. 9) 
     &                      .or. (type(ic) .eq. 9))) .and.    
     &         ((type(jb) .eq. 6) .and. ((type(ja) .eq. 7)
     &                      .or. (type(jc) .eq. 7)))) .or.     
     &        (((type(ib) .eq. 4) .and. ((type(ia) .eq. 7)
     &                      .or. (type(ic) .eq. 7))) .and.    
     &         ((type(jb) .eq. 6) .and. ((type(ja) .eq. 5)
     &                      .or. (type(jc) .eq. 5)))) .or.      
     &        (((type(ib) .eq. 6) .and. ((type(ia) .eq. 7)
     &                      .or. (type(ic) .eq. 7))) .and.    
     &         ((type(jb) .eq. 4) .and. ((type(ja) .eq. 9)
     &                      .or. (type(jc) .eq. 9)))) .or.       
     &        (((type(ib) .eq. 5) .and. ((type(ia) .eq. 6) 
     &                      .or. (type(ic) .eq. 6))) .and.    
     &         ((type(jb) .eq. 7) .and. ((type(ja) .eq. 4)
     &                      .or. (type(jc) .eq. 4)))) .or.    
     &        (((type(ib) .eq. 7) .and. ((type(ia) .eq. 4)
     &                      .or. (type(ic) .eq. 4))) .and.   
     &         ((type(jb) .eq. 5) .and. ((type(ja) .eq. 6) 
     &                    .or. (type(jc) .eq. 6))))) .and.
     & (((ib - jb) .gt. 10) .or. ((jb-ib) .gt. 10)))  then     
c     compute the value of the 1st angle i (alpha) 
              xiab = xia - xib
              yiab = yia - yib
              ziab = zia - zib
              xicb = xic - xib
              yicb = yic - yib
              zicb = zic - zib    
              xip = yiab*zicb - yicb*ziab
              yip = -1*xiab*zicb + xicb*ziab
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
              yjp = -1*xjab*zjcb + xjcb*zjab
              zjp = xjab*yjcb - xjcb*yjab
              rjp = sqrt(xjp*xjp + yjp*yjp + zjp*zjp)
              rjp = max(rjp,0.001d0)
              dotj = xjp*rx + yjp*ry + zjp*rz
              cosinej = dotj / (rjp*rhbond) 
              cosinej = min(1.0d0,max(-1.0d0,cosinej))
              anglej = acos(cosinej)
              anglek = anglei + anglej
              if (anglek .gt. (pi/2) .and. anglek .lt. (3*pi/2)) then
               if (anglek .gt. pi) anglek = 2*pi - anglek
               anglek = anglek*2 - pi
c      potential energy for CG Hbond               
               e = (-1.0d0*ehbmax/2.0d0)*(1-cos(anglek))*
     &                               (r03/(rhbond2*rhbond))
c      first derivatives for angle 
               dedt = (-1.0d0*ehbmax/2.0d0)*sin(anglek)*     
     &                          (r03/(rhbond2*rhbond)) 


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
     &                          ((-1*zicb*ry + yicb*rz)/(rip*rhbond)
     &                          - doti*(yip*-1*zicb + zip*yicb)/
     &                          (rhbond*rip**3))

               dedtyia = dedt * danglekdangle * dangledcosinei * 
     &                          ((zicb*rx - xicb*rz)/(rip*rhbond)
     &                          - doti*(xip*zicb - zip*xicb)/
     &                          (rhbond*rip**3)) 

               dedtzia = dedt * danglekdangle * dangledcosinei *
     &                          ((-1*yicb*rx + xicb*ry)/(rip*rhbond)
     &                          - doti*(-1*xip*yicb + yip*xicb)/
     &                          (rhbond*rip**3))

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
     &                          ((-1*ziab*rx + xiab*rz)/(rip*rhbond)
     &                          - doti*(-1*xip*ziab + zip*xiab)/
     &                          (rhbond*rip**3))

               dedtzic = dedt * danglekdangle * dangledcosinei *
     &                          ((yiab*rx - xiab*ry)/(rip*rhbond)
     &                          - doti*(xip*yiab - yip*xiab)/
     &                          (rhbond*rip**3))

c  ##########  residue j ########

               dedtxja = dedt * danglekdangle * dangledcosinej *
     &                          ((-1*zjcb*ry + yjcb*rz)/(rjp*rhbond)
     &                          - dotj*(-1*yjp*zjcb + zjp*yjcb)/
     &                          (rhbond*rjp**3))

               dedtyja = dedt * danglekdangle * dangledcosinej *
     &                          ((zjcb*rx - xjcb*rz)/(rjp*rhbond)
     &                          - dotj*(xjp*zjcb - zjp*xjcb)/
     &                          (rhbond*rjp**3))

               dedtzja = dedt * danglekdangle * dangledcosinej *
     &                          ((-1*yjcb*rx + xjcb*ry)/(rjp*rhbond)
     &                          - dotj*(-1*xjp*yjcb + yjp*xjcb)/
     &                          (rhbond*rjp**3))

               dedtxjb = dedt * danglekdangle * ( dangledcosinej *
     &                          (((-1*xjp + (zjcb-zjab)*ry + 
     &                          (yjab-yjcb)*rz)/(rjp*rhbond))
     &                          - (dotj*(yjp*(zjcb-zjab) +
     &                          zjp*(yjab-yjcb))/
     &                          (rhbond*rjp**3)) + (dotj * rx/
     &                          (rjp * rhbond**3))) + dangledcosinei *
     &                          ((-1*xip/(rip*rhbond)) + (doti * rx/ 
     &                          (rip*rhbond**3))))

               dedtyjb = dedt * danglekdangle * ( dangledcosinej *
     &                          ((((zjab-zjcb)*rx + -1*yjp + 
     &                          (xjcb-xjab)*rz)/(rjp*rhbond))
     &                          - (dotj*(xjp*(zjab-zjcb) +
     &                          zjp*(xjcb-xjab))/
     &                          (rhbond*rjp**3)) + (dotj * ry/
     &                          (rjp * rhbond**3))) + dangledcosinei *
     &                          ((-1*yip/(rip*rhbond)) + (doti * ry/
     &                          (rip*rhbond**3))))

               dedtzjb = dedt * danglekdangle * ( dangledcosinej *
     &                          ((((yjcb-yjab)*rx + (xjab-xjcb)*ry + 
     &                          -1*zjp)/(rjp*rhbond))
     &                          - (dotj*(xjp*(yjcb-yjab) +
     &                          yjp*(xjab-xjcb))/
     &                          (rhbond*rjp**3)) + (dotj * rz/
     &                          (rjp * rhbond**3))) + dangledcosinei *
     &                          ((-1*zip/(rip*rhbond)) + (doti * rz/
     &                          (rip*rhbond**3))))


               dedtxjc = dedt * danglekdangle * dangledcosinej *
     &                          ((zjab*ry - yjab*rz)/(rjp*rhbond)
     &                          - dotj*(yjp*zjab - zjp*yjab)/
     &                          (rhbond*rjp**3))

               dedtyjc = dedt * danglekdangle * dangledcosinej *
     &                          ((-1*zjab*rx + xjab*rz)/(rjp*rhbond)
     &                          - dotj*(-1*xjp*zjab + zjp*xjab)/
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
c               dedtbix = -1.0d0*dedtbi*xip 
c               dedtbiy = -1.0d0*dedtbi*yip
c               dedtbiz = -1.0d0*dedtbi*zip 
c               end if 
c               dedtbj = -1.0d0*ehbmax/2.0d0*sin(anglek)*
c     &                          (r03/rhbond2/rhbond)/rjp                
c               if (anglej .gt. (pi/2)) then
c               dedtbjx = -1.0d0*dedtbj*xjp 
c               dedtbjy = -1.0d0*dedtbj*yjp
c               dedtbjz = -1.0d0*dedtbj*zjp
c               else
c               dedtbjx = dedtbj*xjp 
c               dedtbjy = dedtbj*yjp
c               dedtbjz = dedtbj*zjp 
c               end if                  

     
c      first derivatives for distance 
               dedr = 1.5d0*ehbmax*(1-cos(anglek))*
     &             (r03/(rhbond2*rhbond2))

               dedrxib = dedr * rx/rhbond  
               dedryib = dedr * ry/rhbond
               dedrzib = dedr * rz/rhbond

               dedrxjb = dedr * -1*rx/rhbond
               dedryjb = dedr * -1*ry/rhbond
               dedrzjb = dedr * -1*rz/rhbond
        
c               write(*,*) "dedr,anglek,rhbond",dedr,anglek,rhbond
c               write(*,*) "dedrxib,rx",dedrxib,rx

c               dedrbix = dedr * rx
c               dedrbiy = dedr * ry
c               dedrbiz = dedr * rz
c               dedrbjx = -1.0d0*dedr * rx
c               dedrbjy = -1.0d0*dedr * ry
c               dedrbjz = -1.0d0*dedr * rz 
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
            end if
           end if 
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
      decghb = decghbo

      return
      end
