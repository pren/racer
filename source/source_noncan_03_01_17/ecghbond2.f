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
      subroutine ecghbond2
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
      use virial
      use hessn
      implicit none
      integer i,ia,ib,ic,id
      integer j,ja,jb,jc,jd,jj
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
      real*8 e
      real*8 dot,cosine
      real*8 anglei,anglej,anglek,fgrp
      real*8 rx,ry,rz
      real*8 xab,yab,zab
      real*8 xcb,ycb,zcb
      real*8 cghbondcutoff, cghbondcutoff2
      real*8 ehbmax,r0,r03
      logical proceed
      real *8 dedr,dedtbi,dedtbj
      real *8 dedtbix,dedtbiy,dedtbiz
      real *8 dedtbjx,dedtbjy,dedtbjz
      real *8 dedrbix,dedrbiy,dedrbiz
      real *8 dedrbjx,dedrbjy,dedrbjz
      real *8 dedbix,dedbiy,dedbiz
      real *8 dedbjx,dedbjy,dedbjz  
      real *8 d2edrdr,d2edtdt  
      real *8 termdr,termrx,termry,termrz
      real *8 termdtbi,termtbix,termtbiy,termtbiz 
      real *8 termdtbj,termtbjx,termtbjy,termtbjz 
      real *8 d2edr(3,3),d2edtib(3,3),d2edtjb(3,3) 
c
c
c      
      cghbondcutoff = 6.0d0
      ehbmax=0.5d0
      r0=2.90d0
      r03=r0*r0*r0
      cghbondcutoff2 = cghbondcutoff*cghbondcutoff
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
c     calculate the bond angle bending energy term
c
      do i = 1, nangle
         ib = iang(2,i)
         ia = iang(1,i)
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
            rx = xjb - xib
            ry = yjb - yib
            rz = zjb - zib               
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
              xab = xia - xib
              yab = yia - yib
              zab = zia - zib
              xcb = xic - xib
              ycb = yic - yib
              zcb = zic - zib    
              xip = ycb*zab - zcb*yab
              yip = zcb*xab - xcb*zab
              zip = xcb*yab - ycb*xab
              rip = sqrt(xip*xip + yip*yip + zip*zip)
              rip = max(rip,0.001d0)
              dot = xip*rx + yip*ry + zip*rz
              cosine = dot / (rip*rhbond) 
              cosine = min(1.0d0,max(-1.0d0,cosine))
              anglei = acos(cosine)
c              anglei = radian * acos(cosine) 
c     computer the second angle j (beta)         
              xab = xja - xjb
              yab = yja - yjb
              zab = zja - zjb
              xcb = xjc - xjb
              ycb = yjc - yjb
              zcb = zjc - zjb    
              xjp = ycb*zab - zcb*yab
              yjp = zcb*xab - xcb*zab
              zjp = xcb*yab - ycb*xab
              rjp = sqrt(xjp*xjp + yjp*yjp + zjp*zjp)
              rjp = max(rjp,0.001d0)
              dot = xjp*rx + yjp*ry + zjp*rz
              cosine = dot / (rjp*rhbond) 
              cosine = min(1.0d0,max(-1.0d0,cosine))
              anglej = acos(cosine)
c              anglej = radjan * acos(cosine)              
              anglek = anglei + anglej
              if (anglek .gt. (pi/2) .and. anglek .lt. (3*pi/2)) then
               if (anglek .gt. pi) anglek = 2*pi - anglek
               anglek = anglek*2 - pi
c      potential energy for CG Hbond               
               e = -1.0d0*ehbmax/2.0d0*(1-cos(anglek))*
     &                               r03/rhbond2/rhbond
c      first derivatives for angle 
             if (anglei .gt. (pi/2)) then
               dedtbi = -1.0d0*ehbmax/2.0d0*sin(anglek)*
     &                            r03/rhbond2/rhbond/rip 
             else 
               dedtbi = ehbmax/2.0d0*sin(anglek)*
     &                            r03/rhbond2/rhbond/rip 
             end if 
               dedtbix = dedtbi*xip 
               dedtbiy = dedtbi*yip
               dedtbiz = dedtbi*zip 
c               
              if (anglej .gt. (pi/2)) then
               dedtbj = ehbmax/2.0d0*sin(anglek)*
     &                          (r03/rhbond2/rhbond)/rjp                
              else
               dedtbj = -1.0d0*ehbmax/2.0d0*sin(anglek)*
     &                          (r03/rhbond2/rhbond)/rjp                
               dedtbjx = dedtbj*xjp 
               dedtbjy = dedtbj*yjp
               dedtbjz = dedtbj*zjp 
               end if                       
c      first derivatives for distance 
               dedr = -1.5d0*ehbmax*(1-cos(anglek))*
     &             r03/rhbond2/rhbond2/rhbond        
               dedrbix = dedr * rx
               dedrbiy = dedr * ry
               dedrbiz = dedr * rz
               dedrbjx = -1.0d0*dedr * rx
               dedrbjy = -1.0d0*dedr * ry
               dedrbjz = -1.0d0*dedr * rz 
c     total derivatives
               dedbix = dedtbix + dedrbix 
               dedbiy = dedtbiy + dedrbiy
               dedbiz = dedtbiz + dedrbiz  
               dedbjx = dedtbjx + dedrbjx 
               dedbjy = dedtbjy + dedrbjy
               dedbjz = dedtbjz + dedrbjz 
c      second derivatives for distance and angle
               d2edrdr = 6.0d0*ehbmax*(1-cos(anglek))*
     &             r03/rhbond2/rhbond2/rhbond 
               d2edtdt =  -1.0d0*ehbmax/2.0d0*r03/rhbond2/
     &                 rhbond*cos(anglek)
c               
c               d2e = de2dr + de2dt
c             dei = sart(dedbix*dedbix+dedbiy*dedbiy+dedbiz*dedbiz)
c             dej = sart(dedbjx*dedbjx+dedbjy*dedbjy+dedbjz*dedbjz)
c      second derivatives for distance       
              if (rhbond2 .eq. 0.0d0) then
               dedr = 0.0d0
               termdr = 0.0d0
              else
              termdr = (d2edrdr - dedr)/rhbond2  
              end if
              termrx = termdr * rx
              termry = termdr * ry
              termrz = termdr * rz
              d2edr(1,1) = termrx*rx + dedr
              d2edr(1,2) = termrx*ry
              d2edr(1,3) = termrx*rz
              d2edr(2,1) = d2edr(1,2)
              d2edr(2,2) = termry*ry + dedr
              d2edr(2,3) = termry*rz
              d2edr(3,1) = d2edr(1,3)
              d2edr(3,2) = d2edr(2,3)
              d2edr(3,3) = termrz*rz + dedr
c     second derivatives for angle
c     for atom bi
              if (rip .eq. 0.0d0) then
               dedtbi = 0.0d0
               termdtbi = 0.0d0
              else
              termdtbi = (d2edtdt - dedtbi)/(rip*rip)
              end if
              termtbix = termdtbi * xip
              termtbiy = termdtbi * yip
              termtbiz = termdtbi * zip
              d2edtib(1,1) = termtbix*xip + dedtbi
              d2edtib(1,2) = termtbix*yip
              d2edtib(1,3) = termtbix*zip
              d2edtib(2,1) = d2edtib(1,2)
              d2edtib(2,2) = termtbiy*yip + dedtbi
              d2edtib(2,3) = termtbiy*zip
              d2edtib(3,1) = d2edtib(1,3)
              d2edtib(3,2) = d2edtib(2,3)
              d2edtib(3,3) = termtbiz*zip + dedtbi
c     second derivatives for angle
c     for atom bj              
              if (rjp .eq. 0.0d0) then
               dedtbj = 0.0d0
               termdtbj = 0.0d0
              else
              termdtbj = (d2edtdt - dedtbj)/(rjp*rjp)
              end if
              termtbjx = termdtbj * xjp
              termtbjy = termdtbj * yjp
              termtbjz = termdtbj * zjp
              d2edtjb(1,1) = termtbjx*xjp + dedtbj
              d2edtjb(1,2) = termtbjx*yjp
              d2edtjb(1,3) = termtbjx*zjp
              d2edtjb(2,1) = d2edtjb(1,2)
              d2edtjb(2,2) = termtbjy*yjp + dedtbj
              d2edtjb(2,3) = termtbjy*zjp
              d2edtjb(3,1) = d2edtjb(1,3)
              d2edtjb(3,2) = d2edtjb(2,3)
              d2edtjb(3,3) = termtbjz*zjp + dedtbj  
c      Add second derivatives for angle and distance togethor
          do jj = 1, 3
           hessx(jj,ib) = hessx(jj,ib) + d2edr(1,jj) + d2edtib(1,jj)
           hessy(jj,ib) = hessy(jj,ib) + d2edr(2,jj) + d2edtib(2,jj)
           hessz(jj,ib) = hessz(jj,ib) + d2edr(3,jj) + d2edtib(3,jj)
           hessx(jj,jb) = hessx(jj,jb) - d2edr(1,jj) + d2edtjb(1,jj)
           hessy(jj,jb) = hessy(jj,jb) - d2edr(2,jj) + d2edtjb(2,jj)
           hessz(jj,jb) = hessz(jj,jb) - d2edr(3,jj) + d2edtjb(3,jj)
          end do
c                                                                                                  
              else 
               e = 0.0d0
               dedbix = 0.0d0 
               dedbiy = 0.0d0
               dedbiz = 0.0d0  
               dedbjx = 0.0d0 
               dedbjy = 0.0d0
               dedbjz = 0.0d0
              end if
              ecghb = ecghb + e
              decghb(1,ib) = decghb(1,ib) + dedbix
              decghb(2,ib) = decghb(2,ib) + dedbiy
              decghb(3,ib) = decghb(3,ib) + dedbiz
              decghb(1,jb) = decghb(1,jb) + dedbjx
              decghb(2,jb) = decghb(2,jb) + dedbjy
              decghb(3,jb) = decghb(3,jb) + dedbjz              
c              write (*,*) "ecghb", e, ib, jb, type(ib), type(jb), 
c     &             radian *anglei, radian *anglej, radian *anglek
             end if
            end if
           end if 
          end do
        end if
      end do
      return
      end
