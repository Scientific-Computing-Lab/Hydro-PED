! Copyright (c) 2017, 
! Eyal Shalev (eyal@gsi.gov.il)
! Vladimir Lyakhovsky
! All rights reserved to Geological Survey of Israel (GSI) 
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!    * Redistributions of source code must retain the above copyright
!      notice, this list of conditions and the following disclaimer.
!    * Redistributions in binary form must reproduce the above copyright
!      notice, this list of conditions and the following disclaimer in the
!      documentation and/or other materials provided with the distribution.
!    * Neither the name of Eyal Shalev or Vladimir Lyakhovsky, nor the
!      names of its contributors may be used to endorse or promote products
!      derived from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
! ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
! WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL Eyal Shalev & Vladimir Lyakhovsky
! BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
! OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND  ON ANY THEORY OF LIABILITY, 
! WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
! OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
! ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

!------------Stress drop -----------------------
!--- F90 damage & porosity ---------------------
      
!> @brief         Calculates stress drop in failed elements
!> @details       Uses explicit time step for iterrative solution
!> of equations of motion. Append data to file "catalog" with 
!> information about failed elements
!>
!> @param[in]     time, time step
!> @param[out]    None

  subroutine drop (time,dt)
 
  use sizes
  use element_data
  use node_data
  use boundary_node_data
  implicit none    

!------------- Variables ------------------------     

  logical:: event,dynamic,space
  integer::n,j,ndrop,mm,nd,nfl,ii,i
  real(kind=8):: time,vvv0,dt,vvv,factor
  real(kind=8):: alp0,dist,si2
  real(kind=8):: dadt,dsf,dsf1,d_al,s(6),rle,rme
  real(kind=8):: d_coup,din,fx,fy,cst,df_el
  real(kind=8):: p_mean,s_mean,pot,sdrop,xx,yy,zz
  real(kind=8):: s11,s22,s33,s12,s13,s23
  real(kind=8):: p11,p22,p33,p12,p13,p23

! final ksi is estimated according to the yield envelope
!-------------------------------------------------
!---------calculate co-seismic velocity related to event
!--- repeat iterrations until boff < boff_min and vvv < vvv0
		vvv0 = 1.e-4
		factor = 1.e-3
	print *,' '
	print *,'	=====  START DROP EVENT ======= '
!--------------------------------------------------------
  !$OMP PARALLEL
  !$OMP DO PRIVATE(n,j)
  do n  = 1,ne
    do j  = 1,6
      stress0(j,n) = stress(j,n)
      strain0(j,n) = strain(j,n) - strainp(j,n)
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

!------- search for epicenter --------------------
!------- create a list of failed elements --------
!------- index(ndrop);  ndrop - number of failed elements ----
  ndrop = 0
  do n  = 1,ne
    if(flag(n).eq.1 .or. flag(n).eq.2 ) then
	ndrop = ndrop + 1
        el_drop(ndrop) = n
        alpha(n) = 1.0_8
  write(6,*)' epicenter: ',el_drop(ndrop),flag(n)
  write(6,*) alpha(n),phi(n),ksi(n),pf_el(n)
  
      if ( flag(n) .eq. 2 ) then		! search for target ksi
            call ksi_final(n)
      end if
    end if
  end do

!-----------------------------------------------------
!	Start main loop for co-seismic displacements
 co_seism: do 
	dynamic = .false.
! reset arrays 
    force=0.0_8
    balance=0.0_8

  !$OMP PARALLEL
  !$OMP DO PRIVATE(n,j)
    do n = 1,ne
! calculate elastic strains for all elements
       do j=1,6
         str_e(j,n) = strain(j,n) - strainp(j,n)
       end do
! update fluid pressure
       pf_el(n) = m_biot(n)*(zi_el(n) - phi(n) - al_biot(n)*i1(n))
   end do
  !$OMP END DO
  !$OMP END PARALLEL

! calculate stress in non-failed elements
 do n = 1,ne
! Stress for SOLID elements
  if ( flag(n) .eq. 0 ) then
! Check to see if anything failed, elastic returns event=.true. 
! and flag(n)=1 or 2 if the damage above critical value
      event=.false.
      call elastic(n,event,s)

     if(event)then
!  Another element failed, add to the list
 write(6,*)' Damage above critical value, element ',n,flag(n),alpha(n),phi(n),pf_el(n)
        ndrop = ndrop + 1
        el_drop(ndrop) = n
        alpha(n) = 1.0_8
	dynamic = .true.
	
! Search for target ksi
        call ksi_final(n)
	
     else
! Add fluid pressure
! Stress for next step
        do j=1,3
          stress(j,n)   = s(j) - al_biot(n)* pf_el(n)
          stress(j+3,n) = s(j+3)
      end do
     end if
  end if
 end do

!---------------------------------------------
! calculate stress in failed elements
 do nd = 1,ndrop
    n = el_drop(nd)
 
! Flow-rules for MODE-II elements

 if(flag(n).eq.2) then

! Mode-2
!-------- stress relaxation in mode-2 ---------------------
! until element is above yield 
!
! Coupling coefficient
     d_coup = coupl(n,1) + coupl(n,3)*num_drop(n) +       & 
              coupl(n,2)*phi(n) + coupl(n,3)*alpha(n)

!---use ksi(n) from previous step ---- 
		event = .true.
   if(ksi(n)*(1.0_8-d_coup*(abs(i1(n))**power)) .le. ksi0(n) ) then
! Back to solid
		event = .false.
      		call elastic(n,event,s)

   else if (event) then 
! Relaxation of elastic strain and porosity change 
                dynamic = .true.

	si2 = dsqrt(i2(n))
   fx = -d_coup*(power+1.0_8)*(abs(i1(n))**power)*si2 +        &
         d_coup*(abs(i1(n))**(power+1.0_8))*ksi(n)/3.0_8 +     &
        2./3.*i1(n)*(ksi(n)-ksi0(n)) + si2*(1.0_8-ksi(n)*ksi(n)/3.0_8)
   
    fy = (ksi(n)-2.0_8*ksi0(n)-ksi(n)*d_coup*(abs(i1(n))**power) )

		df_el  = fx/100.0_8
		phi(n) = phi(n) + df_el
!  print *,' Coup ',d_coup,si2,ksi(n),phi(n),fx
!  print *,' Poro ',n,phi(n),i1(n),ksi(n),ksif(n)
!  print *,' Factor ',factor,fx,fy,df_el,pf_el(n)
!  print *,' Sp ',stress(:,n)
!  print *,' Sn ',s
  
           do j=1,3
    strainp(j,  n) = strainp(j,  n) + df_el/3.0_8
    strainp(j,  n) = strainp(j,  n) + fy*(str_e(j,n)-i1(n)/3.0_8)/100.0_8
    strainp(j+3,n) = strainp(j+3,n) + fy*str_e(j+3,n)/100.0_8
    
    str_e(j,  n) = strain(j,  n) - strainp(j,  n)
    str_e(j+3,n) = strain(j+3,n) - strainp(j+3,n)
           end do

!       invariants of elastic strain
       call strain_invariants(n,si2)
! print *,' Update ',i1(n),phi(n),ksi(n),i2(n)-i1(n)*i1(n)/3.           

! effective moduli
    rle  = lambda(n) - gr(n)/ksif(n)
    rme  = mu(n) + gr(n)*(ksi0(n) - 0.5_8*ksif(n))

      do j=1,3
        s(j) = rle*i1(n) + 2.0_8*rme*str_e(j,n)
        s(j+3) =        2.0_8*rme*str_e(j+3,n)
      end do
   
   end if     

        s_mean = s(1)+ s(2)+ s(3)
        if(s_mean .gt. 0.0_8 ) then
		flag(n) = 1  ! mode changed 
 print *,'  Mode changed (tension): ',n,flag(n),s_mean
        end if
 end if

!-------- stress relaxation in mode-1 ---------------------
  if(flag(n).eq.1.or.flag(n).eq.3) then
! Mode-1
         s = 0.0_8
  end if

! Add fluid pressure
! Save stress for EFDLM and next step
        do j=1,3
   stress(j,n)   = factor*(s(j)-al_biot(n)*pf_el(n)-stress(j,n))+stress(j,n)
   stress(j+3,n) = factor*(s(j+3)-stress(j+3,n))+stress(j+3,n)
        end do
!  print *,ksi(n),ksif,dynamic
  end do

!-------------------------------------------
! ------ steps from EFDLM  ---------------
  do n=1,ne
! calculate and transfer element forces onto nodes
    call force_balance(n)
  end do

! calculate other forces on nodes
  call n_force_balance(dt)

! Estimation max balance-off in system (boff)
  call boff_calc(dt,vvv)
  
  if(vvv .le. 100.0_8*vvv0) then
  	factor = factor*(1.+adp)
  	 if (factor .ge. 1.0) factor = 1.0
  else
         factor = 1.0e-3_8
  end if

!------------------ Move Grid ------------------
        call move_grid (dt)

!-------- change time step ------------------------
        if ( boff .ge. boff_max ) then
                den_scale = den_scale*(1.-adp)
        if ( den_scale .lt. 1. ) den_scale = 1.
                dt = 0.25 / sqrt(vp2 / den_scale)

        else if ( boff .le. boff_min ) then
                den_scale = den_scale*(1.+adp)
                dt = 0.25 / sqrt(vp2 / den_scale)
        end if

!------------------------------------------------
!--- finish or repeate iterration --------------
  if ( .not. dynamic) then
    if ( boff .ge. boff_min ) dynamic = .true.  ! repeate until boff < boff_min
    if ( vvv  .ge. vvv0  ) dynamic = .true.  ! repeate until velocity is big
  end if

!---------------------------------------------------------
  if ( .not. dynamic) then ! test for rupture front (dynamic weakening) 
    do n = 1,ne

          if(flag(n).eq.0) then
		alp0 = alpha(n)

! calculate elastic strains
      do j=1,6
      str_e(j,n) = strain(j,n) - strainp(j,n)
      end do

!       invariants of elastic strain
       call strain_invariants(n,si2)

   if ( ksi(n) .gt. ksi0(n) ) then
!--------------Damage rate correction------------
        dadt = rate(n,1)*i2(n)*(ksi(n)-ksi0(n) )
!----- Inv-distance factor ------------------
        dsf = 0.0_8

                do mm = 1, ndrop
        dist= ( field(n,1) - field(el_drop(mm),1) )**2 +  &
              ( field(n,2) - field(el_drop(mm),2) )**2 +  &
              ( field(n,3) - field(el_drop(mm),3) )**2
	
	dsf1 = 1.0_8/dist
      if(dsf .le. dsf1) dsf = dsf1
                end do
         
        d_al = sqrt(rate(n,4)*dadt*dsf)
        alpha(n) = alp0 + d_al

          event=.false.
	call elastic(n,event,s)

        if(event)then
!  Another element failed, add to the list
   print *,'       Dynamic drop: ',n,flag(n),alpha(n),phi(n),d_al,pf_el(n)
        ndrop = ndrop + 1
        el_drop(ndrop) = n
        alpha(n) = 1.0_8
        dynamic = .true.

        else

        alpha(n) = alp0

        end if
      end if
   end if
     end do
   end if
   if ( .not. dynamic) exit co_seism
  end do co_seism

!---------------------------------------------------------
!--- Only here we really finished iterrations !!!!!!
!
!       New plastic strain and output final result

        open (11,file='catalog',position='append')

                do nd = 1,ndrop
                n = el_drop(nd)

!-------- stress drop occured in this element ------------

                nfl     = flag(n)
        if ( flag(n) .eq. 1 ) then
                flag(n) = 3

       strainp(1,n) = strain(1,n)
       strainp(2,n) = strain(2,n)
       strainp(3,n) = strain(3,n)
       strainp(4,n) = strain(4,n)
       strainp(5,n) = strain(5,n)
       strainp(6,n) = strain(6,n)
        end if


        if ( flag(n) .eq. 2 ) flag(n) = 0
        num_drop(n) = num_drop(n) + 1
     d_coup = coupl(n,1) + coupl(n,3)*num_drop(n) +       & 
              coupl(n,2)*phi(n) + coupl(n,3)*alpha(n)

	if(d_coup .le. 0.0) then
        num_drop(n) = num_drop(n) - 1
                flag(n) = 3
	end if	

!       stress drop in the element

        s11 = stress(1,n) - stress0(1,n)
        s22 = stress(2,n) - stress0(2,n)
        s33 = stress(3,n) - stress0(3,n)
        s12 = stress(4,n) - stress0(4,n)
        s13 = stress(5,n) - stress0(5,n)
        s23 = stress(6,n) - stress0(6,n)

        sdrop = dsqrt((s11**2+s22**2+s33**2)/2.0_8	&
                    + s12**2+s13**2+s23**2)*el_vol(n)

!       potency in the element

        p11 = strain(1,n) - strain0(1,n)
        p22 = strain(2,n) - strain0(2,n)
        p33 = strain(3,n) - strain0(3,n)
        p12 = strain(4,n) - strain0(4,n)
        p13 = strain(5,n) - strain0(5,n)
        p23 = strain(6,n) - strain0(6,n)

        pot = dsqrt((p11**2+p22**2+p33**2)/2.0_8   	&
                  + p12**2+p13**2+p23**2)*el_vol(n)

   print *,' drop ', nd, n, time,' mode-',nfl,alpha(n),phi(n),ksi(n),i2(n),pf_el(n)

!       center of the element
                xx = field(n,1)
                yy = field(n,2)
                zz = field(n,3)

 write (11,100) nd,n,nfl,time,xx,yy,zz,sdrop,pot,pf_el(n)

                end do

100     format(i4,i8,' mode-',i1,1x,g15.8,7(1x,g15.5))
        close(11)
  
  return
end subroutine drop

subroutine ksi_final(n)

 use sizes
 use element_data
 implicit none
 
 integer:: n
 real(kind=8):: fx,fy,d_coup,si2,i1_0,sh_0,factor,csi
 logical:: yld

   call strain_invariants(n,si2)

 if ( i1(n) .ge. 0.0 ) then
   ksif(n) = -dsqrt(3.0_8)*0.1 + ksi0(n)*0.9_8
 
 else

		factor = 1.0e-3_8
		
! Coupling coefficient
     d_coup = coupl(n,1) + coupl(n,3)*num_drop(n) +       & 
              coupl(n,2)*phi(n) + coupl(n,3)*alpha(n)

     sh_0 = dsqrt(i2(n) - i1(n)*i1(n)/3.0_8)
     i1_0 = i1(n)

    if(ksi(n)*(1.0_8-d_coup*(abs(i1(n))**power)) .ge. ksi0(n) ) yld = .true.
     
     if(yld) print *,'  Elelement ',n,' is above yield ',i1(n),ksi(n)

	do while( yld )
 
   fx = -d_coup*(power+1.0_8)*(abs(i1(n))**power)*si2 +        &
         d_coup*(abs(i1(n))**(power+1.0_8))*ksi(n)/3.0_8 +     &
        2./3.*i1(n)*(ksi(n)-ksi0(n)) + si2*(1.0_8-ksi(n)*ksi(n)/3.0_8)
   
    fy = sh_0 * (ksi(n)-2.0_8*ksi0(n)-ksi(n)*d_coup*(abs(i1(n))**power) )

 	i1_0 = i1_0 - fx*factor
 	sh_0 = sh_0 - fy*factor
 	csi = i1_0 / dsqrt(i1_0*i1_0/3.0_8 + sh_0*sh_0)
! print *,yld,i1_0,sh_0,fx,fy,csi,d_coup
! print *,d_coup,power,i1(n),si2,ksi(n),fx
 

    if(csi*(1.0_8-d_coup*(abs(i1_0)**power)) .le. ksi0(n) ) yld = .false.
	
	end do
	
	ksif(n) = csi
  print *,' Target ksi-value ',ksif(n),' expected porosity change ',i1_0-i1(n)

 end if        
 return
end subroutine ksi_final

