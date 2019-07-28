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
! DISCLAIMED. IN NO EVENT SHALL Eyal Shalev @ Vladimir Lyakhovsky
! BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
! OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND  ON ANY THEORY OF LIABILITY, 
! WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
! OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
! ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

!> @brief         Calculates damage evolution
!> @details       Uses explicit time step
!> screen output: 
!> 	min - max rate of damage evolution
!>        min - max value of the damage
!> \f[
!> \frac{d\alpha_D}{dt}=
!> \begin{cases}
!> C_dI_2\left(\xi-\xi_0\right)	                                        & \quad \text{for } \xi>\xi_0 \\
!> C_1\exp\left(\frac{\alpha_D}{C_2}\right)I_2\left(\xi-\xi_0\right)	& \quad \text{for } \xi<\xi_0
!> \end{cases}
!> \f]
!>
!> @param[in]     time step  
!> @param[out]    maximal damage rate accumulation
!> @param[out]    maximal Phi
subroutine damage(dt,damax,dfmax)
  use sizes
  use element_data
  implicit none

  integer:: n,j
  real(kind=8)::dt,damax,dfmax,dfi1,dfi2,fi_eq,fi_n
  real(kind=8)::r1,r2,a1,a2,vol,dlph
  real(kind=8)::f1,f2,rf1,rf2,peff,d_coup

  r1  = 999999.0_8
  r2  = -99999.0_8
  a1  = 999999.0_8
  a2  = -99990.0_8
  f1  = 999999.0_8
  f2  = -99999.0_8
  rf1 = 999999.0_8
  rf2 = -99999.0_8
  
!$OMP PARALLEL
!$OMP DO PRIVATE(n,vol)
  do n = 1,ne

  if ( flag(n).eq.3 ) then
!----- closure of the previously open space         -------
!----- if volume is less than critical (-1.e-1 here) -------

     vol = (strain(1,n) + strain(2,n) + strain(3,n))/3.

     strainp(1,n) = strain(1,n) - vol
     strainp(2,n) = strain(2,n) - vol
     strainp(3,n) = strain(3,n) - vol
     strainp(4,n) = strain(4,n)
     strainp(5,n) = strain(5,n)
     strainp(6,n) = strain(6,n)

    if ( vol .lt. -1.e-1 ) then
              flag(n) = 0
              alpha(n) = 0.0_8
              ksi(n) = -sqrt(3.0_8)
              i1(n) = 3.0_8*vol
              i2(n) = 3.0_8*vol*vol
!------ stress in closed element -------
       do j = 1,3
         stress(j,n) = (lambda(n)-gr(n)/ksi(n))*i1(n) +            &
              (2.0_8*(mu(n)+ksi0(n)*gr(n))- gr(n)*ksi(n))*vol
         stress(j+3,n) = 0.0_8
       end do
 
     print *,' Element',n,' closed, volume =  ',vol,ksi(n),i2(n)
    end if

  end if
  end do

!$OMP END DO
!$OMP END PARALLEL
  
    dalpha = 0.0_8
    dphi   = 0.0_8

!$OMP PARALLEL
!$OMP DO PRIVATE(n,dlph,peff,d_coup,fi_eq,fi_n,dfi1,dfi2)
 do n = 1,ne
   
 if ( flag(n).eq.0 ) then
 
  if ( i1(n) .ge. 0.0 ) then 
!------- damage increase & zero porosity rate ---------------------
      dalpha(n) = rate(n,1) * i2(n)*(ksi(n)-ksi0(n))
      alpha(n)  = alpha(n) + dalpha(n)*dt

  else
   
!  Effective pressure
  peff = -(stress(1,n)+stress(2,n)+stress(3,n))/3.
	if( peff .le. 1. ) peff = 1.
!    
!---------  damage evolution -------------------
!
! Coupling coefficient
     d_coup = coupl(n,1) + coupl(n,3)*num_drop(n) +           & 
              coupl(n,2)*phi(n) + coupl(n,3)*alpha(n)
            
     dalpha(n) = d_coup * ((-i1(n))**(power+1)) * dsqrt(i2(n)) + &
                i2(n)*(ksi(n)-ksi0(n))
                
!   if(n.eq.89298) print *,' Coupl 89298 ',d_coup,alpha(n),phi(n)
    
    if( dalpha(n) .ge. 0.) then 	! damage increase 

!---------  damage increase ----------------------------
	dalpha(n) = rate(n,1) * dalpha(n)
	alpha(n)  = alpha(n) + dalpha(n)*dt

    else
!---------  damage decrease ----------------------------

	dlph = -rate(n,3)*log(1.-rate(n,2)/rate(n,3)*      &
     	        exp(alpha(n)/rate(n,3))*dalpha(n)*dt)
	
	alpha(n)  = alpha(n) + dlph
	dalpha(n) = dlph/dt

!------- alpha >= 0 --------------------------------------
	if(alpha(n) .le. 0. ) alpha(n) = 0.
    endif
!---- end of damage evolution ----------------
!---------------------------------------------
!---------  porosity kinetics ----------------
!--------- pressure-driven compaction ------------------
!------porosity decrease (first term = Athy law) --------------

                fi_eq = phi_eq(n,1)+phi_eq(n,2)*exp(-peff/phi_eq(n,3))
                
                if (phi(n).lt.fi_eq) then
                    dfi1 = 0.0_8
                else
            fi_n = (phi(n)-fi_eq)*exp(-dt*peff*rate(n,5)) + fi_eq
                    dfi1 = fi_n - phi(n)
                end if 

!---------damage-related porosity change (second term) --------------
		dfi2 = 0.
     if ( ksi(n) .ge. ksi0(n) ) then
!---- Regime III - above onset of dilation -------------------
	dfi2 = gr(n)/peff*i2(n)*(ksi(n)-ksi0(n))*dalpha(n)*dt

!----- porosity growth may be reduced by any factor1<1 ----------	
!		dfi2 = dfi2 * factor1	
!   if(n.eq.89298) print *,' 89298 Regime III ',dfi2

     else if ( dalpha(n) .gt. 0. ) then
!---- Regime II - above yield and below onset of dilation ------------
	dfi2 = rate(n,1)*d_coup* ( (-i1(n))**(power+1.0_8) ) *   &
	       dsqrt(i2(n)) * gr(n)/peff *i2(n)*(ksi(n)-ksi0(n)) * dt

!----- compaction rate may be increased by any factor2>1 ----------	
!		dfi2 = dfi2 * factor2	
!   if(n.eq.89298) print *,' 89298 Regime II ',dfi2

     else
!---- Regime I - below yield ------------
	dfi2 = rate(n,2)*exp(alpha(n)/rate(n,3))*d_coup*               &
	       ( (-i1(n))**(power+1.0_8) ) *                           &
	       dsqrt(i2(n)) * gr(n)/peff *i2(n)*(ksi(n)-ksi0(n)) * dt

!----- compaction rate may be increased or decreased by any factor3 ----------	
		dfi2 = 0.       !dfi2 * factor3

!   if(n.eq.89298) print *,' 89298 Regime I ',dfi2

     end if
	
!---------- new porosity and incremental change -------------- 
	 phi(n) = phi(n) + dfi1 + dfi2
	 dphi(n) = dfi1 + dfi2

       if ( phi(n) .le. 0.01 ) then
                phi(n)  = 0.01_8
                dphi(n) = 0.0_8
       end if

!------ change plastic strain ------------------------
	strainp(1,n) = strainp(1,n) + dphi(n) /3.0_8
	strainp(2,n) = strainp(2,n) + dphi(n) /3.0_8
	strainp(3,n) = strainp(3,n) + dphi(n) /3.0_8

  end if
 end if
 end do
!$OMP END DO
!$OMP END PARALLEL
  
  do n = 1,ne
!------- Search for Min-Max -----------------------
	if(dalpha(n) .le. a1  ) a1 = dalpha(n)
	if(dalpha(n) .ge. a2  ) a2 = dalpha(n)

!------ save alpha min-max ----------------------------
	if(alpha(n) .le. r1  ) r1 = alpha(n)
	if(alpha(n) .ge. r2  ) r2 = alpha(n)

!------ save phi min-max ----------------------------
	if(phi(n) .le. rf1  ) rf1 = phi(n)
	if(phi(n) .ge. rf2  ) rf2 = phi(n)

	if( dphi(n) .le. f1  ) f1 = dphi(n)
	if( dphi(n) .ge. f2  ) f2 =  dphi(n)
  end do

	print *,' d alpha min - max: ',a1,a2,' d phi min - max:',f1/dt,f2/dt
	print *,'   alpha min - max: ',r1,r2,'   phi min - max:',rf1,rf2

		a1 = abs(a1) * dt
		a2 = abs(a2) * dt
		
		f1 = abs(f1)
		f2 = abs(f2)

		damax = dmax1(a1,a2)
		dfmax = dfmax + dmax1(f1,f2)

   return   
end subroutine damage
