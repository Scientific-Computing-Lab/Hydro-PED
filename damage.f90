! Copyright (c) 2017, 
! Eyal Shalev (eyal@gsi.gov.il)
! Vladimir Lyakhovsky
! Harel Levin (harellevin@gmail.com)
! Gal Oren (galoren.com@gmail.com)
! All rights reserved to:
! Geological Survey of Israel (GSI) &
! Nuclear Research Center - Negev (NRCN).
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!    * Redistributions of source code must retain the above copyright
!      notice, this list of conditions and the following disclaimer.
!    * Redistributions in binary form must reproduce the above copyright
!      notice, this list of conditions and the following disclaimer in the
!      documentation and/or other materials provided with the distribution.
!    * Neither the name of Eyal Shalev, Vladimir Lyakhovsky, Harel Levin or 
!      Gal Oren, nor the names of its contributors may be used to endorse 
!      or promote products derived from this software without specific prior 
!      written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND
! ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
! WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL Eyal Shalev, Vladimir Lyakhovsky, Harel Levin 
! & Gal Oren BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
! OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
! WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
! OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
! ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#define ALLOC   alloc_if(.true.)
#define FREE    free_if(.true.)
#define RETAIN  free_if(.false.)
#define REUSE   alloc_if(.false.)

! Damage evolution
subroutine damage(dt,dmax)
  use sizes
  use element_data
  implicit none

  integer:: n
  real(kind=8)::dt,dmax,chs
  real(kind=8)::r1,r2,a1,a2,vol,dlph

  	chs = 2.0e-7

  r1 = 999999.0_8
  r2 = -99999.0_8
  a1 = 999999.0_8
  a2 = -99990.0_8
  
  !dir$ offload_transfer target(mic:0) wait(strain_signal) nocopy(strain(1:6,1:phi_length) : REUSE RETAIN)
  !dir$ offload_transfer target(mic:1) wait(strain_signal) nocopy(strain(1:6,phi_length+1:2*phi_length) : REUSE RETAIN)
!  !$OMP PARALLEL
!  !$OMP DO PRIVATE(n,vol)
!          do n = 1,ne

!        if ( flag(n).eq.3 ) then
!----- closure of the previously open space         -------
!----- if volume is less than critical (-1.e-3 here) -------
!
!        vol = (strain(1,n) + strain(2,n) + strain(3,n))/3.
!
!        strainp(1,n) = strain(1,n) - vol
!        strainp(2,n) = strain(2,n) - vol
!        strainp(3,n) = strain(3,n) - vol
!        strainp(4,n) = strain(4,n)
!        strainp(5,n) = strain(5,n)
!        strainp(6,n) = strain(6,n)
!
!        if ( vol .lt. -1.e-3 ) then
!                flag(n) = 0
!                ksi(n) = -sqrt(3.0_8)
!                i2(n) = 3.0_8*vol*vol
!        print *,' Element',n,' closed, volume =  ',vol,ksi(n),i2(n)
!                end if
!
!        end if
!        end do
!
!  !$OMP END DO
!  !$OMP END PARALLEL
  
    dalpha = 0.0_8

  !$OMP PARALLEL
  !$OMP DO PRIVATE(n,dlph)
   do n = 1,ne
   
   if ( flag(n).eq.0 ) then
!    
!---------  damage evolution -------------------
!
        dalpha(n) = i2(n)*(ksi(n)-ksi0(n)) - chs

    if ( dalpha(n) .ge. 0.0_8 ) then    ! damage increase
           dalpha(n) = dalpha(n) * rate(n,1)
	   alpha(n) = alpha(n) + dalpha(n)*dt
	if(alpha(n) .ge. 1.0_8 ) alpha(n) = 1.0_8 ! alpha <= 1

      else ! damage decrease
  dlph = -rate(n,3)*log(1.0_8-rate(n,2)/rate(n,3)*exp(alpha(n)/rate(n,3))*dalpha(n)*dt)
	    alpha(n) = alpha(n) + dlph
	    dalpha(n) = dlph/dt
      end if

  end if
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  
  do n=1,ne
        if(dalpha(n) .le. a1  ) a1 = dalpha(n)
	if(dalpha(n) .ge. a2  ) a2 = dalpha(n)
 ! save alpha min-max 
	if(alpha(n) .le. 0.0_8 ) alpha(n) = 0.0_8 ! alpha >= 0 
	if(alpha(n) .le. r1 .and. flag(n) .eq. 0 ) r1 = alpha(n)
	if(alpha(n) .ge. r2 .and. flag(n) .eq. 0 ) r2 = alpha(n)

  end do

    	dmax = a2*dt
 
  write(6,*)' 	min - max rate  : ',a1,a2
  write(6,*)' 	min - max value : ',r1,r2
      
   return   
end subroutine damage
