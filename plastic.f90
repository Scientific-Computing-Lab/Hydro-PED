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
!    * Neither the name of Harel Levin or Gal Oren, nor the
!      names of its contributors may be used to endorse or promote products
!      derived from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
! ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
! WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL Harel Levin & Gal Oren BE LIABLE FOR ANY
! DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
! (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
! ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

! Plastic strain accumulation
! Effective viscosity proportional to dalpha 
subroutine plastic(dt)
  use sizes
  use element_data
  implicit none
  
  integer n,j
  real(kind=8):: dt,t(6),s_mean,vis
  
  !$OMP PARALLEL
  !$OMP DO PRIVATE(n,s_mean,t,j,vis)
  do n = 1,ne

	vis = 0.0_8
   if( flag(n).eq.0 .and. dalpha(n) .gt. 0.0_8 )then
        vis = ductile(n,1)*dalpha(n)
   end if

! Deviatoric stress
      s_mean = (stress(1,n) + stress(2,n) + stress(3,n))/3.0_8

      do j=1,3
        t(j)   = stress(j,n) - s_mean
        t(j+3) = stress(j+3,n)
      end do

! Damage-related viscosity
! New plastic strain
     do j=1,6
   strainp(j,n)=strainp(j,n) + t(j)*(vis + ductile(n,2))*dt
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  return
end subroutine plastic
