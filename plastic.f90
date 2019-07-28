! Copyright (c) 2017, 
! Eyal Shalev (eyal@gsi.gov.il)
! Vladimir Lyakhovsky
! All rights reserved to:
! Geological Survey of Israel (GSI) &
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
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND
! ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
! WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL Eyal Shalev, Vladimir Lyakhovsky
! BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
! OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
! WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
! OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
! ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

! Plastic strain accumulation
! Effective viscosity proportional to dalpha 
subroutine plastic(dt)
  use sizes
  use element_data
  implicit none
  
  integer n,j
  real(kind=8):: dt,t(6),s_mean
  
  !$OMP PARALLEL
  !$OMP DO PRIVATE(n,s_mean,t,j)
  do n = 1,ne

   if( flag(n).eq.0 )then 
    if( dalpha(n) .gt. 0.0_8 )then

! Deviatoric stress
      s_mean = (stress(1,n) + stress(2,n) + stress(3,n))/3.0_8
       do j=4,6
        t(j) = stress(j,n)
       end do

      do j=1,3
        t(j) = stress(j,n)-s_mean
      end do

! Damage-related viscosity
! New plastic strain
     do j=1,6
   strainp(j,n)=strainp(j,n) + t(j)*ductile(n)*dalpha(n)*dt
     end do
    end if
   end if
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  
!  vsmax=0.0_8
!  do n=1,ne
!    vsmax=max(vsmax,vs1(n))
!  end do
!
!  dtp = 1.0e-6_8 / vsmax
  
  return
end subroutine plastic
