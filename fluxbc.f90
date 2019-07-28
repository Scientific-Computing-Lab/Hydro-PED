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

subroutine fluxbc(time)
  use sizes
!  use element_data
  use node_data
  use diffusion_data
  implicit none

  integer:: i,j,n
  
  real(kind=8):: time,fluxsum,ratiop,ac,wflux

	wflux = 120.0_8	!50.0_8	! 30.0_8
	ac = 0.001_8
		fluxsum=0.
	do i=1,nfs
	n = nface(i)
  fluxsum=fluxsum+q(1,n)*vector(i,1)+q(2,n)*vector(i,2)+q(3,n)*vector(i,3)
	end do

  if (time .gt. 110000.0_8 .and. time .lt. 220000.0_8) then
	over_pres = 0.99*over_pres

  elseif (time .gt. 220000.0_8 .and. time .lt. 230000.0_8) then
	over_pres = 2499.*time - 549770000.
  else

  if (over_pres .lt. 100000.0_8) over_pres =	100000.0_8
  if (over_pres .lt.  10.e+6 )   over_pres = 1.1_8*over_pres
	ratiop=ac*(wflux-fluxsum)/wflux
	over_pres=over_pres*(1+ratiop)

  end if

   do j = 1,njpn
	bc_dfval( njp(j) )=-1000.*9.81*cord(njp(j),3 )+over_pres
   end do

  write(30,*) time,fluxsum,over_pres
  print *,'   Flux  ',time, fluxsum,over_pres
          
 return
end subroutine fluxbc
