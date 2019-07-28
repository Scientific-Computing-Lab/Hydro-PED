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

! Elastic stress-strain test for connectivity
!   flg = 0    -  stable
!   flg = 1    -  1-st type of failure
!   flg = 2    -  2-nd type of failure
!   flg = 3    -  open space

subroutine elastic(n,event,s)
use sizes
use element_data
implicit none

integer:: n,j
logical:: event
real(kind=8)rl,rm,rg,s(6)
real(kind=8)i1,si2,a,b,xmin,p,q,d,cs

 if( flag(n).eq.3 ) then    ! open space 
    do j=1,6
      s(j) = 0.0_8
    end do
 else ! normal cell
    ! invariants of the strain tensor
    call strain_invariants(n,i1,si2)

    ! stress calculation, zero if loss convexity
    ! effective moduli
    rl  = lambda(n)
    rm  = mu(n) + alpha(n)*ksi0(n)*gr(n)
    rg  = alpha(n)*gr(n)
 
   cs = ksi(n)
   
    a =2.0_8*rm-rg*ksi(n)
    p =-(4.0_8*rm+3.0_8*rl-3.0_8*rg*cs)
    q =a**2 + a*(3.0_8*rl-rg*cs) + (rl*rg*cs-rg*rg)*(3.0_8-cs*cs)
    d = 0.25_8*p**2 - q
       
    if ( d .le. 0.0_8 ) then 
      write(6,*)' Discreminant < 0; I do not know how to continue ',d
      stop
    end if
    xmin = -0.5_8*p - sqrt(d)
    if ( xmin .le. 0.0_8 ) then
      flag(n) = 2
      event=.true.
      do j=1,6
        s(j) = 0.0_8
      end do
    else if ( a .le. 0.0_8 ) then
      flag(n) = 1
      event=.true.
      do j=1,6
        s(j) = 0.0_8
      end do
    else
      b=rl*i1 - rg*si2
      do j=1,3
        s(j) = b + a*str_e(j,n)
        s(j+3) =   a*str_e(j+3,n)
      end do
    end if
 end if

  return 
end subroutine elastic
!********************************************************************
subroutine strain_invariants(n,i1,si2)
  use sizes
  use element_data
  implicit none

  integer::n
  real(kind=8)::i1,si2
   
  ! invariants of the strain tensor
  i1    = str_e(1,n) + str_e(2,n) + str_e(3,n)
  i2(n) =      str_e(1,n)**2 + str_e(2,n)**2 + str_e(3,n)**2  &
    +2.0_8*(str_e(4,n)**2 + str_e(5,n)**2 + str_e(6,n)**2)

  si2 = sqrt(i2(n))
  if(si2 .le. 1.0e-10_8 ) then
    ksi(n) = 0.0_8
  else
    ksi(n) = i1/si2
  end if
  
  return
end subroutine strain_invariants
