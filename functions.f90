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

! **************************************************
! minmax
 subroutine minmax(a,n,amin,amax)
  implicit none
  integer:: n,i
  real(kind=8):: a(n),amax,amin

  amin = a(1)
  amax = a(1)
  do i =2,n
    if( a(i) .lt. amin ) amin = a(i) 
    if( a(i) .gt. amax ) amax = a(i)
  end do

  return
 end subroutine minmax
 
 ! **************************************************
 ! Determination abs max velocity
 function  vmax(np,r1,nmax) 
  implicit none
  integer:: np,nmax,nnp
  real(kind=8)::r1(np,3),vm,rr,vmax
    
  vm=0.0_8
  do nnp =1,np
    rr = r1(nnp,1)**2 + r1(nnp,2)**2 + r1(nnp,3)**2
    if(rr.gt.vm) then
      vm = rr
      nmax = nnp
    end if
  end do
  vmax = sqrt ( vm )
  
  return
 end function vmax
 
