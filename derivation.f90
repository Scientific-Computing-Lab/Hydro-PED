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
!    * Neither the name of Harel Levin or Gal Oren, nor the
!      names of its contributors may be used to endorse or promote products
!      derived from this software without specific prior written permission.
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

!--------------- derivations of basic functions --------------
!   element number num  point coordinate inside the element
! dr(k,i);  k=node number (1,2,3,4);  i=coordinate(x,y,z)
!--------------------------------------------------------

subroutine derivation(dr,n)
!dir$ attributes forceinline :: move_grid_pf_inline
!dir$ attributes offload : mic :: derivation
!$omp declare simd(derivation)
  use sizes
  use node_data
  use element_data
  implicit none

  integer, intent(in)::n
  integer::i
  real(kind=8):: xn(4),yn(4),zn(4)
  real(kind=8), intent(out):: dr(4,3)
  real(kind=8):: a,b,c,d,e,vol,v1,v2,v3,v4,onv

  ! coord 4 element nodes
  do i = 1,4
    xn(i) = cord(nop(n,i),1)
    yn(i) = cord(nop(n,i),2)
    zn(i) = cord(nop(n,i),3)
  end do

  ! element volume  (G = 6 * volume )
  !| 1 x1 y1 z1 |
  !| 1 x2 y2 z2 |
  !| 1 x3 y3 z3 |
  !| 1 x4 y4 z4 |

  v1 = xn(2) * (yn(3)*zn(4)-zn(3)*yn(4)) &
     - xn(3) * (yn(2)*zn(4)-zn(2)*yn(4)) &
     + xn(4) * (yn(2)*zn(3)-zn(2)*yn(3))

  v2 = xn(1) * (yn(3)*zn(4)-zn(3)*yn(4)) &
     - xn(3) * (yn(1)*zn(4)-zn(1)*yn(4)) &
     + xn(4) * (yn(1)*zn(3)-zn(1)*yn(3))

  v3 = xn(1) * (yn(2)*zn(4)-zn(2)*yn(4)) &
     - xn(2) * (yn(1)*zn(4)-zn(1)*yn(4)) &
     + xn(4) * (yn(1)*zn(2)-zn(1)*yn(2))

  v4 = xn(1) * (yn(2)*zn(3)-zn(2)*yn(3)) &
     - xn(2) * (yn(1)*zn(3)-zn(1)*yn(3)) &
     + xn(3) * (yn(1)*zn(2)-zn(1)*yn(2))

  vol = v1 - v2 + v3 - v4

  if(vol .ge.0.0_8 )then
!    write(6,*)'	Negative volume ',-vol/6,' element ', n
!    write(6,*)-v1,v2,-v3,v4
!    stop
                flag(n) = 3
                dr = 0.
                return
  end if
  onv= 1.0_8/vol

  !A = | 1 y2 z2 |
  !    | 1 y3 z3 |
  !    | 1 y4 z4 |
  a = yn(3)*zn(4)-zn(3)*yn(4) &
    -(yn(2)*zn(4)-zn(2)*yn(4)) &
    + yn(2)*zn(3)-zn(2)*yn(3)

  !B = | 1 x2 z2 |
  !    | 1 x3 z3 |
  !    | 1 x4 z4 |
  b = xn(3)*zn(4)-zn(3)*xn(4)  &
    -(xn(2)*zn(4)-zn(2)*xn(4)) &
    + xn(2)*zn(3)-zn(2)*xn(3)

  !C = | 1 x2 y2 |
  !    | 1 x3 y3 |
  !    | 1 x4 y4 |
  c = xn(3)*yn(4)-yn(3)*xn(4)  &
    -(xn(2)*yn(4)-yn(2)*xn(4)) &
    + xn(2)*yn(3)-yn(2)*xn(3)

  dr(1,1) = -a * onv
  dr(1,2) =  b * onv
  dr(1,3) = -c * onv

  !A = | 1 y1 z1 |
  !    | 1 y3 z3 |
  !    | 1 y4 z4 |
  a = yn(3)*zn(4)-zn(3)*yn(4)  &
    -(yn(1)*zn(4)-zn(1)*yn(4)) &
    + yn(1)*zn(3)-zn(1)*yn(3)

  !B = | 1 x1 z1 |
  !    | 1 x3 z3 |
  !    | 1 x4 z4 |
  b = xn(3)*zn(4)-zn(3)*xn(4)  &
    -(xn(1)*zn(4)-zn(1)*xn(4)) &
    + xn(1)*zn(3)-zn(1)*xn(3)

  !C = | 1 x1 y1 |
  !    | 1 x3 y3 |
  !    | 1 x4 y4 |
  c = xn(3)*yn(4)-yn(3)*xn(4)  &
    -(xn(1)*yn(4)-yn(1)*xn(4)) &
    + xn(1)*yn(3)-yn(1)*xn(3)

  dr(2,1) =  a * onv
  dr(2,2) = -b * onv
  dr(2,3) =  c * onv

  !A = | 1 y1 z1 |
  !    | 1 y2 z2 |
  !    | 1 y4 z4 |
  a = yn(2)*zn(4)-zn(2)*yn(4)  &
    -(yn(1)*zn(4)-zn(1)*yn(4)) &
    + yn(1)*zn(2)-zn(1)*yn(2)

  !B = | 1 x1 z1 |
  !    | 1 x2 z2 |
  !    | 1 x4 z4 |
  b = xn(2)*zn(4)-zn(2)*xn(4)  &
    -(xn(1)*zn(4)-zn(1)*xn(4)) &
    + xn(1)*zn(2)-zn(1)*xn(2)

  !C = | 1 x1 y1 |
  !    | 1 x2 y2 |
  !    | 1 x4 y4 |
  c = xn(2)*yn(4)-yn(2)*xn(4)  &
    -(xn(1)*yn(4)-yn(1)*xn(4)) &
    + xn(1)*yn(2)-yn(1)*xn(2)

  dr(3,1) = -a * onv
  dr(3,2) =  b * onv
  dr(3,3) = -c * onv

  !A = | 1 y1 z1 |
  !    | 1 y2 z2 |
  !    | 1 y3 z3 |
  a = yn(2)*zn(3)-zn(2)*yn(3)  &
    -(yn(1)*zn(3)-zn(1)*yn(3)) &
    + yn(1)*zn(2)-zn(1)*yn(2)

  !B = | 1 x1 z1 |
  !    | 1 x2 z2 |
  !    | 1 x3 z3 |
  b = xn(2)*zn(3)-zn(2)*xn(3)  &
    -(xn(1)*zn(3)-zn(1)*xn(3)) &
    + xn(1)*zn(2)-zn(1)*xn(2)

  !C = | 1 x1 y1 |
  !    | 1 x2 y2 |
  !    | 1 x3 y3 |
  c = xn(2)*yn(3)-yn(2)*xn(3)  &
    -(xn(1)*yn(3)-yn(1)*xn(3)) &
    + xn(1)*yn(2)-yn(1)*xn(2)

  dr(4,1) =  a * onv
  dr(4,2) = -b * onv
  dr(4,3) =  c * onv

  return
end subroutine derivation
