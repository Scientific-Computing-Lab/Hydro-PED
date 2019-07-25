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

!> @brief         Calculates volumes and masses
!> @details       Calculates volumes for all tetra elements 
!> and node masses (1/4 of the element mass is projected to every node)
!> screen output: 
!> 	   Average element volume
!>        Volume ranges: min - max value
!> <br>\f$
!! V=\frac{1}{6}\text{det} 
!! \begin{vmatrix} 
!! 1 & x_i & y_i & z_i \\ 
!! 1 & x_j & y_j & z_j \\ 
!! 1 & x_m & y_m & z_m \\ 
!! 1 & x_n & y_n & z_n 
!! \end{vmatrix}\f$
!> <br>\f$m=\rho \cdot V\f$
!>
!> @param[in]     None  
!> @param[out]    None
subroutine node_mass
  ! calculate element volumes and masses
  use sizes
  use node_data
  use element_data
  implicit none

  integer::num,i
  real(kind=8)::xn(4),yn(4),zn(4)
  real(kind=8)::vol,vola,volmin,volmax,elmass
  real(kind=8)::v1,v2,v3,v4

  mass=0.0_8 !zero mass array

  vola = 0.0_8
  volmin = 999999999999.0_8
  volmax = 0.0_8

  ! loop for all elements
  do num = 1,ne
    ! coord 4 element nodes
    do i = 1,4
      xn(i) = cord(nop(num,i),1)
      yn(i) = cord(nop(num,i),2)
      zn(i) = cord(nop(num,i),3)
    end do

    !element volume  (G = 6 * volume )
    !| 1 x1 y1 z1 |
    !| 1 x2 y2 z2 |
    !| 1 x3 y3 z3 |
    !| 1 x4 y4 z4 |
    v1= xn(2) * (yn(3)*zn(4)-zn(3)*yn(4)) &
      - xn(3) * (yn(2)*zn(4)-zn(2)*yn(4)) &
      + xn(4) * (yn(2)*zn(3)-zn(2)*yn(3)) 

    v2= xn(1) * (yn(3)*zn(4)-zn(3)*yn(4)) &
      - xn(3) * (yn(1)*zn(4)-zn(1)*yn(4)) &
      + xn(4) * (yn(1)*zn(3)-zn(1)*yn(3))

    v3= xn(1) * (yn(2)*zn(4)-zn(2)*yn(4)) &
      - xn(2) * (yn(1)*zn(4)-zn(1)*yn(4)) &
      + xn(4) * (yn(1)*zn(2)-zn(1)*yn(2))

    v4= xn(1) * (yn(2)*zn(3)-zn(2)*yn(3)) &
      - xn(2) * (yn(1)*zn(3)-zn(1)*yn(3)) &
      + xn(3) * (yn(1)*zn(2)-zn(1)*yn(2))

    vol = v1 - v2 + v3 - v4
    vol = abs(vol)/6.0_8

    vola = vola + vol

    if ( volmin .ge. vol ) volmin = vol
    if ( volmax .le. vol ) volmax = vol

    el_vol(num) = vol
    elmass = dens(num) * vol
    do i = 1,4
      mass(nop(num,i)) = mass(nop(num,i)) + 0.25_8*elmass
    end do
  end do

  vola = vola / ne
  write(6,*)' Average element volume is ',vola
  write(6,*)' Volume ranges from ',volmin,'  to ',volmax

  return
end subroutine node_mass
