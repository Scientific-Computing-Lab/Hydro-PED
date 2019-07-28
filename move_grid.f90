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

!  Move grid along veocity field (Lagrangian) and new total strain tensor  
!> @brief         Moves nodes and calculates strains
!> @details       Uses explicit time step.
!> Calculates new node coordinates, cumulative node displacements, 
!> and incremental displacements
!> Calculates strain rate tensor and integrate it for new total strain
!> \f$\Delta\varepsilon_{ij}=\Delta t \cdot \sum\limits_{k=1}^4{\left(V_I^{(k)}\frac{\partial L_k}{\partial x_j}+V_j^{(k)}\frac{\partial L_k}{\partial x_i}\right)}\f$
!>
!> @param[in]     time step  
!> @param[out]    None
subroutine move_grid(dt)
  
  use sizes
  use node_data
  use element_data

  implicit none
  
  integer:: ii,j,i,n,nn
  real(kind=8)::dt,de(6),dr(4,3)
  
  ! moving of mesh
  !$OMP PARALLEL
  !$OMP DO PRIVATE(i,j)
  do i=1,np
    do j=1,3
      cord(i,j) = cord(i,j) + vel(i,j)*dt
      disp(i,j) = disp(i,j) + vel(i,j)*dt
      dspt(i,j) = dspt(i,j) + vel(i,j)*dt
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL 

  ! New elastic strain tensor
  !$OMP PARALLEL
  !$OMP DO PRIVATE(dr,n,ii,de,j,nn)
  do n = 1,ne
 ! derivation of basic functions
    call derivation(dr,n)
    
    ! Strain rate
    de = 0.0_8
    do ii=1,4
      nn = nop(n,ii)
      de(1)=de(1) + vel(nn,1)*dr(ii,1)
      de(2)=de(2) + vel(nn,2)*dr(ii,2)
      de(3)=de(3) + vel(nn,3)*dr(ii,3)
      de(4)=de(4) + 0.5_8*(vel(nn,1)*dr(ii,2)+vel(nn,2)*dr(ii,1))
      de(5)=de(5) + 0.5_8*(vel(nn,1)*dr(ii,3)+vel(nn,3)*dr(ii,1))
      de(6)=de(6) + 0.5_8*(vel(nn,2)*dr(ii,3)+vel(nn,3)*dr(ii,2))
    end do

 ! New total strains
    do j=1,6
      strain(j,n) = strain(j,n) + de(j)*dt
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  
return
end subroutine move_grid
