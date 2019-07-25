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

!> @brief         Calculates node velocity
!> @details       Uses Lagrangian method for force balance
!> calculations. 
!> <br>\f$F_i=\sum\limits_{faces}{\frac{1}{3}\sigma_{ij}n_j+mg_i}\f$ 
!> <br>\f$\frac{\partial{v_i}}{\partial{t}}=\frac{F_i}{m}\f$ 
!> <br>\f$V_i^{(n)}(t+\Delta t)=V_i^{(n)}(t)+\left[F_i^{(n)}-\chi\left|F_i^{(n)}\right|{sign}\left(V_i^{(n)}\right)\right]\frac{\Delta t}{m_{inert}}\f$
!> @param[in]     time step  
!> @param[out]    event
subroutine efdlm(dt,event)
  use sizes
  use node_data
  use element_data
  implicit none

  logical::event
  integer::n,j,i
  real(kind=8)::dt,vvv,s(6)
  
! Initialization 
  event=.false.
  force=0.0_8
  balance=0.0_8
 
! loop through all elements
  !$OMP PARALLEL
  !$OMP DO PRIVATE(n,s,j)
  do n=1,ne
!calculate elastic strains
    do j=1,6
      str_e(j,n) = strain(j,n) - strainp(j,n)
    end do

! elastic returns event=.true. if the damage above critical value on a new node, calculates stresses
    call elastic(n,event,s)
    if ( alpha(n)-1.0_8 .gt. 1.0e-5_8 ) then
    		event=.true.
    		flag(n) = 2
    end if

! Add fluid pressure
    do j=1,3
     s(j) = s(j) - al_biot(n) * pf_el(n)
    end do

    if ( flag(n) .eq. 0 .or. flag(n) .eq. 3 ) then
! Stress for next step
     do j=1,6
      stress(j,n)=s(j)
    end do
   end if

  end do
  !$OMP END DO
  !$OMP END PARALLEL
  
  if(event) return ! onset of seismic event, move to drop event

  do n=1,ne    
! calculate and transfer element forces onto nodes
    call force_balance(n)
  end do

! calculate other forces on nodes
  call n_force_balance(dt)

! Estimation max balance-off in system (boff)
  call boff_calc(dt,vvv)
  
  return
end subroutine efdlm

!*****************************************************
subroutine force_balance(n)
  use sizes
  use node_data
  use element_data
  implicit none

  logical::event,space
  integer::n,nside,nn,j,ii,n1,n2,n3
  real(kind=8)::x1,x2,y1,y2,z1,z2,xnorm,ynorm,znorm
  real(kind=8)::fs(4,3),f(4,3)

  ! Calculation forces (4 sides, one element)
  ! On the sides (fs(i,j), i-side, j-dimension) fi = sij*nj
  do nside = 1,4
    ! Three nodes of the side
    n1 = nop(n,side(1,nside))
    n2 = nop(n,side(2,nside))
    n3 = nop(n,side(3,nside))
    ! Two vectors of the side
    x1 = cord(n2,1) - cord(n1,1)
    y1 = cord(n2,2) - cord(n1,2)
    z1 = cord(n2,3) - cord(n1,3)

    x2 = cord(n3,1) - cord(n1,1)
    y2 = cord(n3,2) - cord(n1,2)
    z2 = cord(n3,3) - cord(n1,3)

    ! Normal vector to the side (NORM = X1 X X2 / 2.)
    xnorm = 0.5_8*(y1*z2 - z1*y2)
    ynorm = 0.5_8*(z1*x2 - x1*z2)
    znorm = 0.5_8*(x1*y2 - y1*x2)

    fs(nside,1) = stress(1,n)*xnorm + stress(4,n)*ynorm + stress(5,n)*znorm
    fs(nside,2) = stress(4,n)*xnorm + stress(2,n)*ynorm + stress(6,n)*znorm
    fs(nside,3) = stress(5,n)*xnorm + stress(6,n)*ynorm + stress(3,n)*znorm
 end do


  ! In the nodes (f(i,j), i-node,j-dimension) 
  ! each node belongs to three sides (see "side" above)
  do j=1,3
    f(1,j) = (fs(1,j)+fs(2,j)+fs(4,j))/3.0_8
    f(2,j) = (fs(1,j)+fs(3,j)+fs(4,j))/3.0_8
    f(3,j) = (fs(2,j)+fs(3,j)+fs(4,j))/3.0_8
    f(4,j) = (fs(1,j)+fs(2,j)+fs(3,j))/3.0_8
  end do

  ! Insert forces in whole system:
  do ii=1,4
    nn = nop(n,ii)
    do j=1,3
      force(nn,j) = force(nn,j) - f(ii,j)
      balance(nn,j) = balance(nn,j) + abs(f(ii,j))
    end do
  end do
   
  return
end subroutine force_balance

!*****************************************************
subroutine n_force_balance(dt)
  use sizes
  use node_data
  use boundary_node_data
  use element_data
  implicit none

  logical:: space
  integer:: j,nb,i,ii,i3,nrd,node,nn,nd
  real(kind=8):: fs(4,3),f(4,3),fv,vvs,fff,bbb,vvv,vmax,dt
  
  ! Including volumetric force
  !$OMP PARALLEL
  !$OMP DO PRIVATE(i,i3,fv)
  do i = 1, np
    do i3 = 1,3
      fv = mass(i)*g(i3)
      force(i,i3)   = force(i,i3) - fv
      balance(i,i3) = balance(i,i3) + abs(fv)
    end do
  end do
  !$OMP END DO

  ! BOUNDARY CONDITION IN STRESSES
  !$OMP DO PRIVATE(j,nb,nn)
  do j=1,3
    do nb = 1,nbp
      if( force_code(nb,j) )then
        nn = numbn( nb )
        force(nn,j)   = force(nn,j)   + value(nb,j)
        balance(nn,j) = balance(nn,j) + abs (value(nb,j))
      end if
    end do
  end do
  !$OMP END DO  
  
  !$OMP DO PRIVATE(i,j,vvs,fff)
  ! Damping
  do i=1,np
    do j=1,3
      vvs = vel(i,j) + dt*force(i,j)/(mass(i)*den_scale)
      if ( abs(vvs) .gt. 1.0e-20_8) then
        fff = abs(force(i,j))
        force(i,j) = force(i,j) - demf*sign(fff,vvs)
        vel(i,j)=vel(i,j) + dt*force(i,j)/(mass(i)*den_scale)
      end if
    end do
  end do
  !$OMP END DO
  
  ! Boundary condition in velocities:
  !$OMP DO PRIVATE(i,j,node)
  do i=1,nbp 
    do j=1,3 
      if( vel_code(i,j) )then
        node = numbn(i)
        vel(node,j)= value(i,j)
        force(node,j)= 0.0_8
        balance(node,j)= 0.0_8
      end if
    end do
  end do
  !$OMP END DO
  
  ! Put velocity & balance to zero if node is in open space
  !$OMP DO PRIVATE(i,j,space,nrd,nd)
  do i =1,np
    nrd = nos(i,1)
    space = .true.
    do nd = 2,nrd+1
      space = space .and. (flag( nos(i,nd) ) .ne. 0 )
    end do

    if(space)then  	! open space
      do j = 1,3
        vel(i,j) = 0.0_8
        force(i,j) = 0.0_8
        balance(i,j) = 0.0_8
      end do
    end if
  end do
  !$OMP END DO
  !$OMP END PARALLEL
    
  return
end subroutine n_force_balance

!*****************************************************
! Estimation max balance-off in system (boff)
subroutine boff_calc(dt,vvv)
  use sizes
  use node_data
  implicit none
  
!  integer::n,j,i,noff,nmax,nn
  integer::n,j,i,nmax,nn
  real(kind=8)::dt
  real(kind=8)::bbb,vvv,vmax
  
    	bbb = 0.0_8
  do i =1,np
    do j =1,3
      if(bbb.le.balance(i,j)) bbb = balance(i,j)
    end do
  end do


  boff = 0.0_8
  noff = 0

  do i =1,np
    do j =1,3
      balance(i,j)=abs (force(i,j)/bbb)
     if(balance(i,j).gt.boff ) then
        boff = balance(i,j)
        noff = i
     end if
    end do
  end do

  ! Maximum vel
  vvv = vmax (np,vel,nmax)

! Screen output
 write(6,200)dt,vvv,boff,noff
200 format(' dt,vmax,boff:', 3g11.3, i9)
    
  return
end subroutine boff_calc
