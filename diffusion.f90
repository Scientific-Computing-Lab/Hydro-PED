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

      !> @brief         Calculates pore pressure using fluid mass conservation
      !> @details       Uses the Finite Element Method to calculate the
      !> pressure diffusion of fluids in the rock in response to pressure
      !> gradient and rock stress forcing. Uses implicit time step.
      !> \f$\nabla\cdot\left(\frac{\underline{\mathbf{k}}\left(\alpha_D\right)}{\mu_f}\left(\nabla p+\rho g \mathbf{z}\right)\right)=S_\varepsilon\frac{\partial p}{\partial t}+\alpha_B\frac{\partial\varepsilon_{kk}}{\partial t}\f$
      !>
      !> @param[in]     Pore pressure (pfluid) and displacements (disp)  
      !> of previous time step, material properties (D), and timestep (deltat).
      !> @param[out]    Updated pore pressure (pfluid).
      subroutine diffusion(deltat,info)

       use sizes
       use diffusion_data
       use node_data
       use element_data



	implicit double precision (a-h,o-z)
	
	integer info
	real(kind=8) :: tsolve_begin,tsolve_end


!c one time step for diffusion
!c
	f = bc_dfval
	a_matrix = 0.0_8
	xsj = 0.0_8 
!--------- Update permeability ----------------------
		do i = 1,ne

	      D(6,i) = D(14,i) * exp(D(15,i)*alpha(i))
	if (D(6,i) .gt. 1e-10) D(6,i)=1e-10
	      D(7,i) = D(6,i)
	      D(8,i) = D(6,i)
		
		end do

!***************************************************************************
		tweight=0.5_8

!***************************************************************************
!     compute shape functions
!---- set up local arrays
      do  n = 1,ne
         do  i = 1,4
            ii = nop(n,i)
      	       do  j = 1,3
   		  xl(j,i,n) = cord(ii,j)
	       end do
   		  ul(i,n) = pfluid(ii)
                do j = 1,3
    	 	  displ(j,i,n) = disp(ii,j)
	       end do
	   end do
	end do
 
!  ***************************************************************
        do n = 1,ne
         call shape(xl(1,1,n),shpp(1,1,n),xsj(n))
	end do
!  ***************************************************************
!     load vector
        do n = 1,ne
      call tetra4(d(1,n),ul(1,n),displ(1,1,n),xl(1,1,n),	&
      q(1,n),s(1,1,n),p(1,n),shpp(1,1,n),xsj(n),deltat,tweight,6)
	end do
!  ***************************************************************
!     stifness matrix element
        do n = 1,ne
      call tetra4(d(1,n),ul(1,n),displ(1,1,n),xl(1,1,n),	&
      q(1,n),s(1,1,n),p(1,n),shpp(1,1,n),xsj(n),deltat,tweight,3)
	   end do

! Symmetric matrix ... only want lower triangle
 
       do n = 1,ne

           do ka=1, 4
              ii = nop(n,ka)
              f(ii) = f(ii) - p(ka,n)
	   end do

	      a_matrix(ija(1,n))=a_matrix(ija(1,n))-s(1,1,n)
	      a_matrix(ija(2,n))=a_matrix(ija(2,n))-s(1,2,n)
	      a_matrix(ija(3,n))=a_matrix(ija(3,n))-s(1,3,n)
	      a_matrix(ija(4,n))=a_matrix(ija(4,n))-s(1,4,n)
	      a_matrix(ija(5,n))=a_matrix(ija(5,n))-s(2,1,n)
	      a_matrix(ija(6,n))=a_matrix(ija(6,n))-s(2,2,n)
	      a_matrix(ija(7,n))=a_matrix(ija(7,n))-s(2,3,n)
	      a_matrix(ija(8,n))=a_matrix(ija(8,n))-s(2,4,n)
	      a_matrix(ija(9,n))=a_matrix(ija(9,n))-s(3,1,n)
	      a_matrix(ija(10,n))=a_matrix(ija(10,n))-s(3,2,n)
	      a_matrix(ija(11,n))=a_matrix(ija(11,n))-s(3,3,n)
	      a_matrix(ija(12,n))=a_matrix(ija(12,n))-s(3,4,n)
	      a_matrix(ija(13,n))=a_matrix(ija(13,n))-s(4,1,n)
	      a_matrix(ija(14,n))=a_matrix(ija(14,n))-s(4,2,n)
	      a_matrix(ija(15,n))=a_matrix(ija(15,n))-s(4,3,n)
	      a_matrix(ija(16,n))=a_matrix(ija(16,n))-s(4,4,n)
	end do
	
	   f1 = f

	call modify2 

       !open (12,file='a.dat')
       !do ii=1,3414
       !write(12,*)a_matrix(ii)
       !end do
       !close (12)
       !open (12,file='b.dat')
       !do ii=1,509
       !read(12,*)f(ii)
       !end do
       !close (12)
!  ***************************************************************
!     solve
     !call solve87(np,ia,ja,a_matrix,f,order,keep,control,info)
		 
     call cpu_time(tsolve_begin)
     call run_solver_CG(a_matrix, f, info)
     call cpu_time(tsolve_end)
     print *,'DEBUG TRILINOS - Solver took',tsolve_end-tsolve_begin,' seconds.'
!  ***************************************************************
!     fluxes
	
       !open (12,file='x.dat')
       !do ii=1,509
       !write(12,*)f(ii)
       !end do
       !close (12)
       !stop	
	pfluid = f	
       	disp = 0.
      
      	do  n = 1,ne
           do  i = 1,4
   	      ul(i,n) = pfluid(nop(n,i))
	   end do
	end do

         do n = 1,ne
           call flux(d(1,n),ul(1,n),xl(1,1,n),q(1,n),shpp(1,1,n))
	end do

  do n  = 1,ne
!       fluid pressure in elements
                n1 = nop(n,1)
                n2 = nop(n,2)
                n3 = nop(n,3)
                n4 = nop(n,4)

        pf_el(n) = (pfluid(n1)+pfluid(n2)+pfluid(n3)+pfluid(n4))/4.
  end do

	print *, ' 	=== Diffusion step done === ',deltat
      return
      end
