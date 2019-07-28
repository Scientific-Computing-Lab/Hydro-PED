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
! DISCLAIMED. IN NO EVENT SHALL Eyal Shalev, Vladimir Lyakhovsky
! BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
! OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND  ON ANY THEORY OF LIABILITY, 
! WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
! OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
! ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

!  poroelastic element.  tetra element. solves pore pressure at 3 corner nodes.

      !> @brief         Calculates local stiffness matrix and load vector for
      !> a poroelastic element
      !> @details       Uses tetrahedron geometry and linear shape function.
      !> The derivatives of the displacemnts are used as source term for fluid 
      !> with the matrix L Transpose. This subroutine is called first to compute 
      !> the load vector, then it is called for the second time to compute the 
      !> stiffnes matrix. 
      !> <br>\f$\left(S+\Theta H \Delta t\right)p_{n+1}=L^T\Delta u+\left[S-\left(1-\Theta\right)H\Delta t\right]p_n+F\Delta t\f$
      !> <br>\f$\mathbf{S}=\int\limits_\Omega N^TS_\varepsilon N \, d\Omega\f$
      !> <br>\f$\mathbf{H}=\int\limits_\Omega B^T\frac{\underline{\mathbf{k}}}{\mu}B \, d\Omega\f$
      !> <br>\f$\mathbf{L}^T=\int\limits_\Omega N^T\alpha_BB \, d\Omega\f$
      !> <br>\f$\mathbf{F}=-\int\limits_\Gamma N^Tq \, d\Gamma - \int\limits_\Omega B^T\frac{k}{\mu}B\rho g\nabla z \, d\Omega\f$
      !> <br>\f$\Theta\in[0,1]\f$
      !> <br>\f$N^T=\text{transpose of shape functions}\f$
      !> <br>\f$B^T=\text{transpose of shape functions derivatives}\f$
      !>
      !> @param[in]     Pore pressure at 4 nodes (ul) and displacements (displ) 
      !> of previous time step, material properties (D), and timestep (deltat).
      !> @param[out]    Local stiffness matrix and load vector of one element.
      subroutine tetra4(d,ul,displ,xl,ql,s,p,	&
                     shp,xsja,deltat,tweight,isw)
      implicit double precision (a-h,o-z)

      integer a,b

      dimension d(15),ul(4),xl(3,1),s(4,4),p(4),displ(3,4)
      dimension ql(3),shp(4,4),s1(12,4),st(4,12)

!  print *,' d ', d
!  print *,' disp ', displ
!  print *,' xl ', xl
!  print *,' ul ', ul
      


        g=9.81 

!---- compute tangent stiffness matrix and/or internal forces
!C    D(5)=porosity[]
!C    D(6)=k11[m2]
!C    D(7)=k22[m2]
!C    D(8)=k33[m2]
!C    D(9)=viscosity[kg/m yr]
!C    D(10)=density[kg/m3]
!C    D(11)=alpha Biot
!C  Initialize all matrices
	do i=1,4
	   do j=1,4
	      s(j,i)=0.d0
	   end do
!	   p(i)=0.d0
	end do
	 do i=1,4
	      do j=1,3*4
	        s1(j,i)=0.d0
	        st(i,j)=0.d0
	      end do
	 end do

!   compute L
	alpha=d(11)
	 do i=1,4
	   w1=shp(1,i)*xsja
	   w2=shp(2,i)*xsja
	   w3=shp(3,i)*xsja
	      do j=1,4
	        s1(3*i-2,j)=s1(3*i-2,j)+		&
        			(w1*alpha*shp(4,j))
	        s1(3*i-1,j)=s1(3*i-1,j)+		&
        			(w2*alpha*shp(4,j))
	        s1(3*i,j)=s1(3*i,j)+			&
        			(w3*alpha*shp(4,j))
	      end do
	 end do

!   compute L Transpose
	 do i=1,4
	      do j=1,4
	        st(j,3*i-2)=s1(3*i-2,j)
	        st(j,3*i-1)=s1(3*i-1,j)
	        st(j,3*i)=s1(3*i,j)
	      end do
	 end do

!   compute S
	  do i=1,4
	    do j=1,4
	      s(i,j)=s(i,j)+				&
      	             shp(4,i)*d(4)*shp(4,j)*xsja
	    end do
	  end do

!   compute H
	if (isw .eq. 6) then

	do i=1,4
	   p(i)=0.d0
	end do

	 do i=1,4
            w1 = d(6)/d(9)*shp(1,i)*xsja
            w2 = d(7)/d(9)*shp(2,i)*xsja
            w3 = d(8)/d(9)*shp(3,i)*xsja
	     do j=1,4
              s(i,j) = s(i,j)+(w1*shp(1,j) + w2*shp(2,j)		&
                      + w3*shp(3,j))*(-(1-tweight))*deltat 

	     end do	   
	 end do	   


	do i=1,4
	p(i)=s(i,1)*ul(1)+s(i,2)*ul(2)+s(i,3)*ul(3)+			&
      		s(i,4)*ul(4) 
	end do

!   compute p
	   do i=1,4
	      p(i)=p(i)-						&
     	      (shp(3,i)*D(8)/D(9)*D(10)*g*xsja)*deltat
	   end do

!   add LT as source to the diffusion eq.
	do i=1,4
	  p(i)=p(i)-							&
         (st(i,1)*displ(1,1)+st(i,2)*displ(2,1)+st(i,3)*displ(3,1)	&
         +st(i,4)*displ(1,2)+st(i,5)*displ(2,2)+st(i,6)*displ(3,2)	&
         +st(i,7)*displ(1,3)+st(i,8)*displ(2,3)+st(i,9)*displ(3,3)	&
         +st(i,10)*displ(1,4)+st(i,11)*displ(2,4)+st(i,12)*displ(3,4))
	end do

	else

!   compute H

	 do i=1,4
            w1 = d(6)/d(9)*shp(1,i)*xsja
            w2 = d(7)/d(9)*shp(2,i)*xsja
            w3 = d(8)/d(9)*shp(3,i)*xsja
	     do j=1,4
              s(i,j) = s(i,j)+(w1*shp(1,j) + w2*shp(2,j)	&
                      + w3*shp(3,j))*tweight*deltat 

	     end do	   
	 end do	   

!	print *,' TETRA ',s
	end if

      return
      end

!---- compute flux
      !> @brief         Calculates fluid velocity at an element
      !> @details       In case fluid velocities are needed for the transport 
      !> equation, "flux" calculate the velocity after the solver returned the
      !> pore pressure at each node. 
      !>
      !> @param[in]     Pore pressure at 4 nodes (ul) and permeability
      !> @param[out]    Darcy velocity vector
      subroutine flux(d,ul,xl,ql,shp)
      implicit double precision (a-h,o-z)

      integer a,b

      dimension d(15),ul(4),xl(3,1),s(4,4)
      dimension ql(3),shp(4,4)

        g=9.81 

 	      ql(1)=0.d0
	      ql(2)=0.d0
	      ql(3)=0.d0

	   do j=1,4
	      ql(1)=ql(1)-1.0/d(9)*d(6)*shp(1,j)	&
          	        *(ul(j)+d(10)*g*xl(3,j))
	      ql(2)=ql(2)-1.0/d(9)*d(7)*shp(2,j)	&
          	        *(ul(j)+d(10)*g*xl(3,j))
	      ql(3)=ql(3)-1.0/d(9)*d(8)*shp(3,j)	&
          	        *(ul(j)+d(10)*g*xl(3,j))
	   end do

      return
      end

      !> @brief         Calculates linear shape functions and their derivatives 
      !> for a tetrahedron element
      !> @details       Uses coordinates to calculate determinant (volume/6) and 
      !> shape functions
      !> <br>\f$L_k=a_k+b_kx_1+c_kx_2+d_kx_3\f$
      !> <br>\f$N^T=\text{transpose of shape functions}\f$
      !> <br>\f$B^T=\text{transpose of shape functions derivatives}\f$
      !> <br>\f$V=\frac{1}{6}\text{det}
      !> \begin{vmatrix}
      !> 1 & x_i & y_i & z_i \\
      !> 1 & x_j & y_j & z_j \\
      !> 1 & x_m & y_m & z_m \\
      !> 1 & x_n & y_n & z_n
      !> \end{vmatrix}\f$
      !>
      !> @param[in]     Coordinates of the nodes of the tetrahedron (xl) 
      !> @param[out]    Shape functions (shp(4,i)) and their derivatives (shp(1-3,i))

      subroutine shape(xl,shp,xsja)
      implicit double precision (a-h,o-z)

      dimension xl(3,4),shp(4,4)
      dimension a1(4,3)

!       Compute determinants for transformation minors

        do i = 1,4
          j      = mod(i,4) + 1
          k      = mod(j,4) + 1
          l      = mod(k,4) + 1
          a1(i,1) = xl(2,j)*(xl(3,k) - xl(3,l))	&
                + xl(2,k)*(xl(3,l) - xl(3,j))	&
                + xl(2,l)*(xl(3,j) - xl(3,k))

          a1(i,2) = xl(3,j)*(xl(1,k) - xl(1,l))	&
                + xl(3,k)*(xl(1,l) - xl(1,j))	&
                + xl(3,l)*(xl(1,j) - xl(1,k))

          a1(i,3) = xl(1,j)*(xl(2,k) - xl(2,l))	&
                + xl(1,k)*(xl(2,l) - xl(2,j))	&
                + xl(1,l)*(xl(2,j) - xl(2,k))
        end do ! i

!       Correct signs on determinants

        do i = 1,3
          a1(1,i) = -a1(1,i)
          a1(3,i) = -a1(3,i)
        end do ! i

!       Determinant for element volume

        xsja  = (xl(1,1)*a1(1,1) + xl(1,2)*a1(2,1)	&
             + xl(1,3)*a1(3,1) + xl(1,4)*a1(4,1))

        if(xsja.ne.0.0d0) then
          detr = 1.d0/xsja
        else
          print *, ' TETSHP: Determinant =',xsja
	  print *,'  xi  ',0.25
	  print *,'  xl  ',xl
	  stop
          detr = 1.d0
        endif

        do i = 1,4
          shp(1,i) = a1(i,1)*detr
          shp(2,i) = a1(i,2)*detr
          shp(3,i) = a1(i,3)*detr
          shp(4,i) = 0.25d0
        end do ! i

         xsja=xsja/6

      return
      end


