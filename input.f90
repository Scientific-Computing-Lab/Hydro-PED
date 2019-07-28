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
!    * Neither the name of Eyal Shalev, Vladimir Lyakhovsky, Harel Levin or 
!      Gal Oren, nor the names of its contributors may be used to endorse 
!      or promote products derived from this software without specific prior 
!      written permission.
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

#define ALLOC   alloc_if(.true.)
#define FREE    free_if(.true.)
#define RETAIN  free_if(.false.)
#define REUSE   alloc_if(.false.)

subroutine input
  use sizes
  use element_data
  use node_data
  use boundary_node_data
  use diffusion_data
  implicit none
  
  integer:: n,nn,nomat,i,j,ii,jj,k,nbm,ind,nlm,nmat,ne1
  integer:: n1,n2,n3,n4,nj,nf,kf,nbr,kp,nc,nbc(4)
  integer:: ig1,ig2,ig3,ig4,nrd,nd,nside,nfe
  
  logical:: there,outface
  
  real(kind=8):: rl0,rm0,csi0,vp,vs,dist,lmef,muef
  real(kind=8):: rlmax,rmumax,densmin,alp,ran2
  real(kind=8):: s11,s22,s33,s12,s13,s23,x1,x2,y1,y2,z1,z2
  real(kind=8):: xx,yy,zz,dd
  real(kind=8):: q11,q22,fi,xnorm,ynorm,znorm
  real(kind=8):: pres,s_mean,dist1

  real(kind=8),dimension(:),allocatable::rmax
  real(kind=8),dimension(:,:),allocatable::prop
  

 ! read the dimensions of the grid
    inquire(file='grid_size.dat',exist=there)
    if(there)then
      open(1,file="grid_size.dat")
      read(1,*) ne
      read(1,*) np
      read(1,*) nbp
      read(1,*) ke
      read(1,*) mat_size
      close(1)
    else
      write(6,*)' can not find file "grid_size.dat".  Stopping '
      stop
    end if

!-----------------
    write(6,*)' Model size: '
    write(6,*)'  np=',np,'   ne=',ne,'  nbp=',nbp,'   ke=',ke
!-----------------

 ! allocate arrays to the correct sizes
 ! Tetra structure 
      allocate(nop(ne,4),nos(np,ke))

 ! Attributes of elements
  allocate(flag(ne),el_drop(ne), 			     &
    lambda(ne),mu(ne),gr(ne),dens(ne),ductile(ne),           &
    alpha(ne),dalpha(ne),i2(ne),ksi(ne),ksi0(ne),pf_el(ne),  &
    rate(ne,4),el_vol(ne),al_biot(ne),m_biot(ne),            &
    stress0(6,ne),strain0(6,ne),field(ne,3),	             &
    stress(6,ne),strain(6,ne),strainp(6,ne),str_e(6,ne))
           

 ! Attributes of nodes
  allocate(mass(np),force(np,3),pfluid(np),             &
           balance(np,3),cord(np,3),disp(np,3),dspt(np,3),vel(np,3), &
           counter(np) )

 ! Boundary conditions
  allocate(numbn(nbp),vel_code(nbp,3),force_code(nbp,3),value(nbp,3))

 ! Arrays for diffusion
  allocate(id(np),bc_dfval(np),d(15,ne),q(3,ne),a_matrix(mat_size),    &
     order(np),ia(np+1),ja(mat_size),ija(10,ne),ul(4,ne),f(np),f1(np), &
     xl(3,4,ne),shpp(4,4,ne),xsj(ne),s(4,4,ne),p(4,ne),displ(3,4,ne))
 
        over_pres = .0 ! 20.e+6
        print *,'   Over pressure = ',over_pres

!---------- READ dat files --------------------  
!	READ node coordinates 
    inquire(file='coordinates.dat',exist=there)
    if(there)then
      open(1,file='coordinates.dat')
      read(1,*) nn

      if ( nn .ne. np ) then
        write(6,*)' Number of points is not correct ',nn,np
        stop
      end if
      
      do i = 1,np
        read(1,*) cord(i,1),cord(i,2),cord(i,3)
!    initial pressure = hydrostatic SURFACE = 0.
            pfluid(i)=.0 ! -1000.*g(3)*cord(i,3)  ! initial pressure

      end do
      close(1)
      write(6,*) 'end READ node coordinates'  
    else
      write(6,*)'can not find file "coordinates.dat". Stopping'
      stop
    end if

 
                vel  = 0.0_8
 		disp = 0.0_8
                
!-----  Boundary conditions ----
              bc_dfval = 0.0_8
              id    = 0        

    inquire(file='boundary.dat',exist=there)
    if(there)then
        open(1,file='boundary.dat')
        read(1,*) nbm

        if ( nbm .ne. nbp ) then
        print *,'  Number of boundary points is not correct ',nbm,nbp
        stop
        endif

        do i = 1,nbp

        read(1,*) kp, nbc
        
                numbn(i) = kp

!---- type of boundary condition
                do k = 1,3
        vel_code(i,k)   = .false.
        force_code(i,k) = .false.

        if( nbc(k).eq.1 )  vel_code(i,k)   = .true.
        if( nbc(k).eq.2 )  force_code(i,k) = .true.
                value(i,k) =  0.0_8
                end do

              if( nbc(4).eq.1 ) then
        if(force_code(i,1)) bc_dfval(kp)=pfluid(kp) ! +over_pres
        if(  vel_code(i,1)) bc_dfval(kp)=pfluid(kp)

              id(kp)= 1       !  pressure BC

              
              else if ( nbc(4).eq.2 ) then
              bc_dfval(kp)= 0. ! value for BC
              id(kp)= 2       !  flux BC

                end if
        end do

        close(1)
        print *,'end READ BC'
    else
      write(6,*)'can not find file "boundary.dat". Stopping'
      stop
    end if
!-----  end READ BC ----

!-----  READ node connections  ----
    inquire(file='connections.dat',exist=there)
    if(there)then
        open(1,file='connections.dat')
        read(1,*) nn

        if ( nn .ne. np ) then
        print *,'  Number of points is not correct ',nn,np
        stop
        endif

        do i = 1,np
        read(1,*) nc,(nos(i,k+1), k=1,nc)

                nos(i,1) = nc
        if ( nc .gt. ke-1 ) then
        print *,'  ke=number_of_connections+1 is not correct ',i,nc,ke
        stop
        end if
        
        end do

        close (1)
        print *,'end READ node connections'
    else
      write(6,*)'can not find file "connections.dat". Stopping'
      stop
    end if

! READ material properties ------------------
!	Units:
!	density - kg/m^3
!	length  - m
!	stress  - MPa
!	time scale = 1.	!sqrt(1000.) second
  	
! 	tsc = sqrt(1000.0_8)
  	tsc = 1.0_8
  write(6,*)' Numerical time scale: ',tsc
 
  inquire(file='material.dat',exist=there)
  if(there)then
    open(1,file='material.dat')
    read(1,*) nomat
    write(6,*)' number of different materials found =',nomat
  else
    write(6,*)'can not find file "material.dat". Stopping'
    stop
  end if
  
  allocate(prop(nomat,15),rmax(nomat))

  rlmax   = 0.0_8
  rmumax  = 0.0_8
  densmin = 999999.0_8

  do i = 1,nomat
	read(1,*) (prop(i,k), k=1,15)

!	Young nu (poisson) ---> Lambda, Mu (G - rigidity)
	rl0 = prop(i,2)*prop(i,3)/ ( (1.0_8+prop(i,3))*(1.0_8-2.0_8*prop(i,3)) )
	rm0 = prop(i,2) / (2.0_8*(1.0_8+prop(i,3)) )

	prop(i,2) = rl0
	prop(i,3) = rm0

	csi0 = prop(i,5)

!	gamma_r
    rmax(i) = 0.5_8*csi0*((2.0_8*rm0+3.0_8*rl0)/(3.0_8-csi0**2.0_8) + rl0) +         &
      sqrt( ( 0.5_8*csi0*((2.0_8*rm0+3.0_8*rl0)/(3.0_8-csi0**2.0_8) + rl0) )**2 +    &
                2.0_8*rm0*(2.0_8*rm0+3.0_8*rl0)/(3.0_8-csi0**2.0_8) )

  write(6,*) ' -------------------------------------------------'
  write(6,*)'	Material number ',i
  write(6,*)'		Density  = ',prop(i,1)
  write(6,*)'		Lambda   = ',prop(i,2)
  write(6,*)'		Mu       = ',prop(i,3)
  write(6,*)'		M_Biot   = ',prop(i,4)
  write(6,*)'		a_Biot   = ',prop(i,15)
  write(6,*)'		Ksi_0    = ',prop(i,5)
  write(6,*)'		Gamma_r  = ',rmax(i)
  write(6,*)' a1, a2 (permeability) = ',prop(i,13),prop(i,14)
    
!	scale kinetic parameters to time scale
   prop(i, 6) = prop(i, 6)*tsc	! Cd(weakening)
   prop(i, 7) = prop(i, 7)*tsc	! C1 (healing)
   prop(i, 9) = prop(i, 9)/tsc     ! tau_d (dynamic weakening)

!	Pa -- MPa for damage-related viscosity
!	prop(i,12) = prop(i,12)*1.0e+6_8	
  write(6,*)'		R        = ',prop(i,3)*prop(i,12)

!  store maximum Lambda Mu & minimum dens
    if ( rlmax .le. rl0  ) rlmax = rl0
    if ( rmumax .le. rm0 ) rmumax = rm0 
    if ( densmin .ge. prop(i,1) ) densmin = prop(i,1)

      vp = sqrt((rl0+2.*rm0)/prop(i,1))/tsc
	  vs = sqrt(        rm0 /prop(i,1))/tsc
  write(6,*)'		Vp        = ',vp
  write(6,*)'		Vs        = ',vs

  end do
  close (1)
!-----  end READ material properties ------------------
        write(6,*)' end read material properties '

!-----  READ element connections and put properties ----
  inquire(file='elements.dat',exist=there)
  if(there)then
    open(1,file='elements.dat')
    read(1,*) nlm
    if ( nlm .ne. ne ) then
      write(6,*)'  Number of elements is not correct ',nlm,ne
      stop
    end if
  else
    write(6,*)'can not find file "elements.dat". Stopping'
    stop
  end if
    
  do i = 1,ne
   read(1,*) nop(i,1),nop(i,2),nop(i,3),nop(i,4) ! ,nmat (one material type)
	nmat=1
!    if( nmat .gt. nomat )then
!     write(6,*)' Number of material types is not correct ',nmat,nomat
!      stop
!    end if

!  Test for correct volume (det<0)
   call volume_correct_sign(i)

!	Put material properties
        dens   (i) = prop(nmat,1)
        lambda (i) = prop(nmat,2)
        mu     (i) = prop(nmat,3)
        m_biot (i) = prop(nmat,4)
        ksi0   (i) = prop(nmat,5)
        gr     (i) = rmax(nmat)
        rate (i,1) = prop(nmat,6)
        rate (i,2) = prop(nmat,7)
        rate (i,3) = prop(nmat,8)
        rate (i,4) = prop(nmat,9)
!----- Random initial damage ------------        
        call random_number(ran2)
        alpha(i) = prop(nmat,10) + prop(nmat,11)*(ran2-0.5_8) ! rand init damage
        
        ductile(i) = prop(nmat,12)      ! damage related viscosity

        flag(i) = 0


!       alpha Biot 
        al_biot(i) = prop(nmat,15)	!1.0_8

!        put properties for diffusion
!    D(4,i)=Ss [1/Pa]   pressure storage = 1/M = porosity / Kwater
!              Kwater = 2.e+9 Pa
!              M - Biot modulus Pf = M * eps_vol
!    D(5,i)=porosity[]
!    D(6,i)=k11[m2]
!    D(7,i)=k12=k21[m2]
!    D(8,i)=k22[m2]
!    D(9,i)=viscosity [Pa sec]
!    D(10,i)=density[kg/m3]
!    D(11,i)=alpha Biot        assumed const = 1.
!              otherwise change in input.f also


              D(4,i)=1./m_biot(i)
              D(5,i)=0.1

              D(9,i)=1.E-3
              D(10,i)=1000.d0
              D(11,i)=al_biot(i)
              D(12,i)=0.0
              D(13,i)=0.0
              D(14,i)= prop(nmat,13)  !  Values for permeability
              D(15,i)= prop(nmat,14)

              D(6,i)= D(14,i)*exp(D(15,i)*alpha(i))
              D(7,i)= D(6,i)
              D(8,i)= D(6,i)

 end do
!-----  end READ element connections and put properties ----
 close(1)
  write(6,*)' end READ elements and put properties '

! Calculate field points
    field = 0.0_8
  do i = 1,ne
         do j = 1,4
                nj = nop(i,j)
     field(i,1) = field(i,1) + cord(nj,1)/4.0_8 
     field(i,2) = field(i,2) + cord(nj,2)/4.0_8 
     field(i,3) = field(i,3) + cord(nj,3)/4.0_8
                end do  
 
! if (field(i,3) .le. -9.5 .or. field(i,3) .ge. -0.5) rate (i,1) = 0.

 end do

!========================================================       
!-----  search for minimum distance between nodes -----------

   dist = (cord(2,1) - cord(1,1))**2 +     &
          (cord(2,2) - cord(1,2))**2 +     &
          (cord(2,3) - cord(1,3))**2

  do i = 1,ne

    n1 = nop(i,1)
    n2 = nop(i,2)
    n3 = nop(i,3)
    n4 = nop(i,4)

   dist1 = (cord(n1,1) - cord(n2,1))**2 +     &
           (cord(n1,2) - cord(n2,2))**2 +     &
           (cord(n1,3) - cord(n2,3))**2
    if ( dist1 .le. dist ) dist = dist1

   dist1 = (cord(n1,1) - cord(n3,1))**2 +     &
           (cord(n1,2) - cord(n3,2))**2 +     &
           (cord(n1,3) - cord(n3,3))**2
    if ( dist1 .le. dist ) dist = dist1

   dist1 = (cord(n1,1) - cord(n4,1))**2 +     &
           (cord(n1,2) - cord(n4,2))**2 +     &
           (cord(n1,3) - cord(n4,3))**2
    if ( dist1 .le. dist ) dist = dist1

   dist1 = (cord(n2,1) - cord(n3,1))**2 +     &
           (cord(n2,2) - cord(n3,2))**2 +     &
           (cord(n2,3) - cord(n3,3))**2
    if ( dist1 .le. dist ) dist = dist1

   dist1 = (cord(n2,1) - cord(n4,1))**2 +     &
           (cord(n2,2) - cord(n4,2))**2 +     &
           (cord(n2,3) - cord(n4,3))**2
    if ( dist1 .le. dist ) dist = dist1

   dist1 = (cord(n3,1) - cord(n4,1))**2 +     &
           (cord(n3,2) - cord(n4,2))**2 +     &
           (cord(n3,3) - cord(n4,3))**2
    if ( dist1 .le. dist ) dist = dist1

  end do

		ddm = sqrt (dist)
    write(6,*)' min node to node distance : ',ddm ,' [m]'
!       calculate wave-related time scale (vp2=time_scale**2)
        vp2 = ( rlmax + 2.0_8* rmumax ) / densmin / dist

!========================================================
  do n  = 1,ne
!       fluid pressure in elements
                n1 = nop(n,1)
                n2 = nop(n,2)
                n3 = nop(n,3)
                n4 = nop(n,4)

        pf_el(n) = (pfluid(n1)+pfluid(n2)+pfluid(n3)+pfluid(n4))/4.
  end do
!------  Initial stress & strain distribution ---------------- 
  do i = 1,ne
          pres = dens(1)*9.81*(field(i,3)-1000.0)

	s33 = pres
	s11 =  1.5*s33  
	s22 =  0.5*s33 
	s12 = 0.
	s13 = 0.
	s23 = 0.

       s_mean = s11+s22+s33

	stress(1,i) = s11
	stress(2,i) = s22
	stress(3,i) = s33
	stress(4,i) = s12
	stress(5,i) = s13
	stress(6,i) = s23
	
	strain(1,i)=(s11-s_mean*lambda(i)/(3.*lambda(i)+2.*mu(i)))/(2.*mu(i))
	strain(2,i)=(s22-s_mean*lambda(i)/(3.*lambda(i)+2.*mu(i)))/(2.*mu(i))
	strain(3,i)=(s33-s_mean*lambda(i)/(3.*lambda(i)+2.*mu(i)))/(2.*mu(i))
	strain(4,i)= s12/(2.*mu(i))
	strain(5,i) = 0.
	strain(6,i) = 0.
        
  end do

	  strainp = 0.0_8

               
        print *,'   No IJA Array '

  phi_length=4*ne/9
  phi_length=((phi_length+15103)/15104)*15104
  !dir$ offload_transfer target(mic:0) in(cord : ALLOC RETAIN) in(vel : ALLOC RETAIN) in(m_biot : ALLOC RETAIN) in(pf_el : ALLOC RETAIN) in(strain : ALLOC RETAIN) in(nop : ALLOC RETAIN)
  !dir$ offload_transfer target(mic:1) in(cord : ALLOC RETAIN) in(vel : ALLOC RETAIN) in(m_biot : ALLOC RETAIN) in(pf_el : ALLOC RETAIN) in(strain : ALLOC RETAIN) in(nop : ALLOC RETAIN)
  return
 end subroutine input
 !==================================================================
 
subroutine volume_correct_sign(n)
  use sizes
  use element_data
  use node_data
  implicit none
  
  integer:: j,n,itemp
  
  real(kind=8):: x(4),y(4),z(4),det

          do j = 1,4
               x(j) = cord ( abs(nop(n,j)),1 )
               z(j)=  cord ( abs(nop(n,j)),3 )
               y(j) = cord ( abs(nop(n,j)),2 )
	  end do

          det=0
         det=det+(x(2)-x(1))*(y(3)-y(1))*(z(4)-z(1))
         det=det+(y(2)-y(1))*(z(3)-z(1))*(x(4)-x(1))
         det=det+(z(2)-z(1))*(x(3)-x(1))*(y(4)-y(1))

         det=det-(z(2)-z(1))*(y(3)-y(1))*(x(4)-x(1))
         det=det-(y(2)-y(1))*(x(3)-x(1))*(z(4)-z(1))
         det=det-(x(2)-x(1))*(z(3)-z(1))*(y(4)-y(1))


         if( det > 0.0_8 ) then
 print *,'   the tetra vertices must be reordered for a NEGATIVE volume '
 stop
!          itemp=nop(n,1)
!          nop(n,1)=nop(n,2)
!          nop(n,2)=itemp
         end if

 return
end subroutine volume_correct_sign
