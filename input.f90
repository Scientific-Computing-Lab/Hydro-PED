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
!    * Neither the name of Eyal Shalev, Vladimir Lyakhovsky, Harel Levin 
!      or Gal Oren, nor the names of its contributors may be used to endorse 
!      or promote products derived from this software without specific prior 
!      written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
! ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
! WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL Eyal Shalev, Vladimir Lyakhovsky, Harel Levin 
! & Gal Oren BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
! OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND  ON ANY THEORY OF LIABILITY, 
! WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
! OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
! ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

!> @brief         Read and generate all necessary input data and allocate arrays. 
!> @details       Read coordinates, lithology, element connections, and boundary conditions.
!>  boundary conditions. Generate initial conditions. Set material properties to each element
!> gradient and rock stress forcing. Uses implicit time step.
!> using the lithology.
!>
!> @param[in]    
!> @param[out]    All initial values, boundary conditions, and material properties.
subroutine input
  use sizes
  use element_data
  use node_data
  use boundary_node_data
  use diffusion_data
  implicit none
  
  integer:: n,nn,nomat,i,j,k,nbm,nlm,nmat
  integer:: n1,n2,n3,n4,nj,kp,nc,nbc(4)
  integer:: ig1,ig2,ig3,ig4,kk
  
  logical:: there,outface
  character(len=1):: symb
  
  real(kind=8):: rl0,rm0,csi0,vp,vs,lmef,muef
  real(kind=8):: rlmax,rmumax,densmin,ran2
  real(kind=8):: pfl, peff, depth
  real(kind=8):: s11,s22,s33,s12,s13,s23
  real(kind=8):: s_mean,dist,dist1,sx,sy,sz,fi_eq

  real(kind=8),dimension(:),allocatable::rmax
  real(kind=8),dimension(:,:),allocatable::prop
  
!-------------------------------------------------
!       power for pres n=2.0 in kinetic eq.
                power = 2.0
                
!  Define initial stress conditions
		depth = 5000.	! depth 
	sz = -2600.0*9.81*depth ! -rho*g*depth 
	sx = 1.2*sz
	sy = 0.8*sz

	s_mean = (sx+sy+sz)/3.0_8
	
        	pfl = 10.0
	print *,'	Fluid pressure     = ',pfl
		peff = -s_mean - pfl
	print *,'	Effective pressure = ', peff

! read the dimensions of the grid
    inquire(file='grid_size.dat',exist=there)
    if(there)then
      open(1,file="grid_size.dat")
      read(1,*) ne
      read(1,*) np
      read(1,*) nbp
      read(1,*) ke
      read(1,*) mat_size
      read(1,*) nprocs
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
  allocate(flag(ne),el_drop(ne),nop(ne,4),nos(np,ke),num_drop(ne))

 ! Attributes of elements
  allocate(lambda(ne),mu(ne),gr(ne),dens(ne),ductile(ne,2),   &
    el_vol(ne),al_biot(ne),m_biot(ne),phi(ne),dphi(ne),       &
    zi_el(ne),alpha(ne),dalpha(ne),                           &
    i1(ne),i2(ne),ksi(ne),ksi0(ne),ksif(ne),                  &
    pf_el(ne),rate(ne,5),phi_eq(ne,3),coupl(ne,3),            &
    stress(6,ne),strain(6,ne),strainp(6,ne),str_e(6,ne),      &
    stress0(6,ne),strain0(6,ne),field(ne,3))


 ! Attributes of nodes
  allocate(mass(np),force(np,3),pfluid(np),           &
           balance(np,3),cord(np,3),disp(np,3),dspt(np,3),vel(np,3) )

 ! Boundary conditions
  allocate(numbn(nbp),vel_code(nbp,3),force_code(nbp,3),value(nbp,3))

 ! Arrays for diffusion
  allocate(id(np),idt(np),order(np),ia(np+1),ja(mat_size),ija(16,ne),  &
  bc_dfval(np),a_matrix(mat_size),f(np),f1(np),xsj(ne),        &
     d(20,ne),q(3,ne),ul(4,ne),tl(4,ne),p(4,ne),               &
     xl(3,4,ne),shpp(4,4,ne),s(4,4,ne),displ(3,4,ne))
  
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
!    initial pressure = hydrostatic SURFACE 2.7=rock density
            pfluid(i)= pfl	!-1000.*g(3)*cord(i,3)  ! initial pressure
        cord(i,1)=cord(i,1)/100
        cord(i,2)=cord(i,2)/100
        cord(i,3)=cord(i,3)/100
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
        	bc_dfval(kp)=pfluid(kp)

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
  
  allocate(prop(nomat,30),rmax(nomat))

  rlmax   = 0.0_8
  rmumax  = 0.0_8
  densmin = 999999.0_8

!--------- General properties --------------------
        do i = 1,nomat
	read(1,*) (prop(i,k), k=1,4)
	end do

!--------- skip 4 comment lines --------------
	do kk = 1,4
	read(1,'(A1)') symb
		if( symb .ne. '#' ) then
	print *,' Start commemt line in material.dat with # '
		stop
		end if
	end do
		 
!--------- Damage evolution --------------------
        do i = 1,nomat
	read(1,*) (prop(i,k), k=5,12)
	end do

!--------- skip 4 comment lines --------------
	do kk = 1,4
	read(1,'(A1)') symb
		if( symb .ne. '#' ) then
	print *,' Start commemt line in material.dat with # '
		stop
		end if
	end do		 

!--------- Permeability & Porosity  evolution ----------
        do i = 1,nomat
	read(1,*) (prop(i,k), k=13,24)
	end do

!--------- skip 4 comment lines --------------
	do kk = 1,4
	read(1,'(A1)') symb
		if( symb .ne. '#' ) then
	print *,' Start commemt line in material.dat with # '
		stop
		end if
	end do		 
!--------- Heat  evolution ----------
        do i = 1,nomat
	read(1,*) (prop(i,k), k=25,28)
	end do

	close(1)


        do i = 1,nomat

!       Young nu ---> Lambda, Mu
        rl0 = prop(i,2)*prop(i,3)/ ( (1.+prop(i,3))*(1.-2.*prop(i,3)) )
        rm0  = prop(i,2) / (2.*(1.+prop(i,3)) )

        prop(i,2) = rl0
        prop(i,3) = rm0

        csi0 = prop(i,5)

        rmax(i) = 0.5*csi0*((2*rm0+3*rl0)/(3-csi0**2) + rl0) +       &
          sqrt( (0.5*csi0*((2*rm0+3*rl0)/(3-csi0**2) + rl0))**2 +    &
          2*rm0*(2*rm0+3*rl0)/(3-csi0**2) )

        print *,'       Material number ',i
        print *,'               Density = ',prop(i,1)
        print *,'               Lambda  = ',prop(i,2)
        print *,'               Mu      = ',prop(i,3)
        print *,'               M Biot  = ',prop(i,4)
        print *,'               Ksi_0   = ',prop(i,5)
        print *,'               Gamma_r = ',rmax(i)
	print *,' a1, a2 (permeability) = ',prop(i,13),prop(i,14)

        vp = sqrt((rl0+2.*rm0)/prop(i,1))/tsc
        vs = sqrt(        rm0 /prop(i,1))/tsc

        print *,'               Vp      = ',vp
        print *,'               Vs      = ',vs

!       scale kinetic parameters to time scale
        prop(i, 6) = prop(i, 6)*tsc     ! Cd(weakening)
        prop(i, 7) = prop(i, 7)*tsc     ! C1 (healing)
	prop(i, 9) = prop(i, 9)*tsc     ! damage-related viscosity
	prop(i,10) = prop(i,10)/tsc     ! tau_d (dynamic weakening)
	prop(i,17) = prop(i,17)*tsc     ! pressure-driven compaction (exp)
	prop(i,24) = prop(i,24)*tsc     ! fluidity = 1/Maxwell viscosity

!       R-value for damage-related viscosity
        print *,'               R       = ',prop(i,3)*prop(i,9)
           if ( prop(i,3)*prop(i,24) .le. 1.e-15) then
        print *,'         Maxwell time  =   Inf.'
           else
        print *,'         Maxwell time  = ',tsc/(prop(i,3)*prop(i,24))
           end if
!       store maximum Lambda Mu & minimum dens
        if ( rlmax .le. rl0  ) rlmax = rl0
        if ( rmumax .le. rm0 ) rmumax = rm0
        if ( densmin .ge. prop(i,1) ) densmin = prop(i,1)

        end do
!-----	end READ material properties ------------------
        write(6,*)' end read material properties '

!-----  READ element connections and put properties ----

        flag  = 0
        num_drop = 0
!       alpha Biot assumed constant and  = 1.
        al_biot = 1.

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
   read(1,*) nop(i,1),nop(i,2),nop(i,3),nop(i,4),nmat

    if( nmat .gt. nomat )then
      write(6,*)' Number of material types is not correct ',nmat,nomat
      stop
    end if

!  Test for correct volume (det<0)
   call volume_correct_sign(i)

!	Put material properties
        dens   (i)   = prop(nmat,1)
        lambda (i)   = prop(nmat,2)
        mu     (i)   = prop(nmat,3)
	m_biot (i)   = prop(nmat,4)
        ksi0   (i)   = prop(nmat,5)
        gr     (i)   = rmax(nmat)
        rate (i,1)   = prop(nmat,6)      ! Damage Cd
        rate (i,2)   = prop(nmat,7)      ! Damage C1
        rate (i,3)   = prop(nmat,8)      ! Damage C2
        ductile(i,1) = prop(nmat,9)      ! damage related viscosity
        rate (i,4)   = prop(nmat,10)     ! Damage tau_d
!----- Random initial damage ------------        
        call random_number(ran2)
        alpha(i) = prop(nmat,11) + prop(nmat,12)*(ran2-0.5_8) ! rand init damage

!	prop(nmat,13)	! permeability below
!	prop(nmat,14)	! permeability

        call random_number(ran2)
	phi  (i) = prop(nmat,15) + prop(nmat,16)*(ran2-0.5_8) ! rand init porosity

        rate (i,5) = prop(nmat,17)       ! Porosity compaction rate

        coupl (i,1) = prop(nmat,18)      ! coupling D1
        coupl (i,2) = prop(nmat,19)      ! coupling D2
        coupl (i,3) = prop(nmat,20)      ! coupling D3

	phi_eq(i,1)  = prop(nmat,21)
	phi_eq(i,2)  = prop(nmat,22)
	phi_eq(i,3)  = prop(nmat,23)
	ductile(i,2) = prop(nmat,24)	 ! Maxwell viscosity

   
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

	      D(4,i)=1./m_biot(i)	!5.E-9
	      D(5,i)=phi(i)

	      D(9,i)=1.E-3
	      D(10,i)=1000.d0
	      D(11,i)=al_biot(i)	!1.
	      D(12,i)=0.0
	      D(13,i)=0.0
	      D(14,i)= prop(nmat,13)  !  Values for permeability
	      D(15,i)= prop(nmat,14)

	      D(6,i)= D(14,i)*exp(D(15,i)*alpha(i))
	      D(7,i)= D(6,i)
	      D(8,i)= D(6,i)
	      D(16,i)= prop(nmat,25)	!    D(16)=alpha L
	      D(17,i)= prop(nmat,26)	!    D(17)=alpha T
	      D(18,i)=dens(i)		!    D(18)=density solid[kg/m3]
	      D(19,i)= prop(nmat,27)	!    D(19)=heat capacity solid[J/kg C]
	      D(20,i)= prop(nmat,28)	!    D(20)=thermal conductivity solid [W/m/C]


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

!	if ( field(i,3) .le. -4.9 ) rate(i,1) = 0.
!	if ( field(i,3) .ge. -0.1 ) rate(i,1) = 0.

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
!  injection points
!
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

	stress(1,i) = sx
	stress(2,i) = sy
	stress(3,i) = sz
	stress(4,i) = 0.0
	stress(5,i) = 0.0
	stress(6,i) = 0.0

        s11 = sx + al_biot(i)*pf_el(i)  
        s22 = sy + al_biot(i)*pf_el(i)
        s33 = sz + al_biot(i)*pf_el(i)
        s12 =  0.0_8
        s13 =  0.0_8
        s23 =  0.0_8
        s_mean = s11+s22+s33
	
	ksi(i) = ksi0(i)
	lmef = lambda(i) - alpha(i)*gr(i)/ksi(i)
	muef = mu(i) + alpha(i)*gr(i)*ksi0(i) - 0.5_8*alpha(i)*gr(i)*ksi(i)
	
        strain(1,i) = (s11-s_mean*lmef/(3.0_8*lmef+2.0_8*muef))/(2.0_8*muef)
        strain(2,i) = (s22-s_mean*lmef/(3.0_8*lmef+2.0_8*muef))/(2.0_8*muef)
        strain(3,i) = (s33-s_mean*lmef/(3.0_8*lmef+2.0_8*muef))/(2.0_8*muef)
        strain(4,i) =  s12/(2.0_8*muef)
        strain(5,i) =  s13/(2.0_8*muef)
        strain(6,i) =  s23/(2.0_8*muef)

	i1(n) = strain(1,i) + strain(2,i) + strain(3,i)
  zi_el(i) = pf_el(i)/m_biot(i) + i1(n)*al_biot(i) + phi(i)
        
  end do

	  strainp = 0.0_8

!========================================================
        OPEN(12,FILE='ia.in')
        OPEN(14,FILE='ja.in')
        do i=1,np+1
          read(12,*)ia(i)
        end do
  print *,'  read ia '        
        
        do i=1,ia(np+1)-1
          read(14,*)ja(i)
        end do
        close(12)
        close(14)
  print *,'  read ja '        
        do n=1,ne
           ig1=nop(n,1)
           ig2=nop(n,2)
           ig3=nop(n,3)
           ig4=nop(n,4)

	   do j=ia(ig1),ia(ig1+1)-1
                if (ja(j) .eq. ig1) then
                    ija(1,n)=j
                 elseif (ja(j) .eq. ig2) then
                    ija(2,n)=j
                 elseif (ja(j) .eq. ig3) then
                    ija(3,n)=j
                 elseif (ja(j) .eq. ig4) then
                    ija(4,n)=j
                 endif
              end do

	   do j=ia(ig2),ia(ig2+1)-1
                if (ja(j) .eq. ig1) then
                    ija(5,n)=j
                 elseif (ja(j) .eq. ig2) then
                    ija(6,n)=j
                 elseif (ja(j) .eq. ig3) then
                    ija(7,n)=j
                 elseif (ja(j) .eq. ig4) then
                    ija(8,n)=j
                 endif
              end do

	   do j=ia(ig3),ia(ig3+1)-1
                if (ja(j) .eq. ig1) then
                    ija(9,n)=j
                 elseif (ja(j) .eq. ig2) then
                    ija(10,n)=j
                 elseif (ja(j) .eq. ig3) then
                    ija(11,n)=j
                 elseif (ja(j) .eq. ig4) then
                    ija(12,n)=j
                 endif
              end do

	   do j=ia(ig4),ia(ig4+1)-1
                if (ja(j) .eq. ig1) then
                    ija(13,n)=j
                 elseif (ja(j) .eq. ig2) then
                    ija(14,n)=j
                 elseif (ja(j) .eq. ig3) then
                    ija(15,n)=j
                 elseif (ja(j) .eq. ig4) then
                    ija(16,n)=j
                 endif
              end do

        end do
                
        print *,'   Array IJA  is ready '

  return
 end subroutine input
 !==================================================================
 
!> @brief         Check if element node numbers are assined in the correct order 
!> @details       Describing the elements by their nodes must be in the correct 
!>  order to avoid negative values. If the nodes are listed in file "elements.dat"
!>  in the wrong order, the determinant (volume) of the tetrahedron will be negative. 
!>  volume_correct_sign will write check for this order and stop the simulation if 
!>  the order is incorrect
!>
!> @param[in]    element node connectivity
!> @param[out]    
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
! print *,'   the tetra vertices must be reordered for a NEGATIVE volume '
! stop
          itemp=nop(n,1)
          nop(n,1)=nop(n,2)
          nop(n,2)=itemp
         end if

 return
end subroutine volume_correct_sign
