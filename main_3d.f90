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
	
program main_3d
!	Damage 3D + fluid flow converted to F90 + OMP
!	 No double overlap as in original EFDLM
!	Input files:
!		cordinates.dat - node coordimnates and index
!		elememts.dat - element connections andmaterial index
!		material.dat - material properties 
  use sizes
  use element_data
  use node_data
  use boundary_node_data

  implicit none

  logical:: event
  integer:: nloop,nout,n,j
  real(kind=8)::dt,damax,dtp
  real(kind=8):: time,t_start,t_stop,t_delta,tout,	&
  		 timedf0,timedf1,dt_diff

!-------------------------------------------------------------
!       ne   	   - number of elements
!	np   	   - number of nodes
!	nbp        - number of boundary nodes
!	ke   - max number of elenents around the node+1
!
!-------------------------------------------------------------
!
!
!	D	- material properties for diffusion
!    D(4,i)=Ss [1/Pa]   pressure storage 
!    D(5,i)=porosity[]
!    D(6,i)=k11[m2]
!    D(7,i)=k12=k21[m2]
!    D(8,i)=k22[m2]
!    D(9,i)=viscosity [Pa sec]
!    D(10,i)=density[kg/m3]
!    D(11,i)=water bulk mudulus
!-------------------------------------------------------------

           open(4,file='general.dat')
        read(4,*) demf,adp
        read(4,*) boff_min, boff_max
        read(4,*)  g
        read(4,*) t_start,t_stop,t_delta
            close (4)

		 tout = t_start
		 
!------- Remove previous event catalog -------------
        call system('/bin/rm catalog')


!------- input  model geometry,
!------- material properties and
!------- initial conditions -------------------
   print *,'   INPUT  '
      	call input
      	
!!!    call analyze87(np,ia,ja,order,keep,control,info)

!------- time step for elasticity -----------------
		den_scale = 1.e+4
                dt = 0.25_8 / sqrt(vp2 / den_scale)
!		print *,' Sound speed ',sqrt(vp2),'    Delta_t = ',dt
        print *,' Initial time step dt = ',dt,'     Density scale = ',den_scale


!	--------  Massa ----------------
	call node_mass

	    dalpha = 0.0_8
!--------- Initial disp=0 -----------------------
	   dspt = 0.0_8 
           disp = 0.0_8

		dt_diff = 3600.0_8
!!!	 call diffusion(dt_diff,keep,control,info)
!!!  open(30,file='flux.dat')
!!!	 call fluxbc(time)

!dir$ offload_transfer target(mic:0) signal(strain_signal)
!dir$ offload_transfer target(mic:1) signal(strain_signal)
         call output

!****      ------LOOP OF TIME FOR THE INITIAL DISTRIBUTION---------
                time = 0.0_8
		nloop = 0
       do while (boff.ge.0.1*boff_min .or. nloop.le.10)
!       do while (nloop.le.200)
       		nloop = nloop + 1
     
          print 600,nloop,time
600	format(I7,'''s ',g12.4,$)
		time = time + dt*tsc

!      ---------  EFDLM ------------------------- 
	call efdlm (dt,event)

        if ( event ) then
        print *,' ----------------------------------------------'
        print *,'       The initial damage is too large  '
        print *,' Event occures before you get equlibration '
        print *,'     Change initial damage distribution '
        print *,'         and RUN the code again '
        stop
        endif

!------------------ Move Grid ------------------
	call move_grid (dt)

	end do
!      ------ END of loop for the initial distribution ----
	print *,' ------ END of initial loop ---- '

 	call output
	   dspt = 0.0_8 
           disp = 0.0_8
	
		time = 0.
		timedf1 = 10.0_8
		timedf0 = 0.
		nloop = 0
    do while (time .le. t_stop)
         nloop = nloop + 1
	time = time + dt*tsc
	
!      ---------  EFDLM ------------------------- 
          print 600,nloop,time
	call efdlm (dt,event)

!-------------- Stress drop ----------------------
	if ( event ) then

                den_scale = 1.0_8
                dt = 0.25_8 / sqrt(vp2 / den_scale)

	call drop (time,dt)
		dt_diff = 10.0_8
!!!	 call diffusion(dt_diff,keep,control,info)
!	 	 call fluxbc(time)


		timedf0 = time
		timedf1 = time + 10.0_8
	
               event = .false.

	else

!------------------ Move Grid ------------------
	call move_grid (dt)

!------------------ Damage evolution ------------
	call damage(dt,damax)
     
!------------------ Damage-related viscosity --------
 	call plastic(dt)

! ------------------- Diffusion -------------------
	if(time .ge. timedf1 ) then
		dt_diff = (time-timedf0)

!!!	 call diffusion(dt_diff,keep,control,info)
!!!	 	 call fluxbc(time)

		timedf0 = time
		timedf1 = time + 10.0_8
	end if
!-------- change time step ------------------------
        if ( boff .ge. boff_max .or. damax .ge. 0.01_8) then
                den_scale = den_scale*(1.0_8-adp)
        if ( den_scale .lt. 1.0_8 ) den_scale = 1.0_8
                dt = 0.25_8 / sqrt(vp2 / den_scale)

        else if ( boff .le. boff_min ) then
                den_scale = den_scale*(1.0_8+adp)
                dt = 0.25_8 / sqrt(vp2 / den_scale)

        end if

! -------------------  OUTPUT ----------------------
        if ( time .ge. tout   ) then

!dir$ offload_transfer target(mic:0) signal(strain_signal)
!dir$ offload_transfer target(mic:1) signal(strain_signal)
	call output

	tout = tout + t_delta
                endif
   end if
 		end do
!   		close(30)
       stop
 end program main_3d
