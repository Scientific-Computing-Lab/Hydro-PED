! Copyright (c) 2017, 
! Eyal Shalev (eyal@gsi.gov.il)
! Vladimir Lyakhovsky
! All rights reserved to:
! Geological Survey of Israel (GSI) &
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!    * Redistributions of source code must retain the above copyright
!      notice, this list of conditions and the following disclaimer.
!    * Redistributions in binary form must reproduce the above copyright
!      notice, this list of conditions and the following disclaimer in the
!      documentation and/or other materials provided with the distribution.
!    * Neither the name of Eyal Shalev or Vladimir Lyakhovsky, nor the
!      names of its contributors may be used to endorse or promote
!      products
!      derived from this software without specific prior written
!      permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND
! ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
! WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL Eyal Shalev, Vladimir Lyakhovsky
! BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
! OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
! WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
! OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
! ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

!	Convert output data files to TECPLOT ASCII file
!	zoom out the deformation
        parameter ( ne   =      715345,   &
                    np   =      122307,   &
                    nbp  =        8161,   &
		    ke   =          47,   &
		    mat_size =  4)
!-------------------------------------------------------------
!       ne  - number of elements
!	np  - number of nodes
!	nbp - number of boundary nodes
!	ke  - max number of elenents around the node+1
!-------------------------------------------------------------
	real  cord(np,3), ux(np),uy(np),uz(np),pfluid(np)
	integer nop(ne,4), mater(ne),flag(ne)
	real alpha(ne),taue(ne),pe(ne),ksie(ne),per(ne)

	character*14 file1/'    .nodes.dat'/
	character*4 step1
	character*1 st1(4)
	equivalence (step1,file1),(step1,st1(1))

	character*14 file2/'    .tetra.dat'/
	character*4 step2
	equivalence (step2,file2)

	character*16 file3/'    .tecplot.dat'/
	character*4 step3
	equivalence (step3,file3)
!--------------------------------------------------------
	character*4 st_arg
	call getarg(1,st_arg)
		if ( st_arg .eq. ' ' ) then
	print *,'   === Step Number Not Given ==='
		stop
		end if

	read(st_arg,*)  nstep
!--------------------------------------------------------
!-----	READ element connections and put properties ----
	open(1,file='elements.dat')
	read(1,*) nlm

	if ( nlm .ne. ne ) then
	print *,'  Number of elements is not correct ',nlm,ne
	stop
	endif

		do i = 1,ne

	read(1,*) nop(i,1),nop(i,2),nop(i,3),nop(i,4)!, &
!			mater(i)
	mater(i) = 1
		end do

!-----	end READ element connections ----
	close (1)

	print *,'end READ element connections'

!------------------------------------------------
!	generate file-names for input & output

	print 441
441	format('	Step number:  ',$)
!	read(*,*) nstep
	print *,nstep

        write(step1,5) nstep

5       format(i4)

        if (st1(1) .eq. ' ' ) st1(1) = '0'
        if (st1(2) .eq. ' ' ) st1(2) = '0'
        if (st1(3) .eq. ' ' ) st1(3) = '0'

		step2 = step1
		step3 = step1

!------------------------------------------------
! 	input actual node coordinates 
!	node_number  X  Y  Z

        open (3,file=file1)

	do n = 1,np
	read(3,*) ntmp,cord(n,1),cord(n,2),cord(n,3),	&
	ux(n),uy(n),uz(n),pfluid(n)
	end do

		close(3)

! 	node input finished
	print *,'	end READ node coordinates from step ',nstep
!------------------------------------------------
!
!	do n = 1,np
!	
!
!	ux(n) = ux(n) * 1000.
!	uy(n) = uy(n) * 1000.
!	uz(n) = uz(n) * 1000.
!
!	end do
!------------------------------------------------
!	input alpha & alpha_dot

        open(3,file=file2)

	do n = 1,ne
	read(3,*) ntmp,alpha(n),		&
                 s1,s2,s3,s12,s13,s23,		&
                 e1,e2,e3,e12,e13,e23,flag(n)

	pe(n) = -(s1+s2+s3)/3.
	taue(n) = sqrt(((s1+pe(n))**2 + (s2+pe(n))**2 + (s3+pe(n))**2)/2. &
                     + s12**2 + s13**2 + s23**2 )

        per(n) = (e1+e2+e3)/   &
	       sqrt(e1**2+e2**2+e3**2+2.*e12**2+2*e13**2+2*e23**2)

	per(n) = (per(n)+0.7)*(e1**2+e2**2+e3**2+2.*e12**2+2*e13**2+2*e23**2)

	end do

	close(3)

!------------------------------------------------
	
	open(1,file=file3)

	write(1,111) nstep/15.0
	write(1,*) 'VARIABLES = "X", "Y", "Z", "ux", "uy", "uz", ',  &
      '  "P_fluid", "alpha", "Pres", "Shear", "flag" '
	write(1,*) 'ZONE T="Map", N=',np,' E=',ne,	  	&
      ' DATAPACKING=BLOCK ZONETYPE=FETETRAHEDRON',		&
      ' VARLOCATION=([8-12]=CELLCENTERED)'

	write(1,55) ((cord(n,nv), n=1,np), nv=1,3)
	write(1,55) (ux(n), n=1,np)
	write(1,55) (uy(n), n=1,np)
	write(1,55) (uz(n), n=1,np)
	write(1,55) (pfluid(n), n=1,np)
	write(1,55) (alpha(n),  n=1,ne)
	write(1,55) (per(n),  n=1,ne)
	write(1,55) (taue(n), n=1,ne)
	write(1,55) (flag(n), n=1,ne)
	
	
55	format(7(1x,g15.5))
111     format('TITLE =  "',f5.2,'"')

	do n = 1,ne
	write(1,*) nop(n,1), nop(n,2), nop(n,3), nop(n,4)
	end do

	close(1)
!------------------------------------------------

	print *,'	TECPLOT ASCII file is ready '
!------------------------------------------------
       stop
       end

!---------------------------------------------------------

	subroutine comment
!	skip the line starts with # - comment
! if the line starts with another symbol ---> stop

	character*1 sym

	backspace 1

	read(1,1) sym
1	format(a1)

	if ( sym .eq. '#' ) return

	print *,'	Error reading input file '

	stop
	end
