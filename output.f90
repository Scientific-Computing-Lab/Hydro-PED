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

!----------- Output data for VIEWER-----------------

	subroutine output

       use sizes
       use diffusion_data
       use node_data
       use element_data

	character(len=14) file1/'    .nodes.dat'/
	character(len=14) file2/'    .tetra.dat'/
	character(len=4) step1,step2
	character(len=1) st1(4),st2(4)
	equivalence (step1,file1),(step1,st1(1))
	equivalence (step2,file2),(step2,st2(1))

	save nstep
	integer nstep/0/

!	time-step number

        nstep = nstep + 1

!------------------------------------------------
!	generate file-names for output

!..        encode(3,5,step) nstep
        write(step1,'(I4)') nstep
        write(step2,'(I4)') nstep

!       format(i4)

        if (st1(1) .eq. ' ' ) st1(1) = '0'
        if (st1(2) .eq. ' ' ) st1(2) = '0'
        if (st1(3) .eq. ' ' ) st1(3) = '0'

        if (st2(1) .eq. ' ' ) st2(1) = '0'
        if (st2(2) .eq. ' ' ) st2(2) = '0'
        if (st2(3) .eq. ' ' ) st2(3) = '0'

!------------------------------------------------
! 	output node coordinates & displacements
!	node_number  X  Y  Z

        open (3,file=file1)

	do n = 1,np
	write(3,25) n,				&
      		cord(n,1),cord(n,2),cord(n,3),	&
     		dspt(n,1),dspt(n,2),dspt(n,3),	&
     		pfluid(n)
	end do
25	format(i9, 7(1x,g15.7))

		close(3)

! 	node output finished
!------------------------------------------------
! 	output element properties

        open (3,file=file2)

       !dir$ offload_transfer target(mic:0) wait(strain_signal) nocopy(strain(1:6,1:phi_length) : REUSE RETAIN)
       !dir$ offload_transfer target(mic:1) wait(strain_signal) nocopy(strain(1:6,phi_length+1:2*phi_length) : REUSE RETAIN)
	do m = 1,ne

!	convert strain from real*8 to real*4
	exx = strain(1,m)-strainp(1,m)
	eyy = strain(2,m)-strainp(2,m)
	ezz = strain(3,m)-strainp(3,m)
	exy = strain(4,m)-strainp(4,m)
	exz = strain(5,m)-strainp(5,m)
	eyz = strain(6,m)-strainp(6,m)

	write(3,30) m,alpha(m),		&
     	(stress(k,m), k=1,6),		&
     	exx,eyy,ezz,exy,exz,eyz,	&
	flag(m)

!     	d(6,m),q(1,m),q(2,m),q(3,m)

30	format(i9,13(1x,g15.5),i3)
	enddo

	write(3,*) 
	
		close(3)

	return
	end
