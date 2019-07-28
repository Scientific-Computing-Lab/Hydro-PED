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
!  **************************************************************************

      !> @brief         Modify the stiffness matrix to account for boundary conditions
      !> @details       
      !>
      !> @param[in]     Boundary conditions (id) and global stiffnes matrix (a_matrix, ia, ja)  
      !> @param[out]    Updated matrix.
      subroutine modify2

  use sizes
  use diffusion_data
  
      implicit double precision (a-h,o-z)
	integer row

!---- modify rhs to account for essential boundary conditions

      do n = 1,np 
         if (id(n).eq.1) then
	    do i=ia(n),ia(n+1)-1
	       row=ja(i)
		if (n .eq. row) then
		   a_matrix(i)=1.d0
 		else
 		   a_matrix(i)=0.d0
		end if
	    end do
	 else
	    do i=ia(n),ia(n+1)-1
	       row=ja(i)
	       if (id(row) .eq. 1) then
	          f(n)=f(n)-a_matrix(i)*f1(row)
		  a_matrix(i)=0.d0
		end if
	     end do      
	 end if
      end do


      return
      end subroutine modify2

