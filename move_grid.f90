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

!--------------- derivations of basic functions --------------
!   element number num  point coordinate inside the element
! dr(k,i);  k=node number (1,2,3,4);  i=coordinate(x,y,z)
!--------------------------------------------------------

#define ALLOC   alloc_if(.true.)
#define FREE    free_if(.true.)
#define RETAIN  free_if(.false.)
#define REUSE   alloc_if(.false.)

subroutine derivation(dr,n)
!dir$ attributes offload:mic :: derivation
!$omp declare simd(derivation)
  use sizes
  use node_data
  use element_data
  implicit none

  integer, intent(in)::n
  integer::i
  real(kind=8):: xn(4),yn(4),zn(4)
  real(kind=8), intent(out):: dr(4,3)
  real(kind=8):: a,b,c,d,e,vol,v1,v2,v3,v4,onv

  ! coord 4 element nodes
  do i = 1,4
    xn(i) = cord(nop(n,i),1)
    yn(i) = cord(nop(n,i),2)
    zn(i) = cord(nop(n,i),3)
  end do

  ! element volume  (G = 6 * volume )
  !| 1 x1 y1 z1 |
  !| 1 x2 y2 z2 |
  !| 1 x3 y3 z3 |
  !| 1 x4 y4 z4 |

  v1 = xn(2) * (yn(3)*zn(4)-zn(3)*yn(4)) &
     - xn(3) * (yn(2)*zn(4)-zn(2)*yn(4)) &
     + xn(4) * (yn(2)*zn(3)-zn(2)*yn(3))

  v2 = xn(1) * (yn(3)*zn(4)-zn(3)*yn(4)) &
     - xn(3) * (yn(1)*zn(4)-zn(1)*yn(4)) &
     + xn(4) * (yn(1)*zn(3)-zn(1)*yn(3))

  v3 = xn(1) * (yn(2)*zn(4)-zn(2)*yn(4)) &
     - xn(2) * (yn(1)*zn(4)-zn(1)*yn(4)) &
     + xn(4) * (yn(1)*zn(2)-zn(1)*yn(2))

  v4 = xn(1) * (yn(2)*zn(3)-zn(2)*yn(3)) &
     - xn(2) * (yn(1)*zn(3)-zn(1)*yn(3)) &
     + xn(3) * (yn(1)*zn(2)-zn(1)*yn(2))

  vol = v1 - v2 + v3 - v4

  if(vol .ge.0.0_8 )then
!    write(6,*)'	Negative volume ',-vol/6,' element ', n
!    write(6,*)-v1,v2,-v3,v4
!    stop
                flag(n) = 3
                dr = 0.
                return
  end if
  onv= 1.0_8/vol

  !A = | 1 y2 z2 |
  !    | 1 y3 z3 |
  !    | 1 y4 z4 |
  a = yn(3)*zn(4)-zn(3)*yn(4) &
    -(yn(2)*zn(4)-zn(2)*yn(4)) &
    + yn(2)*zn(3)-zn(2)*yn(3)

  !B = | 1 x2 z2 |
  !    | 1 x3 z3 |
  !    | 1 x4 z4 |
  b = xn(3)*zn(4)-zn(3)*xn(4)  &
    -(xn(2)*zn(4)-zn(2)*xn(4)) &
    + xn(2)*zn(3)-zn(2)*xn(3)

  !C = | 1 x2 y2 |
  !    | 1 x3 y3 |
  !    | 1 x4 y4 |
  c = xn(3)*yn(4)-yn(3)*xn(4)  &
    -(xn(2)*yn(4)-yn(2)*xn(4)) &
    + xn(2)*yn(3)-yn(2)*xn(3)

  dr(1,1) = -a * onv
  dr(1,2) =  b * onv
  dr(1,3) = -c * onv

  !A = | 1 y1 z1 |
  !    | 1 y3 z3 |
  !    | 1 y4 z4 |
  a = yn(3)*zn(4)-zn(3)*yn(4)  &
    -(yn(1)*zn(4)-zn(1)*yn(4)) &
    + yn(1)*zn(3)-zn(1)*yn(3)

  !B = | 1 x1 z1 |
  !    | 1 x3 z3 |
  !    | 1 x4 z4 |
  b = xn(3)*zn(4)-zn(3)*xn(4)  &
    -(xn(1)*zn(4)-zn(1)*xn(4)) &
    + xn(1)*zn(3)-zn(1)*xn(3)

  !C = | 1 x1 y1 |
  !    | 1 x3 y3 |
  !    | 1 x4 y4 |
  c = xn(3)*yn(4)-yn(3)*xn(4)  &
    -(xn(1)*yn(4)-yn(1)*xn(4)) &
    + xn(1)*yn(3)-yn(1)*xn(3)

  dr(2,1) =  a * onv
  dr(2,2) = -b * onv
  dr(2,3) =  c * onv

  !A = | 1 y1 z1 |
  !    | 1 y2 z2 |
  !    | 1 y4 z4 |
  a = yn(2)*zn(4)-zn(2)*yn(4)  &
    -(yn(1)*zn(4)-zn(1)*yn(4)) &
    + yn(1)*zn(2)-zn(1)*yn(2)

  !B = | 1 x1 z1 |
  !    | 1 x2 z2 |
  !    | 1 x4 z4 |
  b = xn(2)*zn(4)-zn(2)*xn(4)  &
    -(xn(1)*zn(4)-zn(1)*xn(4)) &
    + xn(1)*zn(2)-zn(1)*xn(2)

  !C = | 1 x1 y1 |
  !    | 1 x2 y2 |
  !    | 1 x4 y4 |
  c = xn(2)*yn(4)-yn(2)*xn(4)  &
    -(xn(1)*yn(4)-yn(1)*xn(4)) &
    + xn(1)*yn(2)-yn(1)*xn(2)

  dr(3,1) = -a * onv
  dr(3,2) =  b * onv
  dr(3,3) = -c * onv

  !A = | 1 y1 z1 |
  !    | 1 y2 z2 |
  !    | 1 y3 z3 |
  a = yn(2)*zn(3)-zn(2)*yn(3)  &
    -(yn(1)*zn(3)-zn(1)*yn(3)) &
    + yn(1)*zn(2)-zn(1)*yn(2)

  !B = | 1 x1 z1 |
  !    | 1 x2 z2 |
  !    | 1 x3 z3 |
  b = xn(2)*zn(3)-zn(2)*xn(3)  &
    -(xn(1)*zn(3)-zn(1)*xn(3)) &
    + xn(1)*zn(2)-zn(1)*xn(2)

  !C = | 1 x1 y1 |
  !    | 1 x2 y2 |
  !    | 1 x3 y3 |
  c = xn(2)*yn(3)-yn(2)*xn(3)  &
    -(xn(1)*yn(3)-yn(1)*xn(3)) &
    + xn(1)*yn(2)-yn(1)*xn(2)

  dr(4,1) =  a * onv
  dr(4,2) = -b * onv
  dr(4,3) =  c * onv

  return
end subroutine derivation
!Move grid along veocity field (Lagrangian) and new total strain tensor  
subroutine move_grid(dt)
  
  use sizes
  use node_data
  use element_data
  use omp_lib

  implicit none
  
  integer:: ii,j,i,n,nn
  real(kind=8)::dt,de(6),dr(4,3)
  integer :: time_start, time_end,time_mic_start,time_mic_end,time_cpu_start,time_cpu_end
 
  call system_clock(time_start) 

  !dir$ offload begin target(mic:0) wait(vel_signal) nocopy(cord : REUSE RETAIN) out(strain(1:6,1:phi_length): REUSE RETAIN) nocopy(vel : REUSE RETAIN) nocopy(m_biot : REUSE RETAIN) nocopy(nop : REUSE RETAIN) nocopy(de,dr,i,j,ii,nn,n) in(dt,np,phi_length) signal(strain_signal)
   call move_grid_inline(1,phi_length)
!dir$ end offload
  
  !dir$ offload begin target(mic:1) wait(vel_signal) nocopy(cord : REUSE RETAIN) out(strain(1:6,phi_length+1:2*phi_length) : REUSE RETAIN) nocopy(vel : REUSE RETAIN) nocopy(m_biot : REUSE RETAIN) nocopy(nop : REUSE RETAIN) nocopy(de,dr,i,j,ii,nn,n) in(dt,np,phi_length) signal(strain_signal)
   call move_grid_inline(phi_length+1,2*phi_length)
!dir$ end offload

  call system_clock(time_cpu_start) 
  ! moving of mesh
  call move_grid_inline(2*phi_length+1,ne)
  call system_clock(time_cpu_end) 

  call system_clock(time_end) 
    
  write (*,*) "n=",n,"time=",(time_end-time_start), "cpu_time=", (time_cpu_end-time_cpu_start)

return
contains
   subroutine move_grid_inline(offset_begin,offset_end)
     !dir$ attributes forceinline :: move_grid_inline
     !dir$ attributes offload:mic :: move_grid_inline
     integer :: offset_begin,offset_end
  
     !$OMP PARALLEL
     do j=1,3
     !$OMP DO SIMD
       do i=1,np
        cord(i,j) = cord(i,j) + vel(i,j)*dt
#ifndef __MIC__
        disp(i,j) = disp(i,j) + vel(i,j)*dt
        dspt(i,j) = dspt(i,j) + vel(i,j)*dt
#endif
       end do
     end do
          
     !$OMP DO SIMD PRIVATE(dr,ii,de,j,nn)
     do n = (offset_begin),(offset_end)
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
     !$OMP END PARALLEL
   end subroutine
end subroutine move_grid

!===============================================================================
!Move grid and update fluid pressure in elements  
subroutine move_grid_pf(dt)
  use sizes
  use node_data
  use element_data
  use omp_lib

  implicit none
  
  integer:: ii,j,i,n,nn
  real(kind=8)::dt,de(6),dr(4,3)
  integer :: time_start, time_end,time_mic_start,time_mic_end,time_cpu_start,time_cpu_end
  
  call system_clock(time_start) 

  !dir$ offload begin target(mic:0) wait(vel_signal) nocopy(cord : REUSE RETAIN) out(strain(1:6,1:phi_length): REUSE RETAIN) out(pf_el(1:phi_length): REUSE RETAIN) nocopy(vel : REUSE RETAIN) nocopy(m_biot : REUSE RETAIN) nocopy(nop : REUSE RETAIN) nocopy(de,dr,i,j,ii,nn,n) in(dt,np,phi_length) signal(strain_signal)
  call move_grid_pf_inline(1,phi_length)
  !dir$ end offload

  !dir$ offload begin target(mic:1) wait(vel_signal) nocopy(cord : REUSE RETAIN) out(strain(1:6,phi_length+1:2*phi_length) : REUSE RETAIN) out(pf_el(phi_length+1:2*phi_length) : REUSE RETAIN) nocopy(vel : REUSE RETAIN) nocopy(m_biot : REUSE RETAIN) nocopy(nop : REUSE RETAIN) nocopy(de,dr,i,j,ii,nn,n) in(dt,np,phi_length) signal(strain_signal)
  call move_grid_pf_inline(phi_length+1,2*phi_length)
  !dir$ end offload
  
! moving of mesh
  call system_clock(time_cpu_start) 
  call move_grid_pf_inline(2*phi_length+1,ne)
  call system_clock(time_cpu_end) 

  call system_clock(time_end) 
 
  write (*,*) "n=",n,"time=",(time_end-time_start), "cpu_time=", (time_cpu_end-time_cpu_start)
  
return
contains
   subroutine move_grid_pf_inline(offset_begin,offset_end)
     !dir$ attributes forceinline :: move_grid_pf_inline
     !dir$ attributes offload:mic :: move_grid_pf_inline
     integer :: offset_begin,offset_end
  
     !$OMP PARALLEL
     do j=1,3
     !$OMP DO SIMD
     do i=1,np
       cord(i,j) = cord(i,j) + vel(i,j)*dt
#ifndef __MIC__
       disp(i,j) = disp(i,j) + vel(i,j)*dt
       dspt(i,j) = dspt(i,j) + vel(i,j)*dt
#endif
     end do
     end do

     ! New elastic strain tensor
     !$OMP DO SIMD PRIVATE(ii,dr,de,j,nn)
     do n=(offset_begin), (offset_end)
!    derivation of basic functions
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

!      New total strains
       do j=1,6
         strain(j,n) = strain(j,n) + de(j)*dt
       end do

!    update fluid pressure
        pf_el(n) = pf_el(n)-m_biot(n)*(de(1)+de(2)+de(3))*dt
     end do
     !$OMP END PARALLEL
   end subroutine
end subroutine move_grid_pf

