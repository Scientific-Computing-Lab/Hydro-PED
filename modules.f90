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

module sizes
! contains dimesnions and commonly used scaling factors 
  integer:: ne,np,nbp,ke,noff, nprocs
  integer, dimension(3,4):: side = reshape([2,1,4, 1,3,4, 3,2,4, 2,3,1],[3,4])
  real(kind=8):: g(3),adp,den_scale,boff,vp2,tsc,ddm,demf,  &
                 boff_min,boff_max,power
  
end module sizes
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

module element_data
  integer,dimension(:),allocatable:: flag,el_drop,num_drop 
  integer,dimension(:,:),allocatable:: nop,nos
  real(kind=8),dimension(:),allocatable::lambda,mu,gr,dens
  real(kind=8),dimension(:),allocatable::el_vol,al_biot,m_biot,phi,dphi,zi_el
  real(kind=8),dimension(:),allocatable::alpha,dalpha,i1,i2,ksi,ksi0,ksif,pf_el
  real(kind=8),dimension(:,:),allocatable::rate,phi_eq,ductile,coupl
  real(kind=8),dimension(:,:),allocatable::stress,strain,strainp,str_e
  real(kind=8),dimension(:,:),allocatable::stress0,strain0,field
    
end module element_data

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

module node_data
  real(kind=8),dimension(:),allocatable:: mass,pfluid
  real(kind=8),dimension(:,:),allocatable:: balance,cord,disp,vel,dspt
  real(kind=8),dimension(:,:),allocatable:: force  
  
end module node_data
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
module boundary_node_data
  logical,dimension(:,:),allocatable:: vel_code,force_code
  integer,dimension(:),allocatable:: numbn
  real(kind=8),dimension(:,:),allocatable:: value
  
end module boundary_node_data
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
module diffusion_data
  use iso_c_binding
  integer:: mat_size
  integer,dimension(:),allocatable:: id,idt,order,ja
  integer(kind=c_size_t),dimension(:),allocatable :: ia
  integer,dimension(:,:),allocatable:: ija
  real(kind=8),dimension(:),allocatable:: bc_dfval,bc_tval,f,f1,ft,ft1,xsj
  real(kind=c_double), dimension(:), allocatable:: a_matrix
  real(kind=8),dimension(:,:),allocatable:: d,q,ul,tl,p
  real(kind=8),dimension(:,:,:),allocatable:: xl,shpp,s,displ
  
end module diffusion_data
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
