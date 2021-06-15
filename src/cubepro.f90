! This file is part of dipro.
! SPDX-Identifier: Apache-2.0
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

!> Implementation of the dimer projection method for extended tight binding methods.
module cubepro
   use mctc_env, only : error_type, fatal_error, get_argument, wp
   use mctc_io, only : structure_type, new
   use mctc_io_convert, only : autoev
   use cubepro_xtb, only : get_calculator
   use tblite_basis_type, only : get_cutoff, basis_type
   use tblite_blas, only : dot, gemv, gemm
   use tblite_context_type, only : context_type
   use tblite_cutoff, only : get_lattice_points
   use tblite_integral_overlap, only : get_overlap
   use tblite_output_format, only : format_string
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction, &
      & get_density_matrix
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_singlepoint, only : xtb_singlepoint
   implicit none
   private

   public :: get_cube, cube_input

   !> Configuration data for calculation
   type :: cube_input
      !> Name of the requested tight binding method
      character(len=:), allocatable :: method
      !> List of fragments, generated if not given here
      !integer, allocatable :: fragment(:)
      !> Threshold for generating fragments from bond orders
      !real(wp) :: thr = 0.1_wp
      !> Calculation accuracy
      !real(wp) :: accuracy = 1.0_wp
      !> Output verbosity
      integer :: verbosity = 2
      !> Electronic temperature in Kelvin
      real(wp) :: etemp = 300.0_wp
   end type cube_input

   !> Conversion factor from temperature to energy (Boltzmann's constant in atomic units)
   real(wp), parameter :: ktoau = 3.166808578545117e-06_wp

contains

!> Entry point for generation of the cube file
subroutine get_cube(input, mol, error)
   !> Input data for the calculation
   type(cube_input), intent(in) :: input
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: spin, charge, stat, unit, ifr, nfrag, nao
   logical :: exist
   !real(wp) :: energy, cutoff, jab, sab, jeff
   !real(wp), allocatable :: gradient(:, :), sigma(:, :)
   type(context_type) :: ctx
   type(xtb_calculator) :: calc, fcalc
   type(structure_type), allocatable :: mfrag(:)
   type(wavefunction_type) :: wfn
   type(wavefunction_type), allocatable :: wfx(:)
   character*(*) fname

   call get_calculator(calc, mol, input%method, error)
   if (allocated(error)) return

   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, &
      & input%etemp * ktoau)

   call xtb_singlepoint(ctx, mol, calc, wfn, input%accuracy, energy, &
      & verbosity=input%verbosity-1)
   if (ctx%failed()) then
      call ctx%get_error(error)
      return
   end if

   fname='density.cube'
   call cube(mol,wfn,fname,basis)

end subroutine get_cube


! Taken from xtb and modified
subroutine cube(mol,wfn,fname,basis)
   implicit none
   type(structure_type), intent(in) :: mol
   type(basis_type), intent(in) :: basis
   type(wavefunction_type), intent(in) :: wfn
   character*(*), intent(in) :: fname

   !real*8, intent ( in ) :: eval(nmo)
   !real*8, intent ( in ) :: occ (nmo)
   !real*8, intent ( in ) :: cmo(nbf,nmo)
   !integer, intent( in ) :: at(n)
   !integer, intent( in ) :: n,nmo,nbf
   integer :: n
   real*8 :: xyz(3,n)

   real*4 ,allocatable  ::cb   (:,:,:)
   real*8 ,allocatable  ::array(:)

   integer i,j,k,m,nm,ii,jj,iii,jjj,npri,nprj,primstart(nbf),nexp
   integer iat,jat,xst,yst,zst,cou,iprimcount,jprimcount,iiii,jjjj
   real*8 thr,thr3,intcut,intcut2
   real*8 w0,w1,t0,t1,r,r2,val,rab,ccc,gridp(3),xyza(3),xyzb(3),dr3
   real*8 step,px,py,pz,xinc,yinc,zinc,nx,ny,nz,vv(1),v,est,nfod
   real*8 f1,f2,dx1,dy1,dz1,dx2,dy2,dz2,r1,dum,r1xy,r2xy
   real*8 dxx1,dyy1,dzz1,dxx2,dyy2,dzz2,ar1,ar2,cc
   integer ifile


   n=mol%nat
   xyz=mol%xyz

!  grid spacing for cube file
   cube_step = 0.4_wp
!  density matrix neglect threshold
   cube_pthr = 0.05_wp
   
   write(*,*)
   write(*,*)'cube file module (SG, 7/16)'
   thr = cube_pthr ! Dmat pre-screen
   step= cube_step ! grid step (Bohr)
   intcut=8.00d0   ! primitive cut
   intcut2=2.0d0*intcut

!  auxiliary array for fast density evaluation
   nexp=100*int(intcut2) ! size of exp array ie arguments 0:nexp
   allocate(array(0:nexp))
   do i=0,nexp
      dum=float(i)/100.0
      array(i)=exp(-dum)
   enddo

   write(*,'('' cube_pthr     : '',f7.3)')cube_pthr
   write(*,'('' cube_step     : '',f7.3)')cube_step
   write(*,'('' non-zero P (%): '',f7.3,''   nmat:'',i8)') &
   & 100.*float(nm)/float(nbf*(nbf+1)/2),nm

   px=maxval(xyz(1,1:n))+3.0
   py=maxval(xyz(2,1:n))+3.0
   pz=maxval(xyz(3,1:n))+3.0
   nx=minval(xyz(1,1:n))-3.0
   ny=minval(xyz(2,1:n))-3.0
   nz=minval(xyz(3,1:n))-3.0
   write(*,*)'Grid Boundaries (x y z) :'
   write(*,*)px,py,pz
   write(*,*)nx,ny,nz

   ! calculate step size and number of steps (step size approx 0.2-0.5)
   ! Alternatively one could choose x=y=z=40-120
   ! x step number
   xst=floor((abs(px)+abs(nx))/step)
   ! x step size
   xinc=(abs(px)+abs(nx))/xst
   ! y step number
   yst=floor((abs(py)+abs(ny))/step)
   ! y step size
   yinc=(abs(py)+abs(ny))/yst
   ! z step number
   zst=floor((abs(pz)+abs(nz))/step)
   ! z step size
   zinc=(abs(pz)+abs(nz))/zst
   write(*,*)'Total # of points', (xst+1)*(yst+1)*(zst+1)
   
   ! Volume for each gridpoint
   dr3=xinc*yinc*zinc

   write(*,*)'writing ',trim(fname)
   call open_file(ifile,fname,'w')
   write(ifile,*)'xtb spin/fo density'
   write(ifile,*)'By S.G.'
   write(ifile,101)n,nx,ny,nz
   write(ifile,101)xst+1,xinc,0.0,0.0
   write(ifile,101)yst+1,0.0,yinc,0.0
   write(ifile,101)zst+1,0.0,0.0,zinc
   do i=1,n
      write(ifile,102)mol%id(i),0.0,xyz(1,i),xyz(2,i),xyz(3,i)
   enddo

   allocate(cb(0:zst,0:yst,0:xst))

   cb  =0

   !     Dmat loop         -----------------------------------
   ! loop over atoms 
   do iat=1, n
   do jat=1, n
      iso = ish_at(iat)
      jso = jsh_at(jat)
      izp=mol%id(iat)
      jzp=mol%id(jat)
      xyza(1:3)=xyz(1:3,iat)
      xyzb(1:3)=xyz(1:3,jat)
      rab=(xyza(1)-xyzb(1))**2  &
      &      +(xyza(2)-xyzb(2))**2  &
      &      +(xyza(3)-xyzb(3))**2
      ! loop over shells
      do ish=1, basis%nsh_at(iat)
      do jsh=1, basis%nsh_at(jat)
         ii = iao_sh(iso + ish)
         jj = iao_sh(jso + jsh)
         icgto = basis%cgto(ish,izp)
         jcgto = basis%cgto(jsh,jzp)
         ! loop over shellblock NAOs
         do iao = 1, nao_sh(iso + ish)
         do jao = 1, nao_sh(jso + jsh)
            pij = wfn%density(ii+iao,jj+jao)
            ! loop over primitives
            pcoef = 0.0_wp
            do npri = 1 , icgto%nprim 
            do nprj = 1 , jcgto%nprim
               pcoef = pcoef+ icgto%coeff(npri) * jcgto%nprim(nprj)
            end do
            end do
            pcoef = pcoef * pij
               ! grid loops
               gridp(1)=nx
               do i=0,xst
                  gridp(1)=nx+(xinc*i)
                  dx1=xyza(1)-gridp(1)
                  dx2=xyzb(1)-gridp(1)
                  dxx1=dx1*dx1
                  dxx2=dx2*dx2
                  gridp(2)=ny
                  do j=0,yst
                     gridp(2)=ny+(yinc*j)
                     dy1=xyza(2)-gridp(2)
                     dy2=xyzb(2)-gridp(2)
                     dyy1=dy1*dy1
                     dyy2=dy2*dy2
                     gridp(3)=nz
                     r1xy=dxx1+dyy1
                     r2xy=dxx2+dyy2
                     do k=0,zst
                        gridp(3)=nz+(zinc*k)
                        dz1=xyza(3)-gridp(3)
                        dz2=xyzb(3)-gridp(3)
                        dzz1=dz1*dz1
                        dzz2=dz2*dz2
                        r1=r1xy+dzz1
                        ! primitive function value, needs to distuinguish between different 
                        !ar1=icgto%alpha(npri)*r1
                        !if(ar1.lt.intcut2)then  ! exp(-16) < 1.d-7 i.e. zero
                           call primvalf(dx1,dy1,dz1,dxx1,dyy1,dzz1,ar1, &
                              &             icgto%ang,nexp,array,f1)
                        r2=r2xy+dzz2
                        !ar2=jcgto%alpha(nrpj)*r2
                        !   if(ar2.lt.intcut2)then
                              call primvalf(dx2,dy2,dz2,dxx2,dyy2,dzz2,ar2, &
                                 &             jcgto%ang,nexp,array,f2)
                              cb(k,j,i)=cb(k,j,i)+pcoef*f1*f2
                           endif
                        endif
                     enddo
                  enddo
            enddo
            endif
         enddo
      enddo
   enddo

   !     Dmat loop end     -----------------------------------
   
   ! write x y z x y z\n
   cou=1
   do i=0,xst
      do j=0,yst
         do k=0,zst
            if (cou.lt.6) then
               write(ifile,'(E14.8,1X)',advance='no')cb(k,j,i)
               cou=cou+1
            else
               write(ifile,'(E14.8)')cb(k,j,i)
               cou=1
            endif
         enddo
      enddo
   enddo
   call close_file(ifile)

   101   format(I5,3F16.6)
   102   format(I5,4F16.6)
 
 end subroutine cube


subroutine primval(dx,dy,dz,dx2,dy2,dz2,alpr2,lao,f)
   implicit none
   integer lao
   real*8 f,dx2,dy2,dz2
   real*8 dx,dy,dz,alpr2

   goto (100,201,202,203,301,302,303,304,305,306) lao

   100  f=exp(-alpr2)
   return
   201  f=exp(-alpr2)*dx
   return
   202  f=exp(-alpr2)*dy
   return
   203  f=exp(-alpr2)*dz
   return
   301  f=exp(-alpr2)*dx2
   return
   302  f=exp(-alpr2)*dy2
   return
   303  f=exp(-alpr2)*dz2
   return
   304  f=exp(-alpr2)*dx*dy
   return
   305  f=exp(-alpr2)*dx*dz
   return
   306  f=exp(-alpr2)*dy*dz
   return

end subroutine primval

! routine with exp table lookup
subroutine primvalf(dx,dy,dz,dx2,dy2,dz2,alpr2,lao,nexp,a,f)
   implicit none
   integer lao,nexp
   real*8 a(0:nexp)
   real*8 f,dx2,dy2,dz2
   real*8 dx,dy,dz,alpr2

   select case (lao)

   case(0)
     f=fastexp(nexp,a,alpr2)
   case(1)  
     f=fastexp(nexp,a,alpr2)*dx
     f=f+fastexp(nexp,a,alpr2)*dy
     f=f+fastexp(nexp,a,alpr2)*dz
   case(2)
     f=fastexp(nexp,a,alpr2)*dx2
     f=f+fastexp(nexp,a,alpr2)*dy2
     f=f+fastexp(nexp,a,alpr2)*dz2
     f=f+fastexp(nexp,a,alpr2)*dx*dy
     f=f+fastexp(nexp,a,alpr2)*dx*dz
     f=f+fastexp(nexp,a,alpr2)*dy*dz
   end select

end subroutine primvalf

real*8 function fastexp(nexp,array,arg)
   integer nexp
   real*8 arg
   real*8 array(0:nexp)

   fastexp=array(int(100.0d0*arg))

end function fastexp

end module cubepro