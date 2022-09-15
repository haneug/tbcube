! This file is part of cubepro.
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
   use tblite_basis_type, only : get_cutoff, basis_type, cgto_type, new_basis
   use tblite_blas, only : dot, gemv, gemm
   use tblite_context_type, only : context_type
   use tblite_cutoff, only : get_lattice_points
   use tblite_integral_overlap, only : get_overlap
   use tblite_output_format, only : format_string
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction, &
      & get_density_matrix
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_singlepoint, only : xtb_singlepoint
   use tblite_data_spin, only : get_spin_constant, read_spin_constants
   use tblite_spin, only : spin_polarization, new_spin_polarization
   use tblite_container, only : container_type, container_cache
   implicit none
   private

   public :: get_cube, cube_input

   !> Configuration data for calculation
   type :: cube_input
      !> Name of the requested tight binding method
      character(len=:), allocatable :: method
      !> Calculation accuracy
      real(wp) :: accuracy = 1.0_wp
      !> Output verbosity
      integer :: verbosity = 2
      !> Electronic temperature in Kelvin
      real(wp) :: etemp = 300.0_wp
      !> Density mode requested
      logical :: sdens = .false.
      !> Number of spin channels
      integer :: nspin = 1 
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

   integer :: spin, charge
   logical :: exist
   type(context_type) :: ctx
   type(xtb_calculator) :: calc, fcalc
   type(structure_type), allocatable :: mfrag(:)
   type(wavefunction_type) :: wfn
   type(wavefunction_type), allocatable :: wfx(:)
   character(len=:), allocatable :: fname
   character(len=:), allocatable :: sp_input
   real(wp) :: energy
   real(wp) :: kt
   integer :: nspin
   logical :: sdens 
   
   !> Spin polarization 
   class(container_type), allocatable :: cont
   type(spin_polarization), allocatable :: spinp
   real(wp), allocatable :: wll(:, :, :)


   allocate(spinp)

   kt=input%etemp * ktoau
   sdens=input%sdens
   nspin=input%nspin

   call get_calculator(calc, mol, input%method, error)
   if (allocated(error)) return

   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, &
      & nspin, kt)

   !> Read external spin constants
   sp_input = "spin_param.txt"
   inquire(file=sp_input, exist=exist)
   if (exist) then
     call read_spin_constants(sp_input)
   end if

   ! if (config%spin_polarized_input) then
   !   call read_spin_constants(config%sp_input)
   ! end if
   call get_spin_constants(wll, mol, calc%bas)
   call new_spin_polarization(spinp, mol, wll, calc%bas%nsh_id)
   call move_alloc(spinp, cont)
   call calc%push_back(cont)

   call xtb_singlepoint(ctx, mol, calc, wfn, input%accuracy, energy, &
      & verbosity=input%verbosity-1)
   if (ctx%failed()) then
      call ctx%get_error(error)
      return
   end if

   if (sdens) then
     fname='spindensity.cube'
   else
     fname='density.cube'
   endif
   
   call cube(mol,wfn,fname,calc%bas,sdens)
   
   fname='homo.cube'
   call mocube(mol,wfn,fname,calc%bas)

end subroutine get_cube


subroutine cube(mol,wfn,fname,basis,sdens)
   implicit none
   type(structure_type), intent(in) :: mol
   type(basis_type), intent(in) :: basis
   type(wavefunction_type), intent(in) :: wfn
   character*(*), intent(in) :: fname
   logical, intent(in) :: sdens

   integer :: n
   real*8 :: xyz(3,mol%nat)
   type(cgto_type) :: icgto, jcgto

   real*4 ,allocatable  ::cb   (:,:,:)
   real*8 ,allocatable  ::array(:)

   integer i,j,k,ii,jj,npri,nprj,nexp
   integer ish,jsh,izp,jzp,iao,jao,iso,jso
   integer iat,jat,xst,yst,zst,cou,iprimcount,jprimcount
   real*8 thr,thr3,intcut,intcut2,dist_cut
   real*8 w0,w1,t0,t1,r,r2,val,rab,ccc,gridp(3),xyza(3),xyzb(3),dr3
   real*8 step,px,py,pz,xinc,yinc,zinc,nx,ny,nz,vv(1),v,est,nfod
   real*8 f1,f2,dx1,dy1,dz1,dx2,dy2,dz2,r1,dum,r1xy,r2xy
   real*8 dxx1,dyy1,dzz1,dxx2,dyy2,dzz2,ar1,ar2
   real*8 pij,pcoef


   real(wp) :: cube_step, cube_pthr
   integer ifile


   n=mol%nat
   xyz=mol%xyz

!  grid spacing for cube file
   cube_step = 0.5_wp
!  density matrix neglect threshold
   cube_pthr = 0.01_wp
!  distance criteria
   dist_cut = 200.0_wp
   
   write(*,*)
   write(*,*)'cube file module (H.N.)'
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
   !xst=floor((abs(px)+abs(nx))/step)
   xst=119
   ! x step size
   xinc=abs(px-nx)/xst
   ! y step number
   !yst=floor((abs(py)+abs(ny))/step)
   yst=119
   ! y step size
   yinc=abs(py-ny)/yst
   ! z step number
   !zst=floor((abs(pz)+abs(nz))/step)
   zst=119
   ! z step size
   zinc=abs(pz-nz)/zst
   write(*,*)'Total # of points', (xst+1)*(yst+1)*(zst+1)
   
   ! Volume for each gridpoint
   dr3=xinc*yinc*zinc

   write(*,*)'writing ',trim(fname)
   open(file=fname, newunit=ifile)
   write(ifile,*)'Cubepro cube file'
   write(ifile,*)'By H.N.'
   write(ifile,101)n,nx,ny,nz
   write(ifile,101)xst+1,xinc,0.0,0.0
   write(ifile,101)yst+1,0.0,yinc,0.0
   write(ifile,101)zst+1,0.0,0.0,zinc
   do i=1,n
      write(ifile,102)mol%num(mol%id(i)),0.0,xyz(1,i),xyz(2,i),xyz(3,i)
   enddo

   allocate(cb(0:zst,0:yst,0:xst))

   cb  =0

   !     Dmat loop         -----------------------------------
   ! loop over atoms 
   do iat=1, n
   do jat=1, n
      ! index shell offset for each atom
      iso = basis%ish_at(iat)
      jso = basis%ish_at(jat)
      ! Atom number
      izp=mol%id(iat)
      jzp=mol%id(jat)
      ! Aufpunkt of the GTOs
      xyza(1:3)=xyz(1:3,iat)
      xyzb(1:3)=xyz(1:3,jat)
      ! Calculte the diastance between Aufpunkt
      rab=(xyza(1)-xyzb(1))**2  &
      &      +(xyza(2)-xyzb(2))**2  &
      &      +(xyza(3)-xyzb(3))**2
      ! loop over shells at each atom
      if(rab.lt.dist_cut)then
         do ish=1, basis%nsh_at(iat)
         do jsh=1, basis%nsh_at(jat)
            ii = basis%iao_sh(iso + ish)
            jj = basis%iao_sh(jso + jsh)
            icgto = basis%cgto(ish,izp)
            jcgto = basis%cgto(jsh,jzp)
            ! loop over shellblock NAOs
            do iao = 1, basis%nao_sh(iso + ish)
            do jao = 1, basis%nao_sh(jso + jsh)
            ! Case: calculate spin polarized spin density
            if (sdens .and. size(wfn%density, DIM = 3).eq. 2) then 
               pij = wfn%density(ii+iao,jj+jao,1)-wfn%density(ii+iao,jj+jao,2)
            ! Case: calculate spin polarized density
            else if (.not.sdens .and. size(wfn%density, DIM = 3).eq. 2) then 
               pij = wfn%density(ii+iao,jj+jao,1)+wfn%density(ii+iao,jj+jao,2)
            ! Case: calculate density
            else
               pij = wfn%density(ii+iao,jj+jao,1)
            endif
               if (abs(pij).gt.thr)then
                  ! loop over primitives
                  do npri = 1 , icgto%nprim 
                  do nprj = 1 , jcgto%nprim
                     pcoef = pij* icgto%coeff(npri) * jcgto%coeff(nprj)
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
                              ar1=icgto%alpha(npri)*r1
                              if(ar1.lt.intcut2)then  ! exp(-16) < 1.d-7 i.e. zero
                                 call primvalf(r1,dx1,dy1,dz1,dxx1,dyy1,dzz1,ar1, &
                                    &             icgto%ang,iao,nexp,array,f1)
                                 r2=r2xy+dzz2
                                 ar2=jcgto%alpha(nprj)*r2
                                 if(ar2.lt.intcut2)then
                                    call primvalf(r2,dx2,dy2,dz2,dxx2,dyy2,dzz2,ar2, &
                                       &             jcgto%ang,jao,nexp,array,f2)
                                    cb(k,j,i)=cb(k,j,i)+pcoef*f1*f2
                                 endif
                              endif
                           enddo
                        enddo
                     enddo
                  enddo
                  enddo
               endif
            enddo
            enddo
         enddo
         enddo
      endif
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
   close(ifile)

   101   format(I5,3F16.6)
   102   format(I5,4F16.6)
 
 end subroutine cube

subroutine mocube(mol,wfn,fname,basis)
   implicit none
   type(structure_type), intent(in) :: mol
   type(basis_type), intent(in) :: basis
   type(wavefunction_type), intent(in) :: wfn
   character*(*), intent(in) :: fname

   integer :: n
   real*8 :: xyz(3,mol%nat)
   type(cgto_type) :: icgto, jcgto

   real*4 ,allocatable  ::cb   (:,:,:)
   real*8 ,allocatable  ::array(:)

   integer i,j,k,ii,jj,npri,nprj,nexp
   integer ish,jsh,izp,jzp,iao,jao,iso,jso
   integer iat,jat,xst,yst,zst,cou,iprimcount,jprimcount
   real*8 thr,thr3,intcut,intcut2,dist_cut
   real*8 w0,w1,t0,t1,r,r2,val,rab,ccc,gridp(3),xyza(3),xyzb(3),dr3
   real*8 step,px,py,pz,xinc,yinc,zinc,nx,ny,nz,vv(1),v,est,nfod
   real*8 f1,f2,dx1,dy1,dz1,dx2,dy2,dz2,r1,dum,r1xy,r2xy
   real*8 dxx1,dyy1,dzz1,dxx2,dyy2,dzz2,ar1,ar2
   real*8 cij,ccoef


   real(wp) :: cube_step, cube_pthr
   integer ifile


   n=mol%nat
   xyz=mol%xyz

!  grid spacing for cube file
   cube_step = 0.5_wp
!  density matrix neglect threshold
   cube_pthr = 0.01_wp
!  distance criteria
   dist_cut = 200.0_wp
   
   write(*,*)
   write(*,*)'cube file module (H.N.)'
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
   !xst=floor((abs(px)+abs(nx))/step)
   xst=119
   ! x step size
   xinc=abs(px-nx)/xst
   ! y step number
   !yst=floor((abs(py)+abs(ny))/step)
   yst=119
   ! y step size
   yinc=abs(py-ny)/yst
   ! z step number
   !zst=floor((abs(pz)+abs(nz))/step)
   zst=119
   ! z step size
   zinc=abs(pz-nz)/zst
   write(*,*)'Total # of points', (xst+1)*(yst+1)*(zst+1)
   
   ! Volume for each gridpoint
   dr3=xinc*yinc*zinc

   write(*,*)'writing ',trim(fname)
   open(file=fname, newunit=ifile)
   write(ifile,*)'Cubepro cube file'
   write(ifile,*)'By H.N.'
   write(ifile,101)n,nx,ny,nz
   write(ifile,101)xst+1,xinc,0.0,0.0
   write(ifile,101)yst+1,0.0,yinc,0.0
   write(ifile,101)zst+1,0.0,0.0,zinc
   do i=1,n
      write(ifile,102)mol%num(mol%id(i)),0.0,xyz(1,i),xyz(2,i),xyz(3,i)
   enddo

   allocate(cb(0:zst,0:yst,0:xst))

   cb  =0

   !     Dmat loop         -----------------------------------
   ! loop over atoms 
   do iat=1, n
      ! index shell offset for each atom
      iso = basis%ish_at(iat)
      ! Atom number
      izp=mol%id(iat)
      ! Aufpunkt of the GTOs
      xyza(1:3)=xyz(1:3,iat)
      ! Calculte the diastance between Aufpunkt
      ! loop over shells at each atom
         do ish=1, basis%nsh_at(iat)
            ii = basis%iao_sh(iso + ish)
            icgto = basis%cgto(ish,izp)
            ! loop over shellblock NAOs
            do iao = 1, basis%nao_sh(iso + ish)
            ! Case: calculate spin polarized spin density
            !if (size(wfn%coeff, DIM = 3).eq. 2) then 
            cij = wfn%coeff(ii+iao,wfn%homo(1),1)
            ! Case: calculate spin polarized density
            !else if (size(wfn%density, DIM = 3).eq. 2) then 
               !pij = wfn%density(ii+iao,jj+jao,1)+wfn%density(ii+iao,jj+jao,2)
            ! Case: calculate density
            !else
               !pij = wfn%density(ii+iao,jj+jao,1)
            !endif
                  ! loop over primitives
                  do npri = 1 , icgto%nprim 
                  ccoef = cij* icgto%coeff(npri)
                     ! grid loops
                     gridp(1)=nx
                     do i=0,xst
                        gridp(1)=nx+(xinc*i)
                        dx1=xyza(1)-gridp(1)
                        dxx1=dx1*dx1
                        gridp(2)=ny
                        do j=0,yst
                           gridp(2)=ny+(yinc*j)
                           dy1=xyza(2)-gridp(2)
                           dyy1=dy1*dy1
                           gridp(3)=nz
                           r1xy=dxx1+dyy1
                           do k=0,zst
                              gridp(3)=nz+(zinc*k)
                              dz1=xyza(3)-gridp(3)
                              dzz1=dz1*dz1
                              r1=r1xy+dzz1
                              ar1=icgto%alpha(npri)*r1
                              if(ar1.lt.intcut2)then  ! exp(-16) < 1.d-7 i.e. zero
                                 call primvalf(r1,dx1,dy1,dz1,dxx1,dyy1,dzz1,ar1, &
                                    &             icgto%ang,iao,nexp,array,f1)
                                    cb(k,j,i)=cb(k,j,i)+ccoef*f1
                              endif
                           enddo
                        enddo
                     enddo
                  enddo
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
   close(ifile)

   101   format(I5,3F16.6)
   102   format(I5,4F16.6)
 
 end subroutine mocube




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
subroutine primvalf(r,dx,dy,dz,dx2,dy2,dz2,alpr2,ang,iao,nexp,a,f)
   implicit none
   integer ang,iao,nexp
   real*8 a(0:nexp)
   real*8 r
   real*8 f,dx2,dy2,dz2
   real*8 dx,dy,dz,alpr2
   real(wp), parameter :: s3 = sqrt(3.0_wp)
   real(wp), parameter :: s3_6 = s3/6.0_wp

   select case (ang)
   ! s-function
   case(0)
     f=fastexp(nexp,a,alpr2)
   ! p-functions
   case(1)
     select case (iao)
     case(1)
       f=fastexp(nexp,a,alpr2)*dy
     case(2)
       f=fastexp(nexp,a,alpr2)*dz
     case(3)
       f=fastexp(nexp,a,alpr2)*dx
     case default
       error stop "Angular momentum not supported"
     end select
   ! d-functions
   case(2)
     select case (iao)
     case(1)
       f=fastexp(nexp,a,alpr2)*dx*dy
     case(2)
       f=fastexp(nexp,a,alpr2)*dy*dz
     case(3)
       f=fastexp(nexp,a,alpr2)*s3_6*(2.0_wp*dz2-dx2-dy2)
     case(4)
       f=fastexp(nexp,a,alpr2)*dx*dz
     case(5)
       f=fastexp(nexp,a,alpr2)*0.5_wp*(dx2-dy2)
     case default
       error stop "Angular momentum not supported"
     end select
   case default
       error stop "Angular momentum not supported"
   end select

end subroutine primvalf

real*8 function fastexp(nexp,array,arg)
   integer nexp
   real*8 arg
   real*8 array(0:nexp)

   fastexp=array(int(100.0d0*arg))

end function fastexp

subroutine get_spin_constants(wll, mol, bas)
   real(wp), allocatable, intent(out) :: wll(:, :, :)
   type(structure_type), intent(in) :: mol
   type(basis_type), intent(in) :: bas

   integer :: izp, ish, jsh, il, jl

   allocate(wll(bas%nsh, bas%nsh, mol%nid), source=0.0_wp)

   do izp = 1, mol%nid
      do ish = 1, bas%nsh_id(izp)
         il = bas%cgto(ish, izp)%ang
         do jsh = 1, bas%nsh_id(izp)
            jl = bas%cgto(jsh, izp)%ang
            wll(jsh, ish, izp) = get_spin_constant(jl, il, mol%num(izp))
         end do
      end do
   end do
end subroutine get_spin_constants



end module cubepro
