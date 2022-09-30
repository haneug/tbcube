! This file is part of tbcube.
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

!> Command-line driver for cube file printout
program main
   use, intrinsic :: iso_fortran_env, only : output_unit, error_unit, input_unit
   use mctc_env, only : error_type, fatal_error, get_argument, wp
   use mctc_io, only : structure_type, read_structure, filetype, get_filetype
   use tbcube, only : get_cube, cube_input
   use tbcube_version, only : get_tbcube_version
   implicit none
   character(len=*), parameter :: prog_name = "tbcube"

   type, extends(cube_input) :: driver_config
      character(len=:), allocatable :: input
      integer, allocatable :: input_format
      integer, allocatable :: charge
      integer, allocatable :: spin
   end type driver_config

   type(driver_config) :: config
   type(structure_type) :: mol
   type(error_type), allocatable :: error
   integer :: spin, charge, stat, unit
   logical :: exist

   call get_arguments(config, error)
   if (allocated(error)) then
      write(error_unit, '(a)') error%message
      error stop
   end if

   if (config%input == "-") then
      if (.not.allocated(config%input_format)) config%input_format = filetype%xyz
      call read_structure(mol, input_unit, config%input_format, error)
   else
      call read_structure(mol, config%input, error, config%input_format)
   end if
   if (allocated(error)) then
      write(error_unit, '(a)') error%message
      error stop
   end if

   if (allocated(config%charge)) then
      mol%charge = config%charge
   else
      inquire(file='.CHRG', exist=exist)
      if (exist) then
         open(file='.CHRG', newunit=unit)
         read(unit, *, iostat=stat) charge
         if (stat == 0) then
            mol%charge = charge
            if (config%verbosity > 0) write(output_unit, '(a)') &
               "[Info] Molecular charge read from .CHRG"
         else
            if (config%verbosity > 0) write(output_unit, '(a)') &
               "[Warn] Could not read molecular charge read from .CHRG"
         end if
         close(unit)
      end if
   end if

   if (allocated(config%spin)) then
      mol%uhf = config%spin
   else
      inquire(file='.UHF', exist=exist)
      if (exist) then
         open(file='.UHF', newunit=unit)
         read(unit, *, iostat=stat) spin
         if (stat == 0) then
            mol%uhf = spin
            if (config%verbosity > 0) write(output_unit, '(a)') &
               "[Info] Molecular spin read from .UHF"
         else
            if (config%verbosity > 0) write(output_unit, '(a)') &
               "[Warn] Could not read molecular spin read from .UHF"
         end if
         close(unit)
      end if
   end if

   if (.not.allocated(config%method)) config%method = "gfn2"


   ! This needs to be modified
   call get_cube(config%cube_input, mol, error)
   if (allocated(error)) then
      write(error_unit, '(a)') error%message
      error stop
   end if

contains


subroutine help(unit)
   integer, intent(in) :: unit

   write(unit, '(a, *(1x, a))') &
      "Usage: "//prog_name//" [options] <input>"

   write(unit, '(a)') &
      "", &
      ""

   write(unit, '(2x, a, t25, a)') &
      "-c, --charge <real>", "Set charge to molecule, overwrites .CHRG file", &
      "    --spin <int>", "Set number of unpaired electrons, overwrites .UHF file", &
      "    --method <name>", "Parametrization of the xTB Hamiltonian to use", &
      "", "Available methods: gfn1, ipea1, gfn2 (Default: gfn2)", &
      "-spGFN  --spin-polarized,    Use spin-polarized xTB Hamiltonian", &
      "    --sdens, Request Spin-Density, enables spin-polarization", &
      "    --dens, Request Density", &
      "    --homo, Request HOMO", &
      "    --lumo, Request LUMO", &
      "    --acc <real>", "Accuracy of the tight-binding calculation", &
      "    --thr <real>", "Threshold for the determination of the fragments", &
      "    --etemp <real>", "Electronic temperature in Kelvin", &
      "-i, --input <format>", "Hint for the format of the input file", &
      "    --version", "Print program version and exit", &
      "    --help", "Show this help message"

   write(unit, '(a)')

end subroutine help


subroutine version(unit)
   integer, intent(in) :: unit
   character(len=:), allocatable :: version_string

   call get_tbcube_version(string=version_string)
   write(unit, '(a, *(1x, a))') &
      & prog_name, "version", version_string

end subroutine version


subroutine get_argument_as_real(iarg, val, error)
   !> Index of command line argument, range [0:command_argument_count()]
   integer, intent(in) :: iarg
   !> Real value
   real(wp), intent(out) :: val
   !> Error handling
   type(error_type), allocatable :: error

   integer :: stat
   character(len=:), allocatable :: arg

   call get_argument(iarg, arg)
   if (.not.allocated(arg)) then
      call fatal_error(error, "Cannot read real value, argument missing")
      return
   end if
   read(arg, *, iostat=stat) val
   if (stat /= 0) then
      call fatal_error(error, "Cannot read real value from '"//arg//"'")
      return
   end if

end subroutine get_argument_as_real


subroutine get_argument_as_int(iarg, val, error)
   !> Index of command line argument, range [0:command_argument_count()]
   integer, intent(in) :: iarg
   !> Real value
   integer, intent(out) :: val
   !> Error handling
   type(error_type), allocatable :: error

   integer :: stat
   character(len=:), allocatable :: arg

   call get_argument(iarg, arg)
   if (.not.allocated(arg)) then
      call fatal_error(error, "Cannot read integer value, argument missing")
      return
   end if
   read(arg, *, iostat=stat) val
   if (stat /= 0) then
      call fatal_error(error, "Cannot read integer value from '"//arg//"'")
      return
   end if

end subroutine get_argument_as_int

subroutine get_arguments(config, error)

   !> Settings for this run
   type(driver_config), intent(out) :: config

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iarg, narg
   logical :: getopts
   character(len=:), allocatable :: arg

   iarg = 0
   getopts = .true.
   narg = command_argument_count()
   do while(iarg < narg)
      iarg = iarg + 1
      call get_argument(iarg, arg)
      if (.not.getopts) then
         if (.not.allocated(config%input)) then
            call move_alloc(arg, config%input)
            cycle
         end if
         call fatal_error(error, "Too many positional arguments present")
         exit
      end if
      select case(arg)
      case("--")
         getopts = .false.
      case("--help")
         call help(output_unit)
         stop
      case("--version")
         call version(output_unit)
         stop
      case default
         if (.not.allocated(config%input)) then
            call move_alloc(arg, config%input)
            cycle
         end if
         call fatal_error(error, "Too many positional arguments present")
         exit
      case("-i", "--input")
         iarg = iarg + 1
         call get_argument(iarg, arg)
         if (.not.allocated(arg)) then
            call fatal_error(error, "Missing argument for input format")
            exit
         end if
         config%input_format = get_filetype("."//arg)
      case("-c", "--charge")
         iarg = iarg + 1
         allocate(config%charge)
         call get_argument_as_int(iarg, config%charge, error)
         if (allocated(error)) exit
      case("--spin")
         iarg = iarg + 1
         allocate(config%spin)
         call get_argument_as_int(iarg, config%spin, error)
         if (allocated(error)) exit
      case("--spin-polarized")
         config%nspin = 2
      case("-spGFN")
         config%nspin = 2
      case("--sdens")
         config%nspin = 2
         config%sdens = .true.
      case("--dens")
         config%dens = .true.
      case("--homo")
         config%homo = .true.
      case("--occ")
         config%occ = .true.
      case("--lumo")
         config%lumo = .true.
      case("--method")
         iarg = iarg + 1
         call get_argument(iarg, config%method)
         if (.not.allocated(config%method)) then
            call fatal_error(error, "Missing argument for method")
            exit
         end if
      case("--acc")
         iarg = iarg + 1
         call get_argument_as_real(iarg, config%accuracy, error)
         if (allocated(error)) exit
      case("--etemp")
         iarg = iarg + 1
         call get_argument_as_real(iarg, config%etemp, error)
         if (allocated(error)) exit
      end select
   end do

   if (.not.(allocated(config%input))) then
      if (.not.allocated(error)) then
         call help(output_unit)
         error stop
      end if
   end if

end subroutine get_arguments


end program main
