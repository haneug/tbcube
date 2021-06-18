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

!> Version information for this project
module cubepro_version
   implicit none
   private

   public :: cubepro_version_string, cubepro_version_compact
   public :: get_cubepro_version


   !> String representation of the cubepro version
   character(len=*), parameter :: cubepro_version_string = "0.1.0"

   !> Numeric representation of the cubepro version
   integer, parameter :: cubepro_version_compact(3) = [0, 1, 0]

contains

!> Getter function to retrieve cubepro version
subroutine get_cubepro_version(major, minor, patch, string)
   !> Major version number of the cubepro version
   integer, intent(out), optional :: major
   !> Minor version number of the cubepro version
   integer, intent(out), optional :: minor
   !> Patch version number of the cubepro version
   integer, intent(out), optional :: patch
   !> String representation of the cubepro version
   character(len=:), allocatable, intent(out), optional :: string

   if (present(major)) then
      major = cubepro_version_compact(1)
   end if
   if (present(minor)) then
      minor = cubepro_version_compact(2)
   end if
   if (present(patch)) then
      patch = cubepro_version_compact(3)
   end if
   if (present(string)) then
      string = cubepro_version_string
   end if

end subroutine get_cubepro_version

end module cubepro_version
