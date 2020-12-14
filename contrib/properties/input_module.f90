!
!MIT License
!
!Copyright (c) 2020 flavianowilliams
!
!Permission is hereby granted, free of charge, to any person obtaining a copy
!of this software and associated documentation files (the "Software"), to deal
!in the Software without restriction, including without limitation the rights
!to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
!copies of the Software, and to permit persons to whom the Software is
!furnished to do so, subject to the following conditions:
!
!THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
!AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
!OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
!SOFTWARE.
!
module input_module

    implicit none

  integer nstp,nvar,nat,nopt,spctot
  real(8) dtime

contains

  subroutine entrada

    implicit none

    real(8) lxf
    character lxc

    !  read(1,*)nstp,dtime,nvar,nat,spctot
    read(2,*)
    read(2,*)
    read(2,*)
    read(2,*)lxc,lxc,lxc,nat,lxf,lxf,dtime,lxf,nstp
    read(2,*)
    read(2,*)

!  write(*,*)'RDF -> 1'
!  write(*,*)'Medias -> 2'

  nopt=1
!  read(*,*)nopt

 end subroutine entrada

end module input_module
