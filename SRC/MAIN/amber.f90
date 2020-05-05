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
module amber

contains

  subroutine amber_vdw(p1,p2,prms)

    implicit none

    real(8) epsi,epsj,ri,rj,prms(2)
    character(2) p1,p2

    select case(p1)
    case('OW')
       epsi=0.006591346d0
       !       ri=1.7683d0
       ri=1.582746d0
    case('HW')
       epsi=0.d0
       ri=0.d0
    case('HO')
       epsi=0.d0
       ri=0.d0
    end select

    select case(p2)
    case('OW')
       epsj=0.006591346d0
       !       rj=1.7683d0
       rj=1.582746d0
    case('HW')
       epsj=0.d0
       rj=0.d0
    case('HO')
       epsj=0.d0
       rj=0.d0
    end select

    prms(1)=sqrt(epsi*epsj)
    prms(2)=ri+rj

    return

  end subroutine amber_vdw

  subroutine amber_bonds(p1,p2,prms)

    implicit none

    integer i
    real(8) prms(2)
    character(2) p1,p2,pa,pb

    prms(1)=0.d0
    prms(2)=0.d0

    pa=p1
    pb=p2

    do i=1,2
       select case(pa)
       case('OW')
          select case(pb)
          case('HW')
             prms(1)=45.9296d0
             prms(2)=1.0120d0
          end select
       end select
       pa=p2
       pb=p1
    end do

    return

  end subroutine amber_bonds

  subroutine amber_bends(p1,p2,p3,prms)

    implicit none

    integer i
    real(8) prms(2)
    character(2) p1,p2,p3,pa,pb,pc

    prms(1)=0.d0
    prms(2)=0.d0

    pa=p1
    pb=p2
    pc=p3

    do i=1,2
       select case(pa)
       case('HW')
          select case(pb)
          case('OW')
             select case(pc)
             case('HW')
                prms(1)=3.2913d0
                prms(2)=104.24d0
             end select
          end select
       end select
       pa=p3
       pc=p1
    end do

    return

  end subroutine amber_bends

end module amber
