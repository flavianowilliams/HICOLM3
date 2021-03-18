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

    integer i,j
    real(4) x1,x2
    real(8) prms(2),epsi(2),ri(2)
    character(2) p1,p2,pa,px

    open(4,file='/tmp/amber/amber_vdw.prm',status='old')

    !-atribuindo valores iniciais

    do i=1,2
       epsi(i)=0.d0
       ri(i)=0.d0
    end do

    px=p1
    do i=1,2
       do j=1,1000
          read(4,*,end=10)pa,x1,x2
          if(pa.eq.px)then
             epsi(i)=dble(x2)
             ri(i)=dble(x1)
          end if
       end do
10     rewind(4)
       px=p2
    end do

    prms(1)=sqrt(epsi(1)*epsi(2))
    prms(2)=ri(1)+ri(2)

    close(4)

    return

  end subroutine amber_vdw

  subroutine amber_bonds(p1,p2,prms)

    implicit none

    integer i
    real(8) prms(2)
    real(4) x1,x2
    character(2) p1,p2,pa,pb

    open(4,file='/tmp/amber/amber_bonds.prm',status='old')

    !-atribuindo valores iniciais

    prms(1)=0.d0
    prms(2)=0.d0

    do i=1,1000
       read(4,*,end=10)pa,pb,x1,x2
       if(pa.eq.p1.and.pb.eq.p2)then
          prms(1)=dble(x1)
          prms(2)=dble(x2)
       elseif(pb.eq.p1.and.pa.eq.p2)then
          prms(1)=dble(x1)
          prms(2)=dble(x2)
       end if
    end do

10  close(4)

    return

  end subroutine amber_bonds

  subroutine amber_bends(p1,p2,p3,prms)

    implicit none

    integer i
    real(8) prms(2)
    real(4) x1,x2
    character(2) p1,p2,p3,pa,pb,pc

    open(4,file='/tmp/amber/amber_angles.prm',status='old')

    !-atribuindo valores iniciais

    prms(1)=0.d0
    prms(2)=0.d0

    do i=1,1000
       read(4,*,end=10)pa,pb,pc,x1,x2
       if(pa.eq.p1.and.pb.eq.p2.and.pc.eq.p3)then
          prms(1)=dble(x1)
          prms(2)=dble(x2)
       elseif(pa.eq.p3.and.pb.eq.p2.and.pc.eq.p1)then
          prms(1)=dble(x1)
          prms(2)=dble(x2)
       end if
    end do

10  close(4)

    return

  end subroutine amber_bends

  subroutine amber_dihedrals(p1,p2,p3,p4,prms)

    implicit none

    real(8) prms(4)
    character(2) p1,p2,p3,p4

    !-atribuindo valores iniciais

    prms(1)=1
    prms(2)=0.d0
    prms(3)=0.d0
    prms(4)=0.d0

    call amber_dihedrals_general(p2,p3,prms)
    call amber_dihedrals_proper(p1,p2,p3,p4,prms)

    return

  end subroutine amber_dihedrals

  subroutine amber_dihedrals_general(p2,p3,prms)

    implicit none

    integer i,x1,x4
    real(8) prms(4)
    real(4) x2,x3
    character(2) p2,p3,pa,pb,pc,pd

    open(4,file='/tmp/amber/amber_dihedrals_general.prm',status='old')

    do i=1,1000
       read(4,*,end=10)pa,pb,pc,pd,x1,x2,x3,x4
       if(pb.eq.p2.and.pc.eq.p3)then
          prms(1)=dfloat(x1)
          prms(2)=dble(x2)
          prms(3)=dble(x3)
          prms(4)=dfloat(x4)
       elseif(pb.eq.p3.and.pc.eq.p2)then
          prms(1)=dfloat(x1)
          prms(2)=dble(x2)
          prms(3)=dble(x3)
          prms(4)=dfloat(x4)
       end if
    end do

10  close(4)

    return

  end subroutine amber_dihedrals_general

  subroutine amber_dihedrals_proper(p1,p2,p3,p4,prms)

    implicit none

    integer i,x1,x4
    real(8) prms(4)
    real(4) x2,x3
    character(2) p1,p2,p3,p4,pa,pb,pc,pd

    open(4,file='/tmp/amber/amber_dihedrals_proper.prm',status='old')

    do i=1,1000
       read(4,*,end=10)pa,pb,pc,pd,x1,x2,x3,x4
       if(pa.eq.p1.and.pb.eq.p2.and.pc.eq.p3.and.pd.eq.p4)then
          prms(1)=dfloat(x1)
          prms(2)=dble(x2)
          prms(3)=dble(x3)
          prms(4)=dfloat(x4)
       elseif(pa.eq.p4.and.pb.eq.p3.and.pc.eq.p2.and.pd.eq.p1)then
          prms(1)=dfloat(x1)
          prms(2)=dble(x2)
          prms(3)=dble(x3)
          prms(4)=dfloat(x4)
       end if
    end do

10  close(4)

    return

  end subroutine amber_dihedrals_proper

end module amber
