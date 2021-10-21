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
module opls_module
  !*******************************************************************************
  !*******************************************************************************

  implicit none

  private
  public :: opls

  type :: opls
     integer                  :: natp
     real(8)                  :: prms_vdw(2)
     real(8)                  :: prms_itors(3)
     real(8)                  :: prms_tors(3)
     real(8)                  :: prms_angles(2)
     real(8)                  :: prms_bonds(2)
     character(6),allocatable :: atp(:)
   contains
     procedure :: set_natp
     procedure :: get_natp
     procedure :: set_oplsbonds
     procedure :: set_oplsangles
     procedure :: set_oplsdihedrals
     procedure :: set_oplsidihedrals
     procedure :: set_oplsvdw
  end type opls

contains

  subroutine set_natp(this,natp)
    class(charmm), intent(inout) :: this
    integer, intent(in)         :: natp
    this%natp=natp
  end subroutine set_natp

  integer function get_natp(this)
    class(charmm), intent(in) :: this
    get_natp=this%natp
  end function get_natp

  subroutine set_oplsbonds(this,p1,p2)
    implicit none
    class(charmm), intent(inout) :: this
    character(6), intent(in)     :: p1,p2
    character(6)                 :: pa,pb
    real(8)                      :: x1,x2
    logical                      :: check
    this%prms_bonds(1)=0.d0
    this%prms_bonds(2)=0.d0
    check=.true.
    do while(check)
       read(12,*,end=1)pa,pb,x1,x2
       if(pa.eq.p1.and.pb.eq.p2)then
          this%prms_bonds(1)=dble(x1)
          this%prms_bonds(2)=dble(x2)
          check=.false.
       elseif(pb.eq.p1.and.pa.eq.p2)then
          this%prms_bonds(1)=dble(x1)
          this%prms_bonds(2)=dble(x2)
          check=.false.
       end if
    end do
1   return
  end subroutine set_oplsbonds

  subroutine set_oplsangles(this,p1,p2,p3)
    implicit none
    class(charmm), intent(inout) :: this
    character(6), intent(in)     :: p1,p2,p3
    character(6)                 :: pa,pb,pc
    real(8)                      :: x1,x2
    logical                      :: check
    this%prms_angles(1)=0.d0
    this%prms_angles(2)=0.d0
    check=.true.
    do while(check)
       read(12,*,end=1)pa,pb,pc,x1,x2
       if(pa.eq.p1.and.pb.eq.p2.and.pc.eq.p3)then
          this%prms_angles(1)=dble(x1)
          this%prms_angles(2)=dble(x2)
          check=.false.
       elseif(pc.eq.p1.and.pb.eq.p2.and.pa.eq.p3)then
          this%prms_angles(1)=dble(x1)
          this%prms_angles(2)=dble(x2)
          check=.false.
       end if
    end do
1   return
  end subroutine set_oplsangles

  subroutine set_oplsdihedrals(this,p1,p2,p3,p4)
    implicit none
    class(charmm), intent(inout) :: this
    character(6), intent(in)     :: p1,p2,p3,p4
    character(6)                 :: pa,pb,pc,pd
    integer                      :: n1
    real(8)                      :: x1,x2
    logical                      :: check
    this%prms_tors(1)=0.d0
    this%prms_tors(2)=0.d0
    this%prms_tors(3)=0.d0
    check=.true.
    do while(check)
       read(12,*,end=1)pa,pb,pc,pd,x1,n1,x2
       if(pa.eq.p1.and.pb.eq.p2.and.pc.eq.p3.and.pd.eq.p4)then
          this%prms_tors(1)=dble(x1)
          this%prms_tors(2)=dble(n1)
          this%prms_tors(3)=dble(x2)
          check=.false.
       elseif(pd.eq.p1.and.pc.eq.p2.and.pb.eq.p3.and.pa.eq.p4)then
          this%prms_tors(1)=dble(x1)
          this%prms_tors(2)=dble(n1)
          this%prms_tors(3)=dble(x2)
          check=.false.
       end if
    end do
1   return
  end subroutine set_oplsdihedrals

  subroutine set_oplsidihedrals(this,p1,p2,p3,p4)
    implicit none
    class(charmm), intent(inout) :: this
    character(6), intent(in)     :: p1,p2,p3,p4
    character(6)                 :: pa,pb,pc,pd
    integer                      :: n1
    real(8)                      :: x1,x2
    logical                      :: check
    this%prms_itors(1)=0.d0
    this%prms_itors(2)=0.d0
    this%prms_itors(3)=0.d0
    check=.true.
    do while(check)
       read(12,*,end=1)pa,pb,pc,pd,x1,n1,x2
       if(pa.eq.p1.and.pb.eq.p2.and.pc.eq.p3.and.pd.eq.p4)then
          this%prms_itors(1)=dble(x1)
          this%prms_itors(2)=dble(n1)
          this%prms_itors(3)=dble(x2)
          check=.false.
       elseif(pd.eq.p1.and.pc.eq.p2.and.pb.eq.p3.and.pa.eq.p4)then
          this%prms_itors(1)=dble(x1)
          this%prms_itors(2)=dble(n1)
          this%prms_itors(3)=dble(x2)
          check=.false.
       end if
    end do
1   return
  end subroutine set_oplsidihedrals

  subroutine set_oplsvdw(this,p1)
    class(charmm), intent(inout) :: this
    character(6), intent(in)     :: p1
    character(6)                 :: pa
    real(8)                      :: x1,x2
    logical                      :: check
    this%prms_vdw(1)=0.d0
    this%prms_vdw(2)=0.d0
    check=.true.
    do while(check)
       read(12,*,end=1)pa,x1,x2
       if(pa.eq.p1)then
          this%prms_vdw(1)=dble(x1)
          this%prms_vdw(2)=dble(x2)
          check=.false.
       end if
    end do
1   return
  end subroutine set_oplsvdw

end module opls_module
