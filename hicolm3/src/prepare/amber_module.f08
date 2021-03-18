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
module amber_module
  !*******************************************************************************************
  !*******************************************************************************************

  implicit none

  integer i

  private
  public :: amber

  type :: amber
     real(8)      :: prms_bonds(2)
     real(8)      :: prms_bends(2)
     real(8)      :: prms_tors(4)
     real(8)      :: prms_vdw(2)
   contains
     procedure :: set_amberbonds
     procedure :: set_amberbends
     procedure :: set_amberdihedrals
     procedure :: set_ambervdw
     generic   :: set_amber => set_amberbonds, set_amberbends, set_ambervdw, &
          set_amberdihedrals
  end type amber

contains

  subroutine set_amberbonds(this,p1,p2)
    class(amber), intent(inout) :: this
    real(4) x1,x2
    character(2) pa,pb,p1,p2
    open(12,file='/tmp/amber/amber_bonds.prm',status='old')
    this%prms_bonds(1)=0.d0
    this%prms_bonds(2)=0.d0
    do i=1,115
       read(12,*,end=1)pa,pb,x1,x2
       if(pa.eq.p1.and.pb.eq.p2.or.pa.eq.p2.and.pb.eq.p1)then
          this%prms_bonds(1)=dble(x1)
          this%prms_bonds(2)=dble(x2)
       end if
    end do
1   close(12)
  end subroutine set_amberbonds

  subroutine set_amberbends(this,p1,p2,p3)
    class(amber), intent(inout) :: this
    real(4) x1,x2
    character(2) pa,pb,pc,p1,p2,p3
    open(12,file='/tmp/amber/amber_angles.prm',status='old')
    this%prms_bends(1)=0.d0
    this%prms_bends(2)=0.d0
    do i=1,281
       read(12,*,end=1)pa,pb,pc,x1,x2
       if(pa.eq.p1.and.pc.eq.p3.or.pa.eq.p3.and.pc.eq.p1)then
          if(pb.eq.p2)then
             this%prms_bends(1)=dble(x1)
             this%prms_bends(2)=dble(x2)
          end if
       end if
    end do
1   close(12)
  end subroutine set_amberbends

  subroutine set_amberdihedrals(this,p1,p2,p3,p4)
    class(amber), intent(inout) :: this
    real(4) x1,x2,x3,x4
    character(2) pa,pb,pc,pd,p1,p2,p3,p4
    open(12,file='/tmp/amber/amber_dihedrals_general.prm',status='old')
    open(13,file='/tmp/amber/amber_dihedrals_proper.prm',status='old')
    this%prms_tors(1)=0.d0
    this%prms_tors(2)=0.d0
    this%prms_tors(3)=0.d0
    this%prms_tors(4)=0.d0
    do i=1,63
       read(12,*,end=1)pa,pb,pc,pd,x1,x2,x3,x4
       if(pb.eq.p2.and.pc.eq.p3.or.pb.eq.p3.and.pc.eq.p2)then
          this%prms_tors(1)=dble(x1)
          this%prms_tors(2)=dble(x2)
          this%prms_tors(3)=dble(x3)
          this%prms_tors(4)=dble(x4)
       end if
    end do
1   do i=1,100
       read(13,*,end=2)pa,pb,pc,pd,x1,x2,x3,x4
       if(pa.eq.p1.and.pb.eq.p2.and.pc.eq.p3.and.pd.eq.p4)then
          this%prms_tors(1)=dble(x1)
          this%prms_tors(2)=dble(x2)
          this%prms_tors(3)=dble(x3)
          this%prms_tors(4)=dble(x4)
       end if
    end do
2   close(12)
    close(13)
  end subroutine set_amberdihedrals

  subroutine set_ambervdw(this,p1)
    class(amber), intent(inout) :: this
    real(4) x1,x2
    character(2) pa,p1
    open(12,file='/tmp/amber/amber_vdw.prm',status='old')
    this%prms_vdw(1)=0.d0
    this%prms_vdw(2)=0.d0
    do i=1,43
       read(12,*,end=1)pa,x1,x2
       if(pa.eq.p1)then
          this%prms_vdw(1)=dble(x1)
          this%prms_vdw(2)=dble(x2)
       end if
    end do
1   close(12)
  end subroutine set_ambervdw

end module amber_module
