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

  private
  public :: amber

  type :: amber
     real(8)                  :: prms_bonds(2)
     real(8)                  :: prms_bends(2)
     real(8)                  :: prms_tors(4)
     real(8)                  :: prms_vdw(2)
     character(2),allocatable :: atp(:)
   contains
     procedure :: set_amberbonds
     procedure :: set_amberbends
     procedure :: set_amberdihedrals
     procedure :: set_ambervdw
     procedure :: set_ambertypes
     generic   :: set_amber => set_amberbonds, set_amberbends, set_ambervdw, &
          set_amberdihedrals
  end type amber

contains

  subroutine set_ambertypes(this)
    class(amber), intent(inout) :: this
    allocate(this%atp(62))
    this%atp(1)='Br'
    this%atp(2)='C'
    this%atp(3)='C*'
    this%atp(4)='C0'
    this%atp(5)='CA'
    this%atp(6)='CB'
    this%atp(7)='CC'
    this%atp(8)='CD'
    this%atp(9)='CK'
    this%atp(10)='Cl'
    this%atp(11)='CM'
    this%atp(12)='CN'
    this%atp(13)='CQ'
    this%atp(14)='CR'
    this%atp(15)='Cs'
    this%atp(16)='CT'
    this%atp(17)='CU'
    this%atp(18)='CV'
    this%atp(19)='CW'
    this%atp(20)='CY'
    this%atp(21)='CZ'
    this%atp(22)='F'
    this%atp(23)='FE'
    this%atp(24)='H'
    this%atp(25)='H1'
    this%atp(26)='H2'
    this%atp(27)='H3'
    this%atp(28)='H4'
    this%atp(29)='H5'
    this%atp(30)='HA'
    this%atp(31)='HC'
    this%atp(32)='HO'
    this%atp(33)='HP'
    this%atp(34)='HS'
    this%atp(35)='HW'
    this%atp(36)='HZ'
    this%atp(37)='I'
    this%atp(38)='IB'
    this%atp(39)='IM'
    this%atp(40)='IP'
    this%atp(41)='K'
    this%atp(42)='Li'
    this%atp(43)='MG'
    this%atp(44)='N'
    this%atp(45)='N*'
    this%atp(46)='N2'
    this%atp(47)='N3'
    this%atp(48)='Na'
    this%atp(49)='NB'
    this%atp(50)='NC'
    this%atp(51)='NT'
    this%atp(52)='NY'
    this%atp(53)='O'
    this%atp(54)='O2'
    this%atp(55)='OH'
    this%atp(56)='OS'
    this%atp(57)='OW'
    this%atp(58)='P'
    this%atp(59)='Rb'
    this%atp(60)='S'
    this%atp(61)='SH'
    this%atp(62)='Zn'
  end subroutine set_ambertypes

  subroutine set_amberbonds(this,p1,p2)
    class(amber), intent(inout) :: this
    integer                     :: i
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
    integer                     :: i
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
    integer                     :: i
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
    integer                     :: i
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
