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
module charmm_module
  !*******************************************************************************
  !*******************************************************************************

  implicit none

  private
  public :: charmm

  type :: charmm
     integer                  :: natp
     real(8), allocatable     :: prms_vdw(:,:)
     real(8), allocatable     :: prms_itors(:,:,:,:,:)
     real(8)                  :: prms_tors(3)
     real(8)                  :: prms_angles(2)
     real(8)                  :: prms_bonds(2)
     character(6),allocatable :: atp(:)
   contains
     procedure :: set_natp
     procedure :: get_natp
!     procedure :: set_charmmtypes
     procedure :: set_charmmbonds
     procedure :: set_charmmangles
     procedure :: set_charmmdihedrals
     procedure :: set_charmmidihedrals
     procedure :: set_charmmvdw
  end type charmm

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

!  subroutine set_charmmtypes(this)
!    class(charmm), intent(inout) :: this
!    allocate(this%atp(this%get_natp()))
!  end subroutine set_charmmtypes

  subroutine set_charmmbonds(this,p1,p2)
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
  end subroutine set_charmmbonds

  subroutine set_charmmangles(this,p1,p2,p3)
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
  end subroutine set_charmmangles

  subroutine set_charmmdihedrals(this,p1,p2,p3,p4)
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
  end subroutine set_charmmdihedrals

  subroutine set_charmmidihedrals(this)
    class(charmm), intent(inout) :: this
    integer                     :: i,j,k,l,m
    allocate(&
         this%prms_itors(this%get_natp(),this%get_natp(),this%get_natp(),this%get_natp(),4))
    do i=1,this%get_natp()
       do j=1,this%get_natp()
          do k=1,this%get_natp()
             do l=1,this%get_natp()
                do m=1,4
                   this%prms_itors(i,j,k,l,m)=0.d0
                end do
             end do
          end do
       end do
    end do
    this%prms_itors(1,28,40,41,1)=120.000d0
    this%prms_itors(1,28,40,41,2)=0.0d0
    this%prms_itors(2,5,42,42,1)=96.000d0
    this%prms_itors(2,5,42,42,2)=0.0d0
    this%prms_itors(2,8,42,42,1)=96.000d0
    this%prms_itors(2,8,42,42,2)=0.0d0
    this%prms_itors(2,19,42,42,1)=96.000d0
    this%prms_itors(2,19,42,42,2)=0.0d0
    this%prms_itors(2,20,42,42,1)=96.000d0
    this%prms_itors(2,20,42,42,2)=0.0d0
    this%prms_itors(4,11,44,38,1)=50.000d0
    this%prms_itors(4,11,44,38,2)=0.0d0
    this%prms_itors(4,28,44,38,1)=50.000d0
    this%prms_itors(4,28,44,38,2)=0.0d0
    this%prms_itors(40,1,7,36,1)=20.000d0
    this%prms_itors(40,1,7,36,2)=0.0d0
    this%prms_itors(40,1,8,36,1)=20.000d0
    this%prms_itors(40,1,8,36,2)=0.0d0
    this%prms_itors(40,1,17,36,1)=20.000d0
    this%prms_itors(40,1,17,36,2)=0.0d0
    this%prms_itors(40,1,28,36,1)=20.000d0
    this%prms_itors(40,1,28,36,2)=0.0d0
    this%prms_itors(43,11,23,3,1)=70.000d0
    this%prms_itors(43,11,23,3,2)=0.0d0
    this%prms_itors(43,28,28,3,1)=70.000d0
    this%prms_itors(43,28,28,3,2)=0.0d0
 end subroutine set_charmmidihedrals

 subroutine set_charmmvdw(this)
   class(charmm), intent(inout) :: this
   integer                     :: i,j
   allocate(this%prms_vdw(this%get_natp(),2))
   do i=1,this%get_natp()
      do j=1,2
         this%prms_vdw(i,j)=0.d0
      end do
   end do
   this%prms_vdw(1,1)=0.1100d0
   this%prms_vdw(1,2)=2.0000d0
   this%prms_vdw(2,1)=0.0700d0
   this%prms_vdw(2,2)=2.0000d0
   this%prms_vdw(3,1)=0.0900d0
   this%prms_vdw(3,2)=2.0000d0
   this%prms_vdw(4,1)=0.0600d0
   this%prms_vdw(4,2)=1.8000d0
   this%prms_vdw(5,1)=0.0320d0
   this%prms_vdw(5,2)=2.0000d0
   this%prms_vdw(6,1)=0.0320d0
   this%prms_vdw(6,2)=2.0000d0
   this%prms_vdw(7,1)=0.0320d0
   this%prms_vdw(7,2)=2.0000d0
   this%prms_vdw(8,1)=0.0320d0
   this%prms_vdw(8,2)=2.0000d0
   this%prms_vdw(9,1)=0.0320d0
   this%prms_vdw(9,2)=2.0000d0
   this%prms_vdw(10,1)=0.0320d0
   this%prms_vdw(10,2)=2.0000d0
   this%prms_vdw(11,1)=0.0320d0
   this%prms_vdw(11,2)=2.0000d0
   this%prms_vdw(12,1)=0.0320d0
   this%prms_vdw(12,2)=2.0000d0
   this%prms_vdw(13,1)=0.0320d0
   this%prms_vdw(13,2)=2.0000d0
   this%prms_vdw(14,1)=0.0320d0
   this%prms_vdw(14,2)=2.0000d0
   this%prms_vdw(15,1)=0.0320d0
   this%prms_vdw(15,2)=2.0000d0
   this%prms_vdw(16,1)=0.0320d0
   this%prms_vdw(16,2)=2.0000d0
   this%prms_vdw(17,1)=0.0320d0
   this%prms_vdw(17,2)=2.0000d0
   this%prms_vdw(18,1)=0.0320d0
   this%prms_vdw(18,2)=2.0000d0
   this%prms_vdw(19,1)=0.0320d0
   this%prms_vdw(19,2)=2.0000d0
   this%prms_vdw(20,1)=0.0560d0
   this%prms_vdw(20,2)=2.0100d0
   this%prms_vdw(21,1)=0.0560d0
   this%prms_vdw(21,2)=2.0100d0
   this%prms_vdw(22,1)=0.0560d0
   this%prms_vdw(22,2)=2.0100d0
   this%prms_vdw(23,1)=0.0560d0
   this%prms_vdw(23,2)=2.0100d0
   this%prms_vdw(24,1)=0.0600d0
   this%prms_vdw(24,2)=2.0200d0
   this%prms_vdw(25,1)=0.0600d0
   this%prms_vdw(25,2)=2.0200d0
   this%prms_vdw(26,1)=0.0560d0
   this%prms_vdw(26,2)=2.0100d0
   this%prms_vdw(27,1)=0.0560d0
   this%prms_vdw(27,2)=2.0100d0
   this%prms_vdw(28,1)=0.0780d0
   this%prms_vdw(28,2)=2.0400d0
   this%prms_vdw(30,1)=0.0450d0
   this%prms_vdw(30,2)=1.3400d0
   this%prms_vdw(31,1)=0.0450d0
   this%prms_vdw(31,2)=1.3400d0
   this%prms_vdw(32,1)=0.0350d0
   this%prms_vdw(32,2)=1.3400d0
   this%prms_vdw(33,1)=0.0350d0
   this%prms_vdw(33,2)=1.3000d0
   this%prms_vdw(34,1)=0.0240d0
   this%prms_vdw(34,2)=1.3400d0
   this%prms_vdw(35,1)=0.0240d0
   this%prms_vdw(35,2)=1.3400d0
   this%prms_vdw(36,1)=0.0460d0
   this%prms_vdw(36,2)=0.2245d0
   this%prms_vdw(37,1)=0.0460d0
   this%prms_vdw(37,2)=0.2245d0
   this%prms_vdw(38,1)=0.0460d0
   this%prms_vdw(38,2)=0.8000d0
   this%prms_vdw(39,1)=0.0460d0
   this%prms_vdw(39,2)=0.2245d0
   this%prms_vdw(40,1)=0.2000d0
   this%prms_vdw(40,2)=1.8500d0
   this%prms_vdw(41,1)=0.1200d0
   this%prms_vdw(41,2)=1.7000d0
   this%prms_vdw(42,1)=0.1200d0
   this%prms_vdw(42,2)=1.7000d0
   this%prms_vdw(43,1)=0.0500d0
   this%prms_vdw(43,2)=1.7000d0
   this%prms_vdw(44,1)=0.1200d0
   this%prms_vdw(44,2)=1.7000d0
   this%prms_vdw(45,1)=0.1200d0
   this%prms_vdw(45,2)=1.7000d0
   this%prms_vdw(46,1)=0.1000d0
   this%prms_vdw(46,2)=1.6500d0
   this%prms_vdw(47,1)=0.1000d0
   this%prms_vdw(47,2)=1.6500d0
   this%prms_vdw(48,1)=0.1000d0
   this%prms_vdw(48,2)=1.6500d0
   this%prms_vdw(49,1)=0.1521d0
   this%prms_vdw(49,2)=1.7700d0
   this%prms_vdw(50,1)=0.1921d0
   this%prms_vdw(50,2)=1.7650d0
   this%prms_vdw(51,1)=0.1921d0
   this%prms_vdw(51,2)=1.7650d0
   this%prms_vdw(52,1)=0.1521d0
   this%prms_vdw(52,2)=1.7700d0
   this%prms_vdw(53,1)=0.1000d0
   this%prms_vdw(53,2)=1.6500d0
   this%prms_vdw(54,1)=0.1000d0
   this%prms_vdw(54,2)=1.6500d0
   this%prms_vdw(55,1)=0.1000d0
   this%prms_vdw(55,2)=1.6500d0
   this%prms_vdw(56,1)=0.1521d0
   this%prms_vdw(56,2)=1.7682d0
   this%prms_vdw(57,1)=0.5850d0
   this%prms_vdw(57,2)=2.1500d0
   this%prms_vdw(58,1)=0.4700d0
   this%prms_vdw(58,2)=2.1000d0
 end subroutine set_charmmvdw

end module charmm_module
