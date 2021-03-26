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
     real(8)                  :: prms_bends(2)
     real(8)                  :: prms_tors(4)
     real(8)                  :: prms_vdw(2)
     real(8), allocatable     :: prms_bonds(:,:,:)
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

  subroutine set_amberbonds(this)
    class(amber), intent(inout) :: this
    allocate(this%prms_bonds(62,62,2))
    this%prms_bonds(1,5,1)=172.0d0
    this%prms_bonds(1,5,2)=1.8900d0
    this%prms_bonds(1,16,1)=159.0d0
    this%prms_bonds(1,16,2)=1.9440d0
    this%prms_bonds(2,2,1)=310.0d0
    this%prms_bonds(2,2,2)=1.5250d0
    this%prms_bonds(2,5,1)=469.0d0
    this%prms_bonds(2,5,2)=1.4090d0
    this%prms_bonds(2,6,1)=447.0d0
    this%prms_bonds(2,6,2)=1.4190d0
    this%prms_bonds(2,11,1)=410.0d0
    this%prms_bonds(2,11,2)=1.4440d0
    this%prms_bonds(2,16,1)=317.0d0
    this%prms_bonds(2,16,2)=1.5220d0
    this%prms_bonds(2,28,1)=367.0d0
    this%prms_bonds(2,28,2)=1.0800d0
    this%prms_bonds(2,29,1)=367.0d0
    this%prms_bonds(2,29,2)=1.0800d0
    this%prms_bonds(2,44,1)=490.0d0
    this%prms_bonds(2,44,2)=1.3350d0
    this%prms_bonds(2,45,1)=424.0d0
    this%prms_bonds(2,45,2)=1.3830d0
    this%prms_bonds(2,50,1)=457.0d0
    this%prms_bonds(2,50,2)=1.3580d0
    this%prms_bonds(2,53,1)=570.0d0
    this%prms_bonds(2,53,2)=1.2290d0
    this%prms_bonds(2,54,1)=656.0d0
    this%prms_bonds(2,54,2)=1.2500d0
    this%prms_bonds(2,55,1)=450.0d0
    this%prms_bonds(2,55,2)=1.3640d0
    this%prms_bonds(2,56,1)=450.0d0
    this%prms_bonds(2,56,2)=1.3230d0
    this%prms_bonds(2,5,1)=418.0d0
    this%prms_bonds(2,5,2)=1.3880d0
    this%prms_bonds(3,6,1)=388.0d0
    this%prms_bonds(3,6,2)=1.4590d0
    this%prms_bonds(3,16,1)=317.0d0
    this%prms_bonds(3,16,2)=1.4950d0
    this%prms_bonds(3,19,1)=546.0d0
    this%prms_bonds(3,19,2)=1.3520d0
    this%prms_bonds(3,31,1)=367.0d0
    this%prms_bonds(3,31,2)=1.0800d0
    this%prms_bonds(5,5,1)=469.0d0
    this%prms_bonds(5,5,2)=1.4000d0
    this%prms_bonds(5,6,1)=469.0d0
    this%prms_bonds(5,6,2)=1.4040d0
    this%prms_bonds(5,11,1)=427.0d0
    this%prms_bonds(5,11,2)=1.4330d0
    this%prms_bonds(5,12,1)=469.0d0
    this%prms_bonds(5,12,2)=1.4000d0
    this%prms_bonds(5,16,1)=317.0d0
    this%prms_bonds(5,16,2)=1.5100d0
    this%prms_bonds(5,28,1)=367.0d0
    this%prms_bonds(5,28,2)=1.0800d0
    this%prms_bonds(5,30,1)=367.0d0
    this%prms_bonds(5,30,2)=1.0800d0
    this%prms_bonds(5,46,1)=481.0d0
    this%prms_bonds(5,46,2)=1.3400d0
    this%prms_bonds(5,50,1)=483.0d0
    this%prms_bonds(5,50,2)=1.3390d0
    this%prms_bonds(5,55,1)=450.0d0
    this%prms_bonds(5,55,2)=1.3640d0
    this%prms_bonds(5,5,1)=427.0d0
    this%prms_bonds(5,5,2)=1.3810d0
    this%prms_bonds(6,6,1)=520.0d0
    this%prms_bonds(6,6,2)=1.3700d0
    this%prms_bonds(6,12,1)=447.0d0
    this%prms_bonds(6,12,2)=1.4190d0
    this%prms_bonds(6,45,1)=436.0d0
    this%prms_bonds(6,45,2)=1.3740d0
    this%prms_bonds(6,49,1)=414.0d0
    this%prms_bonds(6,49,2)=1.3910d0
    this%prms_bonds(6,50,1)=461.0d0
    this%prms_bonds(6,50,2)=1.3540d0
    this%prms_bonds(7,16,1)=317.0d0
    this%prms_bonds(7,16,2)=1.5040d0
    this%prms_bonds(7,18,1)=512.0d0
    this%prms_bonds(7,18,2)=1.3750d0
    this%prms_bonds(7,19,1)=518.0d0
    this%prms_bonds(7,19,2)=1.3710d0
    this%prms_bonds(7,49,1)=410.0d0
    this%prms_bonds(7,49,2)=1.3940d0
    this%prms_bonds(7,5,1)=422.0d0
    this%prms_bonds(7,5,2)=1.3850d0
    this%prms_bonds(8,8,1)=469.0d0
    this%prms_bonds(8,8,2)=1.4000d0
    this%prms_bonds(8,11,1)=549.0d0
    this%prms_bonds(8,11,2)=1.3500d0
    this%prms_bonds(8,16,1)=317.0d0
    this%prms_bonds(8,16,2)=1.5100d0
    this%prms_bonds(8,30,1)=367.0d0
    this%prms_bonds(8,30,2)=1.0800d0
    this%prms_bonds(9,29,1)=367.0d0
    this%prms_bonds(9,29,2)=1.0800d0
    this%prms_bonds(9,45,1)=440.0d0
    this%prms_bonds(9,45,2)=1.3710d0
    this%prms_bonds(9,49,1)=529.0d0
    this%prms_bonds(9,49,2)=1.3040d0
    this%prms_bonds(10,5,1)=193.0d0
    this%prms_bonds(10,5,2)=1.7270d0
    this%prms_bonds(10,16,1)=232.0d0
    this%prms_bonds(10,16,2)=1.7660d0
    this%prms_bonds(11,11,1)=549.0d0
    this%prms_bonds(11,11,2)=1.3500d0
    this%prms_bonds(11,16,1)=317.0d0
    this%prms_bonds(11,16,2)=1.5100d0
    this%prms_bonds(11,28,1)=367.0d0
    this%prms_bonds(11,28,2)=1.0800d0
    this%prms_bonds(11,29,1)=367.0d0
    this%prms_bonds(11,29,2)=1.0800d0
    this%prms_bonds(11,30,1)=367.0d0
    this%prms_bonds(11,30,2)=1.0800d0
    this%prms_bonds(11,45,1)=448.0d0
    this%prms_bonds(11,45,2)=1.3650d0
    this%prms_bonds(11,56,1)=480.0d0
    this%prms_bonds(11,56,2)=1.2400d0
    this%prms_bonds(12,5,1)=428.0d0
    this%prms_bonds(12,5,2)=1.3800d0
    this%prms_bonds(13,29,1)=367.0d0
    this%prms_bonds(13,29,2)=1.0800d0
    this%prms_bonds(13,50,1)=502.0d0
    this%prms_bonds(13,50,2)=1.3240d0
    this%prms_bonds(14,29,1)=367.0d0
    this%prms_bonds(14,29,2)=1.0800d0
    this%prms_bonds(14,49,1)=488.0d0
    this%prms_bonds(14,49,2)=1.3350d0
    this%prms_bonds(14,5,1)=477.0d0
    this%prms_bonds(14,5,2)=1.3430d0
    this%prms_bonds(16,16,1)=310.0d0
    this%prms_bonds(16,16,2)=1.5260d0
    this%prms_bonds(16,20,1)=400.0d0
    this%prms_bonds(16,20,2)=1.4580d0
    this%prms_bonds(16,21,1)=400.0d0
    this%prms_bonds(16,21,2)=1.4590d0
    this%prms_bonds(16,25,1)=340.0d0
    this%prms_bonds(16,25,2)=1.0900d0
    this%prms_bonds(16,26,1)=340.0d0
    this%prms_bonds(16,26,2)=1.0900d0
    this%prms_bonds(16,27,1)=340.0d0
    this%prms_bonds(16,27,2)=1.0900d0
    this%prms_bonds(16,31,1)=340.0d0
    this%prms_bonds(16,31,2)=1.0900d0
    this%prms_bonds(16,33,1)=340.0d0
    this%prms_bonds(16,33,2)=1.0900d0
    this%prms_bonds(16,44,1)=337.0d0
    this%prms_bonds(16,44,2)=1.4490d0
    this%prms_bonds(16,45,1)=337.0d0
    this%prms_bonds(16,45,2)=1.4750d0
    this%prms_bonds(16,46,1)=337.0d0
    this%prms_bonds(16,46,2)=1.4630d0
    this%prms_bonds(16,47,1)=367.0d0
    this%prms_bonds(16,47,2)=1.4710d0
    this%prms_bonds(16,51,1)=367.0d0
    this%prms_bonds(16,51,2)=1.4710d0
    this%prms_bonds(16,55,1)=320.0d0
    this%prms_bonds(16,55,2)=1.4100d0
    this%prms_bonds(16,56,1)=320.0d0
    this%prms_bonds(16,56,2)=1.4100d0
    this%prms_bonds(16,60,1)=227.0d0
    this%prms_bonds(16,60,2)=1.8100d0
    this%prms_bonds(16,61,1)=237.0d0
    this%prms_bonds(16,61,2)=1.8100d0
    this%prms_bonds(18,28,1)=367.0d0
    this%prms_bonds(18,28,2)=1.0800d0
    this%prms_bonds(18,49,1)=410.0d0
    this%prms_bonds(18,49,2)=1.3940d0
    this%prms_bonds(19,28,1)=367.0d0
    this%prms_bonds(19,28,2)=1.0800d0
    this%prms_bonds(19,5,1)=427.0d0
    this%prms_bonds(19,5,2)=1.3810d0
    this%prms_bonds(20,52,1)=600.0d0
    this%prms_bonds(20,52,2)=1.1500d0
    this%prms_bonds(21,21,1)=600.0d0
    this%prms_bonds(21,21,2)=1.2060d0
    this%prms_bonds(21,36,1)=400.0d0
    this%prms_bonds(21,36,2)=1.0560d0
    this%prms_bonds(22,5,1)=386.0d0
    this%prms_bonds(22,5,2)=1.3590d0
    this%prms_bonds(22,16,1)=367.0d0
    this%prms_bonds(22,16,2)=1.3800d0
    this%prms_bonds(24,44,1)=434.0d0
    this%prms_bonds(24,44,2)=1.0100d0
    this%prms_bonds(24,45,1)=434.0d0
    this%prms_bonds(24,45,2)=1.0100d0
    this%prms_bonds(24,46,1)=434.0d0
    this%prms_bonds(24,46,2)=1.0100d0
    this%prms_bonds(24,47,1)=434.0d0
    this%prms_bonds(24,47,2)=1.0100d0
    this%prms_bonds(24,51,1)=434.0d0
    this%prms_bonds(24,51,2)=1.0100d0
    this%prms_bonds(24,5,1)=434.0d0
    this%prms_bonds(24,5,2)=1.0100d0
    this%prms_bonds(32,55,1)=553.0d0
    this%prms_bonds(32,55,2)=0.9600d0
    this%prms_bonds(32,56,1)=553.0d0
    this%prms_bonds(32,56,2)=0.9600d0
    this%prms_bonds(34,61,1)=274.0d0
    this%prms_bonds(34,61,2)=1.3360d0
    this%prms_bonds(37,5,1)=171.0d0
    this%prms_bonds(37,5,2)=2.0750d0
    this%prms_bonds(37,16,1)=148.0d0
    this%prms_bonds(37,16,2)=2.1660d0
    this%prms_bonds(1,47,1)=600.0d0
    this%prms_bonds(1,47,2)=0.2000d0
    this%prms_bonds(1,49,1)=600.0d0
    this%prms_bonds(1,49,2)=0.2000d0
    this%prms_bonds(1,50,1)=600.0d0
    this%prms_bonds(1,50,2)=0.2000d0
    this%prms_bonds(1,51,1)=600.0d0
    this%prms_bonds(1,51,2)=0.2000d0
    this%prms_bonds(1,53,1)=600.0d0
    this%prms_bonds(1,53,2)=0.2000d0
    this%prms_bonds(1,55,1)=600.0d0
    this%prms_bonds(1,55,2)=0.2000d0
    this%prms_bonds(1,56,1)=600.0d0
    this%prms_bonds(1,56,2)=0.2000d0
    this%prms_bonds(1,60,1)=600.0d0
    this%prms_bonds(1,60,2)=0.7000d0
    this%prms_bonds(1,61,1)=600.0d0
    this%prms_bonds(1,61,2)=0.7000d0
    this%prms_bonds(54,58,1)=525.0d0
    this%prms_bonds(54,58,2)=1.4800d0
    this%prms_bonds(55,58,1)=230.0d0
    this%prms_bonds(55,58,2)=1.6100d0
    this%prms_bonds(56,58,1)=230.0d0
    this%prms_bonds(56,58,2)=1.6100d0
    this%prms_bonds(57,35,1)=553.0d0
    this%prms_bonds(57,35,2)=0.9572d0
    this%prms_bonds(60,60,1)=166.0d0
    this%prms_bonds(60,60,2)=2.0380d0
 end subroutine set_amberbonds

!  subroutine set_amberbonds(this,p1,p2)
!    class(amber), intent(inout) :: this
!    integer                     :: i
!    real(4) x1,x2
!    character(2) pa,pb,p1,p2
!    open(12,file='/tmp/amber/amber_bonds.prm',status='old')
!    this%prms_bonds(1)=0.d0
!    this%prms_bonds(2)=0.d0
!    do i=1,115
!       read(12,*,end=1)pa,pb,x1,x2
!       if(pa.eq.p1.and.pb.eq.p2.or.pa.eq.p2.and.pb.eq.p1)then
!          this%prms_bonds(1)=dble(x1)
!          this%prms_bonds(2)=dble(x2)
!       end if
!    end do
!1   close(12)
!  end subroutine set_amberbonds

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
