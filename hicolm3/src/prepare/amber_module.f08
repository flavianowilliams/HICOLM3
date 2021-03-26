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
     integer                  :: natp
     real(8)                  :: prms_tors(4)
     real(8), allocatable     :: prms_vdw(:,:)
     real(8), allocatable     :: prms_angles(:,:,:,:)
     real(8), allocatable     :: prms_bonds(:,:,:)
     character(2),allocatable :: atp(:)
   contains
     procedure :: set_natp
     procedure :: get_natp
     procedure :: set_amberbonds
     procedure :: set_amberangles
     procedure :: set_amberdihedrals
     procedure :: set_ambervdw
     procedure :: set_ambertypes
     generic   :: set_amber => set_amberdihedrals
  end type amber

contains

  subroutine set_natp(this,natp)
    class(amber), intent(inout) :: this
    integer, intent(in)         :: natp
    this%natp=natp
  end subroutine set_natp

  integer function get_natp(this)
    class(amber), intent(inout) :: this
    get_natp=this%natp
  end function get_natp

  subroutine set_ambertypes(this)
    class(amber), intent(inout) :: this
    allocate(this%atp(this%get_natp()))
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
    allocate(this%prms_bonds(this%get_natp(),this%get_natp(),2))
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

  subroutine set_amberangles(this)
    class(amber), intent(inout) :: this
    allocate(this%prms_angles(this%get_natp(),this%get_natp(),this%get_natp(),2))
    this%prms_angles(1,5,5,1)= 70.0d0
    this%prms_angles(1,5,5,2)= 118.8d0
    this%prms_angles(1,16,16,1)= 50.0d0
    this%prms_angles(1,16,16,2)= 108.0d0
    this%prms_angles(1,16,25,1)= 50.0d0
    this%prms_angles(1,16,25,2)= 106.5d0
    this%prms_angles(2,2,53,1)= 80.0d0
    this%prms_angles(2,2,53,2)= 120.0d0
    this%prms_angles(2,2,55,1)= 80.0d0
    this%prms_angles(2,2,55,2)= 120.0d0
    this%prms_angles(2,5,5,1)= 63.0d0
    this%prms_angles(2,5,5,2)= 120.0d0
    this%prms_angles(2,5,30,1)= 50.0d0
    this%prms_angles(2,5,30,2)= 120.0d0
    this%prms_angles(2,6,6,1)= 63.0d0
    this%prms_angles(2,6,6,2)= 119.2d0
    this%prms_angles(2,6,49,1)= 70.0d0
    this%prms_angles(2,6,49,2)= 130.0d0
    this%prms_angles(2,11,11,1)= 63.0d0
    this%prms_angles(2,11,11,2)= 120.7d0
    this%prms_angles(2,11,16,1)= 70.0d0
    this%prms_angles(2,11,16,2)= 119.7d0
    this%prms_angles(2,11,28,1)= 50.0d0
    this%prms_angles(2,11,28,2)= 119.7d0
    this%prms_angles(2,11,30,1)= 50.0d0
    this%prms_angles(2,11,30,2)= 119.7d0
    this%prms_angles(2,16,16,1)= 63.0d0
    this%prms_angles(2,16,16,2)= 111.1d0
    this%prms_angles(2,16,25,1)= 50.0d0
    this%prms_angles(2,16,25,2)= 109.5d0
    this%prms_angles(2,16,31,1)= 50.0d0
    this%prms_angles(2,16,31,2)= 109.5d0
    this%prms_angles(2,16,33,1)= 50.0d0
    this%prms_angles(2,16,33,2)= 109.5d0
    this%prms_angles(2,16,44,1)= 63.0d0
    this%prms_angles(2,16,44,2)= 110.1d0
    this%prms_angles(2,16,47,1)= 80.0d0
    this%prms_angles(2,16,47,2)= 111.2d0
    this%prms_angles(2,16,56,1)= 60.0d0
    this%prms_angles(2,16,56,2)= 109.5d0
    this%prms_angles(2,44,16,1)= 50.0d0
    this%prms_angles(2,44,16,2)= 121.9d0
    this%prms_angles(2,44,24,1)= 50.0d0
    this%prms_angles(2,44,24,2)= 120.0d0
    this%prms_angles(2,45,11,1)= 70.0d0
    this%prms_angles(2,45,11,2)= 121.6d0
    this%prms_angles(2,45,16,1)= 70.0d0
    this%prms_angles(2,45,16,2)= 117.6d0
    this%prms_angles(2,45,24,1)= 50.0d0
    this%prms_angles(2,45,24,2)= 119.2d0
    this%prms_angles(2,50,5,1)= 70.0d0
    this%prms_angles(2,50,5,2)= 120.5d0
    this%prms_angles(2,50,5,1)=150.0d0
    this%prms_angles(2,50,5,2)= 120.0d0
    this%prms_angles(2,53,5,1)=150.0d0
    this%prms_angles(2,53,5,2)= 120.0d0
    this%prms_angles(2,55,32,1)= 50.0d0
    this%prms_angles(2,55,32,2)= 113.0d0
    this%prms_angles(2,55,5,1)=150.0d0
    this%prms_angles(2,55,5,2)= 120.0d0
    this%prms_angles(2,56,16,1)= 60.0d0
    this%prms_angles(2,56,16,2)= 117.0d0
    this%prms_angles(2,56,5,1)=150.0d0
    this%prms_angles(2,56,5,2)= 109.5d0
    this%prms_angles(2,5,2,1)= 70.0d0
    this%prms_angles(2,5,2,2)= 126.4d0
    this%prms_angles(2,5,5,1)= 70.0d0
    this%prms_angles(2,5,5,2)= 125.2d0
    this%prms_angles(2,5,24,1)= 50.0d0
    this%prms_angles(2,5,24,2)= 116.8d0
    this%prms_angles(3,6,5,1)= 63.0d0
    this%prms_angles(3,6,5,2)= 134.9d0
    this%prms_angles(3,6,12,1)= 63.0d0
    this%prms_angles(3,6,12,2)= 108.8d0
    this%prms_angles(3,16,16,1)= 63.0d0
    this%prms_angles(3,16,16,2)= 115.6d0
    this%prms_angles(3,16,31,1)= 50.0d0
    this%prms_angles(3,16,31,2)= 109.5d0
    this%prms_angles(3,19,28,1)= 50.0d0
    this%prms_angles(3,19,28,2)= 120.0d0
    this%prms_angles(3,19,5,1)= 70.0d0
    this%prms_angles(3,19,5,2)= 108.7d0
    this%prms_angles(5,2,5,1)= 63.0d0
    this%prms_angles(5,2,5,2)= 120.0d0
    this%prms_angles(5,2,55,1)= 70.0d0
    this%prms_angles(5,2,55,2)= 120.0d0
    this%prms_angles(5,5,5,1)= 63.0d0
    this%prms_angles(5,5,5,2)= 120.0d0
    this%prms_angles(5,5,6,1)= 63.0d0
    this%prms_angles(5,5,6,2)= 120.0d0
    this%prms_angles(5,5,12,1)= 63.0d0
    this%prms_angles(5,5,12,2)= 120.0d0
    this%prms_angles(5,5,16,1)= 70.0d0
    this%prms_angles(5,5,16,2)= 120.0d0
    this%prms_angles(5,5,28,1)= 50.0d0
    this%prms_angles(5,5,28,2)= 120.0d0
    this%prms_angles(5,5,30,1)= 50.0d0
    this%prms_angles(5,5,30,2)= 120.0d0
    this%prms_angles(5,5,55,1)= 70.0d0
    this%prms_angles(5,5,55,2)= 120.0d0
    this%prms_angles(5,6,6,1)= 63.0d0
    this%prms_angles(5,6,6,2)= 117.3d0
    this%prms_angles(5,6,12,1)= 63.0d0
    this%prms_angles(5,6,12,2)= 116.2d0
    this%prms_angles(5,6,49,1)= 70.0d0
    this%prms_angles(5,6,49,2)= 132.4d0
    this%prms_angles(5,11,11,1)= 63.0d0
    this%prms_angles(5,11,11,2)= 117.0d0
    this%prms_angles(5,11,28,1)= 50.0d0
    this%prms_angles(5,11,28,2)= 123.3d0
    this%prms_angles(5,11,30,1)= 50.0d0
    this%prms_angles(5,11,30,2)= 123.3d0
    this%prms_angles(5,12,6,1)= 63.0d0
    this%prms_angles(5,12,6,2)= 122.7d0
    this%prms_angles(5,12,5,1)= 70.0d0
    this%prms_angles(5,12,5,2)= 132.8d0
    this%prms_angles(5,16,31,1)= 50.0d0
    this%prms_angles(5,16,31,2)= 109.5d0
    this%prms_angles(5,46,16,1)= 50.0d0
    this%prms_angles(5,46,16,2)= 123.2d0
    this%prms_angles(5,46,24,1)= 50.0d0
    this%prms_angles(5,46,24,2)= 120.0d0
    this%prms_angles(5,50,6,1)= 70.0d0
    this%prms_angles(5,50,6,2)= 112.2d0
    this%prms_angles(5,50,13,1)= 70.0d0
    this%prms_angles(5,50,13,2)= 118.6d0
    this%prms_angles(5,50,5,1)=150.0d0
    this%prms_angles(5,50,5,2)= 120.0d0
    this%prms_angles(5,55,32,1)= 50.0d0
    this%prms_angles(5,55,32,2)= 113.0d0
    this%prms_angles(5,5,24,1)= 50.0d0
    this%prms_angles(5,5,24,2)= 118.0d0
    this%prms_angles(6,2,53,1)= 80.0d0
    this%prms_angles(6,2,53,2)= 128.8d0
    this%prms_angles(6,2,5,1)= 70.0d0
    this%prms_angles(6,2,5,2)= 111.3d0
    this%prms_angles(6,3,16,1)= 70.0d0
    this%prms_angles(6,3,16,2)= 128.6d0
    this%prms_angles(6,3,19,1)= 63.0d0
    this%prms_angles(6,3,19,2)= 106.4d0
    this%prms_angles(6,5,28,1)= 50.0d0
    this%prms_angles(6,5,28,2)= 120.0d0
    this%prms_angles(6,5,30,1)= 50.0d0
    this%prms_angles(6,5,30,2)= 120.0d0
    this%prms_angles(6,5,46,1)= 70.0d0
    this%prms_angles(6,5,46,2)= 123.5d0
    this%prms_angles(6,5,50,1)= 70.0d0
    this%prms_angles(6,5,50,2)= 117.3d0
    this%prms_angles(6,6,45,1)= 70.0d0
    this%prms_angles(6,6,45,2)= 106.2d0
    this%prms_angles(6,6,49,1)= 70.0d0
    this%prms_angles(6,6,49,2)= 110.4d0
    this%prms_angles(6,6,50,1)= 70.0d0
    this%prms_angles(6,6,50,2)= 127.7d0
    this%prms_angles(6,12,5,1)= 70.0d0
    this%prms_angles(6,12,5,2)= 104.4d0
    this%prms_angles(6,45,9,1)= 70.0d0
    this%prms_angles(6,45,9,2)= 105.4d0
    this%prms_angles(6,45,16,1)= 70.0d0
    this%prms_angles(6,45,16,2)= 125.8d0
    this%prms_angles(6,45,24,1)= 50.0d0
    this%prms_angles(6,45,24,2)= 125.8d0
    this%prms_angles(6,49,9,1)= 70.0d0
    this%prms_angles(6,49,9,2)= 103.8d0
    this%prms_angles(6,49,5,1)=150.0d0
    this%prms_angles(6,49,5,2)= 126.0d0
    this%prms_angles(6,50,13,1)= 70.0d0
    this%prms_angles(6,50,13,2)= 111.0d0
    this%prms_angles(6,50,5,1)=150.0d0
    this%prms_angles(6,50,5,2)= 120.0d0
    this%prms_angles(7,16,16,1)= 63.0d0
    this%prms_angles(7,16,16,2)= 113.1d0
    this%prms_angles(7,16,31,1)= 50.0d0
    this%prms_angles(7,16,31,2)= 109.5d0
    this%prms_angles(7,18,28,1)= 50.0d0
    this%prms_angles(7,18,28,2)= 120.0d0
    this%prms_angles(7,18,49,1)= 70.0d0
    this%prms_angles(7,18,49,2)= 120.0d0
    this%prms_angles(7,19,28,1)= 50.0d0
    this%prms_angles(7,19,28,2)= 120.0d0
    this%prms_angles(7,19,5,1)= 70.0d0
    this%prms_angles(7,19,5,2)= 120.0d0
    this%prms_angles(7,49,14,1)= 70.0d0
    this%prms_angles(7,49,14,2)= 117.0d0
    this%prms_angles(7,49,5,1)=150.0d0
    this%prms_angles(7,49,5,2)= 126.0d0
    this%prms_angles(7,5,14,1)= 70.0d0
    this%prms_angles(7,5,14,2)= 120.0d0
    this%prms_angles(7,5,24,1)= 50.0d0
    this%prms_angles(7,5,24,2)= 120.0d0
    this%prms_angles(8,8,11,1)= 63.0d0
    this%prms_angles(8,8,11,2)= 120.0d0
    this%prms_angles(8,8,16,1)= 70.0d0
    this%prms_angles(8,8,16,2)= 120.0d0
    this%prms_angles(9,45,16,1)= 70.0d0
    this%prms_angles(9,45,16,2)= 128.8d0
    this%prms_angles(9,45,24,1)= 50.0d0
    this%prms_angles(9,45,24,2)= 128.8d0
    this%prms_angles(9,49,5,1)=150.0d0
    this%prms_angles(9,49,5,2)= 126.0d0
    this%prms_angles(10,5,5,1)= 70.0d0
    this%prms_angles(10,5,5,2)= 118.8d0
    this%prms_angles(10,16,16,1)= 50.0d0
    this%prms_angles(10,16,16,2)= 108.5d0
    this%prms_angles(10,16,25,1)= 50.0d0
    this%prms_angles(10,16,25,2)= 108.5d0
    this%prms_angles(11,2,53,1)= 80.0d0
    this%prms_angles(11,2,53,2)= 125.3d0
    this%prms_angles(11,2,5,1)= 70.0d0
    this%prms_angles(11,2,5,2)= 114.1d0
    this%prms_angles(11,5,46,1)= 70.0d0
    this%prms_angles(11,5,46,2)= 120.1d0
    this%prms_angles(11,5,50,1)= 70.0d0
    this%prms_angles(11,5,50,2)= 121.5d0
    this%prms_angles(11,8,16,1)= 70.0d0
    this%prms_angles(11,8,16,2)= 120.0d0
    this%prms_angles(11,11,16,1)= 70.0d0
    this%prms_angles(11,11,16,2)= 119.7d0
    this%prms_angles(11,11,28,1)= 50.0d0
    this%prms_angles(11,11,28,2)= 119.7d0
    this%prms_angles(11,11,30,1)= 50.0d0
    this%prms_angles(11,11,30,2)= 119.7d0
    this%prms_angles(11,11,45,1)= 70.0d0
    this%prms_angles(11,11,45,2)= 121.2d0
    this%prms_angles(11,11,56,1)= 80.0d0
    this%prms_angles(11,11,56,2)= 125.0d0
    this%prms_angles(11,16,16,1)= 63.0d0
    this%prms_angles(11,16,16,2)= 111.0d0
    this%prms_angles(11,16,56,1)= 50.0d0
    this%prms_angles(11,16,56,2)= 109.5d0
    this%prms_angles(11,45,16,1)= 70.0d0
    this%prms_angles(11,45,16,2)= 121.2d0
    this%prms_angles(11,45,24,1)= 50.0d0
    this%prms_angles(11,45,24,2)= 121.2d0
    this%prms_angles(11,56,16,1)= 60.0d0
    this%prms_angles(11,56,16,2)= 117.0d0
    this%prms_angles(11,56,5,1)=150.0d0
    this%prms_angles(11,56,5,2)= 109.5d0
    this%prms_angles(12,5,30,1)= 50.0d0
    this%prms_angles(12,5,30,2)= 120.0d0
    this%prms_angles(12,5,19,1)= 70.0d0
    this%prms_angles(12,5,19,2)= 111.6d0
    this%prms_angles(12,5,24,1)= 50.0d0
    this%prms_angles(12,5,24,2)= 123.1d0
    this%prms_angles(13,50,5,1)=150.0d0
    this%prms_angles(13,50,5,2)= 120.0d0
    this%prms_angles(14,49,18,1)= 70.0d0
    this%prms_angles(14,49,18,2)= 117.0d0
    this%prms_angles(14,49,5,1)=150.0d0
    this%prms_angles(14,49,5,2)= 126.0d0
    this%prms_angles(14,5,19,1)= 70.0d0
    this%prms_angles(14,5,19,2)= 120.0d0
    this%prms_angles(14,5,24,1)= 50.0d0
    this%prms_angles(14,5,24,2)= 120.0d0
    this%prms_angles(16,2,16,1)= 63.0d0
    this%prms_angles(16,2,16,2)= 117.0d0
    this%prms_angles(16,2,44,1)= 70.0d0
    this%prms_angles(16,2,44,2)= 116.6d0
    this%prms_angles(16,2,53,1)= 80.0d0
    this%prms_angles(16,2,53,2)= 120.4d0
    this%prms_angles(16,2,54,1)= 70.0d0
    this%prms_angles(16,2,54,2)= 117.0d0
    this%prms_angles(16,2,55,1)= 80.0d0
    this%prms_angles(16,2,55,2)= 110.0d0
    this%prms_angles(16,2,56,1)= 80.0d0
    this%prms_angles(16,2,56,2)= 115.0d0
    this%prms_angles(16,3,19,1)= 70.0d0
    this%prms_angles(16,3,19,2)= 125.0d0
    this%prms_angles(16,7,18,1)= 70.0d0
    this%prms_angles(16,7,18,2)= 120.0d0
    this%prms_angles(16,7,19,1)= 70.0d0
    this%prms_angles(16,7,19,2)= 120.0d0
    this%prms_angles(16,7,49,1)= 70.0d0
    this%prms_angles(16,7,49,2)= 120.0d0
    this%prms_angles(16,7,5,1)= 70.0d0
    this%prms_angles(16,7,5,2)= 120.0d0
    this%prms_angles(16,16,5,1)= 63.0d0
    this%prms_angles(16,16,5,2)= 114.0d0
    this%prms_angles(16,16,16,1)= 40.0d0
    this%prms_angles(16,16,16,2)= 109.5d0
    this%prms_angles(16,16,20,1)= 63.0d0
    this%prms_angles(16,16,20,2)= 110.0d0
    this%prms_angles(16,16,21,1)= 63.0d0
    this%prms_angles(16,16,21,2)= 110.0d0
    this%prms_angles(16,16,25,1)= 50.0d0
    this%prms_angles(16,16,25,2)= 109.5d0
    this%prms_angles(16,16,26,1)= 50.0d0
    this%prms_angles(16,16,26,2)= 109.5d0
    this%prms_angles(16,16,31,1)= 50.0d0
    this%prms_angles(16,16,31,2)= 109.5d0
    this%prms_angles(16,16,33,1)= 50.0d0
    this%prms_angles(16,16,33,2)= 109.5d0
    this%prms_angles(16,16,44,1)= 80.0d0
    this%prms_angles(16,16,44,2)= 109.7d0
    this%prms_angles(16,16,45,1)= 50.0d0
    this%prms_angles(16,16,45,2)= 109.5d0
    this%prms_angles(16,16,46,1)= 80.0d0
    this%prms_angles(16,16,46,2)= 111.2d0
    this%prms_angles(16,16,47,1)= 80.0d0
    this%prms_angles(16,16,47,2)= 111.2d0
    this%prms_angles(16,16,51,1)= 80.0d0
    this%prms_angles(16,16,51,2)= 111.2d0
    this%prms_angles(16,16,55,1)= 50.0d0
    this%prms_angles(16,16,55,2)= 109.5d0
    this%prms_angles(16,16,56,1)= 50.0d0
    this%prms_angles(16,16,56,2)= 109.5d0
    this%prms_angles(16,16,60,1)= 50.0d0
    this%prms_angles(16,16,60,2)= 114.7d0
    this%prms_angles(16,16,61,1)= 50.0d0
    this%prms_angles(16,16,61,2)= 108.6d0
    this%prms_angles(16,20,52,1)= 80.0d0
    this%prms_angles(16,20,52,2)= 180.0d0
    this%prms_angles(16,21,21,1)= 80.0d0
    this%prms_angles(16,21,21,2)= 180.0d0
    this%prms_angles(16,44,16,1)= 50.0d0
    this%prms_angles(16,44,16,2)= 118.0d0
    this%prms_angles(16,44,24,1)= 50.0d0
    this%prms_angles(16,44,24,2)= 118.0d0
    this%prms_angles(16,46,24,1)= 50.0d0
    this%prms_angles(16,46,24,2)= 118.4d0
    this%prms_angles(16,47,16,1)= 50.0d0
    this%prms_angles(16,47,16,2)= 109.5d0
    this%prms_angles(16,47,24,1)= 50.0d0
    this%prms_angles(16,47,24,2)= 109.5d0
    this%prms_angles(16,47,5,1)=150.0d0
    this%prms_angles(16,47,5,2)= 109.5d0
    this%prms_angles(16,51,16,1)= 50.0d0
    this%prms_angles(16,51,16,2)= 109.5d0
    this%prms_angles(16,51,24,1)= 50.0d0
    this%prms_angles(16,51,24,2)= 109.5d0
    this%prms_angles(16,51,5,1)=150.0d0
    this%prms_angles(16,51,5,2)= 109.5d0
    this%prms_angles(16,55,32,1)= 55.0d0
    this%prms_angles(16,55,32,2)= 108.5d0
    this%prms_angles(16,55,5,1)=150.0d0
    this%prms_angles(16,55,5,2)= 109.5d0
    this%prms_angles(16,56,16,1)= 60.0d0
    this%prms_angles(16,56,16,2)= 109.5d0
    this%prms_angles(16,56,5,1)=150.0d0
    this%prms_angles(16,56,5,2)= 109.5d0
    this%prms_angles(16,56,58,1)=100.0d0
    this%prms_angles(16,56,58,2)= 120.5d0
    this%prms_angles(16,60,16,1)= 62.0d0
    this%prms_angles(16,60,16,2)=  98.9d0
    this%prms_angles(16,60,5,1)=150.0d0
    this%prms_angles(16,60,5,2)=  90.0d0
    this%prms_angles(16,60,60,1)= 68.0d0
    this%prms_angles(16,60,60,2)= 103.7d0
    this%prms_angles(16,61,34,1)= 43.0d0
    this%prms_angles(16,61,34,2)=  96.0d0
    this%prms_angles(16,61,5,1)=150.0d0
    this%prms_angles(16,61,5,2)=  90.0d0
    this%prms_angles(18,7,5,1)= 70.0d0
    this%prms_angles(18,7,5,2)= 120.0d0
    this%prms_angles(18,49,5,1)=150.0d0
    this%prms_angles(18,49,5,2)= 126.0d0
    this%prms_angles(19,7,49,1)= 70.0d0
    this%prms_angles(19,7,49,2)= 120.0d0
    this%prms_angles(19,7,5,1)= 70.0d0
    this%prms_angles(19,7,5,2)= 120.0d0
    this%prms_angles(19,5,24,1)= 50.0d0
    this%prms_angles(19,5,24,2)= 120.0d0
    this%prms_angles(21,21,36,1)= 50.0d0
    this%prms_angles(21,21,36,2)= 180.0d0
    this%prms_angles(22,5,5,1)= 70.0d0
    this%prms_angles(22,5,5,2)= 121.0d0
    this%prms_angles(22,16,16,1)= 50.0d0
    this%prms_angles(22,16,16,2)= 109.0d0
    this%prms_angles(22,16,22,1)= 77.0d0
    this%prms_angles(22,16,22,2)= 109.1d0
    this%prms_angles(22,16,25,1)= 50.0d0
    this%prms_angles(22,16,25,2)= 109.5d0
    this%prms_angles(22,16,26,1)= 50.0d0
    this%prms_angles(22,16,26,2)= 109.5d0
    this%prms_angles(24,44,24,1)= 35.0d0
    this%prms_angles(24,44,24,2)= 120.0d0
    this%prms_angles(24,46,24,1)= 35.0d0
    this%prms_angles(24,46,24,2)= 120.0d0
    this%prms_angles(24,47,24,1)= 35.0d0
    this%prms_angles(24,47,24,2)= 109.5d0
    this%prms_angles(24,47,5,1)=150.0d0
    this%prms_angles(24,47,5,2)= 109.5d0
    this%prms_angles(24,51,24,1)= 35.0d0
    this%prms_angles(24,51,24,2)= 109.5d0
    this%prms_angles(24,51,5,1)=150.0d0
    this%prms_angles(24,51,5,2)= 109.5d0
    this%prms_angles(25,16,11,1)= 50.0d0
    this%prms_angles(25,16,11,2)= 109.5d0
    this%prms_angles(25,16,20,1)= 50.0d0
    this%prms_angles(25,16,20,2)= 110.0d0
    this%prms_angles(25,16,21,1)= 50.0d0
    this%prms_angles(25,16,21,2)= 110.0d0
    this%prms_angles(25,16,25,1)= 35.0d0
    this%prms_angles(25,16,25,2)= 109.5d0
    this%prms_angles(25,16,44,1)= 50.0d0
    this%prms_angles(25,16,44,2)= 109.5d0
    this%prms_angles(25,16,45,1)= 50.0d0
    this%prms_angles(25,16,45,2)= 109.5d0
    this%prms_angles(25,16,46,1)= 50.0d0
    this%prms_angles(25,16,46,2)= 109.5d0
    this%prms_angles(25,16,51,1)= 50.0d0
    this%prms_angles(25,16,51,2)= 109.5d0
    this%prms_angles(25,16,55,1)= 50.0d0
    this%prms_angles(25,16,55,2)= 109.5d0
    this%prms_angles(25,16,56,1)= 50.0d0
    this%prms_angles(25,16,56,2)= 109.5d0
    this%prms_angles(25,16,60,1)= 50.0d0
    this%prms_angles(25,16,60,2)= 109.5d0
    this%prms_angles(25,16,61,1)= 50.0d0
    this%prms_angles(25,16,61,2)= 109.5d0
    this%prms_angles(26,16,26,1)= 35.0d0
    this%prms_angles(26,16,26,2)= 109.5d0
    this%prms_angles(26,16,45,1)= 50.0d0
    this%prms_angles(26,16,45,2)= 109.5d0
    this%prms_angles(26,16,56,1)= 50.0d0
    this%prms_angles(26,16,56,2)= 109.5d0
    this%prms_angles(28,2,2,1)= 50.0d0
    this%prms_angles(28,2,2,2)= 120.0d0
    this%prms_angles(28,2,11,1)= 50.0d0
    this%prms_angles(28,2,11,2)= 115.0d0
    this%prms_angles(28,2,16,1)= 50.0d0
    this%prms_angles(28,2,16,2)= 115.0d0
    this%prms_angles(28,2,53,1)= 50.0d0
    this%prms_angles(28,2,53,2)= 120.0d0
    this%prms_angles(28,2,55,1)= 50.0d0
    this%prms_angles(28,2,55,2)= 120.0d0
    this%prms_angles(28,11,45,1)= 50.0d0
    this%prms_angles(28,11,45,2)= 119.1d0
    this%prms_angles(28,11,56,1)= 50.0d0
    this%prms_angles(28,11,56,2)= 113.0d0
    this%prms_angles(28,18,49,1)= 50.0d0
    this%prms_angles(28,18,49,2)= 120.0d0
    this%prms_angles(28,19,5,1)= 50.0d0
    this%prms_angles(28,19,5,2)= 120.0d0
    this%prms_angles(29,2,44,1)= 50.0d0
    this%prms_angles(29,2,44,2)= 120.0d0
    this%prms_angles(29,2,53,1)= 50.0d0
    this%prms_angles(29,2,53,2)= 119.0d0
    this%prms_angles(29,2,55,1)= 50.0d0
    this%prms_angles(29,2,55,2)= 107.0d0
    this%prms_angles(29,2,56,1)= 50.0d0
    this%prms_angles(29,2,56,2)= 107.0d0
    this%prms_angles(29,9,45,1)= 50.0d0
    this%prms_angles(29,9,45,2)= 123.0d0
    this%prms_angles(29,9,49,1)= 50.0d0
    this%prms_angles(29,9,49,2)= 123.0d0
    this%prms_angles(29,13,50,1)= 50.0d0
    this%prms_angles(29,13,50,2)= 115.5d0
    this%prms_angles(29,14,49,1)= 50.0d0
    this%prms_angles(29,14,49,2)= 120.0d0
    this%prms_angles(29,14,5,1)= 50.0d0
    this%prms_angles(29,14,5,2)= 120.0d0
    this%prms_angles(30,8,8,1)= 50.0d0
    this%prms_angles(30,8,8,2)= 120.0d0
    this%prms_angles(30,8,11,1)= 50.0d0
    this%prms_angles(30,8,11,2)= 120.0d0
    this%prms_angles(30,8,30,1)= 35.0d0
    this%prms_angles(30,8,30,2)= 119.0d0
    this%prms_angles(30,11,8,1)= 50.0d0
    this%prms_angles(30,11,8,2)= 120.0d0
    this%prms_angles(30,11,16,1)= 50.0d0
    this%prms_angles(30,11,16,2)= 120.0d0
    this%prms_angles(30,11,30,1)= 35.0d0
    this%prms_angles(30,11,30,2)= 120.0d0
    this%prms_angles(31,16,8,1)= 50.0d0
    this%prms_angles(31,16,8,2)= 109.5d0
    this%prms_angles(31,16,11,1)= 50.0d0
    this%prms_angles(31,16,11,2)= 109.5d0
    this%prms_angles(31,16,21,1)= 50.0d0
    this%prms_angles(31,16,21,2)= 110.0d0
    this%prms_angles(31,16,31,1)= 35.0d0
    this%prms_angles(31,16,31,2)= 109.5d0
    this%prms_angles(32,55,5,1)=150.0d0
    this%prms_angles(32,55,5,2)= 109.5d0
    this%prms_angles(32,55,58,1)= 45.0d0
    this%prms_angles(32,55,58,2)= 108.5d0
    this%prms_angles(33,16,33,1)= 35.0d0
    this%prms_angles(33,16,33,2)= 109.5d0
    this%prms_angles(33,16,47,1)= 50.0d0
    this%prms_angles(33,16,47,2)= 109.5d0
    this%prms_angles(34,61,34,1)= 35.0d0
    this%prms_angles(34,61,34,2)=  92.1d0
    this%prms_angles(34,61,5,1)=150.0d0
    this%prms_angles(34,61,5,2)=  90.0d0
    this%prms_angles(35,35,57,1)=  0.0d0
    this%prms_angles(35,35,57,2)= 127.7d0
    this%prms_angles(35,57,35,1)=100.0d0
    this%prms_angles(35,57,35,2)= 104.5d0
    this%prms_angles(37,5,5,1)= 70.0d0
    this%prms_angles(37,5,5,2)= 118.8d0
    this%prms_angles(37,16,16,1)= 50.0d0
    this%prms_angles(37,16,16,2)= 106.0d0
    this%prms_angles(1,53,5,1)=150.0d0
    this%prms_angles(1,53,5,2)= 120.0d0
    this%prms_angles(1,55,5,1)=150.0d0
    this%prms_angles(1,55,5,2)= 109.5d0
    this%prms_angles(1,56,5,1)=150.0d0
    this%prms_angles(1,56,5,2)= 109.5d0
    this%prms_angles(1,60,5,1)=150.0d0
    this%prms_angles(1,60,5,2)= 180.0d0
    this%prms_angles(1,61,5,1)=150.0d0
    this%prms_angles(1,61,5,2)= 180.0d0
    this%prms_angles(44,2,53,1)= 80.0d0
    this%prms_angles(44,2,53,2)= 122.9d0
    this%prms_angles(45,2,50,1)= 70.0d0
    this%prms_angles(45,2,50,2)= 118.6d0
    this%prms_angles(45,2,53,1)= 80.0d0
    this%prms_angles(45,2,53,2)= 120.9d0
    this%prms_angles(45,2,5,1)= 70.0d0
    this%prms_angles(45,2,5,2)= 115.4d0
    this%prms_angles(45,6,50,1)= 70.0d0
    this%prms_angles(45,6,50,2)= 126.2d0
    this%prms_angles(45,9,49,1)= 70.0d0
    this%prms_angles(45,9,49,2)= 113.9d0
    this%prms_angles(46,5,46,1)= 70.0d0
    this%prms_angles(46,5,46,2)= 120.0d0
    this%prms_angles(46,5,50,1)= 70.0d0
    this%prms_angles(46,5,50,2)= 119.3d0
    this%prms_angles(46,5,5,1)= 70.0d0
    this%prms_angles(46,5,5,2)= 116.0d0
    this%prms_angles(50,2,53,1)= 80.0d0
    this%prms_angles(50,2,53,2)= 122.5d0
    this%prms_angles(50,13,50,1)= 70.0d0
    this%prms_angles(50,13,50,2)= 129.1d0
    this%prms_angles(53,2,53,1)= 80.0d0
    this%prms_angles(53,2,53,2)= 126.0d0
    this%prms_angles(53,2,55,1)= 80.0d0
    this%prms_angles(53,2,55,2)= 120.0d0
    this%prms_angles(53,2,56,1)= 80.0d0
    this%prms_angles(53,2,56,2)= 125.0d0
    this%prms_angles(54,2,54,1)= 80.0d0
    this%prms_angles(54,2,54,2)= 126.0d0
    this%prms_angles(54,58,54,1)=140.0d0
    this%prms_angles(54,58,54,2)= 119.9d0
    this%prms_angles(54,58,55,1)= 45.0d0
    this%prms_angles(54,58,55,2)= 108.2d0
    this%prms_angles(54,58,56,1)=100.0d0
    this%prms_angles(54,58,56,2)= 108.2d0
    this%prms_angles(55,58,56,1)= 45.0d0
    this%prms_angles(55,58,56,2)= 102.6d0
    this%prms_angles(56,16,20,1)= 50.0d0
    this%prms_angles(56,16,20,2)= 110.0d0
    this%prms_angles(56,16,21,1)= 50.0d0
    this%prms_angles(56,16,21,2)= 110.0d0
    this%prms_angles(56,16,45,1)= 50.0d0
    this%prms_angles(56,16,45,2)= 109.5d0
    this%prms_angles(56,16,56,1)=160.0d0
    this%prms_angles(56,16,56,2)= 101.0d0
    this%prms_angles(56,58,56,1)= 45.0d0
    this%prms_angles(56,58,56,2)= 102.6d0
    this%prms_angles(58,56,5,1)=150.0d0
    this%prms_angles(58,56,5,2)= 109.5d0
    this%prms_angles(58,56,58,1)=100.0d0
    this%prms_angles(58,56,58,2)= 120.5d0
    this%prms_angles(1,2,53,1)= 80.0d0
    this%prms_angles(1,2,53,2)= 120.6d0
    this%prms_angles(1,5,50,1)= 70.0d0
    this%prms_angles(1,5,50,2)= 123.3d0
    this%prms_angles(1,14,49,1)= 70.0d0
    this%prms_angles(1,14,49,2)= 120.0d0
    this%prms_angles(1,14,5,1)= 70.0d0
    this%prms_angles(1,14,5,2)= 120.0d0
  end subroutine set_amberangles

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

!  subroutine set_amberbends(this,p1,p2,p3)
!    class(amber), intent(inout) :: this
!    integer                     :: i
!    real(4) x1,x2
!    character(2) pa,pb,pc,p1,p2,p3
!    open(12,file='/tmp/amber/amber_angles.prm',status='old')
!    this%prms_bends(1)=0.d0
!    this%prms_bends(2)=0.d0
!    do i=1,281
!       read(12,*,end=1)pa,pb,pc,x1,x2
!       if(pa.eq.p1.and.pc.eq.p3.or.pa.eq.p3.and.pc.eq.p1)then
!          if(pb.eq.p2)then
!             this%prms_bends(1)=dble(x1)
!             this%prms_bends(2)=dble(x2)
!          end if
!       end if
!    end do
!1   close(12)
!  end subroutine set_amberbends

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

  subroutine set_ambervdw(this)
    class(amber), intent(inout) :: this
    allocate(this%prms_vdw(this%get_natp(),2))
    this%prms_vdw(1,1)=2.2200d0
    this%prms_vdw(1,2)=0.3200d0
    this%prms_vdw(2,1)=1.9080d0
    this%prms_vdw(2,2)=0.0860d0
    this%prms_vdw(3,1)=1.9080d0
    this%prms_vdw(3,2)=0.0860d0
    this%prms_vdw(4,1)=1.7131d0
    this%prms_vdw(4,2)=0.4598d0
    this%prms_vdw(10,1)=1.9480d0
    this%prms_vdw(10,2)=0.2650d0
    this%prms_vdw(15,1)=3.3950d0
    this%prms_vdw(15,2)=0.0001d0
    this%prms_vdw(16,1)=1.9080d0
    this%prms_vdw(16,2)=0.1094d0
    this%prms_vdw(22,1)=1.7500d0
    this%prms_vdw(22,2)=0.0610d0
    this%prms_vdw(24,1)=0.6000d0
    this%prms_vdw(24,2)=0.0157d0
    this%prms_vdw(25,1)=1.3870d0
    this%prms_vdw(25,2)=0.0157d0
    this%prms_vdw(26,1)=1.2870d0
    this%prms_vdw(26,2)=0.0157d0
    this%prms_vdw(27,1)=1.1870d0
    this%prms_vdw(27,2)=0.0157d0
    this%prms_vdw(28,1)=1.4090d0
    this%prms_vdw(28,2)=0.0150d0
    this%prms_vdw(29,1)=1.3590d0
    this%prms_vdw(29,2)=0.0150d0
    this%prms_vdw(30,1)=1.4590d0
    this%prms_vdw(30,2)=0.0150d0
    this%prms_vdw(31,1)=1.4870d0
    this%prms_vdw(31,2)=0.0157d0
    this%prms_vdw(32,1)=0.0000d0
    this%prms_vdw(32,2)=0.0000d0
    this%prms_vdw(33,1)=1.1000d0
    this%prms_vdw(33,2)=0.0157d0
    this%prms_vdw(34,1)=0.6000d0
    this%prms_vdw(34,2)=0.0157d0
    this%prms_vdw(35,1)=0.0000d0
    this%prms_vdw(35,2)=0.0000d0
    this%prms_vdw(36,1)=1.4590d0
    this%prms_vdw(36,2)=0.0150d0
    this%prms_vdw(37,1)=2.3500d0
    this%prms_vdw(37,2)=0.4000d0
    this%prms_vdw(38,1)=5.0000d0
    this%prms_vdw(38,2)=0.1000d0
    this%prms_vdw(39,1)=2.4700d0
    this%prms_vdw(39,2)=0.1000d0
    this%prms_vdw(40,1)=1.8680d0
    this%prms_vdw(40,2)=0.0028d0
    this%prms_vdw(41,1)=2.6580d0
    this%prms_vdw(41,2)=0.0003d0
    this%prms_vdw(42,1)=1.1370d0
    this%prms_vdw(42,2)=0.0183d0
    this%prms_vdw(1,1)=0.0000d0
    this%prms_vdw(1,2)=0.0000d0
    this%prms_vdw(43,1)=0.7926d0
    this%prms_vdw(43,2)=0.8947d0
    this%prms_vdw(44,1)=1.8240d0
    this%prms_vdw(44,2)=0.1700d0
    this%prms_vdw(47,1)=1.8240d0
    this%prms_vdw(47,2)=0.1700d0
    this%prms_vdw(48,1)=1.8680d0
    this%prms_vdw(48,2)=0.0028d0
    this%prms_vdw(52,1)=1.8240d0
    this%prms_vdw(52,2)=0.1700d0
    this%prms_vdw(53,1)=1.6612d0
    this%prms_vdw(53,2)=0.2100d0
    this%prms_vdw(54,1)=1.6612d0
    this%prms_vdw(54,2)=0.2100d0
    this%prms_vdw(55,1)=1.7210d0
    this%prms_vdw(55,2)=0.2104d0
    this%prms_vdw(56,1)=1.6837d0
    this%prms_vdw(56,2)=0.1700d0
    this%prms_vdw(57,1)=1.7683d0
    this%prms_vdw(57,2)=0.1520d0
    this%prms_vdw(58,1)=2.1000d0
    this%prms_vdw(58,2)=0.2000d0
    this%prms_vdw(59,1)=2.9560d0
    this%prms_vdw(59,2)=0.0002d0
    this%prms_vdw(60,1)=2.0000d0
    this%prms_vdw(60,2)=0.2500d0
    this%prms_vdw(61,1)=2.0000d0
    this%prms_vdw(61,2)=0.2500d0
    this%prms_vdw(62,1)=1.1000d0
    this%prms_vdw(62,2)=0.0125d0
  end subroutine set_ambervdw

!  subroutine set_ambervdw(this,p1)
!    class(amber), intent(inout) :: this
!    integer                     :: i
!    real(4) x1,x2
!    character(2) pa,p1
!    open(12,file='/tmp/amber/amber_vdw.prm',status='old')
!    this%prms_vdw(1)=0.d0
!    this%prms_vdw(2)=0.d0
!    do i=1,43
!       read(12,*,end=1)pa,x1,x2
!       if(pa.eq.p1)then
!          this%prms_vdw(1)=dble(x1)
!          this%prms_vdw(2)=dble(x2)
!       end if
!    end do
!1   close(12)
!  end subroutine set_ambervdw

end module amber_module
