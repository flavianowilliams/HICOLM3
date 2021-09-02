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
module molecule_module
  !*******************************************************************************************
  !*******************************************************************************************

  use structure_module

  implicit none

  private
  public :: molecule

  type, extends(structure) :: molecule
     integer, allocatable      :: zatmol(:,:)
     real(8), allocatable      :: sf_coul(:)
     real(8), allocatable      :: sf_vdw(:)
     real(8), allocatable      :: qatmol(:,:)
     real(8), allocatable      :: massmol(:,:)
     real(8), allocatable      :: mmolar(:)
     character(6), allocatable :: tpmol(:,:)
   contains
     procedure          :: set_tpmol
     procedure          :: set_zatmol
     procedure          :: set_qatmol
     procedure          :: set_massmol
     procedure          :: set_mmolar
     procedure          :: set_scale_factor
  end type molecule

contains

  subroutine set_tpmol(this,tpmol)
    class(molecule), intent(inout) :: this
    character(6), intent(in)       :: tpmol(:,:)
    integer                        :: i,j
    do i=1,this%get_nmol()
       do j=1,this%nxmol(i)
          this%tpmol(i,j)=tpmol(i,j)
       end do
    end do
  end subroutine set_tpmol

  subroutine set_zatmol(this,zatmol)
    class(molecule), intent(inout) :: this
    integer, intent(in)            :: zatmol(:,:)
    integer                        :: i,j
    do i=1,this%get_nmol()
       do j=1,this%nxmol(i)
          this%zatmol(i,j)=zatmol(i,j)
       end do
    end do
  end subroutine set_zatmol

  subroutine set_qatmol(this,qatmol)
    class(molecule), intent(inout) :: this
    real(8), intent(in)            :: qatmol(:,:)
    integer                        :: i,j
    do i=1,this%get_nmol()
       do j=1,this%nxmol(i)
          this%qatmol(i,j)=qatmol(i,j)
       end do
    end do
  end subroutine set_qatmol

  subroutine set_massmol(this)
    implicit none
    class(molecule), intent(inout) :: this
    integer                        :: i,j
    real(8)                        :: var
    allocate(this%massmol(this%get_nmol(),this%get_natom()))
    do i=1,this%get_nmol()
       do j=1,this%nxmol(i)
          select case(this%zatmol(i,j))
          case(1)
             var=1.007940000d0
          case(2)
             var=4.002602000d0
          case(3)
             var=6.941000000d0
          case(4)
             var=9.012182000d0
          case(5)
             var=10.81100000d0
          case(6)
             var=12.01070000d0
          case(7)
             var=14.00670000d0
          case(8)
             var=15.99940000d0
          case(9)
             var=18.99840300d0
          case(10)
             var=20.17970000d0
          case(11)
             var=22.98976928d0
          case(12)
             var=24.30500000d0
          case(13)
             var=26.98153860d0
          case(14)
             var=28.08550000d0
          case(15)
             var=30.97376200d0
          case(16)
             var=32.06500000d0
          case(17)
             var=35.45300000d0
          case(18)
             var=39.94800000d0
          case(19)
             var=39.09830000d0
          case(20)
             var=40.07800000d0
          case(21)
             var=44.95591200d0
          case(22)
             var=47.86700000d0
          case(23)
             var=50.94150000d0
          case(24)
             var=51.99610000d0
          case(25)
             var=54.93804500d0
          case(26)
             var=55.84500000d0
          case(27)
             var=58.93319500d0
          case(28)
             var=58.69340000d0
          case(29)
             var=58.69340000d0
          case(30)
             var=65.40900000d0
          case(31)
             var=69.72300000d0
          case(32)
             var=72.64000000d0
          case(33)
             var=74.92160000d0
          case(34)
             var=78.96000000d0
          case(35)
             var=79.90400000d0
          case(36)
             var=83.79800000d0
          case(37)
             var=85.46780000d0
          case(38)
             var=87.62000000d0
          case(39)
             var=88.90585000d0
          case(40)
             var=91.22400000d0
          case(41)
             var=92.90638000d0
          case(42)
             var=95.94000000d0
          case(43)
             var=98.00000000d0
          case(44)
             var=101.07000000d0
          case(45)
             var=102.90550000d0
          case(46)
             var=106.42000000d0
          case(47)
             var=107.86820000d0
          case(48)
             var=112.41100000d0
          case(49)
             var=114.81800000d0
          case(50)
             var=118.71000000d0
          case(51)
             var=121.76000000d0
          case(52)
             var=127.60000000d0
          case(53)
             var=126.90447000d0
          case(54)
             var=131.29300000d0
          case(55)
             var=132.90545190d0
          case(56)
             var=137.32700000d0
          case(57)
             var=138.90547000d0
          case(58)
             var=140.11600000d0
          case(59)
             var=140.90765000d0
          case(60)
             var=144.24200000d0
          case(61)
             var=145.00000000d0
          case(62)
             var=150.36000000d0
          case(63)
             var=151.96400000d0
          case(64)
             var=157.25000000d0
          case(65)
             var=158.92535000d0
          case(66)
             var=162.50000000d0
          case(67)
             var=164.93032000d0
          case(68)
             var=167.25900000d0
          case(69)
             var=168.93421000d0
          case(70)
             var=173.04000000d0
          case(71)
             var=174.96700000d0
          case(72)
             var=178.49000000d0
          case(73)
             var=180.94788000d0
          case(74)
             var=183.84000000d0
          case(75)
             var=186.20700000d0
          case(76)
             var=190.23000000d0
          case(77)
             var=192.21700000d0
          case(78)
             var=195.08400000d0
          case(79)
             var=196.96656900d0
          case(80)
             var=200.59000000d0
          case(81)
             var=204.38330000d0
          case(82)
             var=207.20000000d0
          case(83)
             var=208.98040000d0
          case(84)
             var=209.00000000d0
          case(85)
             var=210.00000000d0
          case(86)
             var=222.00000000d0
          case(87)
             var=223.00000000d0
          case(88)
             var=226.00000000d0
          case(89)
             var=227.00000000d0
          case(90)
             var=232.03806000d0
          case(91)
             var=231.03588000d0
          case(92)
             var=238.02891000d0
          case(93)
             var=237.00000000d0
          case(94)
             var=244.00000000d0
          case(95)
             var=243.00000000d0
          case(96)
             var=247.00000000d0
          case(97)
             var=247.00000000d0
          case(98)
             var=251.00000000d0
          case(99)
             var=251.00000000d0
          case(100)
             var=257.00000000d0
          case(101)
             var=258.00000000d0
          case(102)
             var=259.00000000d0
          case(103)
             var=262.00000000d0
          case(104)
             var=261.00000000d0
          case(105)
             var=262.00000000d0
          case(106)
             var=266.00000000d0
          case(107)
             var=264.00000000d0
          case(108)
             var=277.00000000d0
          case(109)
             var=268.00000000d0
          case(110)
             var=281.00000000d0
          case(111)
             var=272.00000000d0
          case(112)
             var=285.00000000d0
          case(113)
             var=284.00000000d0
          case(114)
             var=289.00000000d0
          case(115)
             var=288.00000000d0
          case(116)
             var=292.00000000d0
          case(118)
             var=294.00000000d0
          end select
          this%massmol(i,j)=var
       end do
    end do
  end subroutine set_massmol

  subroutine set_mmolar(this)
    implicit none
    class(molecule), intent(inout) :: this
    integer                        :: i,j
    real(8)                        :: var
    allocate(this%mmolar(this%get_nmol()))
    do i=1,this%get_nmol()
       var=0.d0
       do j=1,this%nxmol(i)
          var=var+this%massmol(i,j)
       end do
       this%mmolar(i)=var
    end do
  end subroutine set_mmolar

  subroutine set_scale_factor(this,sf_coul,sf_vdw)
    implicit none
    class(molecule), intent(inout) :: this
    integer                        :: i
    real(8)                        :: sf_coul
    real(8)                        :: sf_vdw
    allocate(this%sf_coul(this%get_nmol()),this%sf_vdw(this%get_nmol()))
    do i=1,this%get_nmol()
       this%sf_coul(i)=sf_coul
       this%sf_vdw(i)=sf_vdw
    end do
  end subroutine set_scale_factor

end module molecule_module
