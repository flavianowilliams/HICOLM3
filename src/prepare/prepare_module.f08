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
module prepare_module
  !*******************************************************************************************
  ! class to prepare the physical environment for simulation                                 *
  ! HICOLM.sys and HICOLM.top will be create                                                 *
  !*******************************************************************************************

  use forcefield_module

  implicit none

  integer i,j

  private
  public :: prepare

  type, extends(forcefield) :: prepare
   contains
     procedure :: check
     procedure :: print_sys
     procedure :: print_top
     procedure :: print_out
  end type prepare

  interface prepare
     module procedure constructor
  end interface prepare

contains

  type(prepare) function constructor()
    implicit none
    do i=1,3
       do j=1,3
          constructor%v(i,j)=0.d0
       end do
    end do
    constructor%zmatrix_tol=0.5d0
  end function constructor

  subroutine check(this)
    implicit none
    class(prepare), intent(inout) :: this
    if(this%nmol.le.0)then
       write(6,*)'ERROR: The number of types of molecules does not be zero!'
       write(6,*)'Hint: Check the input in the HICOLM.in.'
       stop
    else
       do i=1,this%nmol
          if(this%ntmol(i).le.0)then
             write(6,*)'ERROR: The number of molecules does not be zero!'
             write(6,*)'Hint: Check the input in the HICOLM.in.'
             stop
          end if
          if(this%nxmol(i).le.0)then
             write(6,*)'ERROR: The number of sites at each molecule does not be zero!'
             write(6,*)'Hint: Check the input in the HICOLM.in.'
             stop
          end if
       end do
    end if
    if(this%a.le.1.e-4.or.this%b.le.1.e-4.or.this%c.le.1.e-4)then
       write(6,*)'ERROR: Lattice constant too small!'
       write(6,*)'Hint: Increase the value in the HICOLM.in.'
       stop
    end if
  end subroutine check

  subroutine print_sys(this)
    implicit none
    class(prepare), intent(inout) :: this
    open(10,file='HICOLM.sys',status='unknown')
    do i=1,3
       write(10,'(3f16.8)')(this%v(i,j),j=1,3)
    end do
    do i=1,this%get_natom()
       write(10,'(3f16.8)')this%xa(i),this%ya(i),this%za(i)
    end do
  end subroutine print_sys

  subroutine print_top(this)
    implicit none
    class(prepare), intent(inout) :: this
    open(11,file='HICOLM.top',status='unknown')
    write(11,'(1x,a2)')'MM'
    write(11,'(1x,i2)')this%nmol
    do i=1,this%nmol
       write(11,'(1x,a10,2(1x,i5),2(1x,f8.6))')&
            this%namemol(i),this%ntmol(i),this%nxmol(i),this%sf_coul(i),this%sf_vdw(i)
       write(11,'(15(1x,i2))')(this%zatmol(i,j),j=1,this%nxmol(i))
       write(11,'(15(1x,a2))')(this%tpmol(i,j),j=1,this%nxmol(i))
       write(11,'(15(1x,f8.4))')(this%massmol(i,j),j=1,this%nxmol(i))
       write(11,'(15(1x,f8.4))')(this%qatmol(i,j),j=1,this%nxmol(i))
       write(11,'(1x,a5,1x,i3)')'bonds',this%bondscnt(i)
       do j=1,this%bondscnt(i)
          write(11,'(2(1x,i3),2(1x,f9.4))')this%molbond(i,j,1),this%molbond(i,j,2),&
               this%parbnd(i,j,1),this%parbnd(i,j,2)
       end do
    end do
  end subroutine print_top

  subroutine print_out(this)
    implicit none
    class(prepare), intent(inout) :: this
    write(6,*)('#',i=1,93)
    write(6,*)('SYSTEM ',i=1,13)
    write(6,*)('#',i=1,93)
    write(6,*)
    write(6,'(a18,i5)')'Total of atoms:',this%get_natom()
    write(6,*)
    write(6,*)'Real space:'
    write(6,*)
    write(6,'(a16,3f15.8)')'Lattice constts:',this%a,this%b,this%c
    write(6,'(a16,3f15.8)')' Lattice angles:',this%alpha*this%aconv,this%beta*this%aconv,&
         this%gamma*this%aconv
    write(6,*)
    write(6,'(a16,3f15.8)')'Lattice vectors:',(this%v(1,i),i=1,3)
    write(6,'(16x,3f15.8)')(this%v(2,i),i=1,3)
    write(6,'(16x,3f15.8)')(this%v(3,i),i=1,3)
    write(6,*)
    write(6,'(a14,f14.4)')'       VOLUME:',this%volume
    write(6,*)
    write(6,'(a14,1x,a9)')'     Symmetry:',this%gsym
    write(6,*)
    write(6,*)('#',j=1,93)
    write(6,*)('FORCE FIELD ',j=1,8)
    write(6,*)('#',j=1,93)
    write(6,*)
    write(6,'(39x,a14)')'INTRAMOLECULAR'
    write(6,'(39x,a14)')'=============='
    write(6,*)
    write(6,'(20x,a9)')'Molecules'
    write(6,'(19x,111a1)')('-',i=1,54)
    write(6,'(20x,a4,7x,a3,4x,a6,4(4x,a5))')'Type','Qty','Sites','bonds','bends','dihdl'
    write(6,'(19x,111a1)')('-',i=1,54)
    write(6,'(20x,a6,2x,i5,4(4x,i5))')
    write(6,'(19x,111a1)')('-',i=1,54)
    write(6,'(20x,a6,2x,i5,5(4x,i5))')&
         'Total:'
    write(6,*)
    do i=1,this%nmol
       write(6,'(42x,a6)')this%namemol(i)
       write(6,'(2x,111a1)')('*',j=1,90)
       write(6,*)
       write(6,'(2x,a24,1x,f8.3,1x,a5)')'Molar mass:',this%mmolar(i),'g/mol'
       write(6,'(2x,a24,2x,f8.4)')'1-4 sf (electrostatic):',this%sf_coul(i)
       write(6,'(2x,a24,3x,f7.4)')'1-4 sf (Van der Waals):',this%sf_vdw(i)
       write(6,*)
       if(this%nxmol(i).le.10)then
          write(6,'(7x,a6,10(1x,a2))')'Sites:',(this%tpmol(i,j),j=1,this%nxmol(i))
          write(6,*)
          write(6,'(5x,a8,10(1x,f6.3))')'Charges:',(this%qatmol(i,j),j=1,this%nxmol(i))
       else
          write(6,'(7x,a6,10(1x,a2))')'Sites:',(this%tpmol(i,j),j=1,10)
          write(6,'(13x,10(1x,a2))')(this%tpmol(i,j),j=11,this%nxmol(i))
          write(6,*)
          write(6,'(5x,a8,10(1x,f6.3))')'Charges:',(this%qatmol(i,j),j=1,10)
          write(6,'(13x,10(1x,f6.3))')(this%qatmol(i,j),j=11,this%nxmol(i))
       end if
       write(6,*)
       write(6,'(2x,111a1)')('*',j=1,90)
       write(6,*)
       write(6,*)
    end do
  end subroutine print_out

end module prepare_module
