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
module forcefield_module
  !*******************************************************************************************
  !*******************************************************************************************

  use zmatrix_module
  use amber_module

  implicit none

  integer i,j,k,l

  private
  public :: forcefield

  type, extends(zmatrix) :: forcefield
     type(amber)               :: amber
     real(8), allocatable      :: parbnd(:,:,:)
     character(5), allocatable :: tbonds(:,:)
   contains
     procedure, private :: forcefield_init
     procedure          :: set_potentials
     procedure          :: set_extra_potentials
  end type forcefield

  interface forcefield
     module procedure constructor
  end interface forcefield

contains

  type(forcefield) function constructor()
    implicit none
    call constructor%forcefield_init()
  end function constructor

  subroutine forcefield_init(this)
    implicit none
    class(forcefield), intent(inout) :: this
    allocate(this%parbnd(this%nmol,this%bondmax,2))
    allocate(this%tbonds(this%nmol,this%bondmax))
    do i=1,this%nmol
       do j=1,this%nxmol(i)
          this%tbonds(i,j)='amber'
       end do
    end do
  end subroutine forcefield_init

  subroutine set_potentials(this)
    implicit none
    class(forcefield), intent(inout) :: this
    integer                          :: i1,i2
    call this%forcefield_init()
    call this%amber%set_amber()
    do i=1,this%nmol
       do j=1,this%bondscnt(i)
          i1=this%molbond(i,j,1)
          i2=this%molbond(i,j,2)
          do k=1,115
             if(this%tpmol(i,i1).eq.this%amber%tpam(k,1).and.&
                  this%tpmol(i,i2).eq.this%amber%tpam(k,2))then
                do l=1,2
                   this%parbnd(i,j,l)=this%amber%prms(k,l)
                end do
             end if
          end do
       end do
    end do
  end subroutine set_potentials

  subroutine set_extra_potentials(this)
    implicit none
    class(forcefield), intent(inout) :: this
    integer                          :: i1,i2,i3
    character(12)                    :: key
    character(4)                     :: key2
    character(8)                     :: key3
    character(5)                     :: key4
    character(10)                    :: cvar
    logical                          :: check
    i3=0
1   read(5,*,end=2)key
    if(key.ne.'&FORCE_FIELD')goto 1
    do while (key2.ne.'&END')
       read(5,*)key2
       if(key2.ne.'&END')then
          backspace(5)
          do i=1,this%nmol
             check=.FALSE.
             read(5,*)key3
             if(key3.eq.'molecule')then
                backspace(5)
                read(5,*)key3,cvar
                do j=1,this%nmol
                   if(cvar.eq.this%namemol(j))then
                      check=.TRUE.
                      i3=j
                   end if
                end do
                if(check.eqv..FALSE.)goto 3
                if(i3.ne.0)then
                   read(5,*)
                   read(5,*)
                   read(5,*)
                   read(5,*)key4
                   if(key4.eq.'bonds')then
                      backspace(5)
                      read(5,*)key4,i1
                      if(i1.eq.0)then
                         this%bondscnt(i)=i1
                      elseif(i1.gt.0)then
                         do j=1,i1
                            read(5,*)i2
                            backspace(5)
                            read(5,*)i2,this%molbond(i3,i2,1),this%molbond(i3,i2,2),&
                                 this%tbonds(i3,i2)
                            backspace(5)
                            select case(this%tbonds(i3,i2))
                            case('amber')
                               read(5,*)i2,this%molbond(i3,i2,1),this%molbond(i3,i2,2),&
                                    this%tbonds(i3,i2),(this%parbnd(i3,i2,k),k=1,2)
                            case('harm')
                               read(5,*)i2,this%molbond(i3,i2,1),this%molbond(i3,i2,2),&
                                    this%tbonds(i3,i2),(this%parbnd(i3,i2,k),k=1,2)
                            end select
                         end do
                      end if
                   end if
                end if
             end if
          end do
       end if
    end do
2   rewind(5)
    return
3   write(6,*)'ERROR: There is a molecule that does not belong to the physical system!'
    write(6,*)'Hint: Check the input in the &FORCE_FIELD section.'
    stop
  end subroutine set_extra_potentials

end module forcefield_module
