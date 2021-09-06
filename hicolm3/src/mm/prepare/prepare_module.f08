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
  ! SYSTEM and TOPOLGY files will be create                                                  *
  !*******************************************************************************************

  use forcefield_module

  implicit none

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
    integer :: i,j
    do i=1,3
       do j=1,3
          constructor%v(i,j)=0.d0
       end do
       constructor%sys_shift(i)=0.d0
    end do
    call constructor%set_zmatrix_tol(0.5d0)
    call constructor%set_coulop('fscs')
    call constructor%set_fscsalpha(0.25d0)
    call constructor%set_bondmax(100)
    call constructor%set_bendmax(100)
    call constructor%set_torsmax(100)
    call constructor%set_itorsmax(1000)
    call constructor%charmm%set_natp(58)
    call constructor%charmm%set_charmmtypes()
    call constructor%charmm%set_charmmbonds()
    call constructor%charmm%set_charmmangles()
    call constructor%charmm%set_charmmdihedrals()
    call constructor%charmm%set_charmmidihedrals()
    call constructor%charmm%set_charmmvdw()
  end function constructor

  subroutine check(this)
    implicit none
    class(prepare), intent(inout) :: this
    integer                       :: i,j
    if(this%get_nmol().le.0)then
       write(6,*)'ERROR: The number of types of molecules does not be zero!'
       write(6,*)'Hint: Check the input in the &SYSTEM section.'
       stop
    else
       do i=1,this%get_nmol()
          if(this%ntmol(i).le.0)then
             write(6,*)'ERROR: The number of molecules does not be zero!'
             write(6,*)'Hint: Check the input in the &SYSTEM section.'
             stop
          end if
          if(this%nxmol(i).le.0)then
             write(6,*)'ERROR: The number of sites at each molecule does not be zero!'
             write(6,*)'Hint: Check the input in the &SYSTEM section.'
             stop
          end if
          do j=1,this%nxmol(i)
             if(this%zatmol(i,j).eq.0)then
                write(6,*)'ERROR: The atomic number does not be zero!'
                write(6,*)'Hint: Check the input in the &FORCE_FIELD section.'
                stop
             end if
             if(this%tpmol(i,j).eq.'NA')then
                write(6,*)'ERROR: You must define the type of each site!'
                write(6,*)'Hint: Put the corresponding value in &FORCE_FIELD section.'
                stop
             end if
          end do
       end do
    end if
    if(this%get_a().le.16.d0.or.this%get_b().le.16.d0.or.this%get_c().le.16.d0)then
       write(6,*)'ERROR: Lattice constant too small! It must be higher than 16 angstroms'
       write(6,*)'Hint: Increase the value of the lattice constant.'
       stop
    end if
  end subroutine check

  subroutine print_sys(this)
    implicit none
    class(prepare), intent(inout) :: this
    integer                       :: i,j,k,nx
    open(10,file='SYSTEM',status='unknown')
    do i=1,3
       write(10,'(3f16.8)')(this%v(i,j),j=1,3)
    end do
    write(10,'(1x,2i5)')this%get_nmol(),this%get_natom()
    nx=1
    do i=1,this%get_nmol()
       write(10,'(1x,a10,2(1x,i5))')this%namemol(i),this%ntmol(i),this%nxmol(i)
       do j=1,this%ntmol(i)
          do k=1,this%nxmol(i)
             write(10,'(3f16.8,i5)')&
                  this%xa(nx),this%ya(nx),this%za(nx),this%freeze(nx)
             nx=nx+1
          end do
       end do
    end do
  end subroutine print_sys

  subroutine print_top(this)
    implicit none
    class(prepare), intent(inout) :: this
    integer                       :: i1,i,j,k
    real(8)                       :: f1,f2
    open(11,file='TOPOLOGY',status='unknown')
    write(11,'(1x,a2)')'MM'
    if(this%get_coulop().eq.'fscs')then
       write(11,'(1x,a4,1x,f5.3)')this%get_coulop(),this%get_fscsalpha()
    elseif(this%get_coulop().eq.'coul')then
       write(11,'(1x,a4)')this%get_coulop()
    end if
    write(11,'(1x,i2,4(1x,i3))')this%get_nmol(),this%get_bondmax(),this%get_bendmax(),&
         this%get_torsmax(),this%get_nspcs()
    do i=1,this%get_nmol()
       write(11,'(1x,a10,2(1x,f8.6))')this%namemol(i),this%sf_coul(i),this%sf_vdw(i)
       write(11,'(15(1x,i2))')(this%zatmol(i,j),j=1,this%nxmol(i))
       write(11,'(15(1x,a8))')(this%tpmol(i,j),j=1,this%nxmol(i))
       write(11,'(15(1x,f8.4))')(this%massmol(i,j),j=1,this%nxmol(i))
       write(11,'(15(1x,f8.4))')(this%qatmol(i,j),j=1,this%nxmol(i))
       write(11,'(1x,a5,1x,i3)')'bonds',this%bondscnt(i)
       do j=1,this%bondscnt(i)
          select case(this%tbonds(i,j))
          case('charmm')
             write(11,'(2(1x,i3),1x,a6,2(1x,f9.4))')this%molbond(i,j,1),this%molbond(i,j,2),&
                  this%tbonds(i,j),(this%parbnd(i,j,k),k=1,2)
          case('harm')
             write(11,'(2(1x,i3),1x,a6,2(1x,f9.4))')this%molbond(i,j,1),this%molbond(i,j,2),&
                  this%tbonds(i,j),(this%parbnd(i,j,k),k=1,2)
          end select
       end do
       write(11,'(1x,a5,1x,i3)')'bends',this%bendscnt(i)
       do j=1,this%bendscnt(i)
          select case(this%tbends(i,j))
          case('charmm')
             write(11,'(3(1x,i3),1x,a6,2(1x,f9.4))')(this%molbend(i,j,k),k=1,3),&
                  this%tbends(i,j),(this%parbend(i,j,k),k=1,2)
          case('harm')
             write(11,'(3(1x,i3),1x,a6,2(1x,f9.4))')(this%molbend(i,j,k),k=1,3),&
                  this%tbends(i,j),(this%parbend(i,j,k),k=1,2)
          end select
       end do
       write(11,'(1x,a9,1x,i3)')'dihedrals',(this%torscnt(i)+this%itorscnt(i))
       do j=1,this%torscnt(i)
          select case(this%ttors(i,j))
          case('charmm')
             f1=this%partors(i,j,1)
             i1=nint(this%partors(i,j,2))
             f2=this%partors(i,j,3)
             write(11,'(4(1x,i3),1x,a7,1x,f8.4,1x,i2,1x,f8.4)')&
                  (this%moltors(i,j,k),k=1,4),this%ttors(i,j),f1,i1,f2
          case('icharmm')
             f1=this%partors(i,j,1)
             f2=this%partors(i,j,2)
             write(11,'(4(1x,i3),1x,a7,2(1x,f8.4))')&
                  (this%moltors(i,j,k),k=1,4),this%ttors(i,j),f1,f2
          case('harm')
             f1=this%partors(i,j,2)
             f2=this%partors(i,j,3)
             write(11,'(4(1x,i3),1x,a7,2(1x,f9.4))')(this%moltors(i,j,k),k=1,4),&
                  this%ttors(i,j),this%partors(i,j,1),this%partors(i,j,2)
          end select
       end do
       do j=1,this%itorscnt(i)
          select case(this%titors(i,j))
          case('icharmm')
             f1=this%paritors(i,j,1)
             f2=this%paritors(i,j,2)
             write(11,'(4(1x,i3),1x,a7,2(1x,f8.4))')&
                  (this%molitors(i,j,k),k=1,4),this%titors(i,j),f1,f2
          end select
       end do
    end do
    write(11,'(1x,a3,1x,i3)')'vdw',this%get_nvdw()
    do i=1,this%get_nvdw()
       select case(this%tvdw(i))
       case('charmm')
          write(11,'(2(1x,a6),1x,a6,2(1x,f9.4))')this%spcvdw(i,1),&
               this%spcvdw(i,2),this%tvdw(i),this%parvdw(i,1),this%parvdw(i,2)
       case('lj')
          write(11,'(2(1x,a6),1x,a6,2(1x,f9.4))')this%spcvdw(i,1),&
               this%spcvdw(i,2),this%tvdw(i),this%parvdw(i,1),this%parvdw(i,2)
       end select
    end do
  end subroutine print_top

  subroutine print_out(this)
    implicit none
    class(prepare), intent(inout) :: this
    integer                       :: i1,i,j,k
    real(8)                       :: f1,f2
    write(6,*)('#',i=1,93)
    write(6,*)('SYSTEM ',i=1,13)
    write(6,*)('#',i=1,93)
    write(6,*)
    write(6,'(a18,i5)')'Total of atoms:',this%get_natom()
    write(6,*)
    write(6,*)'Real space:'
    write(6,*)
    write(6,'(a16,3f15.8)')'Lattice constts:',this%get_a(),this%get_b(),this%get_c()
    write(6,'(a16,3f15.8)')' Lattice angles:',this%get_alpha()*this%get_aconv(),&
         this%get_beta()*this%get_aconv(),this%get_gamma()*this%get_aconv()
    write(6,*)
    write(6,'(a16,3f15.8)')'Lattice vectors:',(this%v(1,i),i=1,3)
    write(6,'(16x,3f15.8)')(this%v(2,i),i=1,3)
    write(6,'(16x,3f15.8)')(this%v(3,i),i=1,3)
    write(6,*)
    write(6,'(a14,f14.4)')'       VOLUME:',this%get_volume()
    write(6,*)
    write(6,'(a14,1x,a9)')'     Symmetry:',this%get_gsym()
    write(6,*)
    write(6,*)('#',j=1,93)
    write(6,*)('FORCE FIELD ',j=1,8)
    write(6,*)('#',j=1,93)
    write(6,*)
    write(6,'(39x,a14)')'INTRAMOLECULAR'
    write(6,'(39x,a14)')'=============='
    write(6,*)
    write(6,'(16x,a9)')'Molecules'
    write(6,'(15x,111a1)')('-',i=1,62)
    write(6,'(16x,a4,6x,a3,5x,a6,5(4x,a5))')&
         'Type','Qty','Sites','bonds','bends','dihdl','idihd'
    write(6,'(15x,111a1)')('-',i=1,62)
    do i=1,this%get_nmol()
       write(6,'(17x,a6,1x,i5,5(4x,i5))')this%namemol(i),this%ntmol(i),this%nxmol(i),&
            this%bondscnt(i),this%bendscnt(i),this%torscnt(i),this%itorscnt(i)
    end do
    write(6,'(15x,111a1)')('-',i=1,62)
    write(6,'(17x,a6,1x,i5,5(4x,i5))')'Total:',sum(this%ntmol),sum(this%nxmol*this%ntmol),&
         sum(this%bondscnt*this%ntmol),sum(this%bendscnt*this%ntmol),&
         sum(this%torscnt*this%ntmol),sum(this%itorscnt*this%ntmol)
    write(6,*)
    do i=1,this%get_nmol()
       write(6,'(42x,a6)')this%namemol(i)
       write(6,'(2x,111a1)')('*',j=1,90)
       write(6,*)
       write(6,'(2x,a24,1x,f8.3,1x,a5)')'Molar mass:',this%mmolar(i),'g/mol'
       write(6,'(2x,a24,2x,f8.4)')'1-4 sf (electrostatic):',this%sf_coul(i)
       write(6,'(2x,a24,3x,f7.4)')'1-4 sf (Van der Waals):',this%sf_vdw(i)
       write(6,*)
       if(this%nxmol(i).le.10)then
          write(6,'(7x,a6,10(1x,a6))')'Sites:',(this%tpmol(i,j),j=1,this%nxmol(i))
          write(6,*)
          write(6,'(5x,a8,10(1x,f6.3))')'Charges:',(this%qatmol(i,j),j=1,this%nxmol(i))
       else
          write(6,'(7x,a6,10(1x,a6))')'Sites:',(this%tpmol(i,j),j=1,10)
          write(6,'(13x,10(1x,a6))')(this%tpmol(i,j),j=11,this%nxmol(i))
          write(6,*)
          write(6,'(5x,a8,10(1x,f6.3))')'Charges:',(this%qatmol(i,j),j=1,10)
          write(6,'(13x,10(1x,f6.3))')(this%qatmol(i,j),j=11,this%nxmol(i))
       end if
       write(6,*)
       write(6,'(2x,a6,1x,i3)')'Bonds:',this%bondscnt(i)
       write(6,'(2x,111a1)')('-',j=1,52)
       write(6,'(2x,3(a4,2x),a4,3x,a10)')' i ','Site','Site','Type','Parameters'
       write(6,'(2x,111a1)')('-',j=1,52)
       do j=1,this%bondscnt(i)
          select case(this%tbonds(i,j))
          case('charmm')
             write(6,'(2x,3(i3,3x),a6,2f9.2)')&
                  j,(this%molbond(i,j,k),k=1,2),this%tbonds(i,j),(this%parbnd(i,j,k),k=1,2)
          case('harm')
             write(6,'(2x,3(i3,3x),a4,2f9.2)')&
                  j,(this%molbond(i,j,k),k=1,2),this%tbonds(i,j),(this%parbnd(i,j,k),k=1,2)
          end select
       end do
       write(6,'(2x,111a1)')('-',j=1,52)
       write(6,*)
       write(6,'(2x,a6,1x,i5)')'Bends:',this%bendscnt(i)
       write(6,'(2x,111a1)')('-',j=1,52)
       write(6,'(2x,4(a4,1x),a4,4x,a10)')' i ','Site','Site','Site','Type','Parameters'
       write(6,'(2x,111a1)')('-',j=1,52)
       do j=1,this%bendscnt(i)
          select case(this%tbends(i,j))
          case('charmm')
             write(6,'(2x,4(i3,2x),a6,1x,2f8.1)')&
                  j,(this%molbend(i,j,k),k=1,3),this%tbends(i,j),(this%parbend(i,j,k),k=1,2)
          case('harm')
             write(6,'(2x,4(i3,2x),a4,1x,2f8.1)')&
                  j,(this%molbend(i,j,k),k=1,3),this%tbends(i,j),(this%parbend(i,j,k),k=1,2)
          end select
       end do
       write(6,'(2x,111a1)')('-',j=1,52)
       write(6,*)
       write(6,'(2x,a17,1x,i5)')'Proper dihedrals:',this%torscnt(i)
       write(6,'(2x,111a1)')('-',j=1,52)
       write(6,'(2x,5(a4,1x),a4,4x,a10)')&
            ' i ','Site','Site','Site','Site','Type','Parameters'
       write(6,'(2x,111a1)')('-',j=1,52)
       do j=1,this%torscnt(i)
          select case(this%ttors(i,j))
          case('charmm')
             f1=this%partors(i,j,1)
             i1=nint(this%partors(i,j,2))
             f2=this%partors(i,j,3)
             write(6,'(2x,5(i3,2x),a6,2x,f8.4,1x,i1,1x,f8.4)')&
                  j,(this%moltors(i,j,k),k=1,4),this%ttors(i,j),f1,i1,f2
          case('icharmm')
             f1=this%partors(i,j,1)
             f2=this%partors(i,j,2)
             write(6,'(2x,5(i3,2x),a7,2x,f8.4,1x,f8.4)')&
                  j,(this%moltors(i,j,k),k=1,4),this%ttors(i,j),f1,f2
          case('harm')
             f1=this%partors(i,j,2)
             f2=this%partors(i,j,3)
             write(6,'(2x,5(i3,2x),1x,a4,1x,2f8.1)')j,(this%moltors(i,j,k),k=1,4),&
                  this%ttors(i,j),this%partors(i,j,1),this%partors(i,j,2)
          end select
       end do
       write(6,'(2x,111a1)')('-',j=1,52)
       write(6,*)
       write(6,'(2x,a19,1x,i5)')'Improper dihedrals:',this%itorscnt(i)
       write(6,'(2x,111a1)')('-',j=1,52)
       write(6,'(2x,5(a4,1x),a4,4x,a10)')&
            ' i ','Site','Site','Site','Site','Type','Parameters'
       write(6,'(2x,111a1)')('-',j=1,52)
       do j=1,this%itorscnt(i)
          select case(this%titors(i,j))
          case('icharmm')
             f1=this%paritors(i,j,1)
             f2=this%paritors(i,j,2)
             write(6,'(2x,5(i3,2x),a7,2x,f8.4,1x,f8.4)')&
                  j,(this%molitors(i,j,k),k=1,4),this%titors(i,j),f1,f2
          end select
       end do
       write(6,'(2x,111a1)')('-',j=1,52)
       write(6,*)
       write(6,'(2x,111a1)')('*',j=1,90)
       write(6,*)
    end do
    write(6,'(39x,a14)')'INTERMOLECULAR'
    write(6,'(39x,a14)')'=============='
    write(6,*)
    write(6,'(2x,a14,1x,f7.4)')' Total charge:',this%get_qtotal()
    write(6,*)
    if(this%get_nspcs().le.10)then
       write(6,'(2x,a18,i3,2x,a2,10(1x,a6))')&
            'Total of species:',this%get_nspcs(),'->',(this%spcs(i),i=1,this%get_nspcs())
       write(6,*)
    else
       write(6,'(2x,a18,i3,2x,a2,10(1x,a2))')&
            'Total of species:',this%get_nspcs(),'->',(this%spcs(i),i=1,10)
       write(6,'(27x,10(1x,a2))')(this%spcs(i),i=11,this%get_nspcs())
       write(*,*)
    end if
    select case(this%get_coulop())
    case('coul')
       write(6,'(2x,a53)')'Electrostatic interaction: Direct Coulomb Sum'
       write(6,*)
    case('fscs')
       write(6,'(20x,a41)')'Electrostatic: Force-Shifted Coulomb Sum'
       write(6,'(20x,52a1)')('-',i=1,52)
       write(6,'(20x,a6,1x,f5.3)')'alpha:',this%get_fscsalpha()
       write(6,'(20x,52a1)')('-',i=1,52)
       write(6,*)
    end select
    write(6,'(20x,a15,i5)')'Van der Waals:',this%get_nvdw()
    write(6,'(20x,111a1)')('-',i=1,52)
    write(6,'(22x,a4,4x,a4,5x,a4,6x,a10)')'Site','Site','Type','Parameters'
    write(6,'(20x,111a1)')('-',i=1,52)
    do i=1,this%get_nvdw()
       select case(this%tvdw(i))
       case('charmm')
          write(6,'(21x,a6,2x,a6,4x,a6,3(1x,f9.4))')this%spcvdw(i,1),this%spcvdw(i,2),&
               this%tvdw(i),this%parvdw(i,1),this%parvdw(i,2)
       case('lj')
          write(6,'(21x,a6,2x,a6,4x,a2,3(1x,f9.4))')this%spcvdw(i,1),this%spcvdw(i,2),&
               this%tvdw(i),this%parvdw(i,1),this%parvdw(i,2)
       end select
    end do
    write(6,'(20x,111a1)')('-',i=1,52)
    write(6,*)
  end subroutine print_out

end module prepare_module
