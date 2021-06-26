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
module optimize_module
  !*******************************************************************************************
  !*******************************************************************************************

  use gradientdescent_module

  implicit none

  private
  public :: optimize

  type, extends(gradientdescent) :: optimize
   contains
     procedure :: check
     procedure :: print
     procedure :: convert_units
     procedure :: set_canonicalvariables
     procedure :: print_geometry
  end type optimize

  interface optimize
     module procedure constructor
  end interface optimize

contains

  type(optimize) function constructor()
    implicit none
    call constructor%set_nstep(1000)
    call constructor%set_tolerance(0.001d0)
    call constructor%set_restart('undefine')
  end function constructor

  subroutine check(this)
    implicit none
    class(optimize), intent(inout) :: this
    if(this%get_restart().ne.'undefine'.and.this%get_restart().ne.'position')then
       write(6,*)
       write(6,*)'ERROR: The notifyed restart directive is not an option!'
       write(6,*)'Hint: Check the input in the &OPTIMIZATION section.'
       stop
    end if
  end subroutine check

  subroutine convert_units(this)
    implicit none
    class(optimize), intent(inout) :: this
    integer                        :: i,j
    do i=1,3
       do j=1,3
          this%v(i,j)=this%v(i,j)/this%get_rconv()
       end do
    end do
    call this%set_a(this%get_a()/this%get_rconv())
    call this%set_b(this%get_b()/this%get_rconv())
    call this%set_c(this%get_c()/this%get_rconv())
    do i=1,this%get_natom()
       this%xa(i)=this%xa(i)/this%get_rconv()
       this%ya(i)=this%ya(i)/this%get_rconv()
       this%za(i)=this%za(i)/this%get_rconv()
    end do
    do i=1,this%get_nmol()
       do j=1,this%nxmol(i)
          this%qatmol(i,j)=this%qatmol(i,j)/this%get_elconv()
       end do
    end do
    do i=1,this%get_nmol()
       do j=1,this%bondscnt(i)
          select case(this%tbonds(i,j))
          case('charmm')
             this%parbnd(i,j,1)=this%parbnd(i,j,1)/(this%get_econv()/this%get_rconv()**2)
             this%parbnd(i,j,2)=this%parbnd(i,j,2)/this%get_rconv()
          case('harm')
             this%parbnd(i,j,1)=this%parbnd(i,j,1)/(this%get_econv()/this%get_rconv()**2)
             this%parbnd(i,j,2)=this%parbnd(i,j,2)/this%get_rconv()
          end select
       end do
       do j=1,this%bendscnt(i)
          select case(this%tbends(i,j))
          case('charmm')
             this%parbend(i,j,1)=this%parbend(i,j,1)/this%get_econv()
             this%parbend(i,j,2)=this%parbend(i,j,2)/this%get_aconv()
          case('harm')
             this%parbend(i,j,1)=this%parbend(i,j,1)/this%get_econv()
             this%parbend(i,j,2)=this%parbend(i,j,2)/this%get_aconv()
          end select
       end do
       do j=1,this%torscnt(i)
          select case(this%ttors(i,j))
          case('charmm')
             this%partors(i,j,1)=this%partors(i,j,1)/this%get_econv()
             this%partors(i,j,3)=this%partors(i,j,3)/this%get_aconv()
          case('icharmm')
             this%partors(i,j,1)=this%partors(i,j,1)/this%get_econv()
             this%partors(i,j,2)=this%partors(i,j,2)/this%get_aconv()
          case('harm')
             this%partors(i,j,1)=this%partors(i,j,1)/this%get_econv()
             this%partors(i,j,2)=this%partors(i,j,2)/this%get_aconv()
          end select
       end do
    end do
    do i=1,this%get_nvdw()
       select case(this%tvdw(i))
       case('charmm')
          this%parvdw(i,1)=this%parvdw(i,1)/this%get_econv()
          this%parvdw(i,2)=this%parvdw(i,2)/this%get_rconv()
       case('lj')
          this%parvdw(i,1)=this%parvdw(i,1)/this%get_econv()
          this%parvdw(i,2)=this%parvdw(i,2)/this%get_rconv()
       end select
    end do
    call this%set_fscsalpha(this%get_fscsalpha()/this%get_kconv())
    call this%set_tolerance(this%get_tolerance()/(this%get_econv()/this%get_rconv()))
  end subroutine convert_units

  subroutine set_canonicalvariables(this)
    implicit none
    class(optimize), intent(inout) :: this
    integer                        :: i,j
    select case(this%get_restart())
    case('position')
       open(1,file='hicolm.xsf',status='old')
       do i=1,13
          read(1,*)
       end do
       do i=1,3
          read(1,'(3(3x,f14.8))')(this%v(i,j),j=1,3)
       end do
       read(1,*)
       read(1,*)
       do i=1,this%get_natom()
          read(1,'(5x,3f14.8)')this%xa(i),this%ya(i),this%za(i)
       end do
    end select
  end subroutine set_canonicalvariables

  subroutine print_geometry(this)
    implicit none
    class(optimize), intent(inout) :: this
    integer                        :: i,j
    open(1,file='hicolm.xsf',status='unknown')
    write(1,*)'BEGIN_INFO'
    write(1,*)'  #'
    write(1,*)'  # This is a XCRYSDEN-Structure-File'
    write(1,*)'  #'
    write(1,*)'  # aimed for Visualization of geometry'
    write(1,*)'  #'
    write(1,*)'  # Launch as: xcrysden --xsf HICOLM.XSF'
    write(1,*)'  #'
    write(1,*)'  #'
    write(1,*)'END_INFO'
    write(1,*)'# estrutura final'
    write(1,'(a7)')'CRYSTAL'
    write(1,'(a7)')'PRIMVEC'
    do i=1,3
       write(1,'(3(3x,f14.8))')(this%v(i,j)*this%get_rconv(),j=1,3)
    end do
    write(1,'(a9)')'PRIMCOORD'
    write(1,'(2i5)')this%get_natom(),1
    do i=1,this%get_natom()
       write(1,'(i5,3f14.8)')this%zat(i),&
            this%xa(i)*this%get_rconv(),&
            this%ya(i)*this%get_rconv(),&
            this%za(i)*this%get_rconv()
    end do
    close(1)
    return
  end subroutine print_geometry

  subroutine print(this)
    implicit none
    class(optimize), intent(inout) :: this
    integer                        :: i,j,k,i1
    real(8)                        :: f1,f2
    write(6,*)('#',i=1,93)
    write(6,*)(' MINIMIZATION ',i=1,6)
    write(6,*)('#',i=1,93)
    write(6,*)
    write(6,'(28x,a38)')'Minimization of the energy information'
    write(6,'(28x,38a1)')('-',j=1,38)
    write(6,'(28x,a16,1x,i10)')'Number of steps:',this%get_nstep()
    write(6,'(28x,a16,1x,es10.3,1x,a10)')'      Tolerance:',&
         this%get_tolerance()*(this%get_econv()/this%get_rconv()),'kcal/mol*A'
    write(6,'(28x,38a1)')('-',j=1,38)
    write(6,*)
    write(6,*)('#',i=1,93)
    write(6,*)('SYSTEM ',i=1,13)
    write(6,*)('#',i=1,93)
    write(6,*)
    write(6,'(a18,i5)')'Total of atoms:',this%get_natom()
    write(6,*)
    write(6,*)'Real space:'
    write(6,*)
    write(6,'(a16,3f15.8)')'Lattice constts:',this%get_a()*this%get_rconv(),&
         this%get_b()*this%get_rconv(),this%get_c()*this%get_rconv()
    write(6,'(a16,3f15.8)')' Lattice angles:',this%get_alpha()*this%get_aconv(),&
         this%get_beta()*this%get_aconv(),this%get_gamma()*this%get_aconv()
    write(6,*)
    write(6,'(a16,3f15.8)')'Lattice vectors:',(this%v(1,i)*this%get_rconv(),i=1,3)
    write(6,'(16x,3f15.8)')(this%v(2,i)*this%get_rconv(),i=1,3)
    write(6,'(16x,3f15.8)')(this%v(3,i)*this%get_rconv(),i=1,3)
    write(6,*)
    write(6,'(a14,f14.4)')'       VOLUME:',this%get_volume()*this%get_rconv()**3
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
    write(6,'(21x,a9)')'Molecules'
    write(6,'(20x,111a1)')('-',i=1,53)
    write(6,'(21x,a4,6x,a3,5x,a6,4(4x,a5))')'Type','Qty','Sites','bonds','bends','dihdl'
    write(6,'(20x,111a1)')('-',i=1,53)
    do i=1,this%get_nmol()
       write(6,'(22x,a6,1x,i5,4(4x,i5))')this%namemol(i),this%ntmol(i),this%nxmol(i),&
            this%bondscnt(i),this%bendscnt(i),this%torscnt(i)
    end do
    write(6,'(20x,111a1)')('-',i=1,53)
    write(6,'(22x,a6,1x,i5,5(4x,i5))')'Total:',sum(this%ntmol),sum(this%nxmol*this%ntmol),&
         sum(this%bondscnt*this%ntmol),sum(this%bendscnt*this%ntmol),&
         sum(this%torscnt*this%ntmol)
    write(6,*)
    do i=1,this%get_nmol()
       write(6,'(42x,a6)')this%namemol(i)
       write(6,'(2x,111a1)')('*',j=1,90)
       write(6,*)
       write(6,'(2x,a24,1x,f8.3,1x,a5)')'Molar mass:',this%mmolar(i)*this%get_mconv(),'g/mol'
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
             write(6,'(2x,3(i3,3x),a6,2f9.2)')j,(this%molbond(i,j,k),k=1,2),this%tbonds(i,j),&
                  this%parbnd(i,j,1)*this%get_econv()/this%get_rconv()**2,&
                  this%parbnd(i,j,2)*this%get_rconv()
          case('harm')
             write(6,'(2x,3(i3,3x),a4,2f9.2)')j,(this%molbond(i,j,k),k=1,2),this%tbonds(i,j),&
                  this%parbnd(i,j,1)*this%get_econv()/this%get_rconv()**2,&
                  this%parbnd(i,j,2)*this%get_rconv()
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
             write(6,'(2x,4(i3,2x),a6,1x,2f8.1)')j,(this%molbend(i,j,k),k=1,3),&
                  this%tbends(i,j),&
                  this%parbend(i,j,1)*this%get_econv(),this%parbend(i,j,2)*this%get_aconv()
          case('harm')
             write(6,'(2x,4(i3,2x),a4,1x,2f8.1)')j,(this%molbend(i,j,k),k=1,3),&
                  this%tbends(i,j),&
                  this%parbend(i,j,1)*this%get_econv(),this%parbend(i,j,2)*this%get_aconv()
          end select
       end do
       write(6,'(2x,111a1)')('-',j=1,52)
       write(6,*)
       write(6,'(2x,a10,1x,i5)')'Dihedrals:',this%torscnt(i)
       write(6,'(2x,90a1)')('-',j=1,90)
       write(6,'(2x,5(a4,1x),a4,4x,a10)')&
            ' i ','Site','Site','Site','Site','Type','Parameters'
       write(6,'(2x,90a1)')('-',j=1,90)
       do j=1,this%torscnt(i)
          select case(this%ttors(i,j))
          case('charmm')
             f1=this%partors(i,j,1)*this%get_econv()
             i1=nint(this%partors(i,j,2))
             f2=this%partors(i,j,3)*this%get_aconv()
             write(6,'(2x,5(i3,2x),a6,2x,f8.4,1x,i1,1x,f8.4)')j,&
                  (this%moltors(i,j,k),k=1,4),this%ttors(i,j),f1,i1,f2
          case('icharmm')
             f1=this%partors(i,j,1)*this%get_econv()
             f2=this%partors(i,j,2)*this%get_aconv()
             write(6,'(2x,5(i3,2x),a7,2x,f8.4,1x,f8.4)')j,&
                  (this%moltors(i,j,k),k=1,4),this%ttors(i,j),f1,f2
          case('harm')
             f1=this%partors(i,j,1)*this%get_econv()
             f2=this%partors(i,j,2)*this%get_aconv()
             write(6,'(2x,5(i3,2x),a4,2x,2f8.1)')&
                  j,(this%moltors(i,j,k),k=1,4),this%ttors(i,j),f1,f2
          end select
       end do
       write(6,*)
       write(6,'(2x,90a1)')('*',j=1,90)
       write(6,*)
    end do
    write(6,'(39x,a14)')'INTERMOLECULAR'
    write(6,'(39x,a14)')'=============='
    write(6,*)
    write(6,'(2x,a14,1x,f7.4)')' Total charge:',this%sys%get_qtotal()
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
       write(6,'(20x,a6,1x,f5.3)')'alpha:',this%get_fscsalpha()*this%get_kconv()
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
               this%tvdw(i),this%parvdw(i,1)*this%get_econv(),this%parvdw(i,2)*this%get_rconv()
       case('lj')
          write(6,'(21x,a6,2x,a6,4x,a2,3(1x,f9.4))')this%spcvdw(i,1),this%spcvdw(i,2),&
               this%tvdw(i),this%parvdw(i,1)*this%get_econv(),this%parvdw(i,2)*this%get_rconv()
       end select
    end do
    write(6,'(20x,111a1)')('-',i=1,52)
    write(6,*)
  end subroutine print

end module optimize_module
