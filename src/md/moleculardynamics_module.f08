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
module moleculardynamics_module
  !*******************************************************************************************
  !*******************************************************************************************

  use interaction_module

  implicit none

  integer i,j,k

  private
  public :: moleculardynamics

  type, extends(interaction) :: moleculardynamics
     real(8), private :: time
   contains
     procedure :: print_geometry
     procedure :: read_geometry
     procedure :: print_out
     procedure :: print
     procedure :: set_time
     procedure :: get_time
  end type moleculardynamics

  interface moleculardynamics
     module procedure constructor
  end interface moleculardynamics

contains

  type(moleculardynamics) function constructor()
    implicit none
    call constructor%set_nstep(1)
    call constructor%set_nrelax(1)
    call constructor%set_nframes(1)
    call constructor%set_timestep(0.001d0)
    call constructor%set_press(1.d0)
    call constructor%set_temp(298.d0)
    call constructor%set_rcutoff(8.0d0)
    call constructor%set_drcutoff(0.1d0)
    call constructor%set_ensble('nve')
  end function constructor

  subroutine set_time(this,time)
    implicit none
    class(moleculardynamics), intent(inout) :: this
    real(8), intent(in)                     :: time
    this%time=time
  end subroutine set_time

  double precision function get_time(this)
    implicit none
    class(moleculardynamics), intent(inout) :: this
    get_time=this%time
  end function get_time

  subroutine read_geometry(this)
    implicit none
    class(moleculardynamics), intent(inout) :: this
    open(1,file='HICOLM.XSF',status='old')
    do i=1,13
       read(1,*)
    end do
    do i=1,3
       read(1,'(3(3x,f14.8))')(this%v(i,j),j=1,3)
    end do
    read(1,*)
    read(1,*)
    select case(this%get_restart())
    case(1)
       do i=1,this%get_natom()
          read(1,'(5x,3f14.8,2(2x,3f14.8))')this%xa(i),this%ya(i),this%za(i)
       end do
    case(2)
       do i=1,this%get_natom()
          read(1,'(5x,3f14.8,2(2x,3f14.8))')this%xa(i),this%ya(i),this%za(i),&
               this%fax(i),this%fay(i),this%faz(i)
       end do
    case(3)
       do i=1,this%get_natom()
          read(1,'(5x,3f14.8,2(2x,3f14.8))')this%xa(i),this%ya(i),this%za(i),&
               this%fax(i),this%fay(i),this%faz(i),this%vax(i),this%vay(i),this%vaz(i)
       end do
    end select
    close(1)
  end subroutine read_geometry

 subroutine print_geometry(this,mdstp)
   implicit none
   class(moleculardynamics), intent(inout) :: this
   integer, intent(in)                     :: mdstp
   open(1,file='HICOLM.XSF',status='unknown')
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
      write(1,'(i5,3f14.8,2x,3f14.8,2x,3f14.8)')this%zat(i),&
           this%xa(i)*this%get_rconv(),&
           this%ya(i)*this%get_rconv(),&
           this%za(i)*this%get_rconv(),&
           this%fax(i)*this%get_econv()/this%get_rconv(),&
           this%fay(i)*this%get_econv()/this%get_rconv(),&
           this%faz(i)*this%get_econv()/this%get_rconv(),&
           this%vax(i)*this%get_rconv()/this%get_tconv(),&
           this%vay(i)*this%get_rconv()/this%get_tconv(),&
           this%vaz(i)*this%get_rconv()/this%get_tconv()
   end do
   close(1)
   if(mod(mdstp,2).ne.0)return
   open(2,file='.HICOLM.XSF',status='unknown')
   write(2,*)'BEGIN_INFO'
   write(2,*)'  #'
   write(2,*)'  # This is a XCRYSDEN-Structure-File'
   write(2,*)'  #'
   write(2,*)'  # aimed for Visualization of geometry'
   write(2,*)'  #'
   write(2,*)'  # Launch as: xcrysden --xsf .HICOLM.XSF'
   write(2,*)'  #'
   write(2,*)'  #'
   write(2,*)'END_INFO'
   write(2,*)'# estrutura final'
   write(2,'(a7)')'CRYSTAL'
   write(2,'(a7)')'PRIMVEC'
   do i=1,3
      write(2,'(3(3x,f14.8))')(this%v(i,j)*this%get_rconv(),j=1,3)
   end do
   write(2,'(a9)')'PRIMCOORD'
   write(2,'(2i5)')this%get_natom(),1
   do i=1,this%get_natom()
      write(1,'(i5,3f14.8,2x,3f14.8,2x,3f14.8)')this%zat(i),&
           this%xa(i)*this%get_rconv(),&
           this%ya(i)*this%get_rconv(),&
           this%za(i)*this%get_rconv(),&
           this%fax(i)*this%get_econv()/this%get_rconv(),&
           this%fay(i)*this%get_econv()/this%get_rconv(),&
           this%faz(i)*this%get_econv()/this%get_rconv(),&
           this%vax(i)*this%get_rconv()/this%get_tconv(),&
           this%vay(i)*this%get_rconv()/this%get_tconv(),&
           this%vaz(i)*this%get_rconv()/this%get_tconv()
   end do
   close(2)
    return
  end subroutine print_geometry

  subroutine print_out(this)
    implicit none
    class(moleculardynamics), intent(inout) :: this
    integer                                 :: i1,i2
    real(8)                                 :: f1,f2
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
    write(6,'(20x,a9)')'Molecules'
    write(6,'(19x,111a1)')('-',i=1,54)
    write(6,'(20x,a4,7x,a3,4x,a6,4(4x,a5))')'Type','Qty','Sites','bonds','bends','dihdl'
    write(6,'(19x,111a1)')('-',i=1,54)
    do i=1,this%get_nmol()
       write(6,'(20x,a6,2x,i5,4(4x,i5))')this%namemol(i),this%ntmol(i),this%nxmol(i),&
            this%bondscnt(i),this%bendscnt(i),this%torscnt(i)
    end do
    write(6,'(19x,111a1)')('-',i=1,54)
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
       write(6,'(2x,a6,1x,i3)')'Bonds:',this%bondscnt(i)
       write(6,'(2x,111a1)')('-',j=1,52)
       write(6,'(2x,3(a4,2x),a4,3x,a10)')' i ','Site','Site','Type','Parameters'
       write(6,'(2x,111a1)')('-',j=1,52)
       do j=1,this%bondscnt(i)
          select case(this%tbonds(i,j))
          case('amber')
             write(6,'(2x,3(i3,3x),a5,2f9.2)')j,(this%molbond(i,j,k),k=1,2),this%tbonds(i,j),&
                  this%parbnd(i,j,1)*this%get_econv()/this%get_rconv()**2,&
                  this%parbnd(i,j,2)*this%get_rconv()
          case('harm')
             write(6,'(2x,3(i3,3x),a5,2f9.2)')j,(this%molbond(i,j,k),k=1,2),this%tbonds(i,j),&
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
          case('amber')
             write(6,'(2x,4(i3,2x),a5,1x,2f8.1)')j,(this%molbend(i,j,k),k=1,3),&
                  this%tbends(i,j),&
                  this%parbend(i,j,1)*this%get_econv()/this%get_aconv()**2,&
                  this%parbend(i,j,2)*this%get_aconv()
          case('harm')
             write(6,'(2x,4(i3,2x),a5,1x,2f8.1)')j,(this%molbend(i,j,k),k=1,3),&
                  this%tbends(i,j),&
                  this%parbend(i,j,1)*this%get_econv()/this%get_aconv()**2,&
                  this%parbend(i,j,2)*this%get_aconv()
          end select
       end do
       write(6,'(2x,111a1)')('-',j=1,52)
       write(6,*)
       write(6,'(2x,a10,1x,i5)')'Dihedrals:',this%torscnt(i)
       write(6,'(2x,111a1)')('-',j=1,90)
       write(6,'(2x,5(a4,1x),1x,a4,4x,a10)')&
            ' i ','Site','Site','Site','Site','Type','Parameters'
       write(6,'(2x,111a1)')('-',j=1,90)
       do j=1,this%torscnt(i)
          select case(this%ttors(i,j))
          case('amber')
             i1=nint(this%partors(i,j,1))
             f1=this%partors(i,j,2)*this%get_econv()
             f2=this%partors(i,j,3)*this%get_aconv()
             i2=nint(this%partors(i,j,4))
             write(6,'(2x,5(i3,2x),a5,2x,i2,f8.2,f8.1,1x,i2)')j,&
                  (this%moltors(i,j,k),k=1,4),this%ttors(i,j),i1,f1,f2,i2
          case('harm')
             f1=this%partors(i,j,1)*this%get_econv()/this%get_aconv()**2
             f2=this%partors(i,j,2)*this%get_aconv()
             write(6,'(2x,5(i3,2x),1x,a4,1x,2f8.1)')&
                  j,(this%moltors(i,j,k),k=1,4),this%ttors(i,j),f1,f2
          end select
       end do
       write(6,'(2x,111a1)')('-',j=1,90)
       write(6,*)
       write(6,'(2x,111a1)')('*',j=1,90)
       write(6,*)
    end do
    write(6,'(39x,a14)')'INTERMOLECULAR'
    write(6,'(39x,a14)')'=============='
    write(6,*)
    select case(this%get_coulop())
    case('coul')
       write(6,'(2x,a53)')'Electrostatic interaction: Direct Coulomb Sum'
       write(6,*)
    case('fscs')
       write(6,'(2x,a53)')'Electrostatic interaction: Force-Shifted Coulomb Sum'
       write(6,*)
    case('escl')
       write(6,'(2x,a53)')'Electrostatic interaction: Coulomb scaled potencial'
       write(6,*)
    end select
    write(6,'(2x,a14,1x,f7.4)')' Total charge:',this%sys%get_qtotal()
    write(6,*)
    if(this%get_nspcs().le.10)then
       write(6,'(2x,a18,i3,2x,a2,10(1x,a2))')&
            'Total of species:',this%get_nspcs(),'->',(this%spcs(i),i=1,this%get_nspcs())
       write(6,*)
    else
       write(6,'(2x,a18,i3,2x,a2,10(1x,a2))')&
            'Total of species:',this%get_nspcs(),'->',(this%spcs(i),i=1,10)
       write(6,'(27x,10(1x,a2))')(this%spcs(i),i=11,this%get_nspcs())
       write(*,*)
    end if
    write(6,'(20x,a15,i5)')'Van der Waals:',this%get_nvdw()
    write(6,'(20x,111a1)')('-',i=1,52)
    write(6,'(20x,a4,2x,a4,3x,a4,6x,a10)')'Site','Site','Type','Parameters'
    write(6,'(20x,111a1)')('-',i=1,52)
    do i=1,this%get_nvdw()
       select case(this%tvdw(i))
       case('amber')
          write(6,'(21x,a2,4x,a2,4x,a5,3(1x,f9.4))')this%spcvdw(i,1),this%spcvdw(i,2),&
               this%tvdw(i),this%parvdw(i,1)*this%get_econv(),this%parvdw(i,2)*this%get_rconv()
       case('lj')
          write(6,'(21x,a2,4x,a2,4x,a5,3(1x,f9.4))')this%spcvdw(i,1),this%spcvdw(i,2),&
               this%tvdw(i),this%parvdw(i,1)*this%get_econv(),this%parvdw(i,2)*this%get_rconv()
       end select
    end do
    write(6,'(20x,111a1)')('-',i=1,52)
    write(6,*)
  end subroutine print_out

  subroutine print(this)
    implicit none
    class(moleculardynamics), intent(inout) :: this
    write(6,*)('#',i=1,93)
    write(6,*)('MD RUNNING ',i=1,8)
    write(6,*)('#',i=1,93)
    write(6,*)
    write(6,'(31x,a30)')'Molecular dynamics information'
    write(6,'(28x,36a1)')('-',j=1,36)
    write(6,'(28x,a16,1x,a3,1x,a9)')'Ensemble:',this%get_ensble(),this%get_ensble_mt()
    if(this%get_ensble().eq.'nvt')then
       write(6,'(28x,a16,5x,f8.2)')'Thermostat:',this%get_tstat()*this%get_tconv()
    elseif(this%get_ensble().eq.'npt')then
       write(6,'(28x,a16,5x,f8.2)')'Thermostat:',this%get_tstat()*this%get_tconv()
       if(this%get_ensble_mt().eq.'berendsen')then
          write(6,'(28x,a16,5x,f8.2)')'Barostat:',this%get_pstat()*this%get_tconv()
       elseif(this%get_ensble_mt().eq.'hoover')then
          write(6,'(28x,a16,5x,f8.2)')'Barostat:',this%get_pstat()*this%get_tconv()
       end if
    end if
    write(6,'(28x,a16,5x,f9.3,1x,a1)')'Temperature:',this%get_temp()*this%get_teconv(),'K'
    write(6,'(28x,a16,5x,f9.3,1x,a3)')'Pressure:',this%get_press()*this%get_pconv(),'atm'
    write(6,'(28x,a16,6x,i10)')'Number of steps:',this%get_nstep()
    write(6,'(28x,a16,6x,i10)')'Number of relax:',this%get_nrelax()
    write(6,'(28x,a16,9x,i1)')'Restart:',this%get_restart()
    write(6,'(28x,a16,8x,es10.3,1x,a2)')'Timestep:',this%get_timestep()*this%get_tconv(),'ps'
    write(6,'(28x,a16,8x,2f6.3,1x,a1)')'rcutoff:',this%get_rcutoff()*this%get_rconv(),&
         this%get_drcutoff()*this%get_rconv(),'A'
    write(6,'(28x,36a1)')('-',j=1,36)
    write(6,*)
    write(6,'(2x,a27,f12.6)')'  Correction of VdW energy:',this%get_encorr()*this%get_econv()
    write(6,'(2x,a27,f12.6)')'Correction of VdW pressure:',this%get_vircorr()*this%get_econv()
    write(6,*)
  end subroutine print

end module moleculardynamics_module
