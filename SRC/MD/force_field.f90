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
module force_field
  !******************************************************************************************
  ! Contribuicao de diedros para o campo de força:                                          *
  ! - Energia potencial;                                                                    *
  ! - Forças atômicas;                                                                      *
  ! - Virial;                                                                               *
  ! - Stress.                                                                               *
  !                                                                                         *
  ! Flaviano Williams Fernandes, 05 de dezembro de 2018                                     *
  !******************************************************************************************

  use input
  use alloc_arrays
  use bonds_module
  use bends_module
  use dihedral_module
  use vdw_module
  use coulomb_module
  use neighbour_list

contains

  subroutine ff_prepare
    !***************************************************************************************
    ! Preparando campo de forca                                                            *
    !***************************************************************************************

    implicit none

    integer i,j,k,l,i1,i2
    real(8) chqtot,f1,f2

    !-valores iniciais

    nbondstp=0
    nbendstp=0
    ntorsstp=0

    do i=1,nmolec
       bondsmlc(i)=0
       bendsmlc(i)=0
       torsmlc(i)=0
    end do

    !-unidades de Van der Waals

    call vdw_prepare

    !-unidades de Coulomb

    call coulomb_prepare

    !-alocando arrays (lista de vizinhos de Verlet)

    call neighbour_prepare

    !-unidades de estiramento

    if(nbonds.ne.0)then
       call bonds_alloc
       call bonds_convert
       call bonds_counts
    end if

    !-unidades de deformação angular

    if(nbends.ne.0)then
       call bends_alloc
       call bends_convert
       call bends_counts
    end if

    !-unidades de torção

    if(ntors.ne.0)then
       call tors_alloc
       call tors_convert
       call tors_counts
    end if

    write(6,*)('#',j=1,93)
    write(6,*)('FORCE FIELD ',j=1,8)
    write(6,*)('#',j=1,93)
    WRITE(6,*)

    !-intramolecular description

    write(6,'(39x,a14)')'INTRAMOLECULAR'
    write(6,'(39x,a14)')'=============='
    write(6,*)

    if(nmolec.ne.0)then
       write(6,'(20x,a9)')'Molecules'
       write(6,'(20x,111a1)')('-',i=1,52)
       write(6,'(20x,a4,7x,a3,4x,a6,3(4x,a5))')'Type','Qty','Sites','bonds','bends','dihdl'
       write(6,'(20x,111a1)')('-',i=1,52)
       do i=1,nmolec
          write(6,'(20x,a6,2x,i5,4(4x,i5))')&
               namemol(i),ntmolec(i),nxmolec(i),bondsmlc(i),bendsmlc(i),torsmlc(i)
       end do
       write(6,'(20x,111a1)')('-',i=1,52)
       write(6,'(20x,a6,2x,i5,4(4x,i5))')'Total:',moltot,natom,nbondstp,nbendstp,ntorsstp
       write(6,*)
    end if

    write(6,*)
    write(6,*)

    k=0
    do i=1,nmolec
       write(6,'(42x,a6,1x,a7)')namemol(i),ff_model(i)
       write(6,'(20x,111a1)')('*',j=1,52)
       if(nxmolec(i).le.10)then
          write(6,'(20x,a6,10(1x,a2))')'Sites:',(atsp(j+k),j=1,nxmolec(i))
       else
          write(6,'(20x,a6,10(1x,a2))')'Sites:',(atsp(j+k),j=1,nxmolec(i))
          write(6,'(26x,10(1x,a2))')(atsp(j+k),j=1,nxmolec(i))
       end if
       write(6,*)
       write(6,'(20x,a6,1x,i5)')'Bonds:',bondscnt(i)
       write(6,'(20x,111a1)')('-',j=1,52)
       write(6,'(20x,a4,2x,a4,2x,a4,3x,a10)')'Site','Site','Type','Parameters'
       write(6,'(20x,111a1)')('-',j=1,52)
       do j=1,bondscnt(i)
          select case(bonds(i,j))
          case(1)
             write(6,'(20x,2(i3,3x),a4,1x,3f9.4)')(molbond(i,j,l),l=1,2),&
                  'mors',parbnd(i,j,1)*econv,parbnd(i,j,2)*kconv,parbnd(i,j,3)*rconv
          case(2)
             write(6,'(20x,2(i3,3x),a4,1x,2f9.4)')(molbond(i,j,l),l=1,2),&
                  'harm',parbnd(i,j,1)*(econv/rconv**2.d0),parbnd(i,j,2)*rconv
          case(3)
             write(6,'(20x,2(i3,3x),a5,2f9.4)')(molbond(i,j,l),l=1,2),&
                  'amber',parbnd(i,j,1)*(econv/rconv**2.d0),parbnd(i,j,2)*rconv
          end select
       end do
       write(6,'(20x,111a1)')('-',j=1,52)
       write(6,*)
       write(6,'(20x,a6,1x,i5)')'Bends:',bendscnt(i)
       write(6,'(20x,111a1)')('-',j=1,52)
       write(6,'(20x,3(a4,2x),a4,4x,a10)')'Site','Site','Site','Type','Parameters'
       write(6,'(20x,111a1)')('-',j=1,52)
       do j=1,bendscnt(i)
          select case(bends(i,j))
          case(1)
             write(6,'(20x,3(i3,3x),a4,1x,2f8.3)')(molbend(i,j,l),l=1,3),&
                  'harm',parbend(i,j,1)*econv,parbend(i,j,2)*aconv
          case(2)
             write(6,'(20x,3(i3,3x),a5,2f8.3)')&
                  (molbend(i,j,l),l=1,3),'amber',parbend(i,j,1)*econv,parbend(i,j,2)*aconv
          end select
       end do
       write(6,'(20x,111a1)')('-',j=1,52)
       write(6,*)
       write(6,'(20x,a10,1x,i5)')'Dihedrals:',torscnt(i)
       write(6,'(20x,111a1)')('-',j=1,52)
       write(6,'(20x,4(a4,2x),a4,4x,a10)')'Site','Site','Site','Site','Type','Parameters'
       write(6,'(20x,111a1)')('-',j=1,52)
       do j=1,torscnt(i)
          select case(tors(i,j))
          case(1)
             write(6,'(20x,4(i3,3x),a4,1x,2f8.3)')(moltors(i,j,l),l=1,4),&
                  'harm',partors(i,j,1)*econv,partors(i,j,2)*aconv
          case(4)
             i1=nint(partors(i,j,1))
             f1=partors(i,j,2)*econv
             f2=partors(i,j,3)*aconv
             i2=nint(partors(i,j,4))
             write(6,'(20x,4(i3,3x),a5,2x,i2,f8.5,f8.3,i2)')&
                  (moltors(i,j,l),l=1,4),'amber',i1,f1,f2,i2
          end select
       end do
       write(6,'(20x,111a1)')('-',j=1,52)
       write(6,*)
       write(6,'(20x,111a1)')('*',j=1,52)
       write(6,*)
       write(6,*)
       k=k+nxmolec(i)
    end do
    write(6,*)
    write(6,*)

    !-intermolecular description

    write(6,'(39x,a14)')'INTERMOLECULAR'
    write(6,'(39x,a14)')'=============='
    write(6,*)

    chqtot=0
    do i=1,natom
       chqtot=chqtot+parcoul(atp(i),1)
    end do

    select case(coulop)
    case(1)
       write(6,90)'Electrostatic interaction: Direct Coulomb Sum'
       write(6,*)
    case(2)
       write(6,90)'Electrostatic interaction: Force-Shifted Coulomb Sum'
       write(6,*)
    end select

    if(spctot.le.10)then
       write(6,'(2x,a18,i3,2x,a2,10(1x,a2))')&
            'Total of species:',spctot,'->',(atsp(i),i=1,spctot)
       write(6,*)
       write(6,'(2x,a18,i3,2x,a2,10f7.3)')&
            'Partial charges:',spctot,'->',(parcoul(atnp(i),1),i=1,spctot)
       write(*,*)
       write(6,82)' Total charge:',chqtot
       write(6,*)
       write(6,70)'       Maximum coulomb interaction:',ncoulstp
       write(6,*)
    else
       write(6,'(2x,a18,i3,2x,a2,10(1x,a2))')&
            'Total of species:',spctot,'->',(atsp(i),i=1,10)
       write(6,'(27x,10(1x,a2))')(atsp(i),i=11,spctot)
       write(*,*)
       write(6,'(2x,a18,i3,2x,a2,10f7.3)')&
            'Partial charges:',spctot,'->',(parcoul(atnp(i),1),i=1,10)
       write(6,'(27x,10f7.3)')(parcoul(atnp(i),1),i=11,spctot)
       write(6,*)
       write(6,82)' Total charge:',chqtot
       write(6,*)
       write(6,70)'       Maximum coulomb interaction:',ncoulstp
       write(6,*)
    end if

    write(6,'(20x,a15,i5)')'Van der Waals:',nvdw
    write(6,'(20x,111a1)')('-',i=1,52)
    write(6,'(20x,a4,2x,a4,2x,a4,3x,a10)')'Site','Site','Type','Parameters'
    write(6,'(20x,111a1)')('-',i=1,52)
    do i=1,spctot
       do j=i,spctot
          select case(vdw(i,j))
          case(1)
             write(6,'(20x,i3,3x,i3,2x,a5,2x,3f7.4)')&
                  i,j,'mors',parvdw(i,j,1)*econv,parvdw(i,j,2)*kconv,parvdw(i,j,3)*rconv
          case(2)
             write(6,'(20x,i3,3x,i3,2x,a5,2x,3f7.4)')&
                  i,j,'lj',parvdw(i,j,1)*econv,parvdw(i,j,2)*rconv
          case(3)
!             if(parvdw(i,j,1).ne.0.d0.and.parvdw(i,j,2).ne.0.d0)then
                write(6,'(20x,a3,3x,a3,2x,a5,2x,3f7.4)')&
                     atsp(i),atsp(j),'amber',parvdw(i,j,1)*econv,parvdw(i,j,2)*rconv
!             end if
          end select
       end do
    end do
    write(6,'(20x,111a1)')('-',i=1,52)
    write(6,*)
    write(6,70)' Maximum Van der Waals interaction:',nvdwstp
    write(6,*)

    !-limpando memoria

    deallocate(bondsmlc,bendsmlc,torsmlc)

    return

70  format(2x,a35,i10)
82  format(2x,a18,10f8.4)
90  format(2x,a53)

  end subroutine ff_prepare

  subroutine ff_modules&
       (encoul,enbond,enbend,entors,envdw,virvdw,virbond,virbend,virtors,vircoul)
    !***************************************************************************************
    ! Modulos do campo de forca (intramolecular):                                          *
    ! - Estiramento;                                                                       *
    ! - Deformacao;                                                                        *
    ! - Torsao;                                                                            *
    ! - Van der Waals;                                                                     *
    ! - Eletrostatico.                                                                     *
    !***************************************************************************************

    implicit none

    integer i
    real(8) encoul,enbond,enbend,entors,envdw
    real(8) virvdw,virbond,virbend,virtors,vircoul

    !-valores iniciais

    !-energia

    encoul=0.d0   !coulombiano
    enbond=0.d0   !estiramento
    enbend=0.d0   !deformacao
    entors=0.d0   !torção
    envdw=0.d0    !Van der waals

    !-virial

    virvdw=0.d0   !Van der Waals
    virbond=0.d0  !estiramento
    virbend=0.d0  !deformacao
    virtors=0.d0  !torção
    vircoul=0.d0  !coulombiano

    !-forcas atomicas

    do i=1,natom
       fax(i)=0.d0
       fay(i)=0.d0
       faz(i)=0.d0
    end do

    !-stress

    do i=1,6
       str(i)=0.d0
    end do

    !-calculo da contribuicao de estiramento

    if(nbonds.ne.0)call bonds_calc(enbond,virbond)

    !-calculo da contribuicao de deformacao angular

    if(nbends.ne.0)call bends_calc(enbend,virbend)

    !-calculo da contribuicao de diedros

    if(ntors.ne.0)call tors_calc(entors,virtors)

    !-calculo das contribuicoes intermoleculares (Van der Waals e coulombiano)

    call ff_modules_inter(envdw,encoul,virvdw,vircoul)

    return

  end subroutine ff_modules

  subroutine ff_modules_inter(envdw,encoul,virvdw,vircoul)
    !***************************************************************************************
    ! Modulos do campo de forca:                                                           *
    ! - Van der Waals;                                                                     *
    ! - Eletrostatico.                                                                     *
    !***************************************************************************************

    implicit none

    integer i,j,ni,nj
    real(8) xvz,yvz,zvz
    real(8) envdw,encoul,virvdw,vircoul

    do i=1,natom-1
       do j=1,nlist(i)
          ni=i
          nj=ilist(i,j)
          call mic(ni,nj,xvz,yvz,zvz)
          if(parvdw(atp(ni),atp(nj),1).gt.1.d-8)call vdw_calc(envdw,virvdw,ni,nj,xvz,yvz,zvz)
          if(abs(parcoul(atp(ni),1)).gt.1.e-8)then
             if(abs(parcoul(atp(nj),1)).gt.1.e-8)then
                call coulomb_calc(encoul,vircoul,ni,nj,xvz,yvz,zvz)
             end if
          end if
       end do
    end do

    return

  end subroutine ff_modules_inter

end module force_field
