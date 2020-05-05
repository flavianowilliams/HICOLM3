!
! This file is part of the HICOLM distribution (https://github.com/flavianowilliams/HICOLM).
!
! Copyright (c) 2019 Flaviano Williams Fernandes.
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, version 3.
!
! This program is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
! General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program. If not, see <http://www.gnu.org/licenses/>.
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
  use tersoff_module

contains

  subroutine ff_prepare
    !***************************************************************************************
    ! Preparando campo de forca                                                            *
    !***************************************************************************************

    implicit none

    integer i,j
    real(8) chqtot

    write(6,'(96a1)')('#',i=1,93)
    write(6,*)

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

    if(nvdw.ne.0)then
       call vdw_convert
       call vdw_counts
    end if

    !-unidades de Coulomb

    if(ncoul.ne.0)then
       call coulomb_convert
       call coulomb_counts
    end if

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

    !-unidades de Tersoff

    if(ntrsff.ne.0)call tersoff_convert

    if(nmolec.ne.0)then
       write(6,*)'Molecules'
       write(6,'(111a1)')('-',i=1,52)
       write(6,10)'Type','Qty','Sites','bonds','bends','dihdl'
       write(6,'(111a1)')('-',i=1,52)
       do i=1,nmolec
          write(6,20)namemol(i),ntmolec(i),nxmolec(i),bondsmlc(i),bendsmlc(i),torsmlc(i)
       end do
       write(6,'(111a1)')('-',i=1,52)
       write(6,20)'Total:',moltot,natom,nbondstp,nbendstp,ntorsstp
       write(6,*)
    end if

    chqtot=0
    do i=1,spctot
       chqtot=chqtot+parcoul(atnp(i),1)
    end do

    if(spctot.lt.10)then
       write(6,30)'Total of species:',spctot,'->',(atnp(i),i=1,spctot)
       write(*,*)
       write(6,80)' Partial charges:',spctot,'->',(parcoul(atnp(i),1),i=1,spctot)
       write(*,*)
    else
       write(6,30)'Total of species:',spctot,'->',(atnp(i),i=1,10)
       write(6,31)(atnp(i),i=11,spctot)
       write(*,*)
       write(6,80)' Partial charges:',spctot,'->',(parcoul(atnp(i),1),i=1,10)
       write(6,81)(parcoul(atnp(i),1),i=11,spctot)
       write(6,*)
       write(6,82)' Total charge:',chqtot
       write(6,*)
    end if
    write(6,*)'Intermolecular:'
    write(6,*)

    if(ncoul.ne.0)then
       select case(coulop)
       case('coul')
          write(*,*)'Electrostatic interaction: Force-Shifted Coulomb Sum'
          write(*,*)
       end select
    end if

    if(nvdw.ne.0)then
       write(6,60)'Van der Waals:',nvdw
       write(6,*)
       do i=1,spctot
          do j=i,spctot
             select case(vdw(i,j))
             case(1)
                write(6,50)&
                     i,j,'mors',parvdw(i,j,1)*econv,parvdw(i,j,2)*kconv,parvdw(i,j,3)*rconv
             case(2)
                write(6,50)i,j,'lj',parvdw(i,j,1)*econv,parvdw(i,j,2)*rconv
             end select
          end do
       end do
       write(6,*)
       write(6,70)'       Maximum coulomb interaction:',ncoulstp
       write(6,70)' Maximum Van der Waals interaction:',nvdwstp
       write(6,*)
    end if

    !-limpando memoria

    deallocate(bondsmlc,bendsmlc,torsmlc)

    return

10  format(a4,7x,a3,4x,a6,4x,a5,4x,a5,4x,a5)
20  format(a6,2x,i5,4(4x,i5))
30  format(a18,i3,2x,a2,10i4)
31  format(25x,10i4)
50  format(5x,2i3,a5,5x,3f7.4)
60  format(a14,i5)
70  format(a35,i10)
80  format(a18,i3,2x,a2,10f8.4)
81  format(25x,10f8.4)
82  format(a18,10f8.4)

  end subroutine ff_prepare

  subroutine ff_modules(encoul,&
       enbond,enbend,entors,envdw,entrsff,virvdw,virbond,virbend,virtors,vircoul,virtrsff)
    !***************************************************************************************
    ! Modulos do campo de forca (intramolecular):                                          *
    ! - Estiramento;                                                                       *
    ! - Deformacao;                                                                        *
    ! - Torsao;                                                                            *
    ! - Tersoff;                                                                           *
    ! - Van der Waals;                                                                     *
    ! - Eletrostatico.                                                                     *
    !***************************************************************************************

    implicit none

    integer i
    real(8) encoul,enbond,enbend,entors,envdw,entrsff
    real(8) virvdw,virbond,virbend,virtors,vircoul,virtrsff

    !-valores iniciais

    !-energia

    encoul=0.d0   !coulombiano
    enbond=0.d0   !estiramento
    enbend=0.d0   !deformacao
    entors=0.d0   !torção
    envdw=0.d0    !Van der waals
    entrsff=0.d0  !Tersoff

    !-virial

    virvdw=0.d0   !Van der Waals
    virbond=0.d0  !estiramento
    virbend=0.d0  !deformacao
    virtors=0.d0  !torção
    vircoul=0.d0  !coulombiano
    virtrsff=0.d0 !Tersoff

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

    !-calculo da contribuicao de deformacao angular

    if(ntors.ne.0)call tors_calc(entors,virtors)

    !-calculo da contribuicao de tersoff

    if(ntrsff.ne.0)call tersoff_calc(entrsff,virtrsff)

    !-calculo das contribuicoes intermoleculares (Van der Waals e coulombiano)

    if(nvdw.ne.0.or.ncoul.ne.0)call ff_modules_inter(envdw,encoul,virvdw,vircoul)

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
          if(nvdw.ne.0)then
             if(parvdw(atp(ni),atp(nj),1).ne.0.d0) &
                  call vdw_calc(envdw,virvdw,ni,nj,xvz,yvz,zvz)
          end if
          if(ncoul.ne.0)then
             if(parcoul(atp(ni),1)*parcoul(atp(nj),1).ne.0.d0) &
                  call coulomb_calc(encoul,vircoul,ni,nj,xvz,yvz,zvz)
          end if
       end do
    end do

    return

  end subroutine ff_modules_inter

end module force_field
