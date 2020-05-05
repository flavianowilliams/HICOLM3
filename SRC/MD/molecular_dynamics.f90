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
module molecular_dynamics
  !****************************************************************************************
  ! Dinâmica Molecular                                                                    *
  !                                                                                       *
  ! Flaviano Williams Fernandes, 05 de dezembro de 2018                                   *
  !****************************************************************************************

  use input
  use estrutura
  use force_field
  use ensemble_nve
  use ensemble_nvt
  use ensemble_npt
  use neighbour_list

  real(8) mtot

  save mtot

contains

  subroutine md(t3)
    !***************************************************************************************
    ! Obtencao das variaveis canônicas a cada intervalo de tempo                           *
    !***************************************************************************************

    implicit none

    integer i,ihist
    real(8) t0,t3,time,temp,press
    real(8) enpot,ekinet

    !-contagem do tempo absoluto

    call cpu_time(t0)

    !-preparando Dinâmica molecular

    call md_prepare

    !-imprimindo parametros de entrada

    call md_print

    !-correcao da energia VdW de curto alcance

    if(nvdw.ne.0)call vdw_corr

    !-calculo da temperatura e energia cinetica no instante inicial

    ekinet=0.d0
    do i=1,natom
       ekinet=ekinet+mass(i)*(vax(i)**2+vay(i)**2+vaz(i)**2)
    end do

    ekinet=0.5d0*ekinet

    write(6,'(96a1)')('#',i=1,96)
    write(6,*)
    write(6,'(a27,f12.4)')'  Correction of VdW energy:',envdw_corr*econv
    write(6,'(a27,f12.4)')'Correction of VdW pressure:',virvdw_corr*econv
    write(6,*)

    !-ciclo MD

    write(6,*)

    write(6,'(111a1)')('-',i=1,88)
    write(6,10)'##','Step','Time','VOLUME','TEMPERA','PRESSURE','E(TOTL)'
    write(6,'(111a1)')('-',i=1,88)

    !-relaxação do sistema

    time=dtime
    do i=1,nrelax
       call mdloop(i,temp,press,ekinet,enpot)
       write(6,20)'MD',i,time*tconv,volume*rconv**3,&
            temp*teconv,press*pconv,(ekinet+enpot+envdw_corr)*econv
       time=time+dtime
    end do

    !-sistema em equilibrio termodinâmico

    ihist=1
    do i=nrelax+1,ntrialmax
       call mdloop(i,temp,press,ekinet,enpot)
       write(6,20)'MD',i,time*tconv,volume*rconv**3,&
            temp*teconv,press*pconv,(ekinet+enpot+envdw_corr)*econv
       if(mod(i,nhist).eq.0)call history(ihist)
       time=time+dtime
    end do

    write(6,'(111a1)')('-',i=1,88)
    write(6,*)

    !-contagem do tempo absoluto

    call cpu_time(t3)

    t3=t3-t0

    return

10  format(1x,a2,6x,a4,6x,a5,9x,a6,6x,a8,7x,a8,6x,a7)
20  format(1x,a2,2x,i8,2x,e12.4,2x,e12.4,3(2x,e12.4))

  end subroutine md

  subroutine md_prepare()
    !***************************************************************************************
    ! Preparacao da dinâmica molecular:                                                    *
    ! - Preparando campo de forca;                                                         *
    ! - Preparando lista de vizinhos de Verlet;                                            *
    !***************************************************************************************

    implicit none

    integer i,nwr
    real(8) virvdw,virbond,virbend,virtors,vircoul,virtrsff
    real(8) encoul,enbond,enbend,entors,envdw,entrsff

    !-abrindo ficheiro de dados

    open(2,file='HICOLM.md',status='unknown')

    nwr=6

    write(2,'(i12,e12.4,3i7)')int((ntrialmax-nrelax)/nhist),dtime*tconv,nwr,natom,spctot

    !-preparando Campo de Força

    call ff_prepare

    !-preparando lista de vizinhos de Verlet

    if(nvdw.ne.0.or.ncoul.ne.0)call verlet_list_inter

    if(ntrsff.ne.0)call verlet_list_all

    !-forças iniciais

    if(reuse.eq.0)call ff_modules(encoul,enbond,enbend,entors,envdw,entrsff,virvdw,virbond, &
         virbend,virtors,vircoul,virtrsff)

    !-energia potencial

!    enpot=encoul+enbond+enbend+entors+envdw+entrsff

    !-massa total

    mtot=0.d0
    do i=1,natom
       mtot=mtot+mass(i)
    end do

   return

  end subroutine md_prepare

  subroutine mdloop(mdstp,temp,press,ekinet,enpot)
    !***************************************************************************************
    ! Obtencao das variaveis canonicas;                                                    *
    ! Calculo da energia total;                                                            *
    ! Atualizacao das lista de vizinhos de Verlet                                          *
    ! Impressao das variaveis em ficheiros de saida                                        *
    !***************************************************************************************

    implicit none

    integer mdstp,i,ix
    real(8) temp,press

    real(8) virvdw,virbond,virbend,virtors,vircoul,virtrsff
    real(8) ekinet,encoul,enbond,enbend,entors,envdw,entrsff,enpot

    !-integrando as variáveis canônicas posição e velocidade

    select case(ensble)
    case('nve')
       call nve_vv(temp,press,ekinet,encoul,enbond,enbend,entors,envdw,entrsff,virvdw, &
            virbond,virbend,virtors,vircoul,virtrsff)
    case('nvt')
       call nvt_ber_vv(temp,press,ekinet,encoul,enbond,enbend,entors,envdw,entrsff,virvdw, &
            virbond,virbend,virtors,vircoul,virtrsff)
    case('npt')
       call npt_ber_vv(temp,press,ekinet,encoul,enbond,enbend,entors,envdw,entrsff,virvdw, &
            virbond,virbend,virtors,vircoul,virtrsff)
    end select

    !-calculo da energia potencial

    enpot=encoul+enbond+enbend+entors+envdw+entrsff

    !-imprimindo estrutura

    call geometria

    !-escrevendo variaveis canonicas em arquivo de dados

    if(mdstp.gt.nrelax.and.mod(mdstp-nrelax,nhist).eq.0)then
       ix=(mdstp-nrelax)/nhist
       write(2,10)ix
       write(2,20)(str(i)*econv/rconv**3,i=1,6)
       write(2,20)(ix*dtime)*tconv,volume*rconv**3,temp*teconv,&
            press*pconv,(ekinet+envdw_corr)*econv,(ekinet+enpot+envdw_corr)*econv,&
            (mtot/volume)*mconv/(6.02d-4*rconv**3)
       write(2,*)
       do i=1,3
          write(2,40)v(i,1)*rconv,v(i,2)*rconv,v(i,3)*rconv
       end do
       write(2,*)
       do i=1,natom
          write(2,30)idna(i),atp(i),mass(i)*mconv,parcoul(atp(i),1),xa(i)*rconv, &
               ya(i)*rconv,za(i)*rconv
          write(2,40)vax(i)*rconv/tconv,vay(i)*rconv/tconv,vaz(i)*rconv/tconv
          write(2,40)fax(i)*econv/rconv,fay(i)*econv/rconv,faz(i)*econv/rconv
       end do
    end if

    !-atualizando a lista de vizinhos de Verlet

    if(mod(mdstp,verlchk).eq.0)then
       if(nvdw.ne.0.or.ncoul.ne.0)call verlet_list_inter
    end if

    if(ntrsff.ne.0)call verlet_list_all

    return

10  format(1x,i5)
20  format(6f16.5)
30  format(1x,2i5,5f12.4)
40  format(35x,3f12.4)

  end subroutine mdloop

  subroutine md_print
    !*************************************************************************
    ! Impressão dos dados de entrada                                         *
    !*************************************************************************
    implicit none

    integer j

    write(6,*)('#',j=1,70)
    write(6,*)'$$$',(' ENTRADA',j=1,8),'$$$'
    write(6,*)('#',j=1,70)
    write(6,*)
    write(6,'(5x,a30)')'Molecular mechanics parameters'
    write(6,'(5x,36a1)')('-',j=1,36)
    write(6,'(5x,a16,f7.2)')'Temperature (K):',text*teconv
    write(6,'(7x,a14,e12.4)')'Pressure (atm):',preext*pconv
    write(6,'(12x,a9,a9)')'Ensemble:',ensble
    if(ensble.eq.'nvt')then
       write(6,'(10x,a11,f5.2)')'Thermostat:',tstat*tconv
    elseif(ensble.eq.'npt')then
       write(6,'(10x,a11,f5.2)')'Thermostat:',tstat*tconv
       write(6,'(10x,a11,2f5.2)')'Barostat:',pstat*tconv,bfactor/pconv
    end if
    write(6,'(5x,a16,i10)')'ntrialmax:',ntrialmax
    write(6,'(5x,a16,i10)')'nrelax:',nrelax
    write(6,'(5x,a16,i5)')'reuse:',reuse
    write(6,'(5x,a16,f5.3)')'Timestep (ps):',dtime*tconv
    write(6,'(5x,a16,2f7.3)')'rcutoff:',rcutoff*rconv,drcutoff*rconv
    write(6,'(5x,36a1)')('-',j=1,36)
    write(6,*)

    return

  end subroutine md_print

end module molecular_dynamics
