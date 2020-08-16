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
    ! Calculo das variaveis canônicas a cada ciclo MD                                      *
    !***************************************************************************************

    implicit none

    integer i,ihist,geo_backup
    real(8) t0,t3,time,temp,press,xhi,eta,sigma
    real(8) enpot,ekinet

    !-contagem do tempo absoluto

    call cpu_time(t0)

    !-imprimindo parametros de entrada

    call md_print

    !-preparando Dinâmica molecular

    call md_prepare(xhi,eta,sigma,ekinet,enpot,temp,press)

    write(6,*)('#',i=1,93)
    write(6,*)('MD RUNNING ',i=1,8)
    write(6,*)('#',i=1,93)
    write(6,*)

    write(6,'(2x,a27,f12.4)')'  Correction of VdW energy:',envdw_corr*econv
    write(6,'(2x,a27,f12.4)')'Correction of VdW pressure:',virvdw_corr*econv
    write(6,*)

    !-ciclo MD

    write(6,'(4x,111a1)')('-',i=1,84)
    write(6,10)'##','STEP','TIME','VOLUME','TEMPERATURE' ,'PRESSURE','E(TOTAL)'
    write(6,'(4x,111a1)')('-',i=1,84)

    !-relaxação do sistema

    write(3,5)'i','TIME','VOLUME','TEMPERATURE','PRESSURE','ENERGY'

    geo_backup=-1

    time=dtime
    do i=1,nrelax
       call mdloop(i,geo_backup,xhi,eta,sigma,temp,press,ekinet,enpot)
       write(6,20)'MD',i,time*tconv,volume*rconv**3,&
            temp*teconv,press*pconv,(ekinet+enpot+envdw_corr)*econv
       write(3,30)i,time*tconv,&
            volume*rconv**3,temp*teconv,press*pconv,(ekinet+enpot+envdw_corr)*econv
       time=time+dtime
    end do

    !-sistema em equilibrio termodinâmico

    ihist=1
    do i=nrelax+1,ntrialmax
       call mdloop(i,geo_backup,xhi,eta,sigma,temp,press,ekinet,enpot)
       write(6,20)'MD',i,time*tconv,volume*rconv**3,&
            temp*teconv,press*pconv,(ekinet+enpot+envdw_corr)*econv
       write(3,30)i,time*tconv,&
            volume*rconv**3,temp*teconv,press*pconv,(ekinet+enpot+envdw_corr)*econv
       if(mod(i,nhist).eq.0)call history(ihist)
       time=time+dtime
    end do

    write(6,'(4x,111a1)')('-',i=1,84)
    write(6,*)

    !-contagem do tempo absoluto

    call cpu_time(t3)

    t3=t3-t0

    return

5   format(5x,a8,6x,a4,10x,a6,6x,a11,4x,a8,7x,a6)
10  format(5x,a2,6x,a4,6x,a5,9x,a6,6x,a10,5x,a8,6x,a8)
20  format(5x,a2,2x,i8,2x,es12.4,2x,es12.4,3(2x,es12.4))
30  format(5x,i8,5(2x,es12.4))

  end subroutine md

  subroutine md_prepare(xhi,eta,sigma,ekinet,enpot,temp,press)
    !***************************************************************************************
    ! Preparacao da dinâmica molecular:                                                    *
    ! - Preparando campo de forca;                                                         *
    ! - Preparando lista de vizinhos de Verlet;                                            *
    !***************************************************************************************

    implicit none

    integer i,nwr
    real(8) xhi,eta,sigma,ekinet,temp,press
    real(8) virvdw,virbond,virbend,virtors,vircoul,virtot
    real(8) encoul,enbond,enbend,entors,envdw,enpot

    !-abrindo ficheiro de dados

    open(2,file='HICOLM.md',status='unknown')

    nwr=6

    write(2,'(i12,e12.4,3i7)')int((ntrialmax-nrelax)/nhist),dtime*tconv,nwr,natom,spctot

    !-preparando Campo de Força

    call ff_prepare

    !-preparando lista de vizinhos de Verlet

    call verlet_list_inter

!    if(ntrsff.ne.0)call verlet_list_all

    !-forças iniciais

    if(reuse.eq.0)call ff_modules&
         (encoul,enbond,enbend,entors,envdw,virvdw,virbond,virbend,virtors,vircoul)

    !-energia potencial

    enpot=encoul+enbond+enbend+entors+envdw

    !-correcao da energia de Van der Waals de curto alcance

    call vdw_corr

    !-massa total

    mtot=0.d0
    do i=1,natom
       mtot=mtot+mass(i)
    end do

    !-valores iniciais dos coeficientes de fricção

    xhi=0.d0
    eta=0.d0

    !-calculo energia cinetica desejada

    sigma=0.5d0*nfree*text

    !-valor inicial da energia cinetica

    ekinet=0.d0
    do i=1,natom
       ekinet=ekinet+mass(i)*(vax(i)**2+vay(i)**2+vaz(i)**2)
    end do

    ekinet=0.5d0*ekinet

    !-valor inicial da temperatura

    temp=2.d0*ekinet/nfree

    !-valor inicial da pressão

    virtot=virvdw+virbond+virbend+virtors+vircoul

    press=(2.d0*ekinet+virtot+virvdw_corr)/(3.d0*volume)

    return

  end subroutine md_prepare

  subroutine mdloop(mdstp,geo_backup,xhi,eta,sigma,temp,press,ekinet,enpot)
    !***************************************************************************************
    ! Obtencao das variaveis canonicas;                                                    *
    ! Calculo da energia total;                                                            *
    ! Atualizacao das lista de vizinhos de Verlet                                          *
    ! Impressao das variaveis em ficheiros de saida                                        *
    !***************************************************************************************

    implicit none

    integer mdstp,i,ix,geo_backup
    real(8) temp,press,xhi,eta,sigma

    real(8) virvdw,virbond,virbend,virtors,vircoul
    real(8) ekinet,encoul,enbond,enbend,entors,envdw,enpot

    !-integrando as variáveis canônicas posição e velocidade

    select case(ensble)
    case('nve')
       call nve_vv(temp,press,&
            ekinet,encoul,enbond,enbend,entors,envdw,virvdw,virbond,virbend,virtors,vircoul)
    case('nvt')
       select case(ensble_mt)
       case('berendsen')
          call nvt_ber_vv(sigma,temp,press,ekinet,encoul,enbond,enbend,entors,envdw,&
               virvdw,virbond,virbend,virtors,vircoul)
       case('hoover')
          call nvt_hoover_vv(xhi,sigma,temp,press,ekinet,encoul,enbond,enbend,entors,envdw,&
               virvdw,virbond,virbend,virtors,vircoul)
       end select
    case('npt')
       select case(ensble_mt)
       case('berendsen')
          call npt_ber_vv(temp,press,ekinet,encoul,enbond,enbend,entors,envdw,virvdw, &
               virbond,virbend,virtors,vircoul)
       case('hoover')
          call npt_hoover_vv(xhi,eta,sigma,temp,press,ekinet,&
               encoul,enbond,enbend,entors,envdw,virvdw,virbond,virbend,virtors,vircoul)
       end select
    end select

    !-calculo da energia potencial

    enpot=encoul+enbond+enbend+entors+envdw

    !-imprimindo estrutura

    call geometry(geo_backup)

    !-escrevendo variaveis canonicas em arquivo de dados

    if(mdstp.gt.nrelax.and.mod(mdstp-nrelax,nhist).eq.0)then
       ix=(mdstp-nrelax)/nhist
       write(2,10)ix
       write(2,20)(str(i)*econv/rconv**3,i=1,6)
       write(2,20)(ix*dtime)*nhist*tconv,volume*rconv**3,temp*teconv,&
            press*pconv,(ekinet+envdw_corr)*econv,(ekinet+enpot+envdw_corr)*econv,&
            (mtot/volume)*mconv/(6.02d-4*rconv**3)
       write(2,*)
       do i=1,3
          write(2,40)v(i,1)*rconv,v(i,2)*rconv,v(i,3)*rconv
       end do
       write(2,*)
       do i=1,natom
          write(2,30)idna(i),atp(i),mass(i)*mconv,qat(i),xa(i)*rconv, &
               ya(i)*rconv,za(i)*rconv
          write(2,40)vax(i)*rconv/tconv,vay(i)*rconv/tconv,vaz(i)*rconv/tconv
          write(2,40)fax(i)*econv/rconv,fay(i)*econv/rconv,faz(i)*econv/rconv
       end do
    end if

    !-atualizando a lista de vizinhos de Verlet

    if(mod(mdstp,verlchk).eq.0)call verlet_list_inter

!    if(ntrsff.ne.0)call verlet_list_all

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

    write(6,*)('#',j=1,93)
    write(6,*)('MOLECULAR DYNAMICS ',j=1,5)
    write(6,*)('#',j=1,93)
    write(6,*)
    write(6,'(31x,a30)')'Molecular dynamics information'
    write(6,'(28x,36a1)')('-',j=1,36)
    write(6,'(28x,a12,1x,a3,1x,a9)')'Ensemble:',ensble,ensble_mt
    if(ensble.eq.'nvt')then
       write(6,'(28x,a12,5x,f8.2)')'Thermostat:',tstat*tconv
    elseif(ensble.eq.'npt')then
       write(6,'(28x,a12,5x,f8.2)')'Thermostat:',tstat*tconv
       if(ensble_mt.eq.'berendsen')then
          write(6,'(28x,a12,5x,f8.2,1x,es7.1)')'Barostat:',pstat*tconv,bfactor/pconv
       elseif(ensble_mt.eq.'hoover')then
          write(6,'(28x,a12,5x,f8.2)')'Barostat:',pstat*tconv
       end if
    end if
    write(6,'(28x,a12,5x,f9.3,1x,a1)')'Temperature:',text*teconv,'K'
    write(6,'(28x,a12,5x,f9.3,1x,a3)')'Pressure:',preext*pconv,'atm'
    write(6,'(28x,a12,6x,i10)')'ntrialmax:',ntrialmax
    write(6,'(28x,a12,6x,i10)')'nrelax:',nrelax
    write(6,'(28x,a12,9x,i1)')'reuse:',reuse
    write(6,'(28x,a12,8x,es10.3,1x,a2)')'Timestep:',dtime*tconv,'ps'
    write(6,'(28x,a12,8x,2f6.3,1x,a1)')'rcutoff:',rcutoff*rconv,drcutoff*rconv,'A'
    write(6,'(28x,36a1)')('-',j=1,36)
    write(6,*)

    return

  end subroutine md_print

end module molecular_dynamics
