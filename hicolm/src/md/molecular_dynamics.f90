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

    integer i,ihist
    real(8) t0,t3,time,temp,press,xhi,eta,sigma,mtot
    real(8) enpot,ekinet

    !-contagem do tempo absoluto

    call cpu_time(t0)

    !-imprimindo parametros de entrada

    call md_print

    !-preparando Dinâmica molecular

    call md_prepare(xhi,eta,sigma,mtot,ekinet,enpot,temp,press)

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

    write(3,5)'step',',','time',',','volume',',','temperature',',','pressure',',','ekinet',&
         ',','epotential',',','energy',',','density'

    write(2,40)'step',',','time',',','molecule',',','site',',','Z',',','type',',','mass',',',&
         'charge',',','x',',','y',',','z',',','fx',',','fy',',','fz',',','vx',',','vy',',','vz'

    write(9,50)'step',',','time',',','ax',',','ay',',','az',',','bx',',','by',',','bz',',',&
         'cx',',','cy',',','cz',',','a',',','b',',','c'

    time=dtime
    do i=1,nrelax
       call mdloop(i,xhi,eta,sigma,temp,press,ekinet,enpot)
       if(mod(i,25).eq.0)write(6,20)'MD',i,time*tconv,volume*rconv**3,&
            temp*teconv,press*pconv,(ekinet+enpot+envdw_corr)*econv
       write(3,30)i,',',time*tconv,',',volume*rconv**3,',',temp*teconv,',',press*pconv,',',&
            ekinet*econv,',',(enpot+envdw_corr)*econv,',',(ekinet+enpot+envdw_corr)*econv,&
            ',',(mtot/volume)*mconv/(n0*1.d-24*rconv**3)
       time=time+dtime
    end do

    !-sistema em equilibrio termodinâmico

    ihist=1
    do i=nrelax+1,ntrialmax
       call mdloop(i,xhi,eta,sigma,temp,press,ekinet,enpot)
       if(mod(i,25).eq.0)write(6,20)'MD',i,time*tconv,volume*rconv**3,&
            temp*teconv,press*pconv,(ekinet+enpot+envdw_corr)*econv
       write(3,30)i,',',time*tconv,',',volume*rconv**3,',',temp*teconv,',',press*pconv,',',&
            ekinet*econv,',',(enpot+envdw_corr)*econv,',',(ekinet+enpot+envdw_corr)*econv,&
            ',',(mtot/volume)*mconv/(n0*1.d-24*rconv**3)
       if(mod(i,nhist).eq.0)call history(ihist)
       time=time+dtime
    end do

    write(6,'(4x,111a1)')('-',i=1,84)
    write(6,*)

    !-contagem do tempo absoluto

    call cpu_time(t3)

    t3=t3-t0

    return

5   format(3x,a4,a1,10x,a4,a1,8x,a6,a1,4x,a11,a1,2x,a8,a1,5x,a6,a1,4x,a10,a1,4x,a6,a1,6x,a7)
10  format(5x,a2,6x,a4,6x,a5,9x,a6,6x,a10,5x,a8,6x,a8)
20  format(5x,a2,2x,i8,2x,es12.4,2x,es12.4,3(2x,es12.4))
30  format(1x,i12,a1,8(e12.4,a1))
40  format(3x,a4,a1,10x,a4,a1,7x,a8,a1,2x,a4,a1,a4,a1,a6,a1,4x,a4,a1,7x,a6,a1,8x,3(a1,a1,11x),&
         3(a2,a1,10x),3(a2,a1,10x))
50  format(3x,2(a4,a1,10x),3(a2,a1,10x),3(a2,a1,10x),3(a2,a1,10x),3(a2,a1,10x),2(a2,a1,10x),&
         a2,a1,9x,3(a1,a1,10x))

  end subroutine md

  subroutine md_prepare(xhi,eta,sigma,mtot,ekinet,enpot,temp,press)
    !***************************************************************************************
    ! Preparacao da dinâmica molecular:                                                    *
    ! - Preparando campo de forca;                                                         *
    ! - Preparando lista de vizinhos de Verlet;                                            *
    !***************************************************************************************

    implicit none

    integer i
    real(8) xhi,eta,sigma,ekinet,temp,press,mtot
    real(8) virvdw,virbond,virbend,virtors,vircoul,virtot
    real(8) encoul,enbond,enbend,entors,envdw,enpot

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

  subroutine mdloop(mdstp,xhi,eta,sigma,temp,press,ekinet,enpot)
    !***************************************************************************************
    ! Obtencao das variaveis canonicas;                                                    *
    ! Calculo da energia total;                                                            *
    ! Atualizacao das lista de vizinhos de Verlet                                          *
    ! Impressao das variaveis em ficheiros de saida                                        *
    !***************************************************************************************

    implicit none

    integer mdstp,i,j,k,ix,n1,n2
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

    call geometry(mdstp)

    !-escrevendo variaveis canonicas em arquivo de dados

    if(mdstp.gt.nrelax.and.mod(mdstp-nrelax,nhist).eq.0)then
       ix=(mdstp-nrelax)/nhist
       n1=1
       do i=1,nmolec
          n2=1
          do j=1,ntmolec(i)
             do k =1,nxmolec(i)
                write(2,10)ix,',',ix*dtime*tconv,',',n2,',',k,',',idna(n1),',',&
                     atsp(atp(n1)),',',mass(n1)*mconv,',',qat(n1),',',&
                     xa(n1)*rconv,',',ya(n1)*rconv,',',za(n1)*rconv,',',&
                     fax(n1)*econv/rconv,',',fay(n1)*econv/rconv,',',faz(n1)*econv/rconv,',',&
                     vax(n1)*rconv/tconv,',',vay(n1)*rconv/tconv,',',vaz(n1)*rconv/tconv
                n1=n1+1
             end do
             n2=n2+1
          enddo
       end do
       write(9,20)ix,',',ix*dtime*tconv,',',&
            v(1,1)*rconv,',',v(1,2)*rconv,',',v(1,3)*rconv,',',&
            v(2,1)*rconv,',',v(2,2)*rconv,',',v(2,3)*rconv,',',&
            v(3,1)*rconv,',',v(3,2)*rconv,',',v(3,3)*rconv,',',&
            a*rconv,',',b*rconv,',',c*rconv
    end if

    !-atualizando a lista de vizinhos de Verlet

    if(mod(mdstp,verlchk).eq.0)call verlet_list_inter

!    if(ntrsff.ne.0)call verlet_list_all

    return

10  format(1x,i12,a1,e12.4,a1,2(i8,a1),i5,a1,a5,a1,21(e12.4,a1))
20  format(1x,i12,a1,13(e12.4,a1))

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
