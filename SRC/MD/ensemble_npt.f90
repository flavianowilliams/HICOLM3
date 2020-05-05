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
module ensemble_npt
  !****************************************************************************************
  ! Integradores de movimento do ensemble NPT                                             *
  !                                                                                       *
  ! Flaviano Williams Fernandes, 05 de dezembro de 2018                                   *
  !****************************************************************************************

  use input
  use estrutura
  use force_field

contains

  subroutine npt_ber_vv(temp,press,ekinet,encoul,enbond,enbend,entors,envdw,entrsff,virvdw, &
       virbond,virbend,virtors,vircoul,virtrsff)
    !*************************************************************************************
    ! Controle da temperatura e pressão pelo algoritmo Berendsen                         *
    !*************************************************************************************

    implicit none

    integer i,j
    real(8) temp,press,xhi,eta,ekinet,sigma
    real(8) virvdw,virbond,virbend,virtors,vircoul,virtrsff,virtot
    real(8) encoul,enbond,enbend,entors,envdw,entrsff

    !-preparando parametro de escalonamento isotrópico

    eta=(1.d0-bfactor*dtime*(preext-press)/pstat)**(1.d0/3.d0)

    !-calculando posições no instante t+dt posterior

    do i=1,natom
       vax(i)=vax(i)+fax(i)*(0.5d0*dtime)*fztp(i)/mass(i)
       vay(i)=vay(i)+fay(i)*(0.5d0*dtime)*fztp(i)/mass(i)
       vaz(i)=vaz(i)+faz(i)*(0.5d0*dtime)*fztp(i)/mass(i)
       xa(i)=xa(i)*eta+vax(i)*dtime
       ya(i)=ya(i)*eta+vay(i)*dtime
       za(i)=za(i)*eta+vaz(i)*dtime
    end do

    !-redefinindo volume e vetores de rede

    do i=1,3
       do j=1,3
          v(i,j)=v(i,j)*eta
       end do
    end do

    volume=volume*eta**3

    a=sqrt(v(1,1)**2+v(1,2)**2+v(1,3)**2)
    b=sqrt(v(2,1)**2+v(2,2)**2+v(2,3)**2)
    c=sqrt(v(3,1)**2+v(3,2)**2+v(3,3)**2)

    !-checando dimensoes da caixa

    if((rcutoff+drcutoff).gt.0.5d0*a)stop 'npt_ber_vv: rcutoff exceeds the half-box size'
    if((rcutoff+drcutoff).gt.0.5d0*b)stop 'npt_ber_vv: rcutoff exceeds the half-box size'
    if((rcutoff+drcutoff).gt.0.5d0*c)stop 'npt_ber_vv: rcutoff exceeds the half-box size'

    !-aplicando condicoes de contorno inversa

    call ccp_inv

    !-calculo das forças no instante t+dt posterior

    call ff_modules(encoul,enbond,enbend,entors,envdw,entrsff,virvdw,virbond, &
         virbend,virtors,vircoul,virtrsff)

    !-calculando velocidades no instante t+dt posterior

    do i=1,natom
       vax(i)=vax(i)+fax(i)*(0.5d0*dtime)*fztp(i)/mass(i)
       vay(i)=vay(i)+fay(i)*(0.5d0*dtime)*fztp(i)/mass(i)
       vaz(i)=vaz(i)+faz(i)*(0.5d0*dtime)*fztp(i)/mass(i)
    end do

    !-calculo da energia cinetica total

    ekinet=0.d0
    do i=1,natom
       ekinet=ekinet+mass(i)*(vax(i)**2+vay(i)**2+vaz(i)**2)
    end do

    ekinet=0.5d0*ekinet

    !-preparando parametro de escalonamento

    sigma=0.5d0*nfree*text

    xhi=sqrt(1.d0+dtime*(sigma/ekinet-1.d0)/tstat)

    !-rescalonando velocidades

    do i=1,natom
       vax(i)=vax(i)*xhi
       vay(i)=vay(i)*xhi
       vaz(i)=vaz(i)*xhi
    end do

    !-calculo da energia cinetica

    ekinet=0.d0
    do i=1,natom
       ekinet=ekinet+mass(i)*(vax(i)**2+vay(i)**2+vaz(i)**2)
    end do

    ekinet=0.5d0*ekinet

    !-calculo da temperatura

    temp=2.d0*ekinet/nfree

    !-calculo da pressao

    virtot=virvdw+virbond+virbend+virtors+vircoul+virtrsff

    press=(2.d0*ekinet+virtot+virvdw_corr)/(3.d0*volume)

    !-calculo do stress

    do i=1,natom
       str(1)=str(1)-mass(i)*vax(i)*vax(i)
       str(2)=str(2)-mass(i)*vay(i)*vay(i)
       str(3)=str(3)-mass(i)*vaz(i)*vaz(i)
       str(4)=str(4)-mass(i)*vaz(i)*vay(i)
       str(5)=str(5)-mass(i)*vax(i)*vaz(i)
       str(6)=str(6)-mass(i)*vay(i)*vax(i)
    end do

    do i=1,6
       str(i)=str(i)/volume
    end do

    return

  end subroutine npt_ber_vv

  subroutine npt_hoover_vv(vxhi,veta,sigma,temp,press,ekinet,encoul,enbond,enbend,entors,&
       envdw,entrsff,virvdw,virbond,virbend,virtors,vircoul,virtrsff)
    !*************************************************************************************
    ! Controle da temperatura e pressão pelo algoritmo Berendsen                         *
    !*************************************************************************************

    implicit none

    integer i,j
    real(8) temp,press,vxhi,veta,ekinet,sigma,qmass,pmass,vxhi0,veta0
    real(8) virvdw,virbond,virbend,virtors,vircoul,virtrsff,virtot
    real(8) encoul,enbond,enbend,entors,envdw,entrsff

    !-massa do termostato e barostato

    qmass=2.d0*sigma*tstat**2
    pmass=(nfree+3)*text*pstat**2

    !-parametros de escalonamento isotrópico no instante t+0.5*dt

    vxhi0=vxhi
    veta0=veta

    vxhi=vxhi0+0.5d0*dtime*(2.d0*ekinet+pmass*veta0**2-2.d0*sigma-text)/qmass
    veta=veta0+0.5d0*dtime*(3.d0*volume*(press-preext)/pmass-vxhi0*veta0)

    !-ajustando velocidades através do termostato no instante t

    do i=1,natom
       vax(i)=vax(i)*(1.d0-0.5d0*dtime*vxhi)
       vay(i)=vay(i)*(1.d0-0.5d0*dtime*vxhi)
       vaz(i)=vaz(i)*(1.d0-0.5d0*dtime*vxhi)
    end do

    !-ajustando velocidades através do barostato no instante t

    do i=1,natom
       vax(i)=vax(i)*(1.d0-0.5d0*dtime*veta)
       vay(i)=vay(i)*(1.d0-0.5d0*dtime*veta)
       vaz(i)=vaz(i)*(1.d0-0.5d0*dtime*veta)
    end do

    !-coordenadas no instante t+0.5*dt

    do i=1,natom
       vax(i)=vax(i)+0.5d0*dtime*fax(i)/mass(i)
       vay(i)=vay(i)+0.5d0*dtime*fay(i)/mass(i)
       vaz(i)=vaz(i)+0.5d0*dtime*faz(i)/mass(i)
       xa(i)=xa(i)+vax(i)*dtime*fztp(i)
       ya(i)=ya(i)+vay(i)*dtime*fztp(i)
       za(i)=za(i)+vaz(i)*dtime*fztp(i)
    end do

    !-volume e vetores de rede no instante t+dt

    do i=1,3
       do j=1,3
          v(i,j)=v(i,j)*exp(dtime*veta0)
       end do
    end do

    volume=volume*exp(3.d0*dtime*veta0)

    a=sqrt(v(1,1)**2+v(1,2)**2+v(1,3)**2)
    b=sqrt(v(2,1)**2+v(2,2)**2+v(2,3)**2)
    c=sqrt(v(3,1)**2+v(3,2)**2+v(3,3)**2)

    !-checando dimensoes da caixa

    if((rcutoff+drcutoff).gt.0.5d0*a)stop 'npt_ber_vv: rcutoff exceeds the half-box size'
    if((rcutoff+drcutoff).gt.0.5d0*b)stop 'npt_ber_vv: rcutoff exceeds the half-box size'
    if((rcutoff+drcutoff).gt.0.5d0*c)stop 'npt_ber_vv: rcutoff exceeds the half-box size'

    !-aplicando condicoes de contorno inversa

    call ccp_inv

    !-forças no instante t+dt

    call ff_modules(encoul,enbond,enbend,entors,envdw,entrsff,virvdw,virbond, &
         virbend,virtors,vircoul,virtrsff)

    !-velocidades no instante t+dt

    do i=1,natom
       vax(i)=vax(i)+0.5d0*dtime*fax(i)/mass(i)
       vay(i)=vay(i)+0.5d0*dtime*fay(i)/mass(i)
       vaz(i)=vaz(i)+0.5d0*dtime*faz(i)/mass(i)
    end do

    !-energia cinetica total no instante t+dt

    ekinet=0.d0
    do i=1,natom
       ekinet=ekinet+mass(i)*(vax(i)**2+vay(i)**2+vaz(i)**2)
    end do

    ekinet=0.5d0*ekinet

    !-temperatura no instante t+dt

    temp=2.d0*ekinet/nfree

    !-pressao no instante t+dt

    virtot=virvdw+virbond+virbend+virtors+vircoul+virtrsff
    press=(2.d0*ekinet+virtot+virvdw_corr)/(3.d0*volume)

    !-parametros de escalonamento isotrópico no instante t+dt

    vxhi0=vxhi
    veta0=veta

    vxhi=vxhi0+0.5d0*dtime*(2.d0*ekinet+pmass*veta0**2-2.d0*sigma-text)/qmass
    veta=veta0+0.5d0*dtime*(3.d0*volume*(press-preext)/pmass-vxhi0*veta0)

    vxhi0=vxhi
    veta0=veta

    vxhi=vxhi0+0.5d0*dtime*(2.d0*ekinet+pmass*veta0**2-2.d0*sigma-text)/qmass
    veta=veta0+0.5d0*dtime*(3.d0*volume*(press-preext)/pmass-vxhi0*veta0)

    !-ajustando velocidades através do barostato no instante t+dt

    do i=1,natom
       vax(i)=vax(i)-0.5d0*dtime*veta*vax(i)
       vay(i)=vay(i)-0.5d0*dtime*veta*vay(i)
       vaz(i)=vaz(i)-0.5d0*dtime*veta*vaz(i)
    end do

    !-energia cinetica total no instante t+dt

!    ekinet=0.d0
!    do i=1,natom
!       ekinet=ekinet+mass(i)*(vax(i)**2+vay(i)**2+vaz(i)**2)
!    end do

!    ekinet=0.5d0*ekinet

    !-pressao no instante t+dt

!    virtot=virvdw+virbond+virbend+virtors+vircoul+virtrsff

!    press=(2.d0*ekinet+virtot+virvdw_corr)/(3.d0*volume)

    !-ajustando parametros de escalonamento isotrópico no instante t+dt

    !    vxhi0=vxhi
    !    veta0=veta

    !    vxhi=vxhi0+0.5d0*dtime*(2.d0*ekinet+pmass*veta0**2-2.d0*sigma-text)/qmass
    !    veta=veta0+0.5d0*dtime*(3.d0*volume*(press-preext)/pmass-vxhi0*veta0)

    !-ajustando velocidades através do termostato no instante t+dt

    do i=1,natom
       vax(i)=vax(i)-0.5d0*dtime*vxhi*vax(i)
       vay(i)=vay(i)-0.5d0*dtime*vxhi*vay(i)
       vaz(i)=vaz(i)-0.5d0*dtime*vxhi*vaz(i)
       vax(i)=vax(i)*fztp(i)
       vay(i)=vay(i)*fztp(i)
       vaz(i)=vaz(i)*fztp(i)
    end do

    !-energia cinetica total no instante t+dt

    ekinet=0.d0
    do i=1,natom
       ekinet=ekinet+mass(i)*(vax(i)**2+vay(i)**2+vaz(i)**2)
    end do

    ekinet=0.5d0*ekinet

    !-temperatura no instante t+dt

    temp=2.d0*ekinet/nfree

    !-pressao no instante t+dt

    virtot=virvdw+virbond+virbend+virtors+vircoul+virtrsff

    press=(2.d0*ekinet+virtot+virvdw_corr)/(3.d0*volume)

    !-stress no instante t+dt

    do i=1,natom
       str(1)=str(1)-mass(i)*vax(i)*vax(i)
       str(2)=str(2)-mass(i)*vay(i)*vay(i)
       str(3)=str(3)-mass(i)*vaz(i)*vaz(i)
       str(4)=str(4)-mass(i)*vaz(i)*vay(i)
       str(5)=str(5)-mass(i)*vax(i)*vaz(i)
       str(6)=str(6)-mass(i)*vay(i)*vax(i)
    end do

    do i=1,6
       str(i)=str(i)/volume
    end do

    return

  end subroutine npt_hoover_vv

end module ensemble_npt
