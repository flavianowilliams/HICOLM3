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
module ensemble_nvt
  !****************************************************************************************
  ! Integradores de movimento do ensemble NVT                                             *
  !                                                                                       *
  ! Flaviano Williams Fernandes, 05 de dezembro de 2018                                   *
  !****************************************************************************************

  use input
  use estrutura
  use force_field

contains

  subroutine nvt_ber_vv(sigma,temp,press,ekinet,encoul,enbond,enbend,entors,envdw,&
       virvdw,virbond,virbend,virtors,vircoul)
    !*************************************************************************************
    ! Controle da temperatura pelo algoritmo Berendsen                                   *
    !*************************************************************************************

    implicit none

    integer i
    real(8) temp,press,xhi,ekinet,sigma
    real(8) virvdw,virbond,virbend,virtors,vircoul,virtot
    real(8) encoul,enbond,enbend,entors,envdw

    !-calculando posições no instante t+dt posterior

    do i=1,natom
       vax(i)=vax(i)+fax(i)*(0.5d0*dtime)*fztp(i)/mass(i)
       vay(i)=vay(i)+fay(i)*(0.5d0*dtime)*fztp(i)/mass(i)
       vaz(i)=vaz(i)+faz(i)*(0.5d0*dtime)*fztp(i)/mass(i)
       xa(i)=xa(i)+vax(i)*dtime
       ya(i)=ya(i)+vay(i)*dtime
       za(i)=za(i)+vaz(i)*dtime
    end do

    !-aplicando condicoes de contorno

    call ccp_inv

    !-calculo das forças no instante t+dt posterior

    call ff_modules&
         (encoul,enbond,enbend,entors,envdw,virvdw,virbond,virbend,virtors,vircoul)

    !-calculando velocidades no instante t+dt posterior

    do i=1,natom
       vax(i)=vax(i)+fax(i)*(0.5d0*dtime)*fztp(i)/mass(i)
       vay(i)=vay(i)+fay(i)*(0.5d0*dtime)*fztp(i)/mass(i)
       vaz(i)=vaz(i)+faz(i)*(0.5d0*dtime)*fztp(i)/mass(i)
    end do

    !-calculo da energia cinetica

    ekinet=0.d0
    do i=1,natom
       ekinet=ekinet+mass(i)*(vax(i)**2+vay(i)**2+vaz(i)**2)
    end do

    ekinet=0.5d0*ekinet

    !-preparando parametro de escalonamento da temperatura

    xhi=sqrt(1.d0+dtime*(sigma/ekinet-1.d0)/tstat)

    !-rescalonando velocidades

    do i=1,natom
       vax(i)=vax(i)*xhi
       vay(i)=vay(i)*xhi
       vaz(i)=vaz(i)*xhi
    end do

    !-recalculando temperatura

    temp=2.d0*ekinet/nfree

    !-calculo da pressao

    virtot=virvdw+virbond+virbend+virtors+vircoul

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

  end subroutine nvt_ber_vv

  subroutine nvt_hoover_vv(xhi,sigma,temp,press,ekinet,encoul,enbond,enbend,entors,envdw,&
       virvdw,virbond,virbend,virtors,vircoul)
    !*************************************************************************************
    ! Controle da temperatura pelo algoritmo Berendsen                                   *
    !*************************************************************************************

    implicit none

    integer i
    real(8) temp,press,xhi,ekinet,sigma,qmass
    real(8) virvdw,virbond,virbend,virtors,vircoul,virtot
    real(8) encoul,enbond,enbend,entors,envdw

    !-calculando a massa do termostato

    qmass=2.d0*sigma*tstat**2

    !-calculando o coeficiente de fricção xhi no instante t+0.5dt

    xhi=xhi+0.5d0*dtime*(ekinet-sigma)/qmass

    !-calculando posições no instante t+0.5dt

    do i=1,natom
       vax(i)=vax(i)*(1.d0-0.5d0*dtime*xhi*fztp(i))+fax(i)*(0.5d0*dtime)*fztp(i)/mass(i)
       vay(i)=vay(i)*(1.d0-0.5d0*dtime*xhi*fztp(i))+fay(i)*(0.5d0*dtime)*fztp(i)/mass(i)
       vaz(i)=vaz(i)*(1.d0-0.5d0*dtime*xhi*fztp(i))+faz(i)*(0.5d0*dtime)*fztp(i)/mass(i)
       xa(i)=xa(i)+vax(i)*dtime
       ya(i)=ya(i)+vay(i)*dtime
       za(i)=za(i)+vaz(i)*dtime
    end do

    !-aplicando condicoes de contorno

    call ccp_inv

    !-calculo das forças no instante t+dt posterior

    call ff_modules&
         (encoul,enbond,enbend,entors,envdw,virvdw,virbond,virbend,virtors,vircoul)

    !-calculando velocidades no instante t+dt

    do i=1,natom
       vax(i)=vax(i)+fax(i)*(0.5d0*dtime)*fztp(i)/mass(i)
       vay(i)=vay(i)+fay(i)*(0.5d0*dtime)*fztp(i)/mass(i)
       vaz(i)=vaz(i)+faz(i)*(0.5d0*dtime)*fztp(i)/mass(i)
    end do

    !-calculo da energia cinetica

    ekinet=0.d0
    do i=1,natom
       ekinet=ekinet+mass(i)*(vax(i)**2+vay(i)**2+vaz(i)**2)
    end do

    ekinet=0.5d0*ekinet

    !-calculando o coeficiente de fricção xhi no instante t+dt

    xhi=xhi+0.5d0*dtime*(ekinet-sigma)/qmass

    !-recalculando velocidades no instante t+dt

    do i=1,natom
       vax(i)=vax(i)*(1.d0-0.5d0*dtime*xhi*fztp(i))+fax(i)*(0.5d0*dtime)*fztp(i)/mass(i)
       vay(i)=vay(i)*(1.d0-0.5d0*dtime*xhi*fztp(i))+fay(i)*(0.5d0*dtime)*fztp(i)/mass(i)
       vaz(i)=vaz(i)*(1.d0-0.5d0*dtime*xhi*fztp(i))+faz(i)*(0.5d0*dtime)*fztp(i)/mass(i)
    end do

    !-recalculando a energia cinetica

    ekinet=0.d0
    do i=1,natom
       ekinet=ekinet+mass(i)*(vax(i)**2+vay(i)**2+vaz(i)**2)
    end do

    ekinet=0.5d0*ekinet

    !-calculando propriedades termodinanicas

    !-temperatura

    temp=2.d0*ekinet/nfree

    !-pressao

    virtot=virvdw+virbond+virbend+virtors+vircoul

    press=(2.d0*ekinet+virtot+virvdw_corr)/(3.d0*volume)

    !-stress

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

  end subroutine nvt_hoover_vv

end module ensemble_nvt
