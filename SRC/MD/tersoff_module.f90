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
module tersoff_module
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
  use utils
  use estrutura
  use alloc_arrays
  use neighbour_list

contains

  subroutine tersoff_convert
    !****************************************************************************************
    ! Conversao de unidades de medida:                                                      *
    ! Unidades de entrada ---> a.u.                                                         *
    !****************************************************************************************

    implicit none

    integer i,j

    !-convertendo unidades de medida

    do i=1,spctot
       do j=1,spctot
          select case(tersoff)
          case(1)
             partrsff(i,j,1)=partrsff(i,j,1)/rconv
             partrsff(i,j,2)=partrsff(i,j,2)/rconv
             partrsff(i,j,3)=partrsff(i,j,3)/econv
             partrsff(i,j,4)=partrsff(i,j,4)/kconv
             partrsff(i,j,6)=partrsff(i,j,6)/rconv
          end select
       end do
    end do

    return

  end subroutine tersoff_convert

  subroutine tersoff_calc(entrsff,virtrsff)
    !****************************************************************************************
    ! - Energia potencial;                                                                  *
    ! - Contribuicao para o virial;                                                         *
    ! - Contribuicao para o stress;                                                         *
    ! - Forças atômicas.                                                                    *
    !****************************************************************************************

    implicit none

    integer i,j,ni,nj
    real(8) xvz,yvz,zvz,pot
    real(8) entrsff,virtrsff

    virtrsff=0.d0

    do i=1,natom
       do j=1,nlista(i)
          ni=i
          nj=ilista(i,j)
          call mic(ni,nj,xvz,yvz,zvz)
          call tersoff_force(ni,nj,xvz,yvz,zvz,pot)
          entrsff=entrsff+pot
       end do
    end do

    return

  end subroutine tersoff_calc

  subroutine tersoff_force(i1,i2,xvz,yvz,zvz,pot)
    !****************************************************************************************
    ! - Contribuicao para o virial                                                          *
    ! - Contribuicao para o stress                                                          *
    ! - Forcas atomicas                                                                     *
    !****************************************************************************************

    implicit none

    integer i,ii,i1,i2
    real(8) pot,potp,potc,potg1,potg2,potg,dpotp,potp1,potp2
    real(8) xvz,yvz,zvz,dr,prm(4),fbij(3),fbij1(3),fbij2(3),fcij(3),fra(3)

    !-comecando a zica!!!!!!!!!!!!!<------>:(

    do ii=1,4
       prm(ii)=partrsff(atp(i1),atp(i2),ii+2)
    end do

    dr=sqrt(xvz**2+yvz**2+zvz**2)

    !-atribuindo termo de restricao da ligação

    call tersoff_restrict(i1,i2,dr,xvz,yvz,zvz,potc,fcij)

    !-atribuindo termo de dependencia angular

    call tersoff_bij(i1,i2,xvz,yvz,zvz,potg1,fbij1)
    call tersoff_bij(i2,i1,-xvz,-yvz,-zvz,potg2,fbij2)

    potg=0.5d0*(potg1+potg2)

    do i=1,3
       fbij(i)=0.5d0*(fbij1(i)+fbij2(i))
    end do

    !-VR(r)+bij*VA(r)

    potp1=prm(1)*exp(-sqrt(2.d0*prm(3))*prm(2)*(dr-prm(4)))/(prm(3)-1.d0)
    potp2=prm(1)*prm(3)*exp(-sqrt(2.d0/prm(3))*prm(2)*(dr-prm(4)))/(prm(3)-1.d0)

    potp=potp1-potg*potp2

    !-energia potencial U(r)

    pot=0.5d0*potc*potp

    !-dVR/dr+bij*dVA/dr

    dpotp=-prm(1)*prm(2)*sqrt(2.d0*prm(3))*(exp(-sqrt(2.d0*prm(3))*prm(2)*(dr-prm(4))) &
         -potg*exp(-sqrt(2.d0/prm(3))*prm(2)*(dr-prm(4))))/(prm(3)-1.d0)

    !-(1/r)*Vc(r)*(dVR/dr-bij*dVA/dr)

    dpotp=-potc*dpotp/dr

    fra(1)=dpotp*xvz
    fra(2)=dpotp*yvz
    fra(3)=dpotp*zvz

    !-(VR-bij*VA)*dVc/dr

    fcij(1)=-potp*fcij(1)
    fcij(2)=-potp*fcij(2)
    fcij(3)=-potp*fcij(3)

    !-Vc(r)*VA(r)*dbij/dr

    do i=1,3
       fbij(i)=potc*potp2*fbij(i)
    end do

    !-contribuicao para as forcas atomicas

    fax(i1)=fax(i1)+0.5d0*(-fra(1)-fcij(1)+fbij(1))
    fay(i1)=fay(i1)+0.5d0*(-fra(2)-fcij(2)+fbij(2))
    faz(i1)=faz(i1)+0.5d0*(-fra(3)-fcij(3)+fbij(3))

    return

  end subroutine tersoff_force

  subroutine tersoff_restrict(i,j,dr,xvz,yvz,zvz,potc,fcij)
    !****************************************************************************************
    ! Funcao de corte                                                                       *
    !****************************************************************************************

    implicit none

    integer i,j
    real(8) dr,xvz,yvz,zvz,potc,dpotc,prm(2),fcij(3)

    !-potencial de restricao

    prm(1)=partrsff(atp(i),atp(j),1)
    prm(2)=partrsff(atp(i),atp(j),2)

    if(dr.le.prm(1))then
       potc=1.d0
       dpotc=0.d0
    elseif(dr.gt.prm(1).and.dr.le.prm(2))then
       potc=0.5d0+0.5d0*cos(pi*(dr-prm(1))/(prm(2)-prm(1)))
       dpotc=-0.5d0*pi*sin(pi*(dr-prm(1))/(prm(2)-prm(1)))/(prm(2)-prm(1))
    elseif(dr.gt.prm(2))then
       potc=0.d0
       dpotc=0.d0
    end if

    dpotc=dpotc/dr

    !-gradiente do termo fc(r)

    fcij(1)=dpotc*xvz
    fcij(2)=dpotc*yvz
    fcij(3)=dpotc*zvz

    return

  end subroutine tersoff_restrict

  subroutine tersoff_bij(i1,i2,xvz,yvz,zvz,potg,fbij)
    !****************************************************************************************
    ! Termo de ajuste da atracao dependente das ligacoes entre vizinhos                     *
    !****************************************************************************************

    implicit none

    integer i,j,i1,i2,ii,nj
    real(8) potg,dpotg,potc,drij(3),drik(3),xvz,yvz,zvz,dr1,dr2,prm(4),fbij(3),fcik(3),fbi(3)
    real(8) gtht,sum

    !-valores iniciais

    drij(1)=xvz
    drij(2)=yvz
    drij(3)=zvz

    dr1=sqrt(drij(1)**2+drij(2)**2+drij(3)**2)

    !-derivadas do termo bij

    do ii=1,4
       prm(ii)=partrsff(atp(i1),atp(i2),ii+6)
    end do

    do i=1,3
       fbij(i)=0.d0
    end do

    sum=0.d0
    do i=1,nlista(i1)
       nj=ilista(i1,i)
       if(nj.eq.i2)cycle
       call mic(i1,nj,drik(1),drik(2),drik(3))
       dr2=sqrt(drik(1)**2+drik(2)**2+drik(3)**2)
       call tersoff_restrict(i1,nj,dr2,drik(1),drik(2),drik(3),potc,fcik)
       call tersoff_bij_gtht(prm,dr1,drij,dr2,drik,gtht,fbi)
       do j=1,3
          fbij(j)=fbij(j)-gtht*fcik(j)+potc*fbi(j)
       end do
       sum=sum+potc*gtht
    end do

    !-contribuicao para a energia potencial

    potg=(1.d0+sum)**(-prm(4))

    !-contribuicao para a forca atomica

    dpotg=-prm(4)*potg**(-1)

    return

  end subroutine tersoff_bij

  subroutine tersoff_bij_gtht(prm,dr1,drij,dr2,drik,gtht,fbi)
    !****************************************************************************************
    ! Gradiente do termo de ajuste                                                          *
    !****************************************************************************************

    implicit none

    integer i
    real(8) dr1,dr2,cost,fa,drij(3),drik(3),prm(4),gtht,fbi(3)

    cost=(drij(1)*drik(1)+drij(2)*drik(2)+drij(3)*drik(3))/(dr1*dr2)

    !-g(theta)

    gtht=prm(1)*(1.d0+(prm(2)/prm(3))**2-prm(3)**2/(prm(3)**2+(1.d0+cost)**2))

    !-gradiente de g(theta)

    do i=1,3
       fbi(i)=drik(i)/(dr1*dr2)+drij(i)/(dr1*dr2)-cost*(drij(i)/dr1**2+drik(i)/dr2**2)
    end do

    fa=2.d0*(1.d0+cost)*prm(1)*prm(3)**2/(prm(3)**2+(1.d0+cost)**2)**2

    fbi(1)=fa*fbi(1)
    fbi(2)=fa*fbi(2)
    fbi(3)=fa*fbi(3)

    return

  end subroutine tersoff_bij_gtht

end module tersoff_module
