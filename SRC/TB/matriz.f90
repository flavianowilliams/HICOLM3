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
module matriz

  use input
  use estrutura

  real(8) dens

  real(8), allocatable, dimension(:) :: bnd
  real(8), allocatable, dimension(:,:) :: tpr1, spr1
  real(8), allocatable, dimension(:,:,:) :: onst

  complex(8), allocatable, dimension(:,:) :: ms
  complex(8), allocatable, dimension(:,:) :: mh

contains

  subroutine tb_prepare()
    !****************************************************************************************
    ! Alocacao de arrays globais                                                            *
    !****************************************************************************************
    implicit none

    integer i,j,k,kmx,m1,mlomx,m1mx

    mlomx=0
    kmx=0
    m1mx=0
    do i=1,nlo
       do j=1,mlo(i)
          do m1=1,npril(i,j)
             m1mx=max(m1mx,m1)
          end do
          mlomx=max(mlomx,j)
       end do
       do j=1,i
          do m1=1,lo(j)
             call labels(lo(j),lo(i),m1,k)
             kmx=max(k,kmx)
          end do
       end do
    end do

    allocate(ms(norb,norb),mh(norb,norb))
    allocate(tpr1(nparammax,nprllmmax),spr1(nparammax,nprllmmax))
    allocate(onst(lmax,mmax,nprllmmax))
    allocate(bnd(norb))

    return

  end subroutine tb_prepare

  subroutine neighbourhood(a1,a2,l1,l2,m1,m2,kptt,ix,jx)
    !****************************************************************************************
    ! Elementos da matriz secular                                                           *
    !****************************************************************************************
    implicit none

    integer a1,a2,l1,l2,m1,m2,rx,ry,rz,ix,jx
    real(8) xvz,yvz,zvz,kptt(3),arg
    complex(8) fk,tmn,smn

    !-calculando os elementos H(alpha1,alpha2) e S(alpha1,alpha2)

    mh(ix,jx)=0.d0
    ms(ix,jx)=0.d0
    do rx=-rxmx,rxmx
       do ry=-rymx,rymx
          do rz=-rzmx,rzmx
             call ccp(a1,a2,rx,ry,rz,xvz,yvz,zvz)
             arg=kptt(1)*xvz+kptt(2)*yvz+kptt(3)*zvz
             fk=exp(img*arg)
             call wigner_matriz(a1,a2,l1,l2,m1,m2,xvz,yvz,zvz,tmn,smn)
             mh(ix,jx)=mh(ix,jx)+fk*tmn
             ms(ix,jx)=ms(ix,jx)+fk*smn
          end do
       end do
    end do

    return

  end subroutine neighbourhood

  subroutine wigner_matriz(a1,a2,l1,l2,m1,m2,xvz,yvz,zvz,tmn,smn)
    !***************************************************************************************
    ! Hoppings e overlaps de acordo com a formula de Wigner                                *
    !***************************************************************************************
    implicit none

    integer, INTENT(IN) :: a1,a2,l1,l2,m1,m2
    real(8), INTENT(IN) :: xvz,yvz,zvz

    integer mx,ml,mll,i,la,lb,ma,mb,m1b,m2b
    real(8) dr
    complex(8) func,tmn,smn

    la=l1-1
    lb=l2-1
    m1b=m1-(la+1)
    m2b=m2-(lb+1)

    dr=sqrt(xvz**2+yvz**2+zvz**2)

    !-calculando os termos onsite

    if(a1.eq.a2)then
       if(l1.eq.l2)then
          if(m1.eq.m2)then
             tmn=osenv(a1,l1,m1,dr)
             smn=sprenv(l1,dr)
          elseif(m1.ne.m2)then
             tmn=0.d0
             smn=0.d0
          end if
       elseif(l1.ne.l2)then
          tmn=0.d0
          smn=0.d0
       end if
    end if

    !-calculando os termos offsite

    if(a1.ne.a2)then
       if(m1b.ne.0)then
          ml=1
       else
          ml=0
       end if

       if(m2b.ne.0)then
          mll=1
       else
          mll=0
       end if

       tmn=0.d0
       smn=0.d0
       do ma=-ml,ml,2
          do mb=-mll,mll,2
             do mx=-min(la,lb),min(la,lb)
                call labels(l1,l2,abs(mx)+1,i)
                func=conjg(clmmx(m1b,ma))*clmmx(m2b,mb) &
                     *conjg(dlmmx(la,ma,mx,xvz,yvz,zvz)) &
                     *(dlmmx(lb,mb,mx,xvz,yvz,zvz))
                tmn=tmn+func*skv(i,dr)
                smn=smn+func*sks(i,dr)
             end do
          end do
       end do

       !-aplicando a propriedade H(l,l')=H(l',l)*(-)**(l+l')

       if(la.gt.lb)then
          tmn=tmn*(-1.d0)**(la+lb)
          smn=smn*(-1.d0)**(la+lb)
       end if
    end if

    return

  end subroutine wigner_matriz

  double complex function clmmx(m1,m)
    !***************************************************************************************
    ! Coeficientes de transformacao dos harmonicos esfericos para a base real cartesiana   *
    !***************************************************************************************
    implicit none

    integer, INTENT(IN) :: m,m1

    if(m1.lt.0)then
       if(m.lt.0)then
          clmmx=-0.5d0*sqrt(2.d0)*img
       elseif(m.gt.0)then
          clmmx=-0.5d0*sqrt(2.d0)*img*(-1.d0)**(m1+1)
       end if
    elseif(m1.eq.0)then
       clmmx=1.d0
    elseif(m1.gt.0)then
       if(m.lt.0)then
          clmmx=0.5d0*sqrt(2.d0)
       elseif(m.gt.0)then
          clmmx=0.5d0*sqrt(2.d0)*(-1.d0)**(m1)
       end if
    end if

    return

  end function clmmx

  double complex function dlmmx(l,m,ma,xvz,yvz,zvz)
    !***************************************************************************************
    ! Matriz D de rotacao                                                                  *
    !***************************************************************************************
    implicit none

    integer k,kmax,ma,m,l
    real(8) xvz,yvz,zvz,dr,rho,fx,fx2
    complex(8) func

    kmax=max(0,ma-m,l-ma,l+m)

    dr=sqrt(xvz**2+yvz**2+zvz**2)
    rho=sqrt(xvz**2+yvz**2)

    fx=sqrt(0.5d0*(zvz/dr+1.d0))
    fx2=sqrt(0.5d0*(1.d0-zvz/dr))

    dlmmx=0.0d0
    do k=0,kmax
       func=qklmma(k,l,m,ma)*fx**(2*l+m-ma-2*k) &
            *fx2**(2*k-m+ma)*(-1.d0)**(k-m+ma)
       dlmmx=dlmmx+func
    end do

    dlmmx=dlmmx*(xvz/rho-img*yvz/rho)**m ! Angulo de Euler segundo Eu :)
!    dlmmx=dlmmx*(yvz/rho+img*xvz/rho)**m! Angulo de Euler segundo Sakurai :(

    return

  end function dlmmx

  double precision function qklmma(k,l,m,ma)
    !****************************************************************************************
    !Funcao q(k,l,m,m1) usada pela matriz D de rotacao                                      *
    !****************************************************************************************
    implicit none

    integer k,l,m,ma,i
    real(8) pdt,pdt1

    qklmma=0.d0

    if((l+m-k).lt.0)goto 1
    if((l-ma-k).lt.0)goto 1
    if((k-m+ma).lt.0)goto 1

    pdt=1.d0
    do i=1,(l+m)
       pdt=pdt*i
    end do
    do i=1,(l-m)
       pdt=pdt*i
    end do
    do i=1,(l+ma)
       pdt=pdt*i
    end do
    do i=1,(l-ma)
       pdt=pdt*i
    end do

    pdt1=1.d0
    do i=1,(l+m-k)
       pdt1=pdt1*i
    end do
    do i=1,(l-ma-k)
       pdt1=pdt1*i
    end do
    do i=1,k
       pdt1=pdt1*i
    end do
    do i=1,(k-m+ma)
       pdt1=pdt1*i
    end do

    qklmma=sqrt(pdt)/pdt1

1   return

  end function qklmma

  double precision function osenv(a1,l1,m1,dr)
    !****************************************************************************************
    ! Termos onsite da matriz H                                                             *
    !****************************************************************************************
    implicit none

    integer i

    integer, INTENT(IN) :: l1,a1,m1
    real(8), INTENT(IN) :: dr

    osenv=0.d0

    if(dr.le.err)then
       do i=1,npril(l1,m1)
          osenv=osenv+onst(l1,m1,i)*frho(a1)**(2.5*(i-1))
       end do
    end if

    return

  end function osenv

  double precision function sprenv(l1,dr)
    !****************************************************************************************
    ! Termos onsite da matriz S                                                             *
    !****************************************************************************************
    implicit none

    integer, INTENT(IN) :: l1
    real(8), INTENT(IN) :: dr

    sprenv=0.d0

    if(dr.le.err)sprenv=spr(l1)

    return

  end function sprenv

  double precision function skv(i,dr)
    !****************************************************************************************
    ! Termos offsite da matriz H                                                            *
    !****************************************************************************************
    implicit none

    integer, INTENT(IN) :: i
    real(8), INTENT(IN) :: dr

    tpr1(i,2)=abs(tpr1(i,2))
    tpr1(i,3)=abs(tpr1(i,3))
    tpr1(i,4)=abs(tpr1(i,4))

!    skv=(tpr1(i,1)+tpr1(i,3)*dr)*exp(-abs(tpr1(i,2))*dr)*fcutt(dr)
!    skv=(tpr1(i,1)*dr**(-tpr1(i,2)))*exp(-tpr1(i,3)*dr)*fcutt(dr)
    skv=tpr1(i,1)*dr**(-tpr1(i,4))*exp(-tpr1(i,2)*dr**tpr1(i,3))*fcutt(dr)

  end function skv

  double precision function sks(i,dr)
    !***************************************************************************************
    !Funcao que define os termos onsite da matriz S                                        *
    !***************************************************************************************
    implicit none

    integer, INTENT(IN) :: i
    real(8), INTENT(IN) :: dr

    spr1(i,2)=abs(spr1(i,2))
    spr1(i,3)=abs(spr1(i,3))
    spr1(i,4)=abs(spr1(i,4))

!    sks=(spr1(i,1)+spr1(i,3)*dr)*exp(-abs(spr1(i,2))*dr)*fcutt(dr)
!    sks=(spr1(i,1)*dr**(-spr1(i,2)))*exp(-spr1(i,3)*dr)*fcutt(dr)
    sks=spr1(i,1)*dr**(-spr1(i,4))*exp(-spr1(i,2)*dr**spr1(i,3))*fcutt(dr)

  end function sks

  double precision function frho(a1)
    !***************************************************************************************
    ! Densidade eletronica do atomo a1                                                     *
    !***************************************************************************************
    implicit none

    integer, INTENT(IN) :: a1

    integer i,rx,ry,rz
    real(8) dr,xvz,yvz,zvz

    frho=0.d0
    do i=1,natom
       do rx=-3,3
          do ry=-3,3
             do rz=-3,3
                call ccp(i,a1,rx,ry,rz,xvz,yvz,zvz)
                dr=sqrt(xvz**2+yvz**2+zvz**2)
                frho=frho+exp(-abs(dens)*dr**2)*fcutt(dr)
             end do
          end do
       end do
    end do

    return

  end function frho

  double precision function fcutt(dr)
    !***************************************************************************************
    ! Ocupacao eletronica de Fermi                                                         *
    !***************************************************************************************
    implicit none

    real(8) dr

    fcutt=1.d0/(1.d0+exp((dr-rcutt)/lcutt))

  end function fcutt

end module matriz
