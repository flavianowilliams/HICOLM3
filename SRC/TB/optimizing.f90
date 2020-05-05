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
module optimizing
  !******************************************************************************************
  !Modulo de otimizacao dos parametros Slater-Koster                                        *
  !                                                                                         *
  !1-ss_sigma                                                                               *
  !2-sp_sigma                                                                               *
  !3-pp_sigma                                                                               *
  !4-pp_pi                                                                                  *
  !                                                                                         *
  ! Flaviano Williams Fernandes, 05 de dezembro de 2018                                     *
  !******************************************************************************************
  use sistema
  use input
  use brillouin
  use estrutura
  use fitting
  use diagonalize
  use matriz

  implicit none

  integer kz(nkrlxmax,dkrlxmax),kw,INFOPHM,INFOLMA,npro(lmax),nx
  real(8) kptl(nkrlxmax,dkrlxmax,3),FNORMPHM,FNORMLMA

  save FNORMPHM,FNORMLMA,INFOPHM,INFOLMA

contains

  subroutine optimize
    !****************************************************************************************
    ! Escolha de pontos K                                                                   *
    !****************************************************************************************
    implicit none

    integer i,j,k,l,m1,nkpt,i0,chv,nx
    real(8) kpt(nkmax,3),mkpt(nbandmax)
    real(8) knorm,knorm2,dxpr

    dens=dens0

    do i=1,nlo
       do j=1,mlo(i)
          do k=1,npril(i,j)
             onst(i,j,k)=onst0(i,j,k)
          end do
       end do
    end do

    do i=1,nparamt
       do j=1,nprllm(i)
          tpr1(mskt(i),j)=tpr0(mskt(i),j)
          spr1(mskt(i),j)=spr0(mskt(i),j)
       end do
    end do

    if(nopt.eq.0)goto 3
    if(nopt.ge.1)goto 1

1   call high_symmetry(1,kpt,mkpt,nkpt)
!-------------------------------------
    do j=1,nkrlx
       do l=1,dkrlx(j)
          chv=1
          i0=0
          knorm2=0.d0
          knorm=0.d0
          do i=1,nkpt
             knorm2=knorm2+mkpt(chv)/dband(chv)
             if(krlx(j,l).ge.knorm)then
                if(krlx(j,l).le.knorm2)then
                   do k=1,3
                      kptl(j,l,k)=kpt(i,k)
                   end do
                   kz(j,l)=i
                end if
             end if
             if(i.ge.(i0+dband(chv)))then
                chv=chv+1
                i0=i
             end if
             knorm=knorm2
          end do
       end do
    end do

    nx=1

    write(6,*)('#',j=1,70)
    write(6,*)'$$ ',('FITTING ',j=1,8),'$$$'
    write(6,*)('#',j=1,70)
    write(6,*)
    write(6,*)'    Matrix:',norb
    write(6,*)
    write(6,43)('-',j=1,27)

2   call phm

    dxpr=abs(dens-dens0)

    do i=1,nlo
       do j=1,mlo(i)
          do k=1,npril(i,j)
             dxpr=max(dxpr,abs(onst(i,j,k)-onst0(i,j,k)))
          end do
       end do
    end do

    do i=1,nskpar
       do j=1,nprllm(i)
          if(chk(i).eqv..true.)then
             dxpr=max(dxpr,abs(tpr1(msk(i),j)-tpr0(msk(i),j)))
             dxpr=max(dxpr,abs(spr1(msk(i),j)-spr0(msk(i),j)))
          end if
       end do
    end do

    dens0=dens

    do i=1,nlo
       do j=1,mlo(i)
          do k=1,npril(i,j)
             onst0(i,j,k)=onst(i,j,k)
          end do
       end do
    end do

    do i=1,nparamt
       do j=1,nprllm(i)
          tpr0(mskt(i),j)=tpr1(mskt(i),j)
          spr0(mskt(i),j)=spr1(mskt(i),j)
       end do
    end do

    write(6,42)nx,dxpr,FNORMPHM

    if(dxpr.le.1.d-4)goto 3

    nx=nx+1

    goto 2

    !-escrevendo dados em arquivo de saida

3   write(6,43)('-',j=1,27)
    write(6,*)

    nx=1
    do j=1,nlo
       npro(nx)=lo(j)
       nx=nx+1
    end do

    !-convertendo unidades de medida

    dens0=dens0*(rconv**(-0.5d0))

    do i=1,nlo
       do j=1,mlo(i)
          do k=1,npril(i,j)
             onst0(i,j,k)=onst0(i,j,k)*econv
          end do
       end do
    end do

    do i=1,nparamt
       tpr0(mskt(i),1)=tpr0(mskt(i),1)*econv
       tpr0(mskt(i),2)=tpr0(mskt(i),2)
       tpr0(mskt(i),3)=tpr0(mskt(i),3)*1.d0/rconv
    end do

    do i=1,nparamt
       spr0(mskt(i),2)=spr0(mskt(i),2)
       spr0(mskt(i),3)=spr0(mskt(i),3)*1.d0/rconv
    end do

    !-imprimindo resultados

    nx=0
    do i=1,nparamt
       do j=1,nprllm(i)
          nx=nx+2
       end do
    end do
    do i=1,nlo
       do j=1,mlo(i)
          do k=1,npril(lo(i),mo(i,j))
             nx=nx+1
          end do
       end do
    end do

    if(chkden.eqv..true.)nx=nx+1

    if(nopt.ne.0)then
       WRITE (6,46) FNORMPHM,INFOPHM
       write(6,*)
    end if
    write(6,*)'     Number of parameters:',nx
    write(6,*)'     Number of adjusted parameters:',nparam
    write(6,*)
    do i=1,nlo
       do j=1,mlo(i)
          if(npril(lo(i),mo(i,j)).ge.1)then
             write(6,41)'Density:',dens0,chkden
             goto 4
          end if
       end do
    end do
4   write(6,*)
    write(6,*)'     On-site parameters:'
    write(6,43)('-',j=1,21)
    write(6,44)'l','m','on-site'
    write(6,43)('-',j=1,21)
    do i=1,nlo
       do j=1,mlo(i)
          write(6,47)lo(i)-1,mo(i,j)-lo(i),(onst0(lo(i),mo(i,j),k),k=1,npril(lo(i),mo(i,j)))
       end do
    end do
    write(6,43)('-',j=1,21)
    write(6,*)
    write(6,*)'     Off-site parameters:'
    write(6,43)('-',j=1,41)
    write(6,45)'l','l','m','offsite'
    write(6,43)('-',j=1,41)
    do i=1,nlo
       do j=1,i
          do m1=1,lo(j)
             call labels(lo(j),lo(i),m1,k)
             write(6,48)lo(j)-1,lo(i)-1,m1-1,chk(k),(tpr0(k,l),l=1,nprllm(k))
          end do
       end do
    end do
    write(6,43)('-',j=1,41)
    write(6,*)
    write(6,*)'     Overlap parameters:'
    write(6,43)('-',j=1,41)
    write(6,45)'l','l','m','overlap'
    write(6,43)('-',j=1,41)
    do i=1,nlo
       do j=1,i
          do m1=1,lo(j)
             call labels(lo(j),lo(i),m1,k)
             write(6,48)lo(j)-1,lo(i)-1,m1-1,chk(k),(spr0(k,l),l=1,nprllm(k))
          end do
       end do
    end do
    write(6,43)('-',j=1,41)
    write(6,*)

    return

41  format(5x,a23,f12.4,2x,l)
42  format(5x,i4,2(2x,f8.6))
43  format(5x,43a1)
44  format(5x,2(a2,2x),a8,2x)
45  format(5x,3(a1,2x),3x,a8,2x)
46  FORMAT (5X,31H FINAL L2 NORM OF THE RESIDUALS,D15.7 //   &
         5X,15H EXIT PARAMETER,16X,I10 /)
47  format(5x,2(i2,2x),4(f8.4,2x))
48  format(5x,3(i1,2x),l,2x,10(f8.4,2x))

  end subroutine optimize

  subroutine phm
    !****************************************************************************************
    ! Metodo de minimizacao Powell hibrido fornecido pela biblioteca Minpack                *
    !****************************************************************************************
    implicit none

    integer i,j,k,nx,nmax,MAXFEV,ML,MU,MODE,NPRINT,INFO,LDFJAC,LR,NFEV
    real(8), allocatable, dimension(:) :: WA1,WA2,WA3,WA4
    real(8), allocatable, dimension(:) :: xpr
    real(8), allocatable, dimension(:) :: fvec
    real(8), allocatable, dimension(:) :: diag
    real(8), allocatable, dimension(:) :: rhbrd
    real(8), allocatable, dimension(:) :: qtf
    real(8), allocatable, dimension(:,:) :: fjac
    real(8) XTOL,EPSFCN,FACTOR,ENORM,DPMPAR

    !-alocando arrays

    nmax=nparam

    LR=int(0.5*nmax*(nmax+1))

    allocate(xpr(nmax),fvec(nmax),fjac(nmax,nmax))
    allocate(rhbrd(LR),qtf(nmax),diag(nmax))
    allocate(WA1(nmax),WA2(nmax),WA3(nmax),WA4(nmax))

    XTOL=sqrt(DPMPAR(1))
    MAXFEV = 200*(nmax+1)
    ML = 1
    MU = 1
    EPSFCN = 500.d0
    MODE = 1
    LDFJAC = nmax
    FACTOR = 1.d2
    NPRINT = 1

    !-definindo variaveis para o processo de otimizacao

    nx=1
    if(chkden.eqv..true.)then
       xpr(nx)=dens
       nx=nx+1
    end if

    do i=1,nlopar
       do j=1,mmopar(i)
          do k=1,npril(mlop(i),mmop(i,j))
             xpr(nx)=onst(mlop(i),mmop(i,j),k)
             nx=nx+1
          end do
       end do
    end do

    do i=1,nskpar
       do j=1,nprllm(msk(i))
          xpr(nx)=tpr1(msk(i),j)
          xpr(nx+1)=spr1(msk(i),j)
          nx=nx+2
       end do
    end do

    call hybrd(fcn_PHM,nmax,xpr,fvec,XTOL,MAXFEV,ML,MU,EPSFCN,diag,MODE, &
         FACTOR,NPRINT,INFO,NFEV,fjac,LDFJAC,rhbrd,LR,qtf,WA1,WA2,WA3,WA4)

    !-redefinindo parametros

    nx=1
    if(chkden.eqv..true.)then
       dens=xpr(nx)
       nx=nx+1
    end if

    do i=1,nlopar
       do j=1,mmopar(i)
          do k=1,npril(mlop(i),mmop(i,j))
             onst(mlop(i),mmop(i,j),k)=xpr(nx)
             nx=nx+1
          end do
       end do
    end do

    do i=1,nskpar
       do j=1,nprllm(msk(i))
          tpr1(msk(i),j)=xpr(nx)
          spr1(msk(i),j)=xpr(nx+1)
          nx=nx+2
       end do
    end do

    FNORMPHM=ENORM(nmax,FVEC)

    return

  end subroutine phm

  subroutine fcn_PHM(nmax,xpr,fvec,IFLAG)
    !****************************************************************************************
    ! Rotina requerida pelo metodo Powell hibrido                                           *
    !****************************************************************************************
    implicit none

    integer IFLAG,i,k,l,j,nmax,nx,ll,lx
    real(8) xpr(nmax),fvec(nmax),kptt(3)

    if(IFLAG.ne.0)goto 1

    goto 2

    !-redefinindo parametros

1   nx=1
    if(chkden.eqv..true.)then
       dens=xpr(1)
       nx=nx+1
    end if

    do i=1,nlopar
       do j=1,mmopar(i)
          do k=1,npril(mlop(i),mmop(i,j))
             onst(mlop(i),mmop(i,j),k)=xpr(nx)
             nx=nx+1
          end do
       end do
    end do

    do i=1,nskpar
       do j=1,nprllm(msk(i))
          tpr1(msk(i),j)=xpr(nx)
          spr1(msk(i),j)=xpr(nx+1)
          nx=nx+2
       end do
    end do

    !-checando sistema de equações

   nx=0
    do l=1,nkrlx
       do ll=1,dkrlx(l)
          nx=nx+1
       end do
    end do

    if(nx.ne.nmax)stop 'fcn_PHM: Qty of equations and adjusted variables does not match!'

    !-calculando elemento de matriz H e S

    lx=1
    do l=1,nkrlx
       do ll=1,dkrlx(l)
          do k=1,3
             kptt(k)=kptl(l,ll,k)
          end do

          call bands(kptt)

          fvec(lx)=bndfit(kz(l,ll),nrlx(l))-bnd(ntblx(l))

          lx=lx+1

       end do
    end do

2   return

  end subroutine fcn_PHM

end module optimizing
