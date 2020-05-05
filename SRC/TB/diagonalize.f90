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
module diagonalize
  !*****************************************************************************************
  ! modulo responsavel por resolver o problema de autovalor generalizado                   *
  !*****************************************************************************************
  use input
  use matriz

  implicit none

  integer norbx

  save norbx

contains

  subroutine eigsolve_zhegv()
    !***************************************************************************************
    ! ????????????????????????????                                                         *
    !***************************************************************************************
    implicit none

    integer INFO,LWORK

    real(8), allocatable, dimension(:) :: RWORK

    complex(8), allocatable, dimension(:) :: WORK

    !-alocando arrays

    allocate(WORK(2*norb),RWORK(3*norb-2))

    call zhegv(1,'N','L',norb,mh,norb,ms,norb,bnd,WORK,-1,RWORK,INFO)

    LWORK=min(int(real(WORK(1))),1000)

    call zhegv(1,'N','L',norb,mh,norb,ms,norb,bnd,WORK,LWORK,RWORK,INFO)

    if(INFO.ne.0)then
       write(*,*) 'diagonal: INFO:',INFO,WORK(1)
       stop
    end if

    return

  end subroutine eigsolve_zhegv

  subroutine eigsolve_dgeev()
    !****************************************************************************************
    ! Calculo da equacao matricial                                                          *
    !****************************************************************************************
    implicit none

    real(8), allocatable, dimension(:,:) :: mhr,msr,ml,mc

    norbx=2*norb

    allocate(mhr(norbx,norbx),msr(norbx,norbx))
    allocate(ml(norbx,norbx))
    allocate(mc(norbx,norbx))

    call hermitian(mhr,msr)

    call cholesky(msr,ml)

    call pav(mhr,ml,mc)

    call diagonal(mc)

    return

  end subroutine eigsolve_dgeev

  subroutine hermitian(mhr,msr)
    !****************************************************************************************
    ! Ttransformacao da matriz complexa em real                                             *
    !****************************************************************************************
    implicit none

    integer i,j,ix,jx

    real(8), allocatable, dimension(:,:) :: mhr,msr

    ix=1
    do i=1,norb
       jx=1
       do j=1,norb
          mhr(ix,jx)=real(mh(i,j))
          msr(ix,jx)=real(ms(i,j))
          jx=jx+1
       end do
       do j=1,norb
          mhr(ix,jx)=-aimag(mh(i,j))
          msr(ix,jx)=-aimag(ms(i,j))
          jx=jx+1
       end do
       ix=ix+1
    end do

    do i=1,norb
       jx=1
       do j=1,norb
          mhr(ix,jx)=aimag(mh(i,j))
          msr(ix,jx)=aimag(ms(i,j))
          jx=jx+1
       end do
       do j=1,norb
          mhr(ix,jx)=real(mh(i,j))
          msr(ix,jx)=real(ms(i,j))
          jx=jx+1
       end do
       ix=ix+1
    end do

    return

  end subroutine hermitian

  subroutine cholesky(msr,ml)
    !*****************************************************************************************
    ! Decomposicao L*LT=B usando cholesky                                                    *
    !*****************************************************************************************
    implicit none

    integer i,j,k
    real(8) sum,lx

    real(8), allocatable, dimension(:,:) :: msr,ml

    do i=1,norbx
       do j=1,norbx
          if(abs(msr(i,j)-msr(j,i)).gt.1.d-8)stop 'Cholesky: B is non-symmetric!'
       end do
    end do

    do i=1,norbx
       do j=1,norbx
          sum=msr(i,j)
          do k=1,i-1
             sum=sum-msr(i,k)*ml(j,k)
          end do
          if(i.eq.j)then
             if(sum.le.0.d0)stop 'Cholesky: B is non-positive-definity!'
             lx=sqrt(sum)
             ml(i,i)=msr(i,i)
          else
             ml(j,i)=sum/lx
          end if
       end do
    end do

    return

  end subroutine cholesky

  subroutine pav(mhr,ml,mc)
    !****************************************************************************************
    ! Transformacao do problema generalizado A*x=lambda*B*x na equacao de autovalor         *
    ! (B-1*A)*x=lambda*x                                                                    *
    !****************************************************************************************
    implicit none

    integer i,j,k
    real(8) sum

    real(8), allocatable, dimension(:,:) :: my,ml,mc,mhr

    allocate(my(norbx,norbx))

    do i=1,norbx
       do j=1,norbx
          my(i,j)=mhr(i,j)
          do k=1,j-1
             my(i,j)=my(i,j)-my(i,k)*ml(j,k)
          end do
          my(i,j)=my(i,j)/ml(j,j)
       end do
    end do

    do i=1,norbx
       do j=1,norbx
          sum=my(i,j)
          do k=1,i-1
             sum=sum-ml(i,k)*mc(k,j)
          end do
          mc(i,j)=sum/ml(j,j)
       end do
    end do

    return

  end subroutine pav

  subroutine diagonal(mc)
    !****************************************************************************************
    ! Diagonalizacao da matriz B-1*A-lambda                                                 *
    !****************************************************************************************
    implicit none

    integer LWORK,INFO,i,j,k,chk

    integer, allocatable, dimension(:) :: ix

    real(8) wrmin,wrmax,wx

    real(8), allocatable, dimension(:) :: WI,WR,VL,VR,WORK
    real(8), allocatable, dimension(:,:) :: mc

    LWORK=3*norbx

    allocate(WR(norbx),WI(norbx))
    allocate(VL(norbx),VR(norbx))
    allocate(WORK(LWORK))
    allocate(ix(norbx))

    call dgeev('N','N',norbx,mc,norbx,WR,WI,VL,norbx,VR,norbx,WORK,LWORK,INFO)

    if(INFO.ne.0)then
       write(6,*)'diagonal: INFO:',INFO,WORK(1)
       stop
    end if

    wrmin=0.d0
    do j=1,norbx
       if(WR(j).le.wrmin)then
          wrmin=WR(j)
          ix(1)=j
       end if
    end do

    wrmax=wrmin
    do j=1,norbx
       wrmax=max(wrmax,WR(j))
    end do

    do i=2,norbx
       wx=wrmax
       do j=1,norbx
          if(WR(j).le.wx)then
             chk=1
             do k=1,i-1
                if(j.eq.ix(k))chk=0
             end do
             if(chk.eq.1)then
                ix(i)=j
                wx=WR(j)
             end if
          end if
       end do
    end do

    do i=1,norb
       bnd(i)=WR(ix(2*i))
    end do

    return

    end subroutine diagonal

  subroutine bands(kpt)
    !***************************************************************************************
    ! Calculo dos niveis de energia para cada ponto K                                      *
    !***************************************************************************************
    implicit none

    integer ix,jx,i,j,p,k,m,pp,l1,l2,m1,m2
    real(8) kpt(3)

    ix=1
    do i=1,natom
       do j=1,nlo
          l1=lo(j)
          do p=1,mlo(j)
             m1=mo(j,p)
             jx=1
             do k=1,natom
                do m=1,nlo
                   l2=lo(m)
                   do pp=1,mlo(m)
                      m2=mo(m,pp)
                      call neighbourhood(i,k,l1,l2,m1,m2,kpt,ix,jx)
                      jx=jx+1
                   end do
                end do
             end do
             ix=ix+1
          end do
       end do
    end do

    !-diagonalizando matriz secular

    select case(diagop)
    case(1)
       call eigsolve_zhegv
    case(2)
       call eigsolve_dgeev
    end select

    return

  end subroutine bands

end module diagonalize
