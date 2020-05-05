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
module brillouin
  !*******************************************************************************************
  ! Propriedades referentes ao espaço reciproco:                                             *
  !-Distribuicao de pontos K                                                                 *
  !-Vetores da rede reciproca                                                                *
  !                                                                                          *
  ! Flaviano Williams Fernandes, 05 de dezembro de 2018                                      *
  !*******************************************************************************************

  use sistema
  use input
  use estrutura

  real(8) vr(3,3),vrn(3),vzb(3,3),vzbn(3)
  real(8) kxmax,kymax,kzmax,kxmin,kymin,kzmin

  save vr,vrn,kxmax,kymax,kzmax,kxmin,kymin,kzmin

contains

  subroutine reciprocal
    !****************************************************************************************
    ! Vetores da rede reciproca                                                             *
    !****************************************************************************************

    implicit none

    integer i,j
    real(8) s

    !-produto vetorial ai x aj

    vr(1,1)=(v(2,2)*v(3,3)-v(2,3)*v(3,2))
    vr(1,2)=(v(2,3)*v(3,1)-v(2,1)*v(3,3))
    vr(1,3)=(v(2,1)*v(3,2)-v(2,2)*v(3,1))
    vr(2,1)=(v(3,2)*v(1,3)-v(3,3)*v(1,2))
    vr(2,2)=(v(3,3)*v(1,1)-v(3,1)*v(1,3))
    vr(2,3)=(v(3,1)*v(1,2)-v(3,2)*v(1,1))
    vr(3,1)=(v(1,2)*v(2,3)-v(1,3)*v(2,2))
    vr(3,2)=(v(1,3)*v(2,1)-v(1,1)*v(2,3))
    vr(3,3)=(v(1,1)*v(2,2)-v(1,2)*v(2,1))

    !-volume da celula reciproca a1*(a2 x a3)/2pi

    s=0.d0
    do j=1,3
       s=s+v(1,j)*vr(1,j)
    end do

    s=s/(2.d0*pi)

    !-vetores da rede reciproca

    do i=1,3
       vrn(i)=0.d0
       do j=1,3
          vr(i,j)=vr(i,j)/s
          vrn(i)=vrn(i)+vr(i,j)**2
       end do
       vrn(i)=sqrt(vrn(i))
    end do

    !-limites da zona de Brillouin

    kxmin=min(vr(1,1),min(vr(2,1),vr(3,1)))
    kymin=min(vr(1,2),min(vr(2,2),vr(3,2)))
    kzmin=min(vr(1,3),min(vr(2,3),vr(3,3)))

    kxmax=max(vr(1,1),max(vr(2,1),vr(3,1)))
    kymax=max(vr(1,2),max(vr(2,2),vr(3,2)))
    kzmax=max(vr(1,3),max(vr(2,3),vr(3,3)))

    !-imprimindo informaçoes do espaco reciproco

    write(6,*)
    write(6,*)'Reciprocal space:'
    write(6,*)
    write(6,'(a16,3f15.8)')'Lattice constts:',(vrn(i)*kconv,i=1,3)
    write(6,*)
    write(6,'(a16,3f15.8)')'Lattice vectors:',(vr(1,i)*kconv,i=1,3)
    write(6,'(16x,3f15.8)')(vr(2,i)*kconv,i=1,3)
    write(6,'(16x,3f15.8)')(vr(3,i)*kconv,i=1,3)
    write(6,*)

    return

  end subroutine reciprocal

  subroutine high_symmetry(oprt,lsym,nvec,n)
    !****************************************************************************************
    ! Calculo dos pontos K na linha de simetria                                             *
    !                                                                                       *
    !-bandopt = 1 -> pi/a                                                                   *
    !-bandopt = 2 -> vetores internos                                                       *
    !****************************************************************************************

    implicit none

    integer i,j,l,k,n,sdband,oprt
    real(8) lsym(nkmax,3),vec(nbandmax,3)
    real(8) s,snvec,nvec(nbandmax),bandn(nbandmax,3)

    !-linhas de simetria na zona de brillouin

    select case(bandopt)
    case(1)
       do i=1,nband
          bandn(i,1)=band(i,1)*1.d0*pi/(1.0d0/kconv)
          bandn(i,2)=band(i,2)*1.d0*pi/(1.0d0/kconv)
          bandn(i,3)=band(i,3)*1.d0*pi/(1.0d0/kconv)
       end do
    case(2)
       do k=1,3
          do i=1,nband
             bandn(i,k)=0.d0
             do l=1,3
                bandn(i,k)=bandn(i,k)+band(i,l)*vr(l,k)
             end do
          end do
       end do
    end select

    do j=1,3
       vec(1,j)=0.d0
       do i=2,nband
          vec(i,j)=bandn(i,j)-bandn(i-1,j)
       end do
    end do

    lsym(1,1)=bandn(1,1)
    lsym(1,2)=bandn(1,2)
    lsym(1,3)=bandn(1,3)

    n=2
    do i=2,nband
       do j=1,dband(i)
          do k=1,3
             lsym(n,k)=bandn(i-1,k)+j*vec(i,k)/dband(i)
          end do
          n=n+1
       end do
    end do

    n=n-1

    snvec=0.d0
    sdband=dband(1)
    do i=2,nband
       s=0.d0
       do j=1,3
          s=s+vec(i,j)**2
       end do
       nvec(i)=sqrt(s)
       snvec=snvec+nvec(i)
       sdband=sdband+dband(i)
    end do

    nvec(1)=0.d0

    if(oprt.eq.0)goto 1

    !-imprimindo informacoes do espaço reciproco

    write(6,*)
    write(6,'(a29)')'High-symmetry Brillouin zone!'
    write(6,*)
    write(6,'(a19,i15)')'Number of K-points:',sdband
    write(6,*)
    write(6,'(a25,3f15.8)')'                  length:',snvec
    write(6,'(a25,3f15.8)')'                 Vectors:',(vec(2,j),j=1,3)
    do i=3,nband
       write(6,'(a25,3f15.8)')'                         ',(vec(i,j),j=1,3)
    end do
    write(6,*)

1   return

  end subroutine high_symmetry

end module brillouin
