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
module fitting
  !******************************************************************************************
  ! Bandas de energia do calculo DFT com o uso da plataforma SIESTA                         *
  !                                                                                         *
  ! Flaviano Williams Fernandes, 05 de dezembro de 2018                                     *
  !******************************************************************************************

  use sistema
  use input

  implicit none

  integer ihomo,ilumo,nfit,nkptdft

  real(8), allocatable :: bndfit(:,:)
  real(8), allocatable :: kptfit(:)

  save ihomo,ilumo

contains

  subroutine DFT
    !***************************************************************************************
    ! Valores de entrada no arquivo bands da plataforma SIESTA                             *
    !***************************************************************************************
    implicit none

    integer i,j,k,iband,bandval,bandcnd
    real(8) bottom,top,homo,khomo,gap,lumo,klumo,klen

    data bandval/40/,bandcnd/5/

    open(8,file='TB.dat',status='unknown')
    open(9,file='TB.bands',status='unknown')

    read(9,*)fermi
    read(9,*)
    read(9,*)bottom,top
    read(9,*)nfit,iband,nkptdft

    nfit=nfit*iband

    !-alocando arrays

    allocate(kptfit(nkptdft),bndfit(nkptdft,nfit))

    homo=bottom
    do k=1,nkptdft
       read(9,*)kptfit(k),(bndfit(k,i),i=1,nfit)
       do i=1,nfit
          if(bndfit(k,i).gt.homo.and.bndfit(k,i).le.fermi)then
             khomo=kptfit(k)
             homo=bndfit(k,i)
             ihomo=i
          end if
       end do
    end do

    lumo=top
    do k=1,nkptdft
       do i=1,nfit
          if(bndfit(k,i).lt.lumo.and.bndfit(k,i).ge.fermi)then
             klumo=kptfit(k)
             lumo=bndfit(k,i)
             ilumo=i
          end if
       end do
    end do

    if(homo.lt.fermi)then
       gap=lumo-homo
    else
       gap=0.d0
    end if

    do k=1,nkptdft
       do j=1,nfit
          bndfit(k,j)=bndfit(k,j)/econv
       end do
    end do

    fermi=fermi/econv
    klen=kptfit(nkptdft)

    !-imprimindo informacoes do calculo DFT

    bandval=min(bandval,ihomo)
    bandcnd=min(bandcnd,nfit)

    write(6,*)('#',i=1,70)
    write(6,*)'$$$$ ',('DFT ',i=1,15),'$$$$$'
    write(6,*)('#',i=1,70)
    write(6,*)
    write(6,49)'                (Bands,K):',nfit,klen
    write(6,46)'              Fermi level:',fermi*econv
    write(6,*)
    write(6,45)'            Valence band:',ihomo
    write(6,46)'                (E(k),k):',homo,khomo
    write(6,*)
    write(6,45)'         Conduction band:',ilumo
    write(6,46)'                (E(k),k):',lumo,klumo
    write(6,*)
    write(6,46)'                     GAP:',gap
    write(6,*)

    ! -imprimindo bandas de energia do calculo DFT

    do i=1,nkptdft
       write(8,'(100f8.3)')kptfit(i),&
            ((bndfit(i,j)-fermi)*econv,j=(ihomo-bandval+1),(ihomo+bandcnd))
    end do

    return

45  format(a25,i10)
46  format(a25,2f10.4)
49  format(a25,i10,f10.4)

  end subroutine DFT

  end module fitting
