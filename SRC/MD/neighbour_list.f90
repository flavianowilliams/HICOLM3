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
module neighbour_list
  !*****************************************************************************************
  !Atualizacao da lista de vizinhos segundo o método de Verlet                             *
  !                                                                                        *
  ! Flaviano Williams Fernandes, 05 de dezembro de 2018                                    *
  !*****************************************************************************************

  use input
  use estrutura
  use vdw_module
  use coulomb_module

  integer verlchk,ntrsffstp

  integer, allocatable :: ilist(:,:),ilista(:,:)
  integer, allocatable :: nlist(:),nlista(:)
  integer, allocatable :: trsffi(:)

  save ntrsffstp,verlchk

contains

  subroutine neighbour_prepare
    !****************************************************************************************
    ! Alocacao de variaveis globais                                                         *
    !****************************************************************************************
    implicit none

    integer numb,ierr

    !-valores iniciais

    verlchk=1

    !-alocando arrays (vizinhos intermoleculares)

    if(nvdw.ne.0.or.ncoul.ne.0)then
       numb=int(0.125*max(ncoulstp,nvdwstp))
       allocate(ilist(natom,numb),nlist(natom),stat=ierr)
    end if

    !-alocando arrays (vizinhos intramoleculares)

    if(ntrsff.ne.0)then
       allocate(ilista(natom,natom),nlista(natom))
    end if

    if(ierr.ne.0)stop 'neighbour_prepare: allocation failed'

    return

  end subroutine neighbour_prepare

  subroutine verlet_list_all
    !****************************************************************************************
    ! Atualizacao a lista de vizinhos segundo o método de Verlet                            *
    !****************************************************************************************

    implicit none

    integer i,j,nx
    real(8) dr,xvz,yvz,zvz,drmax,drx,dry,drz

    !-atualizando lista de vizinhos

    do i=1,natom
       nx=1
       do j=1,natom
          if(i.eq.j)cycle
          call mic(i,j,xvz,yvz,zvz)
          dr=sqrt(xvz**2+yvz**2+zvz**2)
          if(dr.le.rcutoff)then
             ilista(i,nx)=j
             nx=nx+1
          end if
          nlista(i)=nx-1
       end do
    end do

    !-atualizando intervalo

    drmax=0.d0
    do i=1,natom
       drx=vax(i)*dtime+0.5d0*fax(i)*dtime**2/mass(i)
       dry=vay(i)*dtime+0.5d0*fay(i)*dtime**2/mass(i)
       drz=vaz(i)*dtime+0.5d0*faz(i)*dtime**2/mass(i)
       dr=sqrt(drx**2+dry**2+drz**2)
       drmax=max(drmax,dr)
    end do

    if(drmax.ne.0.d0)verlchk=min(verlchk,max(1,int(drcutoff/drmax)))

    return

  end subroutine verlet_list_all

  subroutine verlet_list_inter
    !****************************************************************************************
    ! Atualizacao da lista de vizinhos segundo o método de Verlet                           *
    !****************************************************************************************
    implicit none

    integer i,j,ii,jj,nx,ni,nj,nii
    real(8) dr,xvz,yvz,zvz,drmax,drx,dry,drz

    !-atualizando lista de vizinhos

    nii=0
    do i=1,moltot
       do ii=1,nzmolec(i)
          ni=nii+ii
          nj=nii+nzmolec(i)+1
          nx=1
          do j=i+1,moltot
             do jj=1,nzmolec(j)
                call mic(ni,nj,xvz,yvz,zvz)
                dr=sqrt(xvz**2+yvz**2+zvz**2)
                if(dr.le.rcutoff)then
                   ilist(ni,nx)=nj
                   nx=nx+1
                end if
                nj=nj+1
             end do
          end do
          nlist(ni)=nx-1
       end do
       nii=nii+nzmolec(i)
    end do

    !-atualizando intervalo

    drmax=0.d0
    do i=1,natom
       drx=vax(i)*dtime+0.5d0*fax(i)*dtime**2/mass(i)
       dry=vay(i)*dtime+0.5d0*fay(i)*dtime**2/mass(i)
       drz=vaz(i)*dtime+0.5d0*faz(i)*dtime**2/mass(i)
       dr=sqrt(drx**2+dry**2+drz**2)
       drmax=max(drmax,dr)
    end do

    if(drmax.ne.0.d0)verlchk=min(verlchk,max(1,int(drcutoff/drmax)))

    return

  end subroutine verlet_list_inter

end module neighbour_list
