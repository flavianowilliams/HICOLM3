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

  integer verlchk

  integer, allocatable :: ilist(:,:)!,ilista(:,:)
  integer, allocatable :: nlist(:)!,nlista(:)
!  integer, allocatable :: trsffi(:)

  save verlchk

contains

  subroutine neighbour_prepare
    !****************************************************************************************
    ! Alocacao de variaveis globais                                                         *
    !****************************************************************************************
    implicit none

    integer numb,ierr,ierra

    !-valores iniciais

    verlchk=1
    ierr=0
    ierra=0

    !-alocando arrays

!    if(nvdwstp.ne.0.or.ncoulstp.ne.0)then
       numb=max(1,int(0.125*max(ncoulstp,nvdwstp)))
       allocate(ilist(natom,numb),nlist(natom),stat=ierr)
!    end if

    !-alocando arrays (vizinhos intramoleculares)

!    if(ntrsff.ne.0)then
!       allocate(ilista(natom,natom),nlista(natom),stat=ierra)
!    end if

    if(ierr.ne.0)stop 'neighbour_prepare: allocation failed'
    if(ierra.ne.0)stop 'neighbour_prepare: allocation failed'

    return

  end subroutine neighbour_prepare

!  subroutine verlet_list_all
!    !****************************************************************************************
!    ! Atualizacao a lista de vizinhos segundo o método de Verlet                            *
!    !****************************************************************************************
!
!    implicit none
!
!    integer i,j,nx
!    real(8) dr,xvz,yvz,zvz,drmax,drx,dry,drz
!
!    !-atualizando lista de vizinhos
!
!    do i=1,natom
!       nx=1
!       do j=1,natom
!          if(i.eq.j)cycle
!          call mic(i,j,xvz,yvz,zvz)
!          dr=sqrt(xvz**2+yvz**2+zvz**2)
!          if(dr.le.rcutoff)then
!             ilista(i,nx)=j
!             nx=nx+1
!          end if
!          nlista(i)=nx-1
!       end do
!    end do
!
!    !-atualizando intervalo
!
!    drmax=0.d0
!    do i=1,natom
!       drx=vax(i)*dtime+0.5d0*fax(i)*dtime**2/mass(i)
!       dry=vay(i)*dtime+0.5d0*fay(i)*dtime**2/mass(i)
!       drz=vaz(i)*dtime+0.5d0*faz(i)*dtime**2/mass(i)
!       dr=sqrt(drx**2+dry**2+drz**2)
!       drmax=max(drmax,dr)
!    end do
!
!    if(drmax.ne.0.d0)verlchk=min(verlchk,max(1,int(drcutoff/drmax)))
!
!    return
!
!  end subroutine verlet_list_all
!
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
