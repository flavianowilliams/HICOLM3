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
module dihedral_improper_module
  !******************************************************************************************
  ! Contribuicao de diedros para o campo de força:                                          *
  ! - Energia potencial;                                                                    *
  ! - Forças atômicas;                                                                      *
  ! - Virial;                                                                               *
  ! - Stress.                                                                               *
  !                                                                                         *
  ! Flaviano Williams Fernandes, 05 de dezembro de 2018                                     *
  !******************************************************************************************

  use dihedral_module

  integer, allocatable :: itorsijkn(:,:)
  integer, allocatable :: itorsim(:),itorsib(:)

  integer nitorsstp

  save itorsim,itorsib,itorsijkn,nitorsstp

contains

  subroutine itors_alloc
    !****************************************************************************************
    ! Alocacao de variaveis globais                                                         *
    !****************************************************************************************

    implicit none

    integer ierr

    allocate(itorsijkn(4,nitors),itorsim(nitors),itorsib(nitors),stat=ierr)

    return

  end subroutine itors_alloc

  subroutine itors_convert
    !****************************************************************************************
    ! Conversao de unidades de medida:                                                      *
    ! Unidades de entrada ---> a.u.                                                         *
    !****************************************************************************************

    implicit none

    integer i,j

    !-convertendo unidades de medida

    do i=1,nmolec
       do j=1,itorscnt(i)
          select case(itors(i,j))
          case(1)
             paritors(i,j,1)=paritors(i,j,1)/econv
             paritors(i,j,2)=paritors(i,j,2)/aconv
          end select
       end do
    end do

    return

  end subroutine itors_convert

  subroutine itors_counts
    !****************************************************************************************
    ! Contagem do número de diedros                                                         *
    !****************************************************************************************

    implicit none

    integer, allocatable :: chk(:,:)

    integer i,j,k,nkb,nx,nxx,ni,nj,nk,nn,np

    !-alocando array

    allocate(chk(nmolec,nitors))

    !-checando viabilidade de bonds

    nkb=1

    do i=1,nmolec
       do j=1,itorscnt(i)
          select case(itors(i,j))
          case(1)
             nkb=1
          end select
          chk(i,j)=1
          do k=1,nkb
             if(paritors(i,j,k).eq.0.d0)chk(i,j)=0
          end do
       end do
    end do

    !-calculando torsoes intramoleculares

    np=0
    nx=1
    do i=1,nmolec
       nxx=0
       do j=1,ntmolec(i)
          do k=1,itorscnt(i)
             ni=np+molitors(i,k,1)
             nj=np+molitors(i,k,2)
             nk=np+molitors(i,k,3)
             nn=np+molitors(i,k,4)
             if(chk(i,k).eq.1)then
                itorsijkn(1,nx)=ni
                itorsijkn(2,nx)=nj
                itorsijkn(3,nx)=nk
                itorsijkn(4,nx)=nn
                itorsim(nx)=i
                itorsib(nx)=k
                nxx=nxx+1
                nx=nx+1
             end if
          end do
          np=np+nxmolec(i)
       end do
       itorsmlc(i)=nxx
    end do

    nitorsstp=nx-1

    !-limpando memoria

    deallocate(chk)

    return

  end subroutine itors_counts

  subroutine itors_calc(entors,virtors)
    !****************************************************************************************
    ! - Energia potencial;                                                                  *
    ! - Contribuicao para o virial;                                                         *
    ! - Contribuicao para o stress;                                                         *
    ! - Forças atômicas.                                                                    *
    !****************************************************************************************

    implicit none

    integer i,ia,ib,ic,id,im,in
    real(8) pot,fd,drij(3),drjk(3),drkn(3),vc1x,vc1y,vc1z,vc2x,vc2y,vc2z,dvc1,dvc2,phi
    real(8) entors,virtors

    do i=1,nitorsstp

       ia=itorsijkn(1,i)
       ib=itorsijkn(2,i)
       ic=itorsijkn(3,i)
       id=itorsijkn(4,i)
       im=itorsim(i)
       in=itorsib(i)

       call mic(ia,ib,drij(1),drij(2),drij(3))
       call mic(ib,ic,drjk(1),drjk(2),drjk(3))
       call mic(ic,id,drkn(1),drkn(2),drkn(3))

       !-produto vetorial

       vc1x=drij(2)*drjk(3)-drij(3)*drjk(2)
       vc1y=drij(3)*drjk(1)-drij(1)*drjk(3)
       vc1z=drij(1)*drjk(2)-drij(2)*drjk(1)

       vc2x=drjk(2)*drkn(3)-drjk(3)*drkn(2)
       vc2y=drjk(3)*drkn(1)-drjk(1)*drkn(3)
       vc2z=drjk(1)*drkn(2)-drjk(2)*drkn(1)

       dvc1=sqrt(vc1x**2+vc1y**2+vc1z**2) !-|rij x rjk|
       dvc2=sqrt(vc2x**2+vc2y**2+vc2z**2) !-|rjk x rkn|

       !-angulo do diedro

       phi=acos((vc1x*vc2x+vc1y*vc2y+vc1z*vc2z)/(dvc1*dvc2))

       phi=max(1.d-8,phi)

       call itors_flags(im,in,phi,pot,fd)
       call tors_force(ia,ib,ic,id,drij,drjk,drkn,dvc1,dvc2,phi,fd,virtors)

       entors=entors+pot

    end do

    return

  end subroutine itors_calc

  subroutine itors_flags(im,in,phi,pot,fd)
    !****************************************************************************************
    ! Compontente angular                                                                   *
    !****************************************************************************************

    implicit none

    integer im,in
    real(8) phi,fd,pot,p1,p2,p3

    !-valores iniciais

    pot=0.d0
    fd=0.d0

    !-componente dU/dphi

    select case(itors(im,in))
    case(1)
       p1=paritors(im,in,1)
       p2=paritors(im,in,2)
       p3=paritors(im,in,3)
       pot=0.5d0*p1*(1.d0+cos(p3*phi-p2))
       fd=-0.5d0*p3*p1*sin(p3*phi-p2)
    end select

    return

  end subroutine itors_flags

end module dihedral_improper_module
