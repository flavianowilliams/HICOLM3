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
module dihedral_proper_module
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

  integer, allocatable :: torsijkn(:,:)
  integer, allocatable :: torsim(:),torsib(:)

  integer ntorsstp

  save torsim,torsib,torsijkn,ntorsstp

contains

  subroutine tors_alloc
    !****************************************************************************************
    ! Alocacao de variaveis globais                                                         *
    !****************************************************************************************

    implicit none

    integer ierr

    allocate(torsijkn(4,ntors),torsim(ntors),torsib(ntors),stat=ierr)

    return

  end subroutine tors_alloc

  subroutine tors_convert
    !****************************************************************************************
    ! Conversao de unidades de medida:                                                      *
    ! Unidades de entrada ---> a.u.                                                         *
    !****************************************************************************************

    implicit none

    integer i,j

    !-convertendo unidades de medida

    do i=1,nmolec
       do j=1,torscnt(i)
          select case(tors(i,j))
          case(1)
             partors(i,j,1)=partors(i,j,1)/econv
             partors(i,j,2)=partors(i,j,2)/aconv
          case(2)
             partors(i,j,1)=partors(i,j,1)/econv
          case(3)
             partors(i,j,1)=partors(i,j,1)/econv
             partors(i,j,2)=partors(i,j,2)/econv
             partors(i,j,3)=partors(i,j,3)/econv
             partors(i,j,4)=partors(i,j,4)/econv
             partors(i,j,5)=partors(i,j,5)/econv
             partors(i,j,6)=partors(i,j,6)/econv
             partors(i,j,7)=partors(i,j,7)/aconv
          case(4)
             partors(i,j,2)=partors(i,j,2)/econv
             partors(i,j,3)=partors(i,j,3)/aconv
          end select
       end do
    end do

    return

  end subroutine tors_convert

  subroutine tors_counts
    !****************************************************************************************
    ! Contagem do número de diedros                                                         *
    !****************************************************************************************

    implicit none

    integer, allocatable :: chk(:,:)

    integer i,j,k,nkb,nx,nxx,ni,nj,nk,nn,np

    !-alocando array

    allocate(chk(nmolec,ntors))

    !-checando viabilidade de bonds

    nkb=1

    do i=1,nmolec
       do j=1,torscnt(i)
          select case(tors(i,j))
          case(1)
             nkb=1
          case(2)
             nkb=1
          case(3)
             nkb=0
          case(4)
             nkb=2
          end select
          chk(i,j)=1
          do k=1,nkb
             if(partors(i,j,k).eq.0.d0)chk(i,j)=0
          end do
       end do
    end do

    !-calculando torsoes intramoleculares

    np=0
    nx=1
    do i=1,nmolec
       nxx=0
       do j=1,ntmolec(i)
          do k=1,torscnt(i)
             ni=np+moltors(i,k,1)
             nj=np+moltors(i,k,2)
             nk=np+moltors(i,k,3)
             nn=np+moltors(i,k,4)
             if(chk(i,k).eq.1)then
                torsijkn(1,nx)=ni
                torsijkn(2,nx)=nj
                torsijkn(3,nx)=nk
                torsijkn(4,nx)=nn
                torsim(nx)=i
                torsib(nx)=k
                nxx=nxx+1
                nx=nx+1
             end if
          end do
          np=np+nxmolec(i)
       end do
       torsmlc(i)=nxx
    end do

    ntorsstp=nx-1

    !-limpando memoria

    deallocate(chk)

    return

  end subroutine tors_counts

  subroutine tors_calc(entors,virtors)
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

    do i=1,ntorsstp

       ia=torsijkn(1,i)
       ib=torsijkn(2,i)
       ic=torsijkn(3,i)
       id=torsijkn(4,i)
       im=torsim(i)
       in=torsib(i)

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

       call tors_flags(im,in,phi,pot,fd)
       call tors_force(ia,ib,ic,id,drij,drjk,drkn,dvc1,dvc2,phi,fd,virtors)

       entors=entors+pot

    end do

    return

  end subroutine tors_calc

  subroutine tors_flags(im,in,phi,pot,fd)
    !****************************************************************************************
    ! Compontente angular                                                                   *
    !****************************************************************************************

    implicit none

    integer im,in
    real(8) phi,fd,pot,psi,p1,p2,p3,p4

    !-valores iniciais

    pot=0.d0
    fd=0.d0

    !-componente dU/dphi

    select case(tors(im,in))
    case(1)
       pot=0.5d0*partors(im,in,1)*(phi-partors(im,in,2))**2
       fd=partors(im,in,1)*(phi-partors(im,in,2))
    case(2)
       pot=0.5d0*partors(im,in,1)*(cos(phi)-partors(im,in,2))**2
       fd=-sin(phi)*partors(im,in,1)*(cos(phi)-partors(im,in,2))
    case(3)
       psi=phi-partors(im,in,7)
       pot=partors(im,in,1)+partors(im,in,2)*cos(psi)+partors(im,in,3)*cos(psi)**2&
            +partors(im,in,4)*cos(psi)**3+partors(im,in,5)*cos(psi)**4&
            +partors(im,in,6)*cos(psi)**5
       fd=-sin(psi)*(partors(im,in,2)+2.d0*partors(im,in,3)*cos(psi)&
            +3.d0*partors(im,in,4)*cos(psi)**2+4.d0*partors(im,in,5)*cos(psi)**3&
            +5.d0*partors(im,in,6)*cos(psi)**4)
    case(4)
       p1=partors(im,in,1)
       p2=partors(im,in,2)
       p3=partors(im,in,3)
       p4=partors(im,in,4)
       pot=0.5d0*p2*(1.d0+cos(p4*phi-p3))/p1
       fd=-0.5d0*p4*p2*sin(p4*phi-p3)/p1
    end select

    return

  end subroutine tors_flags

end module dihedral_proper_module
