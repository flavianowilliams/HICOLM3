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
module dihedral_module
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
             partors(i,j,1)=partors(i,j,1)/econv
             partors(i,j,2)=partors(i,j,2)/aconv
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
             nkb=1
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

       call ccpmm(ia,ib,drij(1),drij(2),drij(3))
       call ccpmm(ib,ic,drjk(1),drjk(2),drjk(3))
       call ccpmm(ic,id,drkn(1),drkn(2),drkn(3))

       !-produto vetorial

       vc1x=drij(2)*drjk(3)-drij(3)*drjk(2)
       vc1y=drij(3)*drjk(1)-drij(1)*drjk(3)
       vc1z=drij(1)*drjk(2)-drij(2)*drjk(1)

       vc2x=drjk(2)*drkn(3)-drjk(3)*drkn(2)
       vc2y=drjk(3)*drkn(1)-drjk(1)*drkn(3)
       vc2z=drjk(1)*drkn(2)-drjk(2)*drkn(1)

       dvc1=sqrt(vc1x**2+vc1y**2+vc1z**2)
       dvc2=sqrt(vc2x**2+vc2y**2+vc2z**2)

       !-angulo do diedro

       phi=acos((vc1x*vc2x+vc1y*vc2y+vc1z*vc2z)/(dvc1*dvc2))
       phi=max(1.d-4,phi)

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
    real(8) phi,fd,pot,psi

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
       pot=0.5d0*partors(im,in,1)*(phi-partors(im,in,2))**2
       fd=partors(im,in,1)*(phi-partors(im,in,2))
    end select

    return

  end subroutine tors_flags

  subroutine tors_force(i1,i2,i3,i4,drij,drjk,drkn,dvc1,dvc2,phi,fd,virtors)
    !****************************************************************************************
    ! - Contribuicao para o virial                                                          *
    ! - Contribuicao para o stress                                                          *
    ! - Forcas dos atomos i, j, k e n                                                       *
    !****************************************************************************************

    implicit none

    integer i,j,ix(4),i1,i2,i3,i4
    real(8) drij(3),drjk(3),drkn(3),dvc(4,3),fbi(3),fbj(3),fbk(3),fbn(3)
    real(8) dvc1,dvc2,phi,fd,virtors

    !-componente d[(rij x rjk)*(rjk x rkn)/|rij x rjk||rjk x rkn|]

    ix(1)=i1
    ix(2)=i2
    ix(3)=i3
    ix(4)=i4
    !
    do i=1,4
       do j=1,3
          dvc(i,j)=dfunc1(i1,i2,i3,i4,drij,drjk,drkn,ix(i),j)/(dvc1*dvc2) &
               -0.5d0*cos(phi)*(dfunc2(i1,i2,i3,drij,drjk,ix(i),j)/dvc1**2 &
               +dfunc2(i2,i3,i4,drjk,drkn,ix(i),j)/dvc2**2)
       end do
    end do

    !-contribuicao para as forças atomicas

    fbi(1)=fd*dvc(1,1)/sin(phi)
    fbi(2)=fd*dvc(1,2)/sin(phi)
    fbi(3)=fd*dvc(1,3)/sin(phi)

    fbj(1)=fd*dvc(2,1)/sin(phi)
    fbj(2)=fd*dvc(2,2)/sin(phi)
    fbj(3)=fd*dvc(2,3)/sin(phi)

    fbk(1)=fd*dvc(3,1)/sin(phi)
    fbk(2)=fd*dvc(3,2)/sin(phi)
    fbk(3)=fd*dvc(3,3)/sin(phi)

    fbn(1)=fd*dvc(4,1)/sin(phi)
    fbn(2)=fd*dvc(4,2)/sin(phi)
    fbn(3)=fd*dvc(4,3)/sin(phi)

    fax(i1)=fax(i1)+fbi(1)
    fay(i1)=fay(i1)+fbi(2)
    faz(i1)=faz(i1)+fbi(3)

    fax(i2)=fax(i2)+fbj(1)
    fay(i2)=fay(i2)+fbj(2)
    faz(i2)=faz(i2)+fbj(3)

    fax(i3)=fax(i3)+fbk(1)
    fay(i3)=fay(i3)+fbk(2)
    faz(i3)=faz(i3)+fbk(3)

    fax(i4)=fax(i4)+fbn(1)
    fay(i4)=fay(i4)+fbn(2)
    faz(i4)=faz(i4)+fbn(3)

    !-contribuicao para o virial

    virtors=0.d0

    !-calculo da contribuicao para o stress

    str(1)=str(1)-(fbi(1)*drij(1)+fbj(1)*drjk(1)+fbk(1)*drkn(1))
    str(2)=str(2)-(fbi(2)*drij(2)+fbj(2)*drjk(2)+fbk(2)*drkn(2))
    str(3)=str(3)-(fbi(3)*drij(3)+fbj(3)*drjk(3)+fbk(3)*drkn(3))
    str(4)=str(4)-(fbi(3)*drij(2)+fbj(3)*drjk(2)+fbk(3)*drkn(2))
    str(5)=str(5)-(fbi(1)*drij(3)+fbj(1)*drjk(3)+fbk(1)*drkn(3))
    str(6)=str(6)-(fbi(2)*drij(1)+fbj(2)*drjk(1)+fbk(2)*drkn(1))

    return

  end subroutine tors_force

  double precision function dfunc1(i1,i2,i3,i4,drij,drjk,drkn,i,j)
    !****************************************************************************************
    ! d[(rij x rjk)*(rjk x rkn)]                                                            *
    !****************************************************************************************

    implicit none

    integer i,j,i1,i2,i3,i4
    real(8) drij(3),drjk(3),drkn(3)

    dfunc1=drij(j)*acomt(drjk,drjk,j)*(kronij(i,i3)-kronij(i,i4)) &
         +drij(j)*acomt(drjk,drkn,j)*(kronij(i,i3)-kronij(i,i2)) &
         +drjk(j)*acomt(drij,drjk,j)*(kronij(i,i4)-kronij(i,i3)) &
         +drjk(j)*acomt(drjk,drkn,j)*(kronij(i,i2)-kronij(i,i1)) &
         +drkn(j)*acomt(drij,drjk,j)*(kronij(i,i3)-kronij(i,i2)) &
         +drkn(j)*acomt(drjk,drjk,j)*(kronij(i,i1)-kronij(i,i2)) &
         +2.d0*drjk(j)*acomt(drij,drkn,j)*(kronij(i,i2)-kronij(i,i3))

    return

  end function dfunc1

  double precision function dfunc2(i1,i2,i3,dri1,dri2,i,j)
    !****************************************************************************************
    ! d[(rjk x rkn)**2]                                                                     *
    !****************************************************************************************

    implicit none

    integer i,j,i1,i2,i3
    real(8) dri1(3),dri2(3)

    dfunc2=dri1(j)*acomt(dri2,dri2,j)*(kronij(i,i2)-kronij(i,i1)) &
         +dri1(j)*acomt(dri1,dri2,j)*(kronij(i,i2)-kronij(i,i3)) &
         +dri2(j)*acomt(dri1,dri1,j)*(kronij(i,i3)-kronij(i,i2)) &
         +dri2(j)*acomt(dri1,dri2,j)*(kronij(i,i1)-kronij(i,i2))

    dfunc2=2.d0*dfunc2

    return

  end function dfunc2

  double precision function acomt(x1,x2,k)
    !****************************************************************************************
    ! Anticomutador                                                                         *
    !****************************************************************************************

    implicit none

    integer i,k
    real(8) sum,x1(3),x2(3)

    sum=0.d0
    do i=1,3
       sum=sum+(1.d0-kronij(i,k))*x1(i)*x2(i)
    end do

    acomt=sum

    return

  end function acomt

end module dihedral_module
