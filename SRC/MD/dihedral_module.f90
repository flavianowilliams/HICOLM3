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

contains

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

    ix(1)=i1 !-i
    ix(2)=i2 !-j
    ix(3)=i3 !-k
    ix(4)=i4 !-n

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

    virtors=virtors-fbi(1)*xa(i1)
    virtors=virtors-fbi(2)*ya(i1)
    virtors=virtors-fbi(3)*za(i1)

    virtors=virtors-fbj(1)*xa(i2)
    virtors=virtors-fbj(2)*ya(i2)
    virtors=virtors-fbj(3)*za(i2)

    virtors=virtors-fbk(1)*xa(i3)
    virtors=virtors-fbk(2)*ya(i3)
    virtors=virtors-fbk(3)*za(i3)

    virtors=virtors-fbn(1)*xa(i4)
    virtors=virtors-fbn(2)*ya(i4)
    virtors=virtors-fbn(3)*za(i4)

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
    ! d[(rij x rjk)**2]                                                                     *
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
