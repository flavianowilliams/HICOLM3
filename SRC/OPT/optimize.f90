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
module optimize

  use input
  use estrutura
  use force_field

contains

  subroutine opt

    implicit none

    integer i,j
    real(8) dfmax,gax(natom),gay(natom),gaz(natom)
    real(8) encoul,enbond,enbend,entors,envdw,enpot,enpot0
    real(8) virvdw,virbond,virbend,virtors,vircoul

    !-criando ficheiro de saida

    open(3,file='HICOLM_opt.dat',status='unknown')

    !-valores iniciais

    enpot0=0.d0

    do i=1,natom
       gax(i)=0.d0
       gay(i)=0.d0
       gaz(i)=0.d0
    end do

    !-stress

    do i=1,6
       str(i)=0.d0
    end do

    !-preparando Campo de Força

    call ff_prepare

    !-calculando contribuição intramolecular

    write(6,'(93a1)')('#',i=1,93)
    write(6,*)(' OPT',i=1,23)
    write(6,'(93a1)')('#',i=1,93)
    write(6,*)

    do i=1,1000
       encoul=0.d0   !coulombiano
       enbond=0.d0   !estiramento
       enbend=0.d0   !deformacao
       entors=0.d0   !torção
       envdw=0.d0    !Van der waals
       do j=1,natom
          fax(j)=0.d0
          fay(j)=0.d0
          faz(j)=0.d0
       end do
       call ff_modules_intra&
            (enbond,enbend,entors,envdw,encoul,virbond,virbend,virtors,virvdw,vircoul)
       enpot=enbond+enbend+entors+envdw+encoul
       call steepest_descent
       call opt_check(gax,gay,gaz,dfmax)
       call geometria
       if(i.ge.2)write(3,10)i,dfmax*econv/rconv,enpot*econv,(enpot-enpot0)*econv
       do j=1,natom
          gax(j)=fax(j)
          gay(j)=fay(j)
          gaz(j)=faz(j)
       end do
       enpot0=enpot
    end do

    return

10  format(5x,i8,2x,es12.4,2x,es12.4,3(2x,es12.4))

  end subroutine opt

  subroutine steepest_descent

    implicit none

    integer i
    real(8) lambda,df

    lambda=1.e-5

    do i=1,natom
       df=sqrt(fax(i)**2+fay(i)**2+faz(i)**2)
       xa(i)=xa(i)+lambda*fax(i)/df
       ya(i)=ya(i)+lambda*fay(i)/df
       za(i)=za(i)+lambda*faz(i)/df
    end do

    return

  end subroutine steepest_descent

  subroutine opt_check(gax,gay,gaz,dfmax)

    implicit none

    integer i
    real(8) dfmax,gax(natom),gay(natom),gaz(natom),dfx,dfy,dfz

    dfmax=0.d0
    do i=1,natom
       dfx=fax(i)-gax(i)
       dfy=fay(i)-gay(i)
       dfz=faz(i)-gaz(i)
       dfmax=max(dfmax,sqrt(dfx**2+dfy**2+dfz**2))
    end do

    return

  end subroutine opt_check

end module optimize
