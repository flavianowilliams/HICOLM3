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
    real(8) eintra,einter,enpot0,enpot
    real(8) virvdw,virbond,virbend,virtors,vircoul

    !-criando ficheiro de saida

    open(3,file='HICOLM_opt.dat',status='unknown')

    !-valores iniciais

    vircoul=0.d0   !coulombiano
    virbond=0.d0   !estiramento
    virbend=0.d0   !deformacao
    virtors=0.d0   !torção
    virvdw=0.d0    !Van der waals

    do i=1,natom
       gax(i)=0.d0
       gay(i)=0.d0
       gaz(i)=0.d0
    end do

    !-stress

    do i=1,6
       str(i)=0.d0
    end do

    einter=0.d0
    eintra=0.d0
    enpot0=0.d0

    !-preparando Campo de Força

    call ff_prepare

    !-preparando lista de vizinhos de Verlet

    call verlet_list_inter

    !-correcao da energia de Van der Waals de curto alcance

    call vdw_corr

    !-imprimindo informacoes no ficheiro de saida

    call opt_print

    write(3,5)'#','INTRA','INTER','ENERGY','DENERGY','DFORCE'

    write(6,'(4x,111a1)')('-',i=1,84)
    write(6,10)'##','STEP','INTRA','INTER','ENERGY','DENERGY','DFORCE'
    write(6,'(4x,111a1)')('-',i=1,84)

    !-calculando contribuição intramolecular

    enpot0=0.d0

    do i=1,opt_ntrialmax

       do j=1,natom
          fax(j)=0.d0
          fay(j)=0.d0
          faz(j)=0.d0
       end do

       if(i.le.2000)then
          call opt_intra(eintra,einter)
       else
          call opt_inter(einter)
       end if

       enpot=eintra+einter

       call opt_check(gax,gay,gaz,dfmax)

       if(mod(i,50).eq.0)write(6,20)'SD',&
            i,eintra*econv,einter*econv,enpot*econv,abs(enpot-enpot0)*econv,dfmax*econv/rconv
       if(i.ge.2)write(3,30)&
            i,eintra*econv,einter*econv,enpot*econv,abs(enpot-enpot0)*econv,dfmax*econv/rconv

       do j=1,natom
          gax(j)=fax(j)
          gay(j)=fay(j)
          gaz(j)=faz(j)
       end do

       call geometria

       if(dfmax.le.opt_dfmax)exit

       enpot0=enpot

    end do

    write(6,'(4x,111a1)')('-',i=1,84)
    write(6,*)

    return

5   format(5x,a1,14x,a5,9x,a5,8x,a6,8x,a7,7x,a6)
10  format(5x,a2,6x,a4,6x,a5,9x,a5,9x,a6,7x,a7,8x,a6)
20  format(5x,a2,2x,i8,5(2x,es12.4))
30  format(5x,i8,5(2x,es12.4))

  end subroutine opt

  subroutine opt_intra(eintra,einter)

    implicit none

!    integer i,j
!    real(8) dfmax,gax(natom),gay(natom),gaz(natom)
    real(8) encoul,enbond,enbend,entors,envdw,eintra,einter
    real(8) virvdw,virbond,virbend,virtors,vircoul

    encoul=0.d0   !coulombiano
    enbond=0.d0   !estiramento
    enbend=0.d0   !deformacao
    entors=0.d0   !torção
    envdw=0.d0    !Van der waals

    call ff_modules_intra&
         (enbond,enbend,entors,envdw,encoul,virbond,virbend,virtors,virvdw,vircoul)

    call ff_modules_inter(envdw,encoul,virvdw,vircoul)

    eintra=enbond+enbend+entors+envdw+encoul
    einter=envdw+encoul

    call steepest_descent

    return

  end subroutine opt_intra

  subroutine opt_inter(einter)

    implicit none

!    integer i,j
!    real(8) dfmax,gax(natom),gay(natom),gaz(natom)
    real(8) encoul,envdw,einter
    real(8) virvdw,vircoul

    encoul=0.d0   !coulombiano
    envdw=0.d0    !Van der waals

    call ff_modules_inter(envdw,encoul,virvdw,vircoul)

    call steepest_descent_CM

    call verlet_list_inter

    einter=envdw+encoul

    return

  end subroutine opt_inter

  subroutine steepest_descent

    implicit none

    integer i
    real(8) df

    do i=1,natom
       df=sqrt(fax(i)**2+fay(i)**2+faz(i)**2)
       xa(i)=xa(i)+opt_gamma*fax(i)/df
       ya(i)=ya(i)+opt_gamma*fay(i)/df
       za(i)=za(i)+opt_gamma*faz(i)/df
    end do

    return

  end subroutine steepest_descent

  subroutine steepest_descent_CM

    implicit none

    integer i,j,k,nx
    real(8) df,mtotal,fcm(3)

    nx=1
    do i=1,nmolec
       do j=1,ntmolec(i)
          fcm(1)=0.d0
          fcm(2)=0.d0
          fcm(3)=0.d0
          mtotal=0.d0
          do k=1,nxmolec(i)
             fcm(1)=fcm(1)+mass(nx)*fax(nx)
             fcm(2)=fcm(2)+mass(nx)*fay(nx)
             fcm(3)=fcm(3)+mass(nx)*faz(nx)
             mtotal=mtotal+mass(nx)
          end do
          do k=1,nxmolec(i)
             df=sqrt(fcm(1)**2+fcm(2)**2+fcm(3)**2)
             xa(nx)=xa(nx)+opt_gamma*fcm(1)/df
             ya(nx)=ya(nx)+opt_gamma*fcm(2)/df
             za(nx)=za(nx)+opt_gamma*fcm(3)/df
             nx=nx+1
          end do
       end do
    end do

    return

  end subroutine steepest_descent_CM

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

  subroutine opt_print
    !*************************************************************************
    ! Impressão dos dados de entrada                                         *
    !*************************************************************************
    implicit none

    integer i,j

    write(6,'(93a1)')('#',i=1,93)
    write(6,*)(' OPT',i=1,23)
    write(6,'(93a1)')('#',i=1,93)
    write(6,*)
    write(6,'(30x,a30)')'Steepest descent information'
    write(6,'(28x,36a1)')('-',j=1,36)
    write(6,'(28x,a12,6x,i10)')'ntrialmax:',opt_ntrialmax
    write(6,'(28x,a12,8x,es10.3,1x,a4)')'dfmax:',opt_dfmax*econv/rconv,'eV/A'
    write(6,'(28x,a12,8x,es10.3)')'gamma:',opt_gamma
    write(6,'(28x,36a1)')('-',j=1,36)
    write(6,*)

    return

  end subroutine opt_print

end module optimize
