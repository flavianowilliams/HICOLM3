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

    integer i,j,ihist,nx
    real(8) dfmax,gax(natom),gay(natom),gaz(natom)
    real(8) encoul,enbond,enbend,entors,envdw,eintra,einter,enpot0,enpot
    real(8) virvdw,virbond,virbend,virtors,vircoul

    !-criando ficheiro de saida

    open(3,file='HICOLM_opt.dat',status='unknown')

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

    !-calculando valores iniciais

    ihist=1
    nx=1

    do j=1,natom
       fax(j)=0.d0
       fay(j)=0.d0
       faz(j)=0.d0
    end do

    do i=1,natom
       gax(i)=0.d0
       gay(i)=0.d0
       gaz(i)=0.d0
    end do

    encoul=0.d0   !coulombiano
    envdw=0.d0    !Van der waals
    enbond=0.d0   !estiramento
    enbend=0.d0   !deformacao
    entors=0.d0   !torção

    vircoul=0.d0   !coulombiano
    virbond=0.d0   !estiramento
    virbend=0.d0   !deformacao
    virtors=0.d0   !torção
    virvdw=0.d0    !Van der waals

    call ff_modules_intra&
         (enbond,enbend,entors,envdw,encoul,virbond,virbend,virtors,virvdw,vircoul)
    call ff_modules_inter(envdw,encoul,virvdw,vircoul)

    eintra=enbond+enbend+entors+envdw+encoul
    einter=envdw+encoul+envdw_corr
    enpot=eintra+einter

    call opt_check(gax,gay,gaz,dfmax)

    call geometria

    if(mod(i,nhist).eq.0)call history(ihist)

    write(6,20)'SD',&
         1,eintra*econv,einter*econv,enpot*econv,abs(enpot-enpot0)*econv,dfmax*econv/rconv
    write(3,30)&
         1,eintra*econv,einter*econv,enpot*econv,abs(enpot-enpot0)*econv,dfmax*econv/rconv

    enpot0=enpot

    !-calculando contribuição intramolecular

    do i=2,opt_ntotal

       encoul=0.d0   !coulombiano
       envdw=0.d0    !Van der waals
       enbond=0.d0   !estiramento
       enbend=0.d0   !deformacao
       entors=0.d0   !torção

       vircoul=0.d0   !coulombiano
       virbond=0.d0   !estiramento
       virbend=0.d0   !deformacao
       virtors=0.d0   !torção
       virvdw=0.d0    !Van der waals

       call verlet_list_inter

       if(i.le.opt_ninter)then
          call ff_modules_inter(envdw,encoul,virvdw,vircoul)
          call steepest_descent_CM
          einter=envdw+encoul+envdw_corr
       else
          call ff_modules_intra&
               (enbond,enbend,entors,envdw,encoul,virbond,virbend,virtors,virvdw,vircoul)
          call ff_modules_inter(envdw,encoul,virvdw,vircoul)
          call steepest_descent
          eintra=enbond+enbend+entors+envdw+encoul
          einter=envdw+encoul+envdw_corr
       end if

       enpot=eintra+einter

       call opt_check(gax,gay,gaz,dfmax)

       call geometria

       if(mod(i,nhist).eq.0)call history(ihist)

       if(mod(i,50).eq.0)write(6,20)'SD',&
            i,eintra*econv,einter*econv,enpot*econv,abs(enpot-enpot0)*econv,dfmax*econv/rconv
       write(3,30)&
            i,eintra*econv,einter*econv,enpot*econv,abs(enpot-enpot0)*econv,dfmax*econv/rconv

       if(dfmax.le.opt_dfmax)exit

       do j=1,natom
          fax(j)=0.d0
          fay(j)=0.d0
          faz(j)=0.d0
       end do

       enpot0=enpot
       nx=nx+1

    end do

    write(6,'(4x,111a1)')('-',i=1,84)
    write(6,*)

    !-terminando impressao em history

    do i=nx+1,opt_ntotal
       if(mod(i,nhist).eq.0)call history(ihist)
    end do

    return

5   format(5x,a1,14x,a5,9x,a5,8x,a6,8x,a7,7x,a6)
10  format(5x,a2,6x,a4,6x,a5,9x,a5,9x,a6,7x,a7,8x,a6)
20  format(5x,a2,2x,i8,5(2x,es12.4))
30  format(5x,i8,5(2x,es12.4))

  end subroutine opt

  subroutine steepest_descent

    implicit none

    integer i

    do i=1,natom
       xa(i)=xa(i)+opt_gamma*fax(i)
       ya(i)=ya(i)+opt_gamma*fay(i)
       za(i)=za(i)+opt_gamma*faz(i)
       end do

    return

  end subroutine steepest_descent

  subroutine steepest_descent_CM

    implicit none

    integer i,j,k,nx
    real(8) mtotal,fcm(3)

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
             xa(nx)=xa(nx)+opt_gamma*fcm(1)/mtotal
             ya(nx)=ya(nx)+opt_gamma*fcm(2)/mtotal
             za(nx)=za(nx)+opt_gamma*fcm(3)/mtotal
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

    do i=1,natom
       gax(i)=fax(i)
       gay(i)=fay(i)
       gaz(i)=faz(i)
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
    write(6,'(28x,a12,6x,i10)')'n(interm):',opt_ninter
    write(6,'(28x,a12,6x,i10)')' n(total):',opt_ntotal
    write(6,'(28x,a12,6x,i10)')'    nhist:',nhist
    write(6,'(28x,a12,1x,es10.2,1x,a10)')'dfmax:',opt_dfmax*econv/rconv,'kcal/mol*A'
    write(6,'(28x,a12,1x,es10.2,1x,a12)')'gamma:',opt_gamma*rconv**2/econv,'A^2*mol/kcal'
    write(6,'(28x,36a1)')('-',j=1,36)
    write(6,*)

    return

  end subroutine opt_print

end module optimize
