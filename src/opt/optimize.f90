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
          call steepest_descent_CM(gax,gay,gaz)
!          call steepest_descent_rotation
          einter=envdw+encoul+envdw_corr
       else
          if(mod(i,2).eq.0)then
             call ff_modules_inter(envdw,encoul,virvdw,vircoul)
             call steepest_descent_CM(gax,gay,gaz)
             einter=envdw+encoul+envdw_corr
          else
             call ff_modules_intra&
                  (enbond,enbend,entors,envdw,encoul,virbond,virbend,virtors,virvdw,vircoul)
             call steepest_descent(gax,gay,gaz)
             eintra=enbond+enbend+entors+envdw+encoul
          end if
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

    !-imprimindo mensagem de aviso

    if(nx.eq.opt_ntotal)then
       write(6,'(4x,a66)')'Warning: The optimization does not reach the convergence criteria!'
       write(6,*)
    end if

    return

5   format(5x,a1,14x,a5,9x,a5,8x,a6,8x,a7,7x,a6)
10  format(5x,a2,6x,a4,6x,a5,9x,a5,9x,a6,7x,a7,8x,a6)
20  format(5x,a2,2x,i8,5(2x,es12.4))
30  format(5x,i8,5(2x,es12.4))

  end subroutine opt

  subroutine steepest_descent(gax,gay,gaz)

    implicit none

    integer i
    real(8) gax(natom),gay(natom),gaz(natom),dfx,dfy,dfz,dr,df

    do i=1,natom
       dfx=opt_gamma*fax(i)+opt_alpha*(fax(i)-gax(i))
       dfy=opt_gamma*fay(i)+opt_alpha*(fay(i)-gay(i))
       dfz=opt_gamma*faz(i)+opt_alpha*(faz(i)-gaz(i))
       df=sqrt(dfx**2+dfy**2+dfz**2)
       dr=max(1.d-8,min(df,opt_rshift))
       xa(i)=xa(i)+dfx*dr/df
       ya(i)=ya(i)+dfy*dr/df
       za(i)=za(i)+dfz*dr/df
    end do

    return

  end subroutine steepest_descent

  subroutine steepest_descent_CM(gax,gay,gaz)

    implicit none

    integer i,j,k,nx
    real(8) mtotal,fcm(3),dfx,dfy,dfz
    real(8) gax(natom),gay(natom),gaz(natom),dr,df

    nx=1
    do i=1,nmolec
       do j=1,ntmolec(i)
          fcm(1)=0.d0
          fcm(2)=0.d0
          fcm(3)=0.d0
          mtotal=0.d0
          do k=0,(nxmolec(i)-1)
             fcm(1)=fcm(1)+fax(nx+k)
             fcm(2)=fcm(2)+fay(nx+k)
             fcm(3)=fcm(3)+faz(nx+k)
             mtotal=mtotal+mass(nx+k)
          end do
          do k=1,nxmolec(i)
             dfx=opt_gamma*fcm(1)+opt_alpha*(fax(nx)-gax(nx))
             dfy=opt_gamma*fcm(2)+opt_alpha*(fay(nx)-gay(nx))
             dfz=opt_gamma*fcm(3)+opt_alpha*(faz(nx)-gaz(nx))
             df=sqrt(dfx**2+dfy**2+dfz**2)
             dr=max(1.d-8,min(df,opt_rshift))
             xa(nx)=xa(nx)+dfx*dr/df
             ya(nx)=ya(nx)+dfy*dr/df
             za(nx)=za(nx)+dfz*dr/df
             nx=nx+1
          end do
       end do
    end do

    return

  end subroutine steepest_descent_CM

  subroutine steepest_descent_rotation

    implicit none

    integer i,j,k,nx
    real(8) rcm(3),mtotal,mrot(3,3),tx,ty,tz,theta,tr,mi,xl(3),dtx,dty,dtz

    nx=1
    do i=1,nmolec
       do j=1,ntmolec(i)
          rcm(1)=0.d0
          rcm(2)=0.d0
          rcm(3)=0.d0
          mtotal=0.d0
          do k=0,(nxmolec(i)-1)
             rcm(1)=rcm(1)+mass(nx+k)*xa(nx+k)
             rcm(2)=rcm(2)+mass(nx+k)*ya(nx+k)
             rcm(3)=rcm(3)+mass(nx+k)*za(nx+k)
             mtotal=mtotal+mass(nx+k)
          end do
          rcm(1)=rcm(1)/mtotal
          rcm(2)=rcm(2)/mtotal
          rcm(3)=rcm(3)/mtotal
          tx=0.d0
          ty=0.d0
          tz=0.d0
          mi=0.d0
          do k=0,(nxmolec(i)-1)
             tx=tx+((ya(nx+k)-rcm(2))*faz(nx+k)-(za(nx+k)-rcm(3))*fay(nx+k))
             ty=ty+((za(nx+k)-rcm(3))*fax(nx+k)-(xa(nx+k)-rcm(1))*faz(nx+k))
             tz=tz+((xa(nx+k)-rcm(1))*fay(nx+k)-(ya(nx+k)-rcm(2))*fax(nx+k))
             mi=mi+mass(nx+k)*((xa(nx+k)-rcm(1))**2+(ya(nx+k)-rcm(2))**2+(za(nx+k)-rcm(3))**2)
          end do
          tr=sqrt(tx**2+ty**2+tz**2)
          tx=tx/tr
          ty=ty/tr
          tz=tz/tr
          theta=opt_beta!*tr!/mi
          do k=1,nxmolec(i)
             call rotation_matrix(theta,tx,ty,tz,nx,mrot,xl)
             dtx=xl(1)-xa(nx)
             dty=xl(2)-ya(nx)
             dtz=xl(3)-za(nx)
!             dt=sqrt(dtx**2+dty**2+dtz**2)
!             dr=max(1.d-8,min(dt,opt_rshift))
             xa(nx)=xa(nx)+dtx!*dr/dt
             ya(nx)=ya(nx)+dty!*dr/dt
             za(nx)=za(nx)+dtz!*dr/dt
             nx=nx+1
          end do
       end do
    end do

    return

  end subroutine steepest_descent_rotation

  subroutine rotation_matrix(theta,ta,tb,tc,nx,mrot,xl)

    implicit none

    integer nx
    real(8) theta,ta,tb,tc,mrot(3,3),xl(3)

    mrot(1,1)=cos(theta)+(1.d0-cos(theta))*ta**2
    mrot(1,2)=(1.d0-cos(theta))*ta*tb+sin(theta)*tc
    mrot(1,3)=(1.d0-cos(theta))*ta*tc-sin(theta)*tb
    mrot(2,1)=(1.d0-cos(theta))*tb*ta-sin(theta)*tc
    mrot(2,2)=cos(theta)+(1.d0-cos(theta))*tb**2
    mrot(2,3)=(1.d0-cos(theta))*tb*tc+sin(theta)*ta
    mrot(3,1)=(1.d0-cos(theta))*tc*ta+sin(theta)*tb
    mrot(3,2)=(1.d0-cos(theta))*tc*tb-sin(theta)*ta
    mrot(3,3)=cos(theta)+(1.d0-cos(theta))*tc**2

    xl(1)=mrot(1,1)*xa(nx)+mrot(1,2)*ya(nx)+mrot(1,3)*za(nx)
    xl(2)=mrot(2,1)*xa(nx)+mrot(2,2)*ya(nx)+mrot(2,3)*za(nx)
    xl(3)=mrot(3,1)*xa(nx)+mrot(3,2)*ya(nx)+mrot(3,3)*za(nx)

    return

  end subroutine rotation_matrix

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
    write(6,'(28x,a12,1x,es10.2,1x,a12)')'alpha:',opt_alpha*rconv**2/econv,'A^2*mol/kcal'
    write(6,'(28x,a12,8x,2f6.3,1x,a1)')'rcutoff:',rcutoff*rconv,drcutoff*rconv,'A'
    write(6,'(28x,a12,8x,f6.2,1x,a1)')' rshift:',opt_rshift*rconv,'A'
    write(6,'(28x,36a1)')('-',j=1,36)
    write(6,*)

    return

  end subroutine opt_print

end module optimize
