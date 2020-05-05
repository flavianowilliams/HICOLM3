module objetivo

  use input
  use populacao

  real(8) enmax

  save enmax

contains

!  double precision function funcao(xx)
!
!    implicit none
!
!    real(8) xx
!
!    funcao=xx**2
!
!    return
!
!  end function funcao
!
!  subroutine aptidao_bin(xn,fx)
!
!    implicit none
!
!    real(8) fx,fmax,funcx,xn(nparam)
!
!    fmax=funcao(nimax(1))
!    funcx=funcao(xn(1))
!
!    fx=fmax-funcx
!
!    return
!
!  end subroutine aptidao_bin

  subroutine aptidao(xn,en)

    implicit none

    integer i,j,i1,i2,i3,opt,l,k,nx,g,nxx,ii,jj,kk
    real(8) en,diff,fmax,fr,xn(nparam)

    nx=1
    do i=1,nmolec
       do j=1,nprmstr(i)
          if(chkvar(1,i,j).eqv..true.)then
             prmstr(i,j)=xn(nx)
             nx=nx+1
          end if
       end do
       do j=1,nprmscr(i)
          if(chkvar(2,i,j).eqv..true.)then
             prmscr(i,j)=xn(nx)
             nx=nx+1
          end if
       end do
    end do

    do i=1,nvdw
       do j=1,nprmvdw(mvdw(i,1),mvdw(i,2))
          if(chkvdw(mvdw(i,1),mvdw(i,2),j).eqv..true.)then
             prmvdw(mvdw(i,1),mvdw(i,2),j)=xn(nx)
             nx=nx+1
          end if
       end do
    end do

    opt=2

    if(opt.eq.1)goto 1
    if(opt.eq.2)goto 2

1   en=0.d0
    do l=1,mgeo
       nx=0
       do i=1,nmolec
          diff=prmstr(i,1)
!          diff=0.d0
          do j=1,molec(i)
             g=nx+(j-1)*qatom(i)
             do k=1,mnstr(i)
                i1=g+mstr(i,k,1)
                i2=g+mstr(i,k,2)
                diff=diff+bond(l,i,i1,i2)
             end do
             do k=1,mnscr(i)
                i1=g+mscr(i,k,1)
                i2=g+mscr(i,k,2)
                i3=g+mscr(i,k,3)
                diff=diff+bend(l,i,i1,i2,i3)
             end do
          end do
          nx=nx+molec(i)*qatom(i)
       end do
!       write(*,*)diff,func(l)
       en=en+(diff-func(l))**2
    end do

    en=sqrt(en/mgeo)
    en=1.d0/max(en,1.d-8)

    goto 10

    !-zerando forcas

2   do i=1,mgeo
       do j=1,natom
          fx(i,j)=0.d0
          fy(i,j)=0.d0
          fz(i,j)=0.d0
       end do
    end do

    !-calculo das forcas intramoleculares

    en=0.d0
    do l=1,mgeo
       nx=0
       do i=1,nmolec
          do j=1,molec(i)
             g=nx+(j-1)*qatom(i)
             do k=1,mnstr(i)
                i1=g+mstr(i,k,1)
                i2=g+mstr(i,k,2)
                call fbond(l,i,i1,i2)
             end do
             do k=1,mnscr(i)
                i1=g+mscr(i,k,1)
                i2=g+mscr(i,k,2)
                i3=g+mscr(i,k,3)
                call fbend(l,i,i1,i2,i3)
             end do
          end do
          nx=nx+molec(i)*qatom(i)
       end do
!       stop
    end do

    !-calculo das forcas intermoleculares

    do l=1,mgeo
       nx=0
       do i=1,nmolec
          do j=1,molec(i)
             do k=1,qatom(i)
                i1=nx+k+(j-1)*qatom(i)
                nxx=0
                do ii=1,nmolec
                   do jj=j+1,molec(ii)
                      do kk=1,qatom(ii)
                         i2=nxx+kk+(jj-1)*qatom(ii)
                         call fvdw(l,atpx(i1),atpx(i2),i1,i2)
                         call fcoul(l,i1,i2)
                      end do
                   end do
                   nxx=nxx+molec(i)*qatom(ii)
                end do
             end do
          end do
          nx=nx+molec(i)*qatom(i)
       end do
!       stop
    end do

    !-calculo do erro quadratico da forca

    en=0.d0
    fmax=0.d0
    do i=1,mgeo
       do j=1,natom
          fr=(fx(i,j)-dx(i,j))**2+(fy(i,j)-dy(i,j))**2+(fz(i,j)-dz(i,j))**2
          en=en+fr
          fmax=max(fmax,sqrt(fr))
       end do
    end do

    en=sqrt(en/float(mgeo*natom))
    en=1.d0/max(en,1.d-8)

10  return

  end subroutine aptidao

  double precision function bond(l,im,ix1,ix2)

    implicit none

    integer l,ix1,ix2,im,i
    real(8) xx,xn(nprmstr(im)),xvz,yvz,zvz

    xvz=xa(l,ix2)-xa(l,ix1)
    yvz=ya(l,ix2)-ya(l,ix1)
    zvz=za(l,ix2)-za(l,ix1)

    xx=sqrt(xvz**2+yvz**2+zvz**2)

    do i=2,nprmstr(im)
       xn(i)=prmstr(im,i)
    end do

    bond=0.d0

    select case(ntstr(im))
    case(1)
       bond=xn(1)*(exp(-2*xn(2)*(xx-xn(3)))-2*exp(-xn(2)*(xx-xn(3))))
    case(2)
       bond=0.5d0*xn(2)*(xx-xn(3))**2
    end select

    return

  end function bond

  double precision function bend(l,im,ix2,ix1,ix3)

    implicit none

    integer l,ix1,ix2,ix3,im,i
    real(8) xx1sqt,xx2sqt,tht,xx1,yy1,zz1,xx2,yy2,zz2,thtt
    real(8) xn(nprmscr(im))

    xx1=xa(l,ix2)-xa(l,ix1)
    yy1=ya(l,ix2)-ya(l,ix1)
    zz1=za(l,ix2)-za(l,ix1)
    xx2=xa(l,ix3)-xa(l,ix1)
    yy2=ya(l,ix3)-ya(l,ix1)
    zz2=za(l,ix3)-za(l,ix1)

    xx1sqt=sqrt(xx1**2+yy1**2+zz1**2)
    xx2sqt=sqrt(xx2**2+yy2**2+zz2**2)

    tht=acos((xx1*xx2+yy1*yy2+zz1*zz2)/(xx1sqt*xx2sqt))

    do i=1,nprmscr(im)
       xn(i)=prmscr(im,i)
    end do

    bend=0.d0

    select case(ntscr(im))
    case(1)
       !       thtt=180.d0*tht/acos(-1.d0)
       !       bend=0.5d0*xn(1)*(thtt-xn(2))**2
       thtt=xn(2)*acos(-1.d0)/180.d0
       bend=0.5d0*xn(1)*(tht-thtt)**2
    end select

    return

  end function bend

  double precision function vdw(op,l,ix1,ix2,xn)

    implicit none

    integer op,l,ix1,ix2
    real(8) xx,xn(nparam)

    xx=sqrt((xa(l,ix2)-xa(l,ix1))**2+(ya(l,ix2)-ya(l,ix1))**2+(za(l,ix2)-za(l,ix1))**2)

    vdw=0.d0

    select case(op)
    case(1)
       vdw=xn(6)*(exp(-2*xn(7)*(xx-xn(8)))-2*exp(-xn(7)*(xx-xn(8))))
    case(2)
       vdw=4.d0*xn(6)*((xn(7)/xx)**12-(xn(7)/xx)**6)
    end select

    return

  end function vdw

  double precision function coul(l,ix1,ix2,i,j,ii,jj)

    implicit none

    integer l,ix1,ix2,i,j,ii,jj
    real(8) xx

    xx=sqrt((xa(l,ix2)-xa(l,ix1))**2+(ya(l,ix2)-ya(l,ix1))**2+(za(l,ix2)-za(l,ix1))**2)

    coul=-kelect*chq(i,j)*chq(ii,jj)/xx

    return

  end function coul

  subroutine fbond(l,im,ix1,ix2)

    implicit none

    integer l,ix1,ix2,im,i
    real(8) xx,xvz,yvz,zvz,fr,xn(nprmstr(im))

    xvz=xa(l,ix2)-xa(l,ix1)
    yvz=ya(l,ix2)-ya(l,ix1)
    zvz=za(l,ix2)-za(l,ix1)

    xx=sqrt(xvz**2+yvz**2+zvz**2)

    do i=1,nprmstr(im)
       xn(i)=prmstr(im,i)
    end do

    fr=0.d0

    select case(ntstr(im))
    case(1)
       fr=xn(1)*(exp(-2*xn(2)*(xx-xn(3)))-2*exp(-xn(2)*(xx-xn(3))))
    case(2)
       fr=xn(2)*(xx-xn(3))
    end select

    fr=-fr/xx

    fx(l,ix2)=fx(l,ix2)+fr*xvz
    fy(l,ix2)=fy(l,ix2)+fr*yvz
    fz(l,ix2)=fz(l,ix2)+fr*zvz

    fx(l,ix1)=fx(l,ix1)-fr*xvz
    fy(l,ix1)=fy(l,ix1)-fr*yvz
    fz(l,ix1)=fz(l,ix1)-fr*zvz

    return

  end subroutine fbond

  subroutine fbend(l,im,ix2,ix1,ix3)

    implicit none

    integer l,ix1,ix2,ix3,ix(3),i,j,im
    real(8) tht,thtt,fa,derij(3,3),drij(3),drik(3),dr1,dr2
    real(8) xn(nprmscr(im))

    drij(1)=xa(l,ix2)-xa(l,ix1)
    drij(2)=ya(l,ix2)-ya(l,ix1)
    drij(3)=za(l,ix2)-za(l,ix1)
    drik(1)=xa(l,ix3)-xa(l,ix1)
    drik(2)=ya(l,ix3)-ya(l,ix1)
    drik(3)=za(l,ix3)-za(l,ix1)

    dr1=sqrt(drij(1)**2+drij(2)**2+drij(3)**2)
    dr2=sqrt(drik(1)**2+drik(2)**2+drik(3)**2)

    tht=acos((drij(1)*drik(1)+drij(2)*drik(2)+drij(3)*drik(3))/(dr1*dr2))

    do i=1,nprmscr(im)
       xn(i)=prmscr(im,i)
    end do

    fa=0.d0

    select case(ntscr(im))
    case(1)
       thtt=xn(2)*acos(-1.d0)/180.d0
       fa=xn(1)*(tht-thtt)
    end select

    !-calculo das forças atomicas  i, j e k

    ix(1)=ix1
    ix(2)=ix2
    ix(3)=ix3

    do j=1,3
       do i=1,3
          derij(i,j)=(kronij(ix(i),ix(2))-kronij(ix(i),ix(1)))*drik(j)/(dr1*dr2) &
               +(kronij(ix(i),ix(3))-kronij(ix(i),ix(1)))*drij(j)/(dr1*dr2) &
               -cos(tht)*((kronij(ix(i),ix(2))-kronij(ix(i),ix(1)))*drij(j)/dr1**2 &
               +(kronij(ix(i),ix(3))-kronij(ix(i),ix(1)))*drik(j)/dr2**2)
       end do
    end do

    fx(l,ix1)=fx(l,ix1)+fa*derij(1,1)/sin(tht)
    fy(l,ix1)=fy(l,ix1)+fa*derij(1,2)/sin(tht)
    fz(l,ix1)=fz(l,ix1)+fa*derij(1,3)/sin(tht)

    fx(l,ix2)=fx(l,ix2)+fa*derij(2,1)/sin(tht)
    fy(l,ix2)=fy(l,ix2)+fa*derij(2,2)/sin(tht)
    fz(l,ix2)=fz(l,ix2)+fa*derij(2,3)/sin(tht)

    fx(l,ix3)=fx(l,ix3)+fa*derij(3,1)/sin(tht)
    fy(l,ix3)=fy(l,ix3)+fa*derij(3,2)/sin(tht)
    fz(l,ix3)=fz(l,ix3)+fa*derij(3,3)/sin(tht)

    return

  end subroutine fbend

  subroutine fvdw(l,im,in,ix1,ix2)

    implicit none

    integer l,ix1,ix2,im,in,i
    real(8) xx,xvz,yvz,zvz,fr,xn(nparam)

    xvz=xa(l,ix2)-xa(l,ix1)
    yvz=ya(l,ix2)-ya(l,ix1)
    zvz=za(l,ix2)-za(l,ix1)

    xx=sqrt(xvz**2+yvz**2+zvz**2)

    do i=1,nprmvdw(im,in)
       xn(i)=prmvdw(im,in,i)
    end do

    fr=0.d0

    select case(ntvdw(im,in))
    case(1)
       fr=-24.d0*xn(1)*(2.d0*(xn(2)/xx)**12-(xn(2)/xx)**6)/xx
    end select

    fr=-fr/xx

!    write(*,*)im,in,xn(1),xn(2),nprmvdw(im,in),fr,xx

    fx(l,ix2)=fx(l,ix2)+fr*xvz
    fy(l,ix2)=fy(l,ix2)+fr*yvz
    fz(l,ix2)=fz(l,ix2)+fr*zvz

    fx(l,ix1)=fx(l,ix1)-fr*xvz
    fy(l,ix1)=fy(l,ix1)-fr*yvz
    fz(l,ix1)=fz(l,ix1)-fr*zvz

    return

  end subroutine fvdw

  subroutine fcoul(l,ix1,ix2)

    implicit none

    integer l,ix1,ix2
    real(8) xx,xvz,yvz,zvz,fr

    xvz=xa(l,ix2)-xa(l,ix1)
    yvz=ya(l,ix2)-ya(l,ix1)
    zvz=za(l,ix2)-za(l,ix1)

    xx=sqrt(xvz**2+yvz**2+zvz**2)

    fr=-kelect*chqx(ix1)*chqx(ix2)/xx**2

    fr=-fr/xx

    fx(l,ix2)=fx(l,ix2)+fr*xvz
    fy(l,ix2)=fy(l,ix2)+fr*yvz
    fz(l,ix2)=fz(l,ix2)+fr*zvz

    fx(l,ix1)=fx(l,ix1)-fr*xvz
    fy(l,ix1)=fy(l,ix1)-fr*yvz
    fz(l,ix1)=fz(l,ix1)-fr*zvz

    return

  end subroutine fcoul

  integer function kronij(i,j)
    !************************************************************************************
    ! Delta de Kronecker                                                                *
    !************************************************************************************

    implicit none

    integer i,j

    kronij=int((float(i+j)-abs(i-j))/(float(i+j)+abs(i-j)))

    return

  end function kronij

 subroutine funcaomax

    implicit none

    integer i,j
    real(8) en,xn(nparam)

    enmax=0.d0

    do i=1,npop
       do j=1,nparam
          xn(j)=ind(i,j)
       end do
!       if(opcod.eq.1)call aptidao_bin(xn,en)
       call aptidao(xn,en)
       enmax=max(abs(en),enmax)
    end do

    return

  end subroutine funcaomax

end module objetivo
