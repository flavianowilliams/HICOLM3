program teste_pes

  implicit none

  integer npes,i,j,natom
  integer bnd(2,2),bend(1,3),nbnd,nbend,qatom,molec
  real(8) x(100,3),y(100,3),z(100,3),dx(100,3),dy(100,3),dz(100,3)
  real(8) pes(100),xx1,yy1,zz1,xx2,yy2,zz2,rsq,xx1sqt,xx2sqt,tht,thtt
  real(8) eps,prmbnd(2),prmbend(2)

  common /params/nbnd,nbend,eps,prmbnd,prmbend,bend,bnd
  common /paratm/natom,npes,qatom,molec
  common /parcoor/x,y,z

  npes=30
  natom=3
  qatom=1

  open(1,file='structure-old.pes',status='old')
  open(2,file='structure.pes',status='unknown')

  !-lendo PES em arquivo de entrada

  do i=1,npes
     do j=1,natom
        read(1,*)x(i,j),y(i,j),z(i,j),dx(i,j),dy(i,j),dz(i,j)
     end do
     read(1,*)
     read(1,*)
  end do

  !-convertendo unidades medida

  do i=1,npes
     do j=1,natom
        x(i,j)=x(i,j)*0.529177249d0
        y(i,j)=y(i,j)*0.529177249d0
        z(i,j)=z(i,j)*0.529177249d0
        !        dx(i,j)=dx(i,j)*27.2114d0/0.529177d0
        !        dy(i,j)=dy(i,j)*27.2114d0/0.529177d0
        !        dz(i,j)=dz(i,j)*27.2114d0/0.529177d0
     end do
  end do

  !-calculando PES

  nbnd=2

  bnd(1,1)=1
  bnd(1,2)=2
  bnd(2,1)=1
  bnd(2,2)=3

  eps=-2700d0

  prmbnd(1)=45.d0
  prmbnd(2)=1.0d0

  nbend=1

  bend(1,1)=1
  bend(1,2)=2
  bend(1,3)=3

  prmbend(1)=3.d0
  prmbend(2)=110.d0

  call energia(pes)

  call forca(dx,dy,dz)

  !-convertendo unidades medida

  do i=1,npes
     do j=1,natom
        x(i,j)=x(i,j)/0.529177249d0
        y(i,j)=y(i,j)/0.529177249d0
        z(i,j)=z(i,j)/0.529177249d0
        dx(i,j)=dx(i,j)*0.529177d0/27.2114d0
        dy(i,j)=dy(i,j)*0.529177d0/27.2114d0
        dz(i,j)=dz(i,j)*0.529177d0/27.2114d0
     end do
     pes(i)=pes(i)/27.2114d0
  end do

  !-imprimindo PES em arquivo de saida

  do i=1,npes
     do j=1,natom
        write(2,'(6f15.8)')x(i,j),y(i,j),z(i,j),dx(i,j),dy(i,j),dz(i,j)
     end do
     write(2,'(f21.15)')pes(i)
     write(2,*)
  end do

end program teste_pes

subroutine energia(pes)

  implicit none

  integer i,j
  integer natom,npes
  integer bnd(2,2),bend(1,3),nbnd,nbend
  real(8) x(100,3),y(100,3),z(100,3),pes(100)
  real(8) xx1,yy1,zz1,rsq,xx1sqt,xx2sqt,tht,thtt,xx2,yy2,zz2
  real(8) eps,prmbnd(2),prmbend(2)

  common /params/nbnd,nbend,eps,prmbnd,prmbend,bend,bnd
  common /paratm/natom,npes
  common /parcoor/x,y,z
  
  do i=1,npes
     pes(i)=eps
     do j=1,nbnd
        xx1=x(i,bnd(j,2))-x(i,bnd(j,1))
        yy1=y(i,bnd(j,2))-y(i,bnd(j,1))
        zz1=z(i,bnd(j,2))-z(i,bnd(j,1))
        rsq=sqrt(xx1**2+yy1**2+zz1**2)
        pes(i)=pes(i)+0.5d0*prmbnd(1)*(rsq-prmbnd(2))**2
     end do
     do j=1,nbend
        xx1=x(i,bend(j,2))-x(i,bend(j,1))
        yy1=y(i,bend(j,2))-y(i,bend(j,1))
        zz1=z(i,bend(j,2))-z(i,bend(j,1))
        xx2=x(i,bend(j,3))-x(i,bend(j,1))
        yy2=y(i,bend(j,3))-y(i,bend(j,1))
        zz2=z(i,bend(j,3))-z(i,bend(j,1))
        xx1sqt=sqrt(xx1**2+yy1**2+zz1**2)
        xx2sqt=sqrt(xx2**2+yy2**2+zz2**2)
        tht=acos((xx1*xx2+yy1*yy2+zz1*zz2)/(xx1sqt*xx2sqt))
        thtt=prmbend(2)*acos(-1.d0)/180.d0
        pes(i)=pes(i)+0.5d0*prmbend(1)*(tht-thtt)**2
     end do
  end do

  return

end subroutine energia

subroutine forca(dx,dy,dz)

  implicit none

  integer i,j,g,nx
  integer natom,npes
  integer bnd(2,2),bend(1,3),nbnd,nbend,qatom,molec
  real(8) x(100,3),y(100,3),z(100,3),dx(100,3),dy(100,3),dz(100,3)
  real(8) xx1,yy1,zz1,rsq,xx1sqt,xx2sqt,tht,thtt,xx2,yy2,zz2
  real(8) eps,prmbnd(2),prmbend(2),xn(2)

  common /params/nbnd,nbend,eps,prmbnd,prmbend,bend,bnd
  common /paratm/natom,npes,qatom,molec
  common /parcoor/x,y,z

  do i=1,npes
     do j=1,natom
        dx(i,j)=0.d0
        dy(i,j)=0.d0
        dz(i,j)=0.d0
     end do
  end do

  do i=1,2
     xn(i)=prmbnd(i)
  end do

  do l=1,npes
     nx=0
     do i=1,molec
        g=nx+(i-1)*qatom
        do j=1,nbnd
           call fbond(l,bnd(j,g+1),bnd(j,g+2),xn,dx,dy,dz)
        end do
     end do
     nx=nx+molec*qatom
  end do

  do i=1,2
     xn(i)=prmbend(i)
  end do

  do i=1,npes
     do j=1,nbend
        call fbend(i,bend(j,1),bend(j,2),bend(j,3),xn,dx,dy,dz)
     end do
  end do

  return

end subroutine forca

subroutine fbond(l,ix1,ix2,xn,dx,dy,dz)

  implicit none

  integer l,ix1,ix2,im,i
  real(8) xx,xvz,yvz,zvz,fr,xn(2)
  real(8) x(100,3),y(100,3),z(100,3),dx(100,3),dy(100,3),dz(100,3)

  common /parcoor/x,y,z

  xvz=x(l,ix2)-x(l,ix1)
  yvz=y(l,ix2)-y(l,ix1)
  zvz=z(l,ix2)-z(l,ix1)

  xx=sqrt(xvz**2+yvz**2+zvz**2)

  fr=xn(1)*(xx-xn(2))

  fr=-fr/xx

  dx(l,ix2)=dx(l,ix2)+fr*xvz
  dy(l,ix2)=dy(l,ix2)+fr*yvz
  dz(l,ix2)=dz(l,ix2)+fr*zvz

  dx(l,ix1)=dx(l,ix1)-fr*xvz
  dy(l,ix1)=dy(l,ix1)-fr*yvz
  dz(l,ix1)=dz(l,ix1)-fr*zvz

  return

end subroutine fbond

subroutine fbend(l,ix1,ix2,ix3,xn,dx,dy,dz)

  implicit none

  integer l,ix1,ix2,ix3,ix(3),i,j,im,kronij

  external kronij

  real(8) tht,thtt,fa,derij(3,3),drij(3),drik(3),dr1,dr2
  real(8) xn(2),x(100,3),y(100,3),z(100,3),dx(100,3),dy(100,3),dz(100,3)

    common /parcoor/x,y,z

    drij(1)=x(l,ix2)-x(l,ix1)
    drij(2)=y(l,ix2)-y(l,ix1)
    drij(3)=z(l,ix2)-z(l,ix1)
    drik(1)=x(l,ix3)-x(l,ix1)
    drik(2)=y(l,ix3)-y(l,ix1)
    drik(3)=z(l,ix3)-z(l,ix1)

    dr1=sqrt(drij(1)**2+drij(2)**2+drij(3)**2)
    dr2=sqrt(drik(1)**2+drik(2)**2+drik(3)**2)

    tht=acos((drij(1)*drik(1)+drij(2)*drik(2)+drij(3)*drik(3))/(dr1*dr2))
    thtt=xn(2)*acos(-1.d0)/180.d0

    fa=xn(1)*(tht-thtt)

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

    dx(l,ix1)=dx(l,ix1)+fa*derij(1,1)/sin(tht)
    dy(l,ix1)=dy(l,ix1)+fa*derij(1,2)/sin(tht)
    dz(l,ix1)=dz(l,ix1)+fa*derij(1,3)/sin(tht)

    dx(l,ix2)=dx(l,ix2)+fa*derij(2,1)/sin(tht)
    dy(l,ix2)=dy(l,ix2)+fa*derij(2,2)/sin(tht)
    dz(l,ix2)=dz(l,ix2)+fa*derij(2,3)/sin(tht)

    dx(l,ix3)=dx(l,ix3)+fa*derij(3,1)/sin(tht)
    dy(l,ix3)=dy(l,ix3)+fa*derij(3,2)/sin(tht)
    dz(l,ix3)=dz(l,ix3)+fa*derij(3,3)/sin(tht)

    return

  end subroutine fbend

  integer function kronij(i,j)
    !************************************************************************************
    ! Delta de Kronecker                                                                *
    !************************************************************************************

    implicit none

    integer i,j

    kronij=int((float(i+j)-abs(i-j))/(float(i+j)+abs(i-j)))

    return

  end function kronij
