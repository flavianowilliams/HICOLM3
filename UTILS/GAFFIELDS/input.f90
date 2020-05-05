module input

  integer mgeo,natom,nparam,npop,ncron,prec,nstp,nvdw,nscr,nstr,opcr
  integer tvdw(10),svdw(10,2)
  integer nmolec,molec(10),qatom(10),moltot,atp(10,100),atpx(10*100),idax(10*100),ida(10,100)
!  integer nbnd

  integer mnstr(10),mnscr(10),ntstr(10),ntscr(10),ntvdw(10,10),mstr(10,10,2),mscr(10,10,3)
  integer mvdw(20,2),nprmvdw(10,10)

  integer nxbnd,bndtp(10)

  real(8) txgr,pmut,tol,dind
  real(8) ind0(500),nimin(500),nimax(500)
  real(8) chq(10,100),chqx(10*100)
  real(8) kelect

  integer nprmstr(10),nprmscr(10)
  real(8) prmstr(10,3),prmscr(10,2),prmvdw(10,10,3)

  real(8), allocatable :: dx(:,:),dy(:,:),dz(:,:)
  real(8), allocatable :: fx(:,:),fy(:,:),fz(:,:)
  real(8), allocatable :: xa(:,:),ya(:,:),za(:,:)
  real(8), allocatable :: func(:)

  character(6) idmolec(10)

  logical chkvar(2,10,3),chkvdw(10,10,3)

  save txgr,nimin,nimax,mgeo,nparam,prec,pmut,tol,ind0,dind
  save tvdw,svdw,nstp,nvdw,nscr,nstr,nprmstr,nprmscr,ntvdw,nprmvdw,mvdw,chkvar,chkvdw
  save nmolec,molec,qatom,chq,kelect,moltot,idmolec,bndtp,nxbnd,atpx

contains

  subroutine alloc_inicial

    implicit none

    !-alocando arrays

    allocate(xa(mgeo,natom),ya(mgeo,natom),za(mgeo,natom),func(mgeo))
    allocate(fx(mgeo,natom),fy(mgeo,natom),fz(mgeo,natom))
    allocate(dx(mgeo,natom),dy(mgeo,natom),dz(mgeo,natom))

    return

  end subroutine alloc_inicial

  subroutine parametros

    implicit none

    integer i,j,k

    open(1,file='gaffields.in',status='old')
    open(2,file='structure.pes',status='old')

    read(1,*)
    read(1,*)mgeo
    read(1,*)
    read(1,*)npop
    read(1,*)
    read(1,*)nstp
    read(1,*)
    read(1,*)dind
    read(1,*)
    read(1,*)opcr
    read(1,*)
    read(1,*)txgr,pmut
    read(1,*)
    read(1,*)tol
    read(1,*)

    do i=1,10
       do j=1,10
          do k=1,3
             prmvdw(i,j,k)=0.d0
          end do
          nprmvdw(i,j)=2
       end do
       do k=1,2
          prmstr(i,k)=0.d0
       end do
       do k=1,2
          prmscr(i,k)=0.d0
       end do
    end do

    call input_molecs

    !-definindo cromossomo inicial

    call input_cromoss

    write(*,*)

    kelect=14.36559d0

    return

  end subroutine parametros

  subroutine input_molecs

    implicit none

    integer i,j,k

    read(1,*)nmolec

    do i=1,nmolec
       read(1,*)idmolec(i)
       read(1,*)molec(i),qatom(i)
       read(1,*)(atp(i,j),j=1,qatom(i))
       read(1,*)(ida(i,j),j=1,qatom(i))
       read(1,*)(chq(i,j),j=1,qatom(i))
       read(1,*)
       read(1,*)mnstr(i),ntstr(i)
       call param_check(1,ntstr(i),nprmstr(i))
       read(1,*)(prmstr(i,j),j=1,nprmstr(i))
       read(1,*)(chkvar(1,i,j),j=1,nprmstr(i))
       do j=1,mnstr(i)
          read(1,*)(mstr(i,j,k),k=1,2)
       end do
       read(1,*)
       read(1,*)mnscr(i),ntscr(i)
       call param_check(2,ntscr(i),nprmscr(i))
       read(1,*)(prmscr(i,j),j=1,nprmscr(i))
       read(1,*)(chkvar(2,i,j),j=1,nprmscr(i))
       do j=1,mnscr(i)
          read(1,*)(mscr(i,j,k),k=1,3)
       end do
    end do

    read(1,*)
    read(1,*)nvdw
    do i=1,nvdw
       read(1,*)mvdw(i,1),mvdw(i,2),k
       backspace(1)
       call vdwchk(k,nprmvdw(mvdw(i,1),mvdw(i,2)))
       read(1,*)mvdw(i,1),mvdw(i,2),ntvdw(mvdw(i,1),mvdw(i,2)),&
            (prmvdw(mvdw(i,1),mvdw(i,2),j),j=1,nprmvdw(mvdw(i,1),mvdw(i,2)))
       read(1,*)(chkvdw(mvdw(i,1),mvdw(i,2),j),j=1,nprmvdw(mvdw(i,1),mvdw(i,2)))
    end do

    natom=1
    moltot=0
    nstr=0
    nscr=0
    do i=1,nmolec
       do j=1,molec(i)
          do k=1,qatom(i)
             atpx(natom)=atp(i,k)
             chqx(natom)=chq(i,k)
             idax(natom)=ida(i,k)
             natom=natom+1
          end do
          moltot=moltot+1
       end do
       nstr=nstr+mnstr(i)
       nscr=nscr+mnscr(i)
    end do

    natom=natom-1

    !-imprimindo informacoes em output

    write(*,*)'Estrutura:'
    write(*,*)

    write(*,'(4x,111a1)')('-',i=1,52)
    write(*,10)'Tipo','Qde','Sites','bonds','bends','dihdl'
    write(*,'(4x,111a1)')('-',i=1,52)

    do i=1,nmolec
       write(6,20)idmolec(i),molec(i),qatom(i),mnstr(i),mnscr(i)
    end do

    write(*,'(4x,111a1)')('-',i=1,52)
    write(*,20)'Total:',moltot,natom,nstr,nscr
    write(*,*)

    write(*,*)'Parametros intermoleculares:'
    write(*,*)
    write(*,*)'Van der waals:',nvdw
    write(*,'(4x,111a1)')('-',i=1,52)
    write(*,30)'Sites','Sites','Tipo'
    write(*,'(4x,111a1)')('-',i=1,52)

    do i=1,nvdw
       write(*,40)mvdw(i,1),mvdw(i,2)
    end do

    write(*,'(4x,111a1)')('-',i=1,52)

    return

10  format(a4,7x,a3,4x,a6,4x,a5,4x,a5,4x,a5)
20  format(4x,a6,2x,i5,4(4x,i5))
30  format(4x,a3,7x,a5,4x,a5)
40  format(4x,4(4x,i5))

  end subroutine input_molecs

  subroutine param_check(ffop,opt,nx)

    implicit none

    integer opt,ffop,nx

    if(ffop.eq.1)goto 1
    if(ffop.eq.2)goto 2

    !-contagem de parametros para ajuste

1   select case(opt)
    case(1)
       nx=3
    case(2)
       nx=3
    end select

    return

2   select case(opt)
    case(1)
       nx=2
    end select

    return

  end subroutine param_check

  subroutine vdwchk(i,op)

    implicit none

    integer i,op

    select case(i)
    case(1)
       op=2
    end select

    return

  end subroutine vdwchk

  subroutine input_cromoss

    implicit none

    integer i,j

    nparam=1
    do i=1,nmolec
       do j=1,nprmstr(i)
          if(chkvar(1,i,j).eqv..true.)then
             ind0(nparam)=prmstr(i,j)
             nparam=nparam+1
          end if
       end do
       do j=1,nprmscr(i)
          if(chkvar(2,i,j).eqv..true.)then
             ind0(nparam)=prmscr(i,j)
             nparam=nparam+1
          end if
       end do
    end do

    do i=1,nvdw
       do j=1,nprmvdw(mvdw(i,1),mvdw(i,2))
          if(chkvdw(mvdw(i,1),mvdw(i,2),j).eqv..true.)then
             ind0(nparam)=prmvdw(mvdw(i,1),mvdw(i,2),j)
             nparam=nparam+1
          end if
       end do
    end do

    nparam=nparam-1

    do i=1,nparam
       nimin(i)=ind0(i)-abs(dind*ind0(i))
       nimax(i)=ind0(i)+abs(dind*ind0(i))
    end do

!    if(opcod.eq.1)ncron=int(log2((nimax(1)-nimin(1))*10**prec))+1
    ncron=nparam

    write(*,*)'Cromossomo:'
    write(*,*)
    write(*,'(4x,a27)')'Codificacao ponto flutuante'
!    write(*,'(4x,a9,i5)')'Precisao:',prec
!    write(*,'(4x,a22,i5)')'Tamanho do cromossomo:',ncron
!    write(*,*)
    write(*,'(4x,a36,i5)')'Quantidade de parametros a otimizar:',nparam
    write(*,*)
    write(*,'(4x,a16,15f14.4)')'Valores minimos:',(nimin(i),i=1,nparam)
    write(*,'(4x,a16,15f14.4)')'Valores maximos:',(nimax(i),i=1,nparam)
    write(*,*)

    return

  end subroutine input_cromoss

  subroutine coords

    implicit none

    integer i,j

    do i=1,mgeo
       do j=1,natom
          read(2,*)xa(i,j),ya(i,j),za(i,j),dx(i,j),dy(i,j),dz(i,j)
       end do
       read(2,*)func(i)
       read(2,*)
    end do

    write(*,*)'PES:'
    write(*,*)
    write(*,'(4x,a21,i5)')'Quantidade de passos:',mgeo
    write(*,'(4x,a11,f7.4)')'Tolerancia:',tol
    write(*,*)

    return

  end subroutine coords

  subroutine convert

    implicit none

    integer i,j

    !-convertendo unidades de medida

    do i=1,mgeo
       do j=1,natom
          xa(i,j)=xa(i,j)*0.529177d0
          ya(i,j)=ya(i,j)*0.529177d0
          za(i,j)=za(i,j)*0.529177d0
          dx(i,j)=dx(i,j)*(-27.2114d0/0.529177d0)
          dy(i,j)=dy(i,j)*(-27.2114d0/0.529177d0)
          dz(i,j)=dz(i,j)*(-27.2114d0/0.529177d0)
       end do
       func(i)=func(i)*27.2114d0
    end do

!    do i=1,nmolec
!       prmscr(i,2)=acos(-1.d0)*prmscr(i,2)/180.d0
!    end do

    return

  end subroutine convert

  double precision function log2(x)

    implicit none

    real(8) x

    log2=log(x)/log(2.d0)

    return

  end function log2

end module input
