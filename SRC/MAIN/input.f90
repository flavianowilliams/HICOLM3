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
module input
  !****************************************************************************************
  ! Leitura dos dados de entrada                                                          *
  ! Flaviano Williams Fernandes, 05 de dezembro de 2018                                   *
  !****************************************************************************************

  use sistema
  use elements
  use amber

  implicit none
  !----------------------------------------------------------------------------
  !-constantes
  !
  real(8) pi,rconv,econv,kconv,pconv,aconv,mconv,tconv,teconv,hcconv
  real(8) elconv,keconv,n0,kb,kelect
  !
  complex(8) img
  !
  save rconv,econv,pconv,kconv,aconv,mconv,tconv,teconv,hcconv,elconv,keconv
  save kb,pi,img,n0
  !----------------------------------------------------------------------------
  !-variaveis de estrutura
  integer rxmx,rymx,rzmx,gsymopt,natom,reuse,idna(natmax),atp(natmax),fztp(natmax)
  !
  real(8) a,b,c,xa(natmax),ya(natmax),za(natmax),v(3,3),volume
  real(8) vax(natmax),vay(natmax),vaz(natmax),fax(natmax),fay(natmax),faz(natmax)
  !
  save rxmx,rymx,rzmx,gsymopt,natom,reuse,idna,atp,fztp
  !----------------------------------------------------------------------------
  !-variaveis da mecanica molecular
  integer nhist,ntrialmax,nrelax,bendscnt(molecmax),bondscnt(molecmax),nzmolec(molecmax)
  integer ato(natmax),nrl,atnp(ntpmax),natnp(ntpmax),nxmolec(molecmax),ntmolec(molecmax)
  integer vdw(ntpmax,ntpmax),bonds(molecmax,bondmax),bends(molecmax,bendmax),tersoff,coulop
  integer tors(molecmax,torsmax),molbend(molecmax,bendmax,3),torscnt(molecmax)
  integer dstp,ndstp,xstp,nmolec,moltot,nfree,spctot,molbond(molecmax,bondmax,2)
  integer nvdw,ncoul,nbonds,nbends,ntors,moltors(molecmax,torsmax,4)
  !
  real(8) dtime,drmax,fmstp,text,tstat,preext,pstat,bfactor,lrmax,zmatrix_tol
  real(8) parbnd(molecmax,bondmax,5),parvdw(ntpmax,ntpmax,3),fzstr(6)
  real(8) parbend(molecmax,bendmax,4),parcoul(ntpmax,1)
  real(8) partors(molecmax,torsmax,7)
  real(8) mass(natmax),massmin,massmax,rcutoff,drcutoff,lambdain,lambdafi,sf_vdw,sf_coul
  !
  character(2) att
  character(2) atsp(ntpmax)
  character(7) prop,ensble
  character(9) ensble_mt
  character(10) method,namemol(molecmax)
  character(7) ff_model(molecmax)
  !
  save dtime,drmax,nhist,ntrialmax,nrelax,fmstp,dstp,xstp,ndstp,drcutoff,zmatrix_tol
  save mass,massmin,massmax,prop,namemol,nmolec,moltot,nfree,spctot,ensble,ensble_mt,rcutoff
  save att,ato,nrl,atnp,natnp,text,preext,pstat,bfactor,tstat,lrmax,fzstr,sf_vdw,sf_coul
  save parbnd,parvdw,parbend,parcoul,lambdain,lambdafi,partors
  save vdw,bonds,bends,nbonds,nbends,ntors,coulop,tersoff,tors,atsp
  save bendscnt,molbend,bondscnt,molbond,torscnt,ntmolec,nzmolec,nxmolec,ff_model
  !
  !----------------------------------------------------------------------------
  !-variaveis da mecanica molecular
  integer opt_ntotal,opt_ninter
  real(8) opt_dfmax,opt_gamma
  !
  save opt_ntotal,opt_ninter,opt_dfmax,opt_gamma
  !
contains

  subroutine entrada(t0)
    !***************************************************************************************
    ! Leitura dos parametros de entrada                                                    *
    !***************************************************************************************
    implicit none

    real(8) t0,tf
    character(10) in
    logical lval

    open(5,file='HICOLM.in',status='old')
    open(10,file='HICOLM.str',status='old')

    !-contagem de tempo

    call cpu_time(t0)

    !-atribuindo valores default

    call default

    !-lendo estrutura

    call sys_input

    !-lendo parametros do metodo escolhido

    rewind(5)

    lval=.false.
    do while (lval.eqv..false.)

       read(5,*,end=11)in

       if(in.eq.'@MDPREPARE')then
          call opt_input
          lval=.true.
       elseif(in.eq.'@MDRUNNING')then
          call md_input
          lval=.true.
       end if

    end do

    !-lendo informacoes do campo de força

    call ff_input

    !-contagem de tempo

    call cpu_time(tf)

    t0=tf-t0

    return

11  stop 'entrada: Method does not found'

  end subroutine entrada

  subroutine sys_input
    !***************************************************************************************
    ! Leitura dos dados do sistema                                                         *
    !***************************************************************************************
    implicit none

    integer i,j,k,natfx,ival(20)
    character(7) in
    character(9) key

    rewind(5)

    !-base

1   read(5,*,end=11)in

    if(in.ne.'&STRUCT')goto 1

    do i=1,100
       read(5,*)key
       if(key.eq.'&END')exit
       if(key.eq.'cell')then
          do j=1,3
             read(5,*)(v(j,k),k=1,3)
          end do
          a=sqrt(v(1,1)**2+v(1,2)**2+v(1,3)**2)
          b=sqrt(v(2,1)**2+v(2,2)**2+v(2,3)**2)
          c=sqrt(v(3,1)**2+v(3,2)**2+v(3,3)**2)
       end if
       if(key.eq.'molecs')then
          backspace(5)
          read(5,*)key,nmolec
          do j=1,nmolec
             read(5,*)namemol(j),ntmolec(j),nxmolec(j)
          end do
       end if
       if(key.eq.'reuse')then
          backspace(5)
          read(5,*)key,ival(1)
          reuse=ival(1)
       end if
    end do

    !-calculando moleculas e atomos totais

11  moltot=1
    natom=0
    do i=1,nmolec
       do j=1,ntmolec(i)
          do k=1,nxmolec(i)
             natom=natom+1
          end do
          nzmolec(moltot)=nxmolec(i)
          moltot=moltot+1
       end do
    end do

    moltot=moltot-1

    !-lendo estrutura inicial

    do i=1,natom
       read(10,*)idna(i),xa(i),ya(i),za(i),atp(i),fztp(i)
    end do

    !-definindo especies e total de especies

!    nx=1
!    nxx=0
!    do i=1,nmolec
!       do j=1,ntmolec(i)
!          do k=1,nxmolec(i)
!             atp(nx)=k+nxx
!             nx=nx+1
!          end do
!       end do
!       nxx=nxx+nxmolec(i)
!    end do

    atnp(1)=atp(1)

    spctot=1
    do i=2,natom
       do j=1,spctot
          if(atp(i).eq.atnp(j))goto 553
       end do
       atnp(spctot+1)=atp(i)
       spctot=spctot+1
553    continue
    end do

    do j=1,spctot
       natnp(j)=0
    end do

    do i=1,natom
       do j=1,spctot
          if(atp(i).eq.atnp(j))natnp(j)=natnp(j)+1
       end do
    end do

    !-calculando graus de liberdade

    natfx=0
    do i=1,natom
       if(fztp(i).eq.0)natfx=natfx+1
    end do

    nfree=3*(natom-natfx-1)

    !-deslocando sistema para o centro da supercelula

    do i=1,natom
       xa(i)=xa(i)+0.5d0*(v(1,1)+v(2,1)+v(3,1))
       ya(i)=ya(i)+0.5d0*(v(1,2)+v(2,2)+v(3,2))
       za(i)=za(i)+0.5d0*(v(1,3)+v(2,3)+v(3,3))
    end do

    !-definindo massa atomica

    call atomic_mass

    return

  end subroutine sys_input

  subroutine zmatrix(imol)

    implicit none

    integer i,j,k,kk,ii,nx,imol,ia,ib
    real(8) dr,rca,rcb

    nx=0
    do i=2,imol
       nx=nx+nxmolec(i)*ntmolec(i)
    end do

    !-calculo das ligacoes quimicas

    bondscnt(imol)=1
    do i=1,nxmolec(imol)
       ia=i+nx
       do j=i+1,nxmolec(imol)
          ib=j+nx
          dr=sqrt((xa(ia)-xa(ib))**2+(ya(ia)-ya(ib))**2+(za(ia)-za(ib))**2)
          call covalent_radius(idna(ia),rca)
          call covalent_radius(idna(ib),rcb)
          if(dr.gt.(rca+rcb-zmatrix_tol).and.dr.le.(rca+rcb+zmatrix_tol))then
             molbond(imol,bondscnt(imol),1)=i
             molbond(imol,bondscnt(imol),2)=j
             bondscnt(imol)=bondscnt(imol)+1
          end if
       end do
    end do

    bondscnt(imol)=bondscnt(imol)-1

    !-calculo das deformacoes angulares

    bendscnt(imol)=1
    do i=1,bondscnt(imol)
       do j=i+1,bondscnt(imol)
          do k=1,2
             do kk=1,2
                if(molbond(imol,i,k).eq.molbond(imol,j,kk))then
                   molbend(imol,bendscnt(imol),1)=molbond(imol,i,3-k)
                   molbend(imol,bendscnt(imol),2)=molbond(imol,i,k)
                   molbend(imol,bendscnt(imol),3)=molbond(imol,j,3-kk)
                   bendscnt(imol)=bendscnt(imol)+1
                end if
             end do
          end do
       end do
    end do

    bendscnt(imol)=bendscnt(imol)-1

    !-calcudos dos diedros

    torscnt(imol)=1
    do i=1,bondscnt(imol)
       do j=1,i-1
          do k=1,j-1
             do ii=1,2
                do kk=1,2
                   if(molbond(imol,j,1).eq.molbond(imol,i,3-ii))then
                      if(molbond(imol,j,2).eq.molbond(imol,k,kk))then
                         moltors(imol,torscnt(imol),1)=molbond(imol,i,ii)
                         moltors(imol,torscnt(imol),2)=molbond(imol,j,1)
                         moltors(imol,torscnt(imol),3)=molbond(imol,j,2)
                         moltors(imol,torscnt(imol),4)=molbond(imol,k,3-kk)
                         torscnt(imol)=torscnt(imol)+1
                      end if
                   end if
                end do
             end do
          end do
          do k=j+1,bondscnt(imol)
             do ii=1,2
                do kk=1,2
                   if(molbond(imol,j,1).eq.molbond(imol,i,3-ii))then
                      if(molbond(imol,j,2).eq.molbond(imol,k,kk))then
                         moltors(imol,torscnt(imol),1)=molbond(imol,i,ii)
                         moltors(imol,torscnt(imol),2)=molbond(imol,j,1)
                         moltors(imol,torscnt(imol),3)=molbond(imol,j,2)
                         moltors(imol,torscnt(imol),4)=molbond(imol,k,3-kk)
                         torscnt(imol)=torscnt(imol)+1
                      end if
                   end if
                end do
             end do
          end do
       end do
       do j=i+1,bondscnt(imol)
          do k=1,j-1
             do ii=1,2
                do kk=1,2
                   if(molbond(imol,j,1).eq.molbond(imol,i,3-ii))then
                      if(molbond(imol,j,2).eq.molbond(imol,k,kk))then
                         moltors(imol,torscnt(imol),1)=molbond(imol,i,ii)
                         moltors(imol,torscnt(imol),2)=molbond(imol,j,1)
                         moltors(imol,torscnt(imol),3)=molbond(imol,j,2)
                         moltors(imol,torscnt(imol),4)=molbond(imol,k,3-kk)
                         torscnt(imol)=torscnt(imol)+1
                      end if
                   end if
                end do
             end do
          end do
          do k=j+1,bondscnt(imol)
             do ii=1,2
                do kk=1,2
                   if(molbond(imol,j,1).eq.molbond(imol,i,3-ii))then
                      if(molbond(imol,j,2).eq.molbond(imol,k,kk))then
                         moltors(imol,torscnt(imol),1)=molbond(imol,i,ii)
                         moltors(imol,torscnt(imol),2)=molbond(imol,j,1)
                         moltors(imol,torscnt(imol),3)=molbond(imol,j,2)
                         moltors(imol,torscnt(imol),4)=molbond(imol,k,3-kk)
                         torscnt(imol)=torscnt(imol)+1
                      end if
                   end if
                end do
             end do
          end do
       end do
    end do

    torscnt(imol)=torscnt(imol)-1

    return

  end subroutine zmatrix

  subroutine opt_input

    implicit none

    integer i,ival(20)
    real(8) val(20)
    character(7) in
    character(10) key

    rewind(5)

1   read(5,*,end=11)in

    if(in.ne.'&OPT')goto 1

    do i=1,100
       read(5,*)key
       if(key.eq.'&END')exit
       if(key.eq.'ntrialmax')then
          backspace(5)
          read(5,*)key,ival(1),ival(2)
          opt_ninter=ival(1)
          opt_ntotal=ival(2)
       end if
       if(key.eq.'gamma')then
          backspace(5)
          read(5,*)key,val(1)
          opt_gamma=val(1)
       end if
       if(key.eq.'dfmax')then
          backspace(5)
          read(5,*)key,val(1)
          opt_dfmax=val(1)
       end if
       if(key.eq.'rcutoff')then
          backspace(5)
          read(5,*)key,val(1),val(2)
          rcutoff=val(1)
          drcutoff=val(2)
       end if
    end do

    !-escolha do metodo

11  method='@MDPREPARE'

    return

  end subroutine opt_input

  subroutine md_input
    !***************************************************************************************
    ! Leitura dos parametros de entrada da Dinamica molecular                              *
    !***************************************************************************************
    implicit none

    integer i,j,ival(20)
    real(8) val(20)
    character(7) in,char
    character(9) char2
    character(10) key

    rewind(5)

1   read(5,*,end=11)in

    if(in.ne.'&MD')goto 1

    do i=1,100
       read(5,*)key
       if(key.eq.'&END')exit
       if(key.eq.'nhist')then
          backspace(5)
          read(5,*)key,ival(1)
          nhist=ival(1)
       end if
       if(key.eq.'ntrialmax')then
          backspace(5)
          read(5,*)key,ival(1)
          ntrialmax=ival(1)
       end if
       if(key.eq.'nrelax')then
          backspace(5)
          read(5,*)key,ival(1)
          nrelax=ival(1)
       end if
       if(key.eq.'rcutoff')then
          backspace(5)
          read(5,*)key,val(1),val(2)
          rcutoff=val(1)
          drcutoff=val(2)
       end if
       if(key.eq.'timestep')then
          backspace(5)
          read(5,*)key,val(1)
          dtime=val(1)
       end if
       if(key.eq.'drmax')then
          backspace(5)
          read(5,*)key,val(1),val(2),ival(1),ival(2),ival(3)
          drmax=val(1)
          fmstp=val(2)
          dstp=ival(1)
          ndstp=ival(2)
          xstp=ival(3)
       end if
       if(key.eq.'lambda')then
          backspace(5)
          read(5,*)key,(val(j),j=1,2)
          lambdain=val(1)
          lambdafi=val(2)
       end if
       if(key.eq.'tstat')then
          backspace(5)
          read(5,*)key,val(1)
          tstat=val(1)
       end if
       if(key.eq.'text')then
          backspace(5)
          read(5,*)key,val(1)
          text=val(1)
       end if
       if(key.eq.'preext')then
          backspace(5)
          read(5,*)key,val(1)
          preext=val(1)
       end if
       if(key.eq.'pstat')then
          backspace(5)
          read(5,*)key,val(1),val(2)
          pstat=val(1)
          bfactor=val(2)
       end if
       if(key.eq.'ensemble')then
          backspace(5)
          read(5,*)key,char
          if(char.eq.'nve')then
             ensble=char
          elseif(char.eq.'nvt')then
             backspace(5)
             read(5,*)key,char,char2
             if(char2.eq.'berendsen')then
                backspace(5)
                read(5,*)key,char,char2,val(1)
                ensble=char
                ensble_mt=char2
                tstat=val(1)
             elseif(char2.eq.'hoover')then
                backspace(5)
                read(5,*)key,char,char2,val(1)
                ensble=char
                ensble_mt=char2
                tstat=val(1)
             end if
          elseif(char.eq.'npt')then
             backspace(5)
             read(5,*)key,char,char2
             if(char2.eq.'berendsen')then
                backspace(5)
                read(5,*)key,char,char2,val(1),val(2),val(3)
                ensble=char
                ensble_mt=char2
                tstat=val(1)
                pstat=val(2)
                bfactor=val(3)
             elseif(char2.eq.'hoover')then
                backspace(5)
                read(5,*)key,char,char2,val(1),val(2)
                ensble=char
                ensble_mt=char2
                tstat=val(1)
                pstat=val(2)
             end if
          end if
       end if
       if(key.eq.'fzstr')then
          backspace(5)
          read(5,*)key,(val(j),j=1,6)
          do j=1,6
             fzstr(j)=val(j)
          end do
       end if
    end do

    !-escolha do metodo

11  method='@MDRUNNING'

    return

  end subroutine md_input

  subroutine ff_input

    implicit none

    integer i,ii,iii,iv,v,j,jj,k,g,p,m,numt,nx,ival(20),spctt
    real(8) val(20)
    character(7) in,char
    character(10) lxmol,key

    !-Campo de Forca

    rewind(5)

1   read(5,*,end=11)in

    if(in.ne.'&FORCE')goto 1

    do i=1,100
       read(5,*)key
       if(key.eq.'&END')exit
       if(key.eq.'$AMBER')then
          read(5,*)key
          backspace(5)
          if(key.eq.'zmatrix')then
             read(5,*)key,zmatrix_tol
          end if
          spctt=0
          do j=1,nmolec
             if(key.eq.'$END')goto 433
             read(5,*)key,lxmol
             do g=1,nmolec
                if(lxmol.eq.namemol(g))then
                   read(5,*)(atsp(k+spctt),k=1,nxmolec(g))
                   read(5,*)(parcoul(k+spctt,1),k=1,nxmolec(g))
                   ff_model(g)='(AMBER)'
                   spctt=spctt+nxmolec(g)
                   nx=g
                end if
             end do
             call zmatrix(nx)
             jj=spctt-nxmolec(nx)
             do k=1,bondscnt(nx)
                ii=jj+molbond(nx,k,1)
                iii=jj+molbond(nx,k,2)
                call amber_bonds(atsp(ii),atsp(iii),val)
                do p=1,2
                   parbnd(nx,k,p)=val(p)
                end do
                bonds(nx,k)=3
             end do
             do k=1,bendscnt(nx)
                ii=jj+molbend(nx,k,1)
                iii=jj+molbend(nx,k,2)
                iv=jj+molbend(nx,k,3)
                call amber_bends(atsp(ii),atsp(iii),atsp(iv),val)
                do p=1,2
                   parbend(nx,k,p)=val(p)
                end do
                bends(nx,k)=2
             end do
             do k=1,torscnt(nx)
                ii=jj+moltors(nx,k,1)
                iii=jj+moltors(nx,k,2)
                iv=jj+moltors(nx,k,3)
                v=jj+moltors(nx,k,4)
                call amber_dihedrals(atsp(ii),atsp(iii),atsp(iv),atsp(v),val)
                do p=1,4
                   partors(nx,k,p)=val(p)
                end do
                tors(nx,k)=4
             end do
          end do
433       do j=1,spctt
             do jj=j,spctt
                call amber_vdw(atsp(j),atsp(jj),val)
                do k=1,2
                   parvdw(j,jj,k)=val(k)
                   parvdw(jj,j,k)=val(k)
                end do
                vdw(j,jj)=3
                vdw(jj,j)=3
             end do
          end do
          coulop=2
       elseif(key.eq.'$INTRA')then
          do j=1,nmolec
             read(5,*)key,lxmol
             do g=1,molecmax
                if(lxmol.eq.namemol(g))nx=g
             end do
             do g=1,100
                read(5,*)key
                backspace(5)
                if(key.eq.'molecule')exit
                if(key.eq.'$END')goto 523
                if(key.eq.'bends#')then
                   read(5,*)key,bendscnt(nx)
                   do k=1,bendscnt(nx)
                      read(5,*)ival(1),ival(2),ival(3),char
                      call bend_opt(char,m,numt)
                      backspace(5)
                      read(5,*)ival(1),ival(2),ival(3),char,(val(p),p=1,numt)
                      do p=1,numt
                         parbend(nx,k,p)=val(p)
                      end do
                      bends(nx,k)=m
                      molbend(nx,k,1)=ival(1)
                      molbend(nx,k,2)=ival(2)
                      molbend(nx,k,3)=ival(3)
                   end do
                elseif(key.eq.'bends*')then
                   read(5,*)key,ii
                   read(5,*)ival(1),ival(2),ival(3),char
                   call bend_opt(char,m,numt)
                   backspace(5)
                   read(5,*)ival(1),ival(2),ival(3),char,(val(p),p=1,numt)
                   do p=1,numt
                      parbend(nx,ii,p)=val(p)
                   end do
                   bends(nx,ii)=m
                   molbend(nx,ii,1)=ival(1)
                   molbend(nx,ii,2)=ival(2)
                   molbend(nx,ii,3)=ival(3)
                elseif(key.eq.'bonds#')then
                   read(5,*)key,bondscnt(nx)
                   do k=1,bondscnt(nx)
                      read(5,*)ival(1),ival(2),char
                      call bonds_opt(char,m,numt)
                      backspace(5)
                      read(5,*)ival(1),ival(2),char,(val(p),p=1,numt)
                      do p=1,numt
                         parbnd(nx,k,p)=val(p)
                      end do
                      bonds(nx,k)=m
                      molbond(nx,k,1)=ival(1)
                      molbond(nx,k,2)=ival(2)
                   end do
                elseif(key.eq.'bonds*')then
                   read(5,*)key,ii
                   read(5,*)ival(1),ival(2),char
                   call bonds_opt(char,m,numt)
                   backspace(5)
                   read(5,*)ival(1),ival(2),char,(val(p),p=1,numt)
                   do p=1,numt
                      parbnd(nx,ii,p)=val(p)
                   end do
                   bonds(nx,ii)=m
                   molbond(nx,ii,1)=ival(1)
                   molbond(nx,ii,2)=ival(2)
                elseif(key.eq.'dihedrals#')then
                   read(5,*)key,torscnt(nx)
                   do k=1,torscnt(nx)
                      read(5,*)ival(1),ival(2),ival(3),ival(4),char
                      call dihedral_opt(char,m,numt)
                      backspace(5)
                      read(5,*)ival(1),ival(2),ival(3),ival(4),char,(val(p),p=1,numt)
                      do p=1,numt
                         partors(nx,k,p)=val(p)
                      end do
                      tors(nx,k)=m
                      moltors(nx,k,1)=ival(1)
                      moltors(nx,k,2)=ival(2)
                      moltors(nx,k,3)=ival(3)
                      moltors(nx,k,4)=ival(4)
                   end do
                elseif(key.eq.'dihedrals*')then
                   read(5,*)key,ii
                   read(5,*)ival(1),ival(2),ival(3),ival(4),char
                   call dihedral_opt(char,m,numt)
                   backspace(5)
                   read(5,*)ival(1),ival(2),ival(3),ival(4),char,(val(p),p=1,numt)
                   do p=1,numt
                      partors(nx,ii,p)=val(p)
                   end do
                   tors(nx,ii)=m
                   moltors(nx,ii,1)=ival(1)
                   moltors(nx,ii,2)=ival(2)
                   moltors(nx,ii,3)=ival(3)
                   moltors(nx,ii,4)=ival(4)
                end if
             end do
          end do
       elseif(key.eq.'$INTER')then
          do ii=1,100
             read(5,*)key
             backspace(5)
             if(key.eq.'$END')goto 523
             if(key.eq.'vdw')then
                read(5,*)key,nx
                do j=1,nx
                   read(5,*)ival(1),ival(2),char
                   call vdw_opt(char,m,numt)
                   backspace(5)
                   read(5,*)ival(1),ival(2),char,(val(k),k=1,numt)
                   do k=1,numt
                      parvdw(ival(1),ival(2),k)=val(k)
                      parvdw(ival(2),ival(1),k)=val(k)
                   end do
                   vdw(ival(1),ival(2))=m
                   vdw(ival(2),ival(1))=m
                end do
             elseif(key.eq.'elect')then
                read(5,*)key,ncoul,m
                do j=1,ncoul
                   call coul_opt(m,numt)
                   read(5,*)ival(1),(val(k),k=1,numt)
                   do k=1,numt
                      parcoul(ival(1),k)=val(k)
                   end do
                end do
                coulop=m
             end if
          end do
       end if
523    continue
    end do

    !-calculando parametros intramoleculares totais

11  nbends=0
    nbonds=0
    ntors=0
    do i=1,nmolec
       do j=1,ntmolec(i)
          nbends=nbends+bendscnt(i)
          nbonds=nbonds+bondscnt(i)
          ntors=ntors+torscnt(i)
       end do
    end do

    return

  end subroutine ff_input

  subroutine default
    !***************************************************************************************
    ! Valores padrão                                                                       *

    integer i,j,k

    rcutoff=10.d0
    drcutoff=0.1d0

    reuse=0

    gsymopt=0

    nbonds=0
    nbends=0
    ntors=0

    coulop=1

    nmolec=0
    natom=0

    rxmx=1
    rymx=1
    rzmx=1

    tersoff=0

    do i=1,molecmax
       ff_model(i)=''
       do j=1,torsmax
          do k=1,7
             partors(i,j,k)=0.d0
          end do
          tors(i,j)=0
          moltors(i,j,1)=0
          moltors(i,j,2)=0
          moltors(i,j,3)=0
          moltors(i,j,4)=0
       end do
       do j=1,bendmax
          do k=1,4
             parbend(i,j,k)=0.d0
          end do
          bends(i,j)=0
          molbend(i,j,1)=0
          molbend(i,j,2)=0
          molbend(i,j,3)=0
       end do
       do j=1,bondmax
          do k=1,5
             parbnd(i,j,k)=0.d0
          end do
          bonds(i,j)=0
          molbond(i,j,1)=0
          molbond(i,j,2)=0
       end do
       bendscnt(i)=0
       bondscnt(i)=0
       torscnt(i)=0
    end do

    do i=1,ntpmax
       do j=1,ntpmax
          do k=1,3
             parvdw(i,j,k)=0.d0
          end do
       end do
       parcoul(i,1)=0.d0
       atsp(i)='X'
    end do

    a=0.d0
    b=0.d0
    c=0.d0
    att='C'

    do i=1,natmax
       xa(i)=0.d0
       ya(i)=0.d0
       za(i)=0.d0
       atp(i)=0
       fztp(i)=1
    end do

    nfree=0

    do i=1,6
       fzstr(i)=1.d0
    end do

    !-parametros da Dinâmica Molecular

    dtime=0.0d0
    zmatrix_tol=0.5d0

    sf_vdw=1.d0/2.d0
    sf_coul=1.d0/1.2d0

    !-parametros de otimizacao

    opt_ninter=95000
    opt_ntotal=100000
    opt_dfmax=1.e-4
    opt_gamma=1.e-7

    !-variaveis canonicas

    do i=1,natmax
       vax(i)=0.d0
       vay(i)=0.d0
       vaz(i)=0.d0
       fax(i)=0.d0
       fay(i)=0.d0
       faz(i)=0.d0
    end do

    return

  end subroutine default

  subroutine convert
    !**************************************************************************************
    ! Conversão para o sistema padrão u.a. de unidades                                    *
    !**************************************************************************************

    implicit none

    integer i,j

    !-unidades fundamentais (SI)

    mconv=9.10938291d-31                      ! massa (kg)
    aconv=180.d0/acos(-1.d0)                  ! ângulo (radianos)
    hcconv=1.054571726d-34                    ! Constante de Planck (J*s).
    elconv=1.602176565d-19                    ! carga eletrica (C)
    keconv=8.9875517873681d+9                 ! constante eletrostática (N*m^2/C^2)

    !-unidades derivadas (SI)

    econv=mconv*elconv**4*keconv**2/hcconv**2 ! energia (J)
    rconv=hcconv**2/(mconv*keconv*elconv**2)  ! comprimento (m)
    pconv=econv/rconv**3                      ! pressao (N/m^2)
    tconv=hcconv/econv                        ! tempo (s)
    teconv=econv/kb                           ! temperatura (K)

    !-convertendo unidades SI para unidades de entrada

    mconv=mconv*(1.d+3*n0) !OK
    rconv=rconv*(1.d+10) !OK
    kconv=1.d0/rconv !OK
    econv=econv*(6.241506363094d+18) !OK
    pconv=pconv*(0.9872d-5) !OK
    tconv=tconv*(1.d+12) !OK
    elconv=1.d0 !OK
    teconv=teconv*(6.241506363094d+18) !OK

    !-vetores de rede e coordenadas atomicas

    do i=1,3
       do j=1,3
          v(i,j)=v(i,j)/rconv
       end do
    end do

    a=sqrt(v(1,1)**2+v(1,2)**2+v(1,3)**2)
    b=sqrt(v(2,1)**2+v(2,2)**2+v(2,3)**2)
    c=sqrt(v(3,1)**2+v(3,2)**2+v(3,3)**2)

    do i=1,natom
       xa(i)=xa(i)/rconv
       ya(i)=ya(i)/rconv
       za(i)=za(i)/rconv
    end do

    !-velocidades

    do i=1,natom
       vax(i)=vax(i)/(rconv/tconv)
       vay(i)=vay(i)/(rconv/tconv)
       vaz(i)=vaz(i)/(rconv/tconv)
    end do

    !-forcas

    do i=1,natom
       fax(i)=fax(i)/(econv/rconv)
       fay(i)=fay(i)/(econv/rconv)
       faz(i)=faz(i)/(econv/rconv)
    end do

    !-parametros de otimizacao

    opt_dfmax=opt_dfmax/(econv/rconv)
    opt_gamma=opt_gamma/(rconv**2/econv)

    !-parametros da Dinâmica Molecular

    drmax=drmax/rconv
    preext=preext/pconv
    text=text/teconv
    fmstp=fmstp/(econv/rconv)
    kb=kb/econv

    dtime=dtime/tconv
    tstat=tstat/tconv
    pstat=pstat/tconv
    bfactor=bfactor*pconv
    rcutoff=rcutoff/rconv
    drcutoff=drcutoff/rconv

    !-parametros do Campo de Forca

    kelect=kelect/keconv

    !-unidades de massa atomica

    do i=1,natom
       mass(i)=mass(i)/mconv
    end do
    massmin=massmin/mconv
    massmax=massmax/mconv

    return

  end subroutine convert

  subroutine cte
    !***************************************************************************************
    ! Constantes essenciais                                                                *
    !***************************************************************************************
    implicit none

    img=dcmplx(0.d0,1.d0)      ! numero imaginario
    pi=acos(-1.d0)             ! Pi
    kelect=1.d0                ! constante eletrostática
    n0=6.022d+23               ! numero de Avogadro
    kb=8.6173324d-5            ! constante de Boltzmann

    return

  end subroutine cte

  subroutine vdw_opt(char,m,nprvdw)
    !***************************************************************************************
    ! Flags do potencial de Van der Waals                                                  *
    !***************************************************************************************
    implicit none

    integer nprvdw,m
    character(7) char

    select case(char)
    case('mors')
       m=1
       nprvdw=3
    case('lj')
       m=2
       nprvdw=2
    end select

    return

  end subroutine vdw_opt

  subroutine coul_opt(m,nprcoul)
    !***************************************************************************************
    ! Flags do potencial eletrostatico                                                     *
    !***************************************************************************************
    implicit none

    integer nprcoul,m

    select case(m)
    case(1)
       nprcoul=1
    case(2)
       nprcoul=1
    case(3)
       nprcoul=1
    end select

    return

  end subroutine coul_opt

  subroutine dihedral_opt(char,m,nprtors)
    !***************************************************************************************
    ! Flags do potencial intramolecular de torsao                                          *
    !***************************************************************************************
    implicit none

    integer nprtors,m
    character(7) char

    select case(char)
    case('harm')
       m=1
       nprtors=2
    case('hrm2')
       m=2
      nprtors=2
    case('ryck')
       m=3
       nprtors=7
    end select

    return

  end subroutine dihedral_opt

  subroutine bend_opt(char,m,nprbend)
    !***************************************************************************************
    ! Flags do potencial intramolecular de deformacao                                      *
    !***************************************************************************************
    implicit none

    integer nprbend,m
    character(7) char

    select case(char)
    case('harm')
       m=1
       nprbend=2
    case('amber')
       m=2
       nprbend=2
    end select

    return

  end subroutine bend_opt

  subroutine bonds_opt(char,m,nprpar)
    !***************************************************************************************
    ! Flags do potencial intramolecular de estiramento                                     *
    !***************************************************************************************
    implicit none

    integer nprpar,m
    character(7) char

    select case(char)
    case('mors')
       m=1
       nprpar=3
    case('harm')
       m=2
       nprpar=2
    end select

    return

  end subroutine bonds_opt

  subroutine atomic_mass
    !***************************************************************************************
    ! Massas atomicas                                                                      *
    !***************************************************************************************
    implicit none

    integer i

    !-definindo massa atomica

    do i=1,natom
       select case(idna(i))
       case(1)
          mass(i)=1.0079400d0
       case(6)
          mass(i)=12.010700d0
       case(7)
          mass(i)=14.006700d0
       case(8)
          mass(i)=15.999400d0
       case(15)
          mass(i)=30.973762d0
       case(18)
          mass(i)=39.948000d0
       case(80)
          mass(i)=200.59000d0
       end select
    end do

    !-definindo menor e maior massa atomica

    massmax=0.d0
    do i=1,natom
       massmax=max(massmax,mass(i))
    end do

    massmin=massmax
    do i=1,natom
       massmin=min(massmin,mass(i))
    end do

    return

  end subroutine atomic_mass

end module input
