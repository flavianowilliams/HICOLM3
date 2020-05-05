!
! This file is part of the HICOLM distribution (https://github.com/flavianowilliams/HICOLM).
!
! Copyright (c) 2019 Flaviano Williams Fernandes.
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, version 3.
!
! This program is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
! General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program. If not, see <http://www.gnu.org/licenses/>.
!
module input
  !****************************************************************************************
  ! Leitura dos dados de entrada                                                          *
  ! Flaviano Williams Fernandes, 05 de dezembro de 2018                                   *
  !****************************************************************************************

  use sistema

  implicit none
  !----------------------------------------------------------------------------
  !-constantes
  !
  real(8) pi,rconv,econv,kconv,pconv,aconv,mconv,tconv,teconv,hcconv
  real(8) elconv,keconv,n0,kb,kelect
  !
  complex(8) img
  character(7) method
  !
  save rconv,econv,pconv,kconv,aconv,mconv,tconv,teconv,hcconv,elconv,keconv
  save kb,pi,img,n0
  !----------------------------------------------------------------------------
  !-variaveis de estrutura
  integer rxmx,rymx,rzmx,gsymopt,natom,reuse,idna(natmax)
  !
  real(8) a,b,c,xa(natmax),ya(natmax),za(natmax),v(3,3),volume
  real(8) vax(natmax),vay(natmax),vaz(natmax),fax(natmax),fay(natmax),faz(natmax)
  !
  save rxmx,rymx,rzmx,gsymopt,natom,reuse,idna
  !----------------------------------------------------------------------------
  !-variaveis da mecanica molecular
  integer nhist,ntrialmax,nrelax,nxmolec(molecmax),bendscnt(molecmax),bondscnt(molecmax)
  integer atp(natmax),ato(natmax),fztp(natmax),nrl,atnp(ntpmax),natnp(ntpmax)
  integer vdw(ntpmax,ntpmax),bonds(molecmax,bondmax),bends(molecmax,bendmax),tersoff
  integer tors(molecmax,torsmax),molbend(molecmax,bendmax,3),torscnt(molecmax)
  integer dstp,ndstp,xstp,nmolec,moltot,nfree,spctot,molbond(molecmax,bondmax,2)
  integer nvdw,nbonds,nbends,ncoul,ntrsff,ntors,moltors(molecmax,torsmax,4)
  !
  real(8) dtime,drmax,fmstp,text,tstat,preext,pstat,bfactor,lrmax,fzstr(6)
  real(8) parbnd(molecmax,bondmax,5),parvdw(ntpmax,ntpmax,3)
  real(8) parbend(molecmax,bendmax,4),parcoul(ntpmax,1),partrsff(ntpmax,ntpmax,16)
  real(8) partors(molecmax,torsmax,2)
  real(8) mass(natmax),massmin,massmax,rcutoff,drcutoff,lambdain,lambdafi
  !
  character(2) att
  character(7) prop,ensble
  character(7) coulop
  character(10) namemol(molecmax)
  !
  save dtime,drmax,nhist,ntrialmax,nrelax,nxmolec,fmstp,dstp,xstp,ndstp,drcutoff
  save mass,massmin,massmax,prop,namemol,nmolec,moltot,nfree,spctot,ensble,rcutoff
  save att,atp,ato,nrl,fztp,atnp,natnp,text,preext,pstat,bfactor,tstat,lrmax,fzstr
  save parbnd,parvdw,parbend,parcoul,partrsff,lambdain,lambdafi,partors
  save vdw,bonds,bends,nvdw,nbonds,nbends,ncoul,ntrsff,ntors,coulop,tersoff,tors
  save bendscnt,molbend,bondscnt,molbond,torscnt
  !----------------------------------------------------------------------------
  !-variaveis do método Tight-binding
  integer ntype,ntmolec(molecmax),nzmolec(molecmax),mmop(lmax,mmax)
  integer nskpar,nlopar,mlop(lmax),mesh,nband,diagop,bandopt,msk(nparammax),mskt(nparammax)
  integer nprllm(nparammax),npril(lmax,mmax),mpack(3),ebndtot,dband(5),mmopar(lmax)
  integer nlo,lo(lmax),mo(lmax,mmax),mlo(mmax),norb,nonst,nparam,nparamt
  integer dkrlx(nkrlxmax),nrlx(nkrlxmax),ntblx(nkrlxmax),nkrlx,nopt,lfermi
  !
  real(8) band(nbandmax,3),lcutt,rcutt,err,spr(lmax),onst0(lmax,mmax,nprllmmax)
  real(8) fermi,ksh(3),tpr0(nparammax,nprllmmax),spr0(nparammax,nprllmmax)
  real(8) dens0,kptt0(3),krlx(nkrlxmax,dkrlxmax)
  !
  logical chk(nparammax),chko(lmax,mmax),chkden
  !
  save ntype,ntmolec,nzmolec
  save nskpar,nlopar,mlop,mesh,nband,diagop,bandopt,chk,chko,chkden,lfermi
  save nprllm,npril,mpack,ebndtot,nlo,lo,mo,mlo,norb,nonst,nparam,nparamt,msk
  save band,lcutt,rcutt,err,spr,ksh,kptt0,nopt,dkrlx,nrlx,ntblx,nkrlx,krlx
  !---------------------------------------------------------------------------
contains

  subroutine entrada(t0)
    !***************************************************************************************
    ! Leitura dos parametros de entrada                                                    *
    !***************************************************************************************
    implicit none

    integer i,j
    real(8) val(20),t0,tf
    character(7) in,char
    character(9) key

    open(5,file='HICOLM.in',status='old')

    !-contagem de tempo

    call cpu_time(t0)

    !-atribuindo valores default

    call default

    !-Escolha do metodo

1   read(5,*,end=11)in

    if(in.ne.'&INIT')goto 1

    do i=1,100
       read(5,*)key
       if(key.eq.'&END')exit
       backspace(5)
       read(5,*)key,char
       if(key.eq.'method')method=char
    end do

11  rewind(5)

2   read(5,*,end=22)in

    !-Parametros de rede

    if(in.ne.'&CELL')goto 2

    do i=1,3
       read(5,*)(v(i,j),j=1,3)
    end do

    !-base

22  rewind(5)

3   read(5,*,end=33)in

    if(in.ne.'&COORD')goto 3

    do i=1,100
       read(5,*)key
       if(key.eq.'&END')exit
       backspace(5)
       read(5,*)key,val(1)
       if(key.eq.'reuse')reuse=int(val(1))
       if(key.eq.'gsymopt')gsymopt=int(val(1))
       if(key.eq.'err')err=val(1)
       if(key.eq.'replc')then
          backspace(5)
          read(5,*)key,rxmx,rymx,rzmx
       end if
       if(key.eq.'natom')then
          backspace(5)
          read(5,*)key,val(1)
          natom=int(val(1))
          do j=1,natom
             read(5,*)idna(j),xa(j),ya(j),za(j),atp(j),fztp(j)
          end do
       end if
    end do

33  rewind(5)

    select case(method)
    case('TB')
       call tb_input
!    rewind(5)
    case('MD')
       call md_input
    end select

    !-convertendo unidades de medida

    call convert

    !-contagem de tempo

    call cpu_time(tf)

    t0=tf-t0

    return

  end subroutine entrada

  subroutine tb_input
    !***************************************************************************************
    ! Leitura dos parametros de entrada do metodo Tight-binding                            *
    !***************************************************************************************
    implicit none

    integer i,j,k,m,l1,l2,m1,nsk
    real(8) val(20)
    character(7) in,char
    character(9) key
    logical loc

    !-parametros para o calculo tight-binding

1   read(5,*,end=11)in

    if(in.ne.'&TBIND')goto 1

    do i=1,100
       read(5,*)key
       if(key.eq.'&END')exit
       backspace(5)
       read(5,*)key,val(1)
       if(key.eq.'diagop')diagop=int(val(1))
       if(key.eq.'rcutt')rcutt=val(1)
       if(key.eq.'lcutt')lcutt=val(1)
       if(key.eq.'dens')then
          backspace(5)
          read(5,*)key,val(1),loc
          dens0=val(1)
          chkden=loc
       end if
       if(key.eq.'nonst')then
          backspace(5)
          read(5,*)key,val(1)
          nonst=int(val(1))
          do j=1,nonst
             read(5,*)l1,m1,loc,val(1)
             m1=m1+(l1+1)
             l1=l1+1
             chko(l1,m1)=loc
             npril(l1,m1)=int(val(1))
             backspace(5)
             read(5,*)val(1),val(1),loc,val(1),(onst0(l1,m1,k),k=1,npril(l1,m1))
          end do
       end if
       if(key.eq.'nsk')then
          backspace(5)
          read(5,*)key,val(1)
          nsk=int(val(1))
          do j=1,nsk
             read(5,*)l1,l2,m1,val(1),loc
             call labels(l1+1,l2+1,m1+1,m)
             nprllm(m)=int(val(1))
             read(5,*)(tpr0(m,k),k=1,nprllm(m))
             read(5,*)(spr0(m,k),k=1,nprllm(m))
             chk(m)=loc
          end do
       end if
    end do

    !-definindo parametros para o fitting

11  rewind(5)

2   read(5,*,end=22)in

    if(in.ne.'&OPTIM')goto 2

    do i=1,100
       read(5,*)key
       if(key.eq.'&END')exit
       backspace(5)
       read(5,*)key,val(1)
       if(key.eq.'nopt')nopt=int(val(1))
       if(key.eq.'nkrlx')then
          backspace(5)
          read(5,*)key,val(1)
          nkrlx=int(val(1))
          do j=1,nkrlx
             read(5,*)ntblx(j),nrlx(j),dkrlx(j)
             backspace(5)
             read(5,*)ntblx(j),nrlx(j),dkrlx(j),(krlx(j,k),k=1,dkrlx(j))
          end do
       end if
    end do

    !-momento angular e orbitais atomicos

22  rewind(5)

3   read(5,*,end=33)in

    if(in.ne.'&ORBIT')goto 3

    do i=1,100
       read(5,*)key
       if(key.eq.'&END')exit
       if(key.eq.'nlo')then
          backspace(5)
          read(5,*)key,val(1)
          nlo=int(val(1))
          do k=1,nlo
             read(5,*)lo(k),mlo(k)
             backspace(5)
             read(5,*)val(1),val(2),(mo(k,j),j=1,mlo(k))
          end do
       end if
    end do

    !-propriedades do metodo tight-binding

33  rewind(5)

4   read(5,*,end=44)in

    if(in.ne.'&PROPER')goto 4

    do i=1,100
       read(5,*)key
       if(key.eq.'&END')exit
       if(key.eq.'prop')then
          backspace(5)
          read(5,*)key,char
          prop=char
       end if
       if(key.eq.'bandopt')then
          backspace(5)
          read(5,*)key,val(1)
          bandopt=int(val(1))
       end if
       if(key.eq.'nband')then
          backspace(5)
          read(5,*)key,val(1)
          nband=int(val(1))
          do j=1,nband
             read(5,*)dband(j),(band(j,k),k=1,3)
          end do
       end if
       if(key.eq.'mesh')then
          backspace(5)
          read(5,*)key,val(1)
          mesh=int(val(1))
       end if
       if(key.eq.'fermi')then
          backspace(5)
          read(5,*)key,val(1)
          fermi=int(val(1))
       end if
       if(key.eq.'lfermi')then
          backspace(5)
          read(5,*)key,val(1)
          lfermi=int(val(1))
       end if
       if(key.eq.'kptt0')then
          backspace(5)
          read(5,*)key,(kptt0(j),j=1,3)
       end if
    end do

    !-calculando qtdes de parametros e orbitais atomicos

44  norb=0
    do i=1,nlo
       do j=1,mlo(i)
          norb=norb+1
       end do
    end do

    norb=norb*natom

    do i=1,nlo
       do j=1,mlo(i)
          mo(i,j)=mo(i,j)+(lo(i)+1)
       end do
       lo(i)=lo(i)+1
    end do

    nlopar=0
    do i=1,nlo
       mmopar(nlopar+1)=0
       do j=1,mlo(i)
          if(chko(lo(i),mo(i,j)).eqv..true..and.chko(lo(i),mo(i,j)).neqv..false.)&
               then
             mlop(nlopar+1)=lo(i)
             mmop(nlopar+1,mmopar(nlopar+1)+1)=mo(i,j)
             mmopar(nlopar+1)=mmopar(nlopar+1)+1
          end if
       end do
       nlopar=nlopar+1
    end do

    nskpar=0
    nparamt=0
    do i=1,nlo
       do j=i,nlo
          l1=min(lo(i),lo(j))
          l2=max(lo(i),lo(j))
          do k=1,l1
             call labels(l1,l2,k,m)
             if(m.ne.0)then
                if(chk(m).eqv..true..and.chk(m).neqv..false.)then
                   msk(nskpar+1)=m
                   nskpar=nskpar+1
                end if
                mskt(nparamt+1)=m
                nparamt=nparamt+1
             end if
          end do
       end do
    end do

    nparam=0

    if(chkden.eqv..true.)nparam=nparam+1

    do i=1,nlopar
       do j=1,mmopar(i)
          nparam=nparam+npril(mlop(i),mmop(i,j))
       end do
    end do
    do i=1,nskpar
       nparam=nparam+2*nprllm(msk(i))
    end do

    return

  end subroutine tb_input

  subroutine md_input
    !***************************************************************************************
    ! Leitura dos parametros de entrada da Dinamica molecular                              *
    !***************************************************************************************
    implicit none

    integer i,j,k,g,p,m,numt,natfx
    real(8) val(20)
    character(7) in,char
    character(9) key

    !-parametros de dinâmica molecular

1   read(5,*,end=11)in

    if(in.ne.'&MD')goto 1

    do i=1,100
       read(5,*)key
       if(key.eq.'&END')exit
       if(key.eq.'nhist')then
          backspace(5)
          read(5,*)key,val(1)
          nhist=int(val(1))
       end if
       if(key.eq.'ntrialmax')then
          backspace(5)
          read(5,*)key,val(1)
          ntrialmax=int(val(1))
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
          read(5,*)key,(val(j),j=1,5)
          drmax=val(1)
          fmstp=val(2)
          dstp=int(val(3))
          ndstp=int(val(4))
          xstp=int(val(5))
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
             read(5,*)key,char,val(1)
             tstat=val(1)
             ensble=char
          elseif(char.eq.'npt')then
             backspace(5)
             read(5,*)key,char,val(1),val(2),val(3)
             tstat=val(1)
             pstat=val(2)
             bfactor=val(3)
             ensble=char
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

    !-Campo de Forca

11  rewind(5)

2   read(5,*,end=22)in

    if(in.ne.'&FORCE')goto 2

    do i=1,100
       read(5,*)key
       if(key.eq.'&END')exit
       if(key.eq.'ebndtot')then
          backspace(5)
          read(5,*)key,val(1)
          ebndtot=int(val(1))
       end if
       if(key.eq.'vdw')then
          backspace(5)
          read(5,*)key,nvdw
          do j=1,nvdw
             read(5,*)val(1),val(2),char
             call vdw_opt(char,m,numt)
             backspace(5)
             read(5,*)val(1),val(2),char,(val(k),k=3,2+numt)
             do k=1,numt
                parvdw(int(val(1)),int(val(2)),k)=val(k+2)
                parvdw(int(val(2)),int(val(1)),k)=val(k+2)
             end do
             vdw(int(val(1)),int(val(2)))=m
             vdw(int(val(2)),int(val(1)))=m
          end do
       end if
       if(key.eq.'elect')then
          backspace(5)
          read(5,*)key,char,ncoul
          do j=1,ncoul
             call coul_opt(char,numt)
             read(5,*)(val(k),k=1,numt+1)
             do k=1,numt
                parcoul(int(val(1)),k)=val(k+1)
             end do
          end do
          coulop=char
       end if
       if(key.eq.'tersoff')then
          backspace(5)
          read(5,*)key,char,ntrsff
          call trsff_opt(char,m,numt)
          do j=1,ntrsff
             read(5,*)(val(k),k=1,numt+2)
             do k=1,numt
                partrsff(int(val(1)),int(val(2)),k)=val(k+2)
                partrsff(int(val(2)),int(val(1)),k)=val(k+2)
             end do
             tersoff=m
          end do
       end if
    end do

    !-Descricao das moleculas por especie

22  rewind(5)

3   read(5,*,end=33)in

    if(in.ne.'&MOLEC')goto 3

    do i=1,100
       read(5,*)key
       if(key.eq.'&END')exit
       if(key.eq.'molecs')then
          backspace(5)
          read(5,*)key,nmolec
          do j=1,nmolec
             read(5,*)namemol(j),ntmolec(j),nxmolec(j)
             do g=1,3
                read(5,*)key
                if(key.eq.'&END')goto 33
                if(key.eq.'bends')then
                   backspace(5)
                   read(5,*)key,bendscnt(j)
                   do k=1,bendscnt(j)
                      read(5,*)val(1),val(2),val(3),char
                      call bend_opt(char,m,numt)
                      backspace(5)
                      read(5,*)val(1),val(2),val(3),char,(val(p),p=4,3+numt)
                      do p=1,numt
                         parbend(j,k,p)=val(p+3)
                      end do
                      bends(j,k)=m
                      molbend(j,k,1)=int(val(1))
                      molbend(j,k,2)=int(val(2))
                      molbend(j,k,3)=int(val(3))
                   end do
                elseif(key.eq.'bonds')then
                   backspace(5)
                   read(5,*)key,bondscnt(j)
                   do k=1,bondscnt(j)
                      read(5,*)val(1),val(2),char
                      call bonds_opt(char,m,numt)
                      backspace(5)
                      read(5,*)val(1),val(2),char,(val(p),p=3,2+numt)
                      do p=1,numt
                         parbnd(j,k,p)=val(p+2)
                      end do
                      bonds(j,k)=m
                      molbond(j,k,1)=int(val(1))
                      molbond(j,k,2)=int(val(2))
                   end do
                elseif(key.eq.'dihedral')then
                   backspace(5)
                   read(5,*)key,torscnt(i)
                   do k=1,torscnt(i)
                      read(5,*)val(1),val(2),val(3),val(4),char
                      call dihedral_opt(char,m,numt)
                      backspace(5)
                      read(5,*)val(1),val(2),val(3),val(4),char,(val(p),p=5,4+numt)
                      do p=1,numt
                         partors(j,k,p)=val(p+4)
                      end do
                      tors(j,k)=m
                      moltors(j,k,1)=int(val(1))
                      moltors(j,k,2)=int(val(2))
                      moltors(j,k,3)=int(val(3))
                      moltors(j,k,4)=int(val(4))
                   end do
                end if
             end do
             backspace(5)
          end do
       end if
    end do

    !-calculando moleculas totais

33  moltot=1
    do i=1,nmolec
       do j=1,ntmolec(i)
          nzmolec(moltot)=nxmolec(i)
          moltot=moltot+1
       end do
    end do

    moltot=moltot-1

    !-calculando parametros intramoleculares totais

    nbends=0
    nbonds=0
    ntors=0
    do i=1,nmolec
       do j=1,ntmolec(i)
          nbends=nbends+bendscnt(i)
          nbonds=nbonds+bondscnt(i)
          ntors=ntors+torscnt(i)
       end do
    end do

    !-calculando total de especies

    atnp(1)=atp(1)

    spctot=1
    do i=2,natom
       do j=1,spctot
          if(atp(i).eq.atnp(j))goto 553
       end do
      atnp(spctot+1)=atp(i)
      spctot=spctot+1
553   continue
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

    !-redefinindo parametros de rede e posicoes atomicas

    if(reuse.eq.1)call frame

    !-definindo massa atomica

    call atomic_mass

  return

  end subroutine md_input

  subroutine default
    !***************************************************************************************
    ! Valores padrão                                                                       *
    !***************************************************************************************

    implicit none

    integer i,j,k,l1,l2

    prop='NONE'
    mesh=1

    nhist=1
    text=273.d0
    tstat=0.5d0
    pstat=0.5d0
    bfactor=4.679d-5
    preext=1.d0
    ntrialmax=1
    nrelax=1
    drmax=0.d0
    fmstp=0.d0
    ndstp=0
    dstp=0
    xstp=0
    lrmax=0.d0
    nrl=2

    chkden=.false.

    lambdain=0.d0
    lambdafi=1.d0

    ensble='nve'

    rcutoff=0.d0
    drcutoff=0.d0

    reuse=0

    gsymopt=0

    nbonds=0
    nbends=0
    nvdw=0
    ncoul=0
    ntrsff=0
    ntors=0

    ebndtot=0

    coulop='coul'

    nmolec=0
    natom=0

    rxmx=1
    rymx=1
    rzmx=1

    tersoff=0

    do i=1,molecmax
       do j=1,torsmax
          do k=1,2
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
          do k=1,16
             partrsff(i,j,k)=0.d0
          end do
          do k=1,3
             parvdw(i,j,k)=0.d0
          end do
       end do
       parcoul(i,1)=0.d0
    end do

    do i=1,3
       kptt0(i)=0.d0
       mpack(i)=8
       ksh(i)=0.d0
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

    fermi=0.d0

    lfermi=1

    err=0.3d0

    bandopt=1

    diagop=1

    ntype=1

    dens0=0.d0

    lcutt=0.15d0

    rcutt=4.d0

    nlo=3
    do i=1,nlo
       lo(i)=i
       do j=1,lo(i)
          mo(i,j)=j
       end do
    end do

    norb=0
    do i=1,nlo
       do j=1,lo(i)
          norb=norb+2*j-1
       end do
    end do

    norb=norb*natmax

    nskpar=0
    do i=1,nlo
       do j=i,nlo
          nskpar=nskpar+max(lo(i),lo(j))
       end do
    end do

    do i=1,nlo
       do j=i+1,nlo
          l1=min(lo(i),lo(j))
          l2=max(lo(i),lo(j))
          l1=l1-1
          l2=l2-1
          if(l1.eq.0.and.l2.eq.1)nskpar=nskpar-1
          if(l1.eq.0.and.l2.eq.2)nskpar=nskpar-2
          if(l1.eq.1.and.l2.eq.2)nskpar=nskpar-1
       end do
       do j=1,2*i+1
          spr(i)=1.d0
          chko(i,j)=.false.
          npril(i,j)=1
          do k=1,npril(i,j)
             onst0(i,j,k)=0.d0
          end do
       end do
       k=k+2*i-1
    end do

    nparam=2*nskpar+nlo+1

    do i=1,nskpar
       do j=1,nprllmmax
          tpr0(i,j)=0.d0
          spr0(i,j)=0.d0
       end do
       nprllm(i)=1
       chk(i)=.false.
    end do

    do i=1,6
       fzstr(i)=1.d0
    end do

    !-parametros da Dinâmica Molecular

    dtime=0.0d0

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

  subroutine labels(l1,l2,m1,i)
    !***************************************************************************************
    ! Parametros Slater-Koster                                                             *
    !***************************************************************************************

    implicit none

    integer l1,l2,la,lb,m1,i,ma

    i=0

    la=min(l1,l2)
    lb=max(l1,l2)

    la=la-1
    lb=lb-1
    ma=m1-1

    if(la.eq.0)then
       if(lb.eq.0)then
          if(ma.eq.0)i=1
       elseif(lb.eq.1)then
          if(ma.eq.0)i=2
       elseif(lb.eq.2)then
          if(ma.eq.0)i=3
       end if
    elseif(la.eq.1)then
       if(lb.eq.1)then
          if(ma.eq.0)i=4
          if(ma.eq.1)i=5
       elseif(lb.eq.2)then
          if(ma.eq.0)i=6
          if(ma.eq.1)i=7
       end if
    elseif(la.eq.2)then
       if(lb.eq.2)then
          if(ma.eq.0)i=8
          if(ma.eq.1)i=9
          if(ma.eq.2)i=10
       end if
    end if

    return

  end subroutine labels

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

  subroutine frame
    !***************************************************************************************
    ! Leitura do ficheiro de entrada:                                                      *
    ! - Coordenadas espaciais;                                                             *
    ! - Velocidades;                                                                       *
    ! - Forca.                                                                             *
    !***************************************************************************************
    implicit none

    integer i,j

    open(1,file='HICOLM.XSF',status='old')

    do i=1,13
       read(1,*)
    end do

    do i=1,3
       read(1,'(3(3x,f14.8))')(v(i,j),j=1,3)
    end do

    read(1,*)
    read(1,*)natom

    do i=1,natom
       read(1,'(i5,3f14.8,2(2x,3f14.8))') &
            idna(i),xa(i),ya(i),za(i),fax(i),fay(i),faz(i),vax(i),vay(i),vaz(i)
    end do

    !-fechando arquivo XSF com as coordenadas atomicas

    close(1)

    return

  end subroutine frame

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

  subroutine trsff_opt(char,m,nprtrsff)
    !***************************************************************************************
    ! Flags do potencial Tersoff                                                           *
    !***************************************************************************************
    implicit none

    integer nprtrsff,m
    character(7) char

    select case(char)
    case('ters')
       nprtrsff=10
       m=1
    end select

    return

  end subroutine trsff_opt

  subroutine coul_opt(char,nprcoul)
    !***************************************************************************************
    ! Flags do potencial eletrostatico                                                     *
    !***************************************************************************************
    implicit none

    integer nprcoul
    character(7) char

    select case(char)
    case('coul')
       nprcoul=1
    case('escl')
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
          mass(i)=1.00794d0
       case(6)
          mass(i)=12.0107d0
       case(7)
          mass(i)=14.0067d0
       case(8)
          mass(i)=15.9994d0
       case(18)
          mass(i)=39.9480d0
       case(80)
          mass(i)=200.5900d0
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
