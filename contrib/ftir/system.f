!     
!     MIT License
!     
!     Copyright (c) 2020 flavianowilliams
!     
!     Permission is hereby granted, free of charge, to any person obtaining a copy
!     of this software and associated documentation files (the "Software"), to deal
!     in the Software without restriction, including without limitation the rights
!     to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
!     copies of the Software, and to permit persons to whom the Software is
!     furnished to do so, subject to the following conditions:
!     
!     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!     IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!     FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
!     AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!     LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
!     OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
!     SOFTWARE.
!     
      module system
**********************************************************************
*     Leitura das variaveis atomicas e parametros de entrada         *
*     Flaviano Williams Fernandes, 20 de agosto de 2014              *
**********************************************************************
      use sizes
      use error
c
      implicit none
c
      integer qmol,nmolec(molmax),qatom(molmax)
      integer opm,opx,ddw,d0w,dtw,bnd,flx,opsmth
      integer atmax,tmax,timestep(fmax),keyi(5)
      integer oph(2,10,3),natk(2,10,3),opz(2)
      real(kind=8) strcut,qat(nmax),rmin,rmax,mmol(nmax)
      real(kind=8) r(fmax,nmax,iz),v(fmax,nmax,iz),a(fmax,nmax,iz)
      real(kind=8) l(fmax,iz),keyf(8),dtime,dsmth
      character*5 at(nmax),att
c
      save qmol,nmolec,qatom
      save atmax,tmax,r,v,a,at,l
      save dtime,timestep,keyf,keyi
      save opm,oph,bnd,flx,strcut,natk,opz,qat
      save d0w,dtw,ddw
      save opx,rmin,rmax,att
      save opsmth,dsmth
c
      contains
c
      subroutine input
c
      implicit none
c
      integer i
c
      write(*,*)'Quantidade de especies'
      read(*,*)qmol
c      qmol=1

      write(*,*)'Quantidade de atomos e moleculas por especie'
c
      do i=1,qmol
         read(*,*)qatom(i),nmolec(i)
      end do
c      qatom(1)=3
c      nmolec(1)=309
c
      atmax=0
      do i=1,qmol
         atmax=atmax+qatom(i)*nmolec(i)
      end do
c
      call opinit()
c
      write(iwrt,'(10x,14a3)')('ITA',i=1,14)
      write(iwrt,'(7x,a43)')'Calculo do sistema solvente-substrato'
      write(iwrt,'(10x,14a3)')('ITA',i=1,14)
      write(iwrt,*)
      write(iwrt,*)'Parametros de entrada:'
      write(iwrt,*)
      write(iwrt,*)'Quantidade de especies:',qmol
      write(iwrt,*)
      write(iwrt,'(a7,a15,a15,a15)')
     1     'Especie','atom/molec','moleculas'
      write(iwrt,*)('-',i=1,35)

      do i=1,qmol
         write(iwrt,'(i4,i13,i18,f17.4)')i,qatom(i),nmolec(i)
      end do

      write(iwrt,*)('-',i=1,35)
c
      return

      end subroutine input
c
      subroutine coord(check,w)
c
      implicit none
c
      integer w,i,j,k,h,p,s,g,check
      real(8) dtt
      character*5 lixc
c
      read(ird,100,end=1)timestep(w),dtt,(l(w,i),i=1,3)
c
      do i=1,atmax
         read(ird2,200)at(i),mmol(i),qat(i),(r(w,i,j),j=1,iz),
     1        (a(w,i,j),j=1,iz),(v(w,i,j),j=1,iz)
      end do
c
      if(w.eq.1)then
         dtime=dtt
      elseif(w.eq.2)then
         dtime=dtt-dtime
      end if
c
      tmax=w
c
      goto 2
c
 1    check=1

 2    return
c
 100  format(1x,i12,1x,e12.4,119x,3(e12.4,1x))
 200  format(51x,a5,1x,11(e12.4,1x))
c
      end subroutine coord
c
      subroutine opinit()
c
      implicit none
c
      integer i,j
c--------------------------------------------------------
c-atribuindo valores iniciais
      att='OW'
      dsmth=0.
c--------------------------------------------------------
      write(*,*)'Calculo ACF!'
      write(*,*)'Molecula de referencia:'
      read(*,*)opm
c      opm=1
c
      write(*,*)'Controle do calculo ACF:'
      write(*,*)'Inicial, amostras, intervalo amostral:'
      read(*,*)d0w,dtw,ddw
c      d0w=1
c      dtw=1
c      ddw=200
c
      write(*,*)'Range:'
      write(*,*)'None?         -> 1'
      write(*,*)'Eixo z?       -> 2'
      write(*,*)'Pesonalizado? -> 3'
      read(*,*)opx
c      opx=1
c
      if(opx.eq.2)then
         write(*,*)'Alcance minimo e maximo:'
         read(*,*)rmin,rmax
      end if
c
      if(opx.eq.3)then
         write(*,*)'Atomo de referencia e alcance min. e max.:'
         read(*,*)att,rmin,rmax
      end if
c
      write(*,*)'Estiramento pelo metodo TCF!'
      write(*,*)'Quantidade a ser analisada'
      read(*,*)bnd
c      bnd=1


      write(*,*)'Intermolecular? -> 1'
      write(*,*)'Intramolecular? -> 2'
      read(*,*)opz(1)
c      opz(1)=2
c
      write(*,*)'Moleculas e atomos da coordenada relativa:'
c
      do i=1,bnd
         write(*,*)'->',i
         if(opz(1).eq.1)then
            do j=1,2
               read(*,*)oph(1,i,j),natk(1,i,j)
            end do
         else
      read(*,*)(natk(1,i,j),j=1,2)
c            natk(1,i,1)=1
c            natk(1,i,2)=2
            do j=1,2
               oph(1,i,j)=opm
            end do
         end if
      end do
c
      write(*,*)'Deformacao pelo metodo TCF!'
      write(*,*)'Quantidade a ser analisada'
      read(*,*)flx
c      flx=1
c
      write(*,*)'Intermolecular? -> 1'
      write(*,*)'Intramolecular? -> 2'
      read(*,*)opz(2)
c      opz(2)=2
c
      do i=1,flx
         write(*,*)'->',i
         if(opz(2).eq.1)then
            write(*,*)'Moleculas e atomos da coordenada relativa'
            do j=1,3
               read(*,*)oph(2,i,j),natk(2,i,j)
            end do
         else
            write(*,*)'Atomos da coordenada relativa'
      read(*,*)(natk(2,i,j),j=1,3)
c            natk(2,i,1)=3
c            natk(2,i,2)=1
c            natk(2,i,3)=2
            do j=1,3
               oph(2,i,j)=opm
            end do
         end if
      end do
c
      write(*,*)'Distancia maxima de estiramento:'
      read(*,*)strcut
c      strcut=2.0
c
      write(*,*)'Escolha da funcao janela para o metodo de suavizacao'
      write(*,*)'1 -> retangular'
      write(*,*)'2 -> Poisson'
      write(*,*)'3 -> Hann-Poisson'
      write(*,*)'4 -> Cauchy'
      write(*,*)'5 -> Gaussiana'
      read(*,*)opsmth
c      opsmth=1
c
      if(opsmth.ne.1)then
         write(*,*)'Parametro de decaimento'
         read(*,*)dsmth
      end if
c
      return
c
      end subroutine opinit
c
      end module system
