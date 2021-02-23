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
      module printfile
**********************************************************************
*     Imprimindo probabilidades                                      *
*     Flaviano Williams Fernandes, 20 de agosto de 2014              *
**********************************************************************
c
      use dipmol
c
      contains
c
      subroutine prbwrt
c
      implicit none
c
      integer maxn,i,strmax,nmax,flxmax
c
      parameter (maxn=500)
c
      real(kind=8) mi,strmx,dstr,strv,dstr0,dflx0,dflx,flxmx,flxv
      real(kind=8) str(maxn),flx(maxn)
c
!      common/dipavecalc/ mimax,dmi,nn
      common/tcfstrdata/ strmax,strmx,dstr0,str
      common/tcfflexdata/ flxmax,flxmx,dflx0,flx
c-----------------------------------------------------
      nmax=min(maxn,min(mimax,min(strmax,flxmax)))
c
      dstr=(strmx-dstr0)/nmax
      dflx=(flxmx-dflx0)/nmax
c
      write(iwrx,'(1x,a1,16x,a6,17x,a7,19x,a4)')
     1     '#','dipole','stretch','bend'
c     
      mi=0.
      strv=0.+dstr0
      flxv=0.+dflx0
      do i=1,nmax
         write(iwrx,'(6f12.4)')
     1        mi,nn(i),strv,str(i),flxv*180.d0/3.141593,flx(i)
         mi=mi+dmi
         strv=strv+dstr
         flxv=flxv+dflx
      end do
c-----------------------------------------------------
      return
c
      end subroutine prbwrt
c
      subroutine output(strmm,flexmm,int)

      implicit none
c
      integer i
      real(kind=8) strmm,flexmm,int
c
      write(iwrt,*)
      write(iwrt,*)
      write(iwrt,*)'Valores de saida'
      write(iwrt,*)'Molecula',opm
      write(iwrt,*)('-',i=1,41)
      write(iwrt,'(4x,a8,7x,a7,7x,a5)')'Grandeza','Unidade','Valor'
      write(iwrt,*)('-',i=1,41)
      write(iwrt,'(a13,4x,a8,7x,f5.2)')'Dipolo','Debye',int
      write(iwrt,'(a13,4x,a9,6x,f5.2)')'Estiramento','Angstrom',strmm
      write(iwrt,'(a13,4x,a8,6x,f6.2)')'Angulo flexao','Graus',flexmm
      write(iwrt,*)('-',i=1,41)
c
      return
c
      end subroutine output
c
      end module printfile
