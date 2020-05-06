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
      subroutine output(vacf,strmm,flexmm,int)

      implicit none

      integer ird,irdd,iwrx,iwrt,iwtt,iwrz,iwrh
      integer opm,ddw,bnd,flx,i
      integer oph(2,10,3),natk(2,10,3),opz(2)
      real(kind=4) dni,strcut,vacf,strmm,flexmm,int

      common/acfdata/ opm,opz,oph,ddw,bnd,flx,dni,strcut,natk
      common/units/ ird,irdd,iwrh,iwrx,iwtt,iwrz,iwrt

      write(iwrt,*)
      write(iwrt,*)
      write(iwrt,*)'Valores de saida'
      write(iwrt,*)'Molecula',opm
      write(iwrt,*)('-',i=1,41)
      write(iwrt,'(4x,a8,7x,a7,7x,a5)')'Grandeza','Unidade','Valor'
      write(iwrt,*)('-',i=1,41)
      write(iwrt,'(a13,4x,a11,4x,f5.2)')'Difusao','x10-5 cm2/s',vacf
      write(iwrt,'(a13,4x,a8,7x,f5.2)')'Dipolo','Debye',int
      write(iwrt,'(a13,4x,a9,6x,f5.2)')'Estiramento','Angstrom',strmm
      write(iwrt,'(a13,4x,a8,6x,f6.2)')'Angulo flexao','Graus',flexmm
      write(iwrt,*)('-',i=1,41)

      return

      end subroutine output
