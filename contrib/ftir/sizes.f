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
      module sizes
**********************************************************************
*     Dimensionando arrays                                           *
*     Flaviano Williams Fernandes, 20 de agosto de 2014              *
**********************************************************************
      public
c********************************************
c-fmax:   limite de frames                  * 
c-nmax:   limite de átomos                  *  
c-molmax: limite da quantidade de espécies  *
c-nmmax:  limite de moléculas               *
c-iz:     dimensão do espaço vetorial       *
c-atmmax: limite de átomos por molécula     *
c********************************************
      integer fmax,nmax,iz,molmax,nmmax,atmmax
c
      parameter (fmax=2500)
      parameter (nmax=6000)
      parameter (molmax=4)
      parameter (nmmax=2000)
      parameter (atmmax=200)
      parameter (iz=3)
c
      integer ird,irdd,iwrx,iwrt,iwtt,iwrz,iwrh,ird2
      parameter (ird=1,irdd=2,iwrx=3,iwtt=4,iwrz=8,iwrh=6,iwrt=7,ird2=9)
c     
      end module sizes
