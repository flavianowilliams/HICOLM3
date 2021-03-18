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
      subroutine prbwrt
c
      implicit none
c
      integer maxn,mimax,i,strmax,nmax,flxmax
c
      parameter (maxn=500)
c
      integer ird,irdd,iwrx,iwrt,iwtt,iwrz,iwrh
      real(kind=4) dmi,mi,strmx,dstr,strv,dstr0,dflx0,dflx,flxmx,flxv
      real(kind=4) nn(maxn),str(maxn),flx(maxn)
c
      common/units/ ird,irdd,iwrh,iwrx,iwtt,iwrz,iwrt
      common/dipavecalc/ mimax,dmi,nn
      common/tcfstrdata/ strmax,strmx,dstr0,str
      common/tcfflexdata/ flxmax,flxmx,dflx0,flx
c-----------------------------------------------------
      nmax=min(maxn,min(mimax,min(strmax,flxmax)))
c
      dstr=(strmx-dstr0)/nmax
      dflx=(flxmx-dflx0)/nmax
c
      mi=0.
      strv=0.+dstr0
      flxv=0.+dflx0
      do i=1,nmax
         write(iwrx,'(3(a4,4x,f7.4,f12.4))')
     1        '1',mi,nn(i),'2',strv,str(i),'3',flxv,flx(i)
         mi=mi+dmi
         strv=strv+dstr
         flxv=flxv+dflx
      end do
c-----------------------------------------------------
      return
c
      end subroutine prbwrt
