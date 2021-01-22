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
      module utilsrange
**********************************************************************
*     Teste dos vínculos impostos para o cálculo VDOS                *
*     Flaviano Williams Fernandes, 20 de agosto de 2014              *
**********************************************************************
      use sizes
      use system
c
      contains
c
      subroutine testrange(wop,rc,ans)
c
      implicit none
c
      integer ans,wop,j,i
      real(kind=4) rint,rend,rx,rc(iz)
c
      ans=0
      rx=0.
      select case(opx)
      case(1)
         ans=1
      case(2)
         call range(wop,rint,rend)
         if(rc(3).ge.rint.and.rc(3).lt.rend)ans=1
      case(3)
         do j=1,atmax
            if(at(j).eq.att)then
               rx=0.
               do i=1,iz
                  rx=rx+(rc(i)-r(wop,j,i))**2
               end do
               rx=sqrt(rx)
               if(rx.ge.rmin.and.rx.lt.rmax)ans=1
            end if
         end do
      end select
c
      return
c
      end subroutine testrange
c
      subroutine range(wop,rint,rend)
c
      implicit none
c
      integer c,i,wop
      real(kind=4) rint,rend,az
      character*5 ats
c
      ats='OW'
      az=-0.5*l(wop,3)
      c=0
      do i=1,atmax
         if(at(i).eq.ats)az=0.
      end do
c
      do i=1,atmax
         if(at(i).eq.ats.and.r(wop,i,3).lt.0.)then
            az=az+r(wop,i,3)
            c=c+1
         end if
      end do
c
      az=az/float(max(1,c))
      rint=rmin+az
      rend=rmax+az
c
      return
c
      end subroutine range
c
      end module utilsrange
