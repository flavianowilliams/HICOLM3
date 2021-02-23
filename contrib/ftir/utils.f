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
      module utils
**********************************************************************
*     Subrotinas úteis para o pacote FTIR-Class                      *
*     Flaviano Williams Fernandes, 20 de agosto de 2014              *
**********************************************************************
      use sizes
      use system
c
      contains
c
      subroutine smooth(www,dni,nimax,smed)
c
      implicit none
c
      integer i,op,www
      real(kind=8) smed(fmax),ni,dni,nimax,pi,phi
      real(kind=8) d
c
      op=opsmth
      d=dsmth
c
      pi=acos(-1.d0)
c
      ni=0.+dni
      select case(op)
      case(1)                   !retangular
         do i=2,www
            smed(i)=smed(i)*1.0
         end do
      case(2)                   !Poisson
         do i=2,www
            smed(i)=smed(i)*exp(-2.*d*abs(ni)/nimax)
            ni=ni+dni
         end do
      case(3)                   !Hann-Poisson
         do i=2,www
            phi=pi*ni/nimax
            smed(i)=smed(i)*.5*(1+cos(2.*phi))*exp(-2.*d*abs(ni)/nimax)
         end do
      case(4)                   !Cauchy
         do i=2,www
            smed(i)=smed(i)/(1+(2.*d*ni/nimax)**2)
         end do
      case(5)                   !Gaussiana
         do i=2,www
            smed(i)=smed(i)*exp(-.5*(2.*d*ni/nimax)**2)
         end do
      end select
c
      return
c
      end subroutine smooth
c
      subroutine periodicidade(w,img)

      implicit none

      integer w
      real(kind=8) aaa,bbb,ccc,img(iz)

      aaa=1.d0/l(w,1)
      bbb=1.d0/l(w,2)
      ccc=1.d0/l(w,3)
c     
      img(1)=img(1)-l(w,1)*nint(aaa*img(1))
      img(2)=img(2)-l(w,2)*nint(bbb*img(2))
      img(3)=img(3)-l(w,3)*nint(ccc*img(3))
c     
      return

      end subroutine periodicidade
c
      subroutine desvio(nmax,valy,dsv,md)
      
      implicit none

      integer ymax,nmax,i

      parameter (ymax=10000000)

      real (kind=8) valy(ymax)
      real (kind=8) sum,md,dsv

      md=0.
      do i=1,nmax
         md=md+valy(i)
      end do

      md=md/float(nmax)

      sum=0.
      do i=1,nmax
         sum=sum+(valy(i)-md)**2
      end do

      dsv=sqrt(sum/float(nmax-1))

      return
      
      end subroutine desvio
c
      subroutine info(t0)

      implicit none

      real(kind=8) t0,tf
      character(10) host,time
      character(8) date

      call date_and_time(Date=date)
      call date_and_time(TIME=time)
c      call hostnm(host)
      call cpu_time(tf)

      host='undefined'

      write(iwrt,*)
      write(iwrt,*)('Host: '),host
      write(iwrt,'(''Date: '',2x,a8)')date
      write(iwrt,'(''Time: '',2x,a10)')time
      write(iwrt,'(''Tempo estimado = '',f8.2,'' segundos'')')tf-t0

      return

      end subroutine info
c
      end module utils
