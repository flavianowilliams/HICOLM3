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
      integer ats
      real(kind=4) rint,rend,az
c      character*8 ats
c
      ats=1
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
