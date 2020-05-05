      module corr
**********************************************************************
*     Cálculo de autocorrelacao                                      *
*     Flaviano Williams Fernandes, 20 de agosto de 2014              *
**********************************************************************
      use sizes
      use utils
      use utilsrange
      use system
      use cmmolec
c
      contains
c
      subroutine acfcms(nm,smed,smedi)
c
      implicit none
c
      integer atm(nmmax)
      integer ww,w,n,www
      real(kind=4) vs,med,smed(fmax),smedi,dni,nimax
      real(kind=4) vm0(nmmax,iz),nm(fmax,molmax,nmmax,iz)
c
      smedi=0.
      vs=0.
c
      www=1
      do ww=d0w,ddw,dtw
         smed(www)=0.
         www=www+1
      end do
c
      w=0
 1    www=2
      call racfcms(d0w+w,d0w+w,n,atm,nm,vm0,med)
      smedi=smedi+med
      vs=vs+sqrt(n*med)
      do ww=d0w+1,ddw,dtw
         call racfcms(ww+w,d0w+w,n,atm,nm,vm0,med)
         smed(www)=smed(www)+med
         if((ww+w).ge.tmax)goto 2
         www=www+1
      end do
      w=w+1
      goto 1
c
 2    continue
c
      smed(1)=smedi
      do ww=1,www
         smed(ww)=smed(ww)/smedi
      end do
c
      smedi=smedi/(w+1)
c
      dni=dtw
      nimax=ddw
c
      call smooth(www,dni,nimax,smed)
c
      return
c
      end subroutine acfcms
c
      subroutine racfcms(w,wop,n,atm,nm,vm0,med)

      implicit none

      integer i,j,n,w,wop,ans,g,atm(nmmax)
      real(kind=4) vm0(nmmax,iz),med
      real(kind=4) nm(fmax,molmax,nmmax,iz),rc(iz)
c
      med=0.
c
      if(w.eq.wop)then
         n=0
         do j=1,nmolec(opm)
            do g=1,iz
               rc(g)=rm(w,opm,j,g)
            end do
            call testrange(w,rc,ans)
            if(ans.eq.1)then
               do i=1,iz
                  vm0(j,i)=nm(w,opm,j,i)
               end do
               n=n+1
               atm(n)=j
            end if
         end do
         do j=1,n
            do i=1,iz
               med=med+nm(w,opm,atm(j),i)*vm0(atm(j),i)
            end do
         end do
         med=med/float(n)
      else
         do j=1,n
            do i=1,iz
               med=med+nm(w,opm,atm(j),i)*vm0(atm(j),i)
            end do
         end do
         med=med/float(n)
      end if

      return

      end subroutine racfcms
c
      subroutine acf(atx,smed,smedi)
c
      implicit none
c
      integer ww,www,w,n,atm(nmmax)
      integer atx
      real(kind=4) vs,med,smedi
      real(kind=4) v0(nmax,iz),smed(fmax)
c      character*8 atx
c
      smedi=0.
      vs=0.
c
      www=1
      do ww=d0w,ddw,dtw
         smed(www)=0.
         www=www+1
      end do
c
      w=0
 1    www=2
      call racf(d0w+w,d0w+w,atx,n,atm,v0,med)
      smedi=smedi+med
      vs=vs+sqrt(n*med)
      do ww=d0w+1,ddw,dtw
         call racf(ww+w,d0w+w,atx,n,atm,v0,med)
         smed(www)=smed(www)+med
         if((ww+w).ge.tmax)goto 2
         www=www+1
      end do
      w=w+1
      goto 1
c
 2    smed(1)=smedi
      do ww=1,www
         smed(ww)=smed(ww)/smedi
      end do
c
      smedi=smedi/(w+1)
c
      return
c
      end subroutine acf
c
      subroutine racf(w,wop,atx,n,atm,v0,med)
c
      implicit none
c
      integer i,j,n,w,wop
      integer atx
      integer atm(nmmax),ans,g
      real(kind=4) med,rc(iz),v0(nmax,iz)
c      character*8 atx
c
      med=0.
c
      if(w.eq.wop)then
         n=0
         do j=1,atmax
            if(at(j).eq.atx)then
               do g=1,iz
                  rc(g)=r(w,j,g)
               end do
               call testrange(w,rc,ans)
               if(ans.eq.1)then
                  do i=1,iz
                     v0(j,i)=v(w,j,i)
                  end do
                  n=n+1
                  atm(n)=j
               end if
            end if
         end do
         do j=1,n
            do i=1,iz
               med=med+v(w,atm(j),i)*v0(atm(j),i)
            end do
         end do
         med=med/n
      else
         do j=1,n
            do i=1,iz
               med=med+v(w,atm(j),i)*v0(atm(j),i)
            end do
         end do
         med=med/n
      end if
c
      return
c
      end subroutine racf
c
      subroutine tcfacf(nm,nt,smed,smedi)
c
      implicit none
c
      integer dblemol
c
      parameter (dblemol=2*nmmax)
c
      integer ww,w,n
      integer atm(dblemol),nt(fmax)
      real(kind=4) v0(dblemol),smed(fmax),nm(fmax,dblemol)
      real(kind=4) med,smedi
c
      smedi=0.
c
      do ww=1,ddw
         smed(ww)=0.
      end do
c
      w=0
 1    call tcfracf(1+w,1+w,nm,nt,n,atm,v0,med)
      smedi=smedi+med
      do ww=2,ddw
         call tcfracf(ww+w,1+w,nm,nt,n,atm,v0,med)
         smed(ww)=smed(ww)+med
         if((ww+w).ge.tmax)goto 2
      end do
      w=w+1
      goto 1
c
 2    smed(1)=smedi
      do ww=1,ddw
         smed(ww)=smed(ww)/smedi
      end do
c
      smedi=smedi/(w+1)
c
      return
c
      end subroutine tcfacf
c
      subroutine tcfracf(w,wop,nm,nt,n,atm,v0,med)
c
      implicit none
c
      integer j,n,w,wop,dblemol,nt(fmax)
c
      parameter (dblemol=2*nmmax)
c
      integer atm(dblemol)
      real(kind=4) nm(fmax,dblemol),med,v0(dblemol)
c
      med=0.
c
      if(w.eq.wop)then
         n=0
         do j=1,nt(w)
            v0(j)=nm(w,j)
            n=n+1
            atm(n)=j
         end do
         do j=1,nt(w)
            med=med+nm(w,atm(j))*v0(atm(j))
         end do
         med=med/float(n)
      else
         do j=1,nt(wop)
            med=med+nm(w,atm(j))*v0(atm(j))
         end do
         med=med/float(n)
      end if
c
      return
c
      end subroutine tcfracf
c
      end module corr
