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
      module ftir
**********************************************************************
*     Cálculo do espectro vibracional usando transformada de Fourier *
*     Flaviano Williams Fernandes, 20 de agosto de 2014              *
**********************************************************************
      use sizes
      use utils
      use system
      use cmmolec
      use dipmol
      use corr
      use tcf
c
      contains
c
      subroutine infrared(op,strmm,flexmm)
c
      implicit none
c
      integer pmax,dblemol
c
      parameter (pmax=4000)
      parameter (dblemol=2*nmmax)
c
      integer nimx,w,www,k,i,j,p,op,nx,nt(fmax)
      real(kind=4) str(fmax,dblemol)
      real(kind=4) smedi,smed(fmax),dm(fmax,molmax,nmmax,iz)
      real(kind=4) vel(atmmax,fmax),fint(3,pmax),integ(3,atmmax+3,pmax)
      real(kind=4) dni,ni,nimax,vol,temp,cte,light,dt,strmm,flexmm
      real(kind=4) tcfmed(fmax),flexmed(fmax),dip(fmax)
      character*5 atx(atmmax)
c
      data light/3.e-2/,cte/4.553024e5/,dni/1.e-2/
c
      www=ifix((1+(ddw-d0w))/float(dtw))
c
      vol=l(1,1)*l(1,2)*l(1,3)
      temp=133
c
      nimax=1/((timestep(2)-timestep(1))*2*dtime*light)
      nimx=int(nimax/dni)
c
      if(nimx.gt.pmax)then
         nimx=pmax
         dni=nimax/nimx
      end if
c---------------------------------------------------------
c     -calculo FT-IR pelo método ACF do momento de dipolo
c     -OBS.: constantes no SI
c---------------------------------------------------------
      write(*,*)'-> FT-IR pelo método ACF do momento de dipolo'
c
      call dipole(dm)
      call acfcms(dm,dip,smedi)
      call fourier(op,dip,fint)
c
      cte=cte*smedi/(vol*temp)
c
      ni=0.
      do k=1,nimx
         do i=1,3
            integ(i,1,k)=fint(i,k) !*cte*ni**2
            ni=ni+dni
         end do
      end do
c---------------------------------------------------------
c-calculo FT-IR por sítio molecular
c---------------------------------------------------------
      write(*,*)'-> FT-IR por sitio molecular'
c
      nx=1
      atx(1)=atmolec(opm,1)
      do i=2,qatom(opm)
         do j=1,nx
            if(atmolec(opm,i).eq.atx(j))goto 2
         end do
         nx=nx+1
         atx(nx)=atmolec(opm,i)
 2       continue
      end do
c
      if(nx.gt.12)nx=12
c
      do p=1,nx
         call acf(atx(p),smed,smedi)
         call fourier(op,smed,fint)
c
         do w=1,www
            vel(p,w)=smed(w)
         end do
c
         do k=1,nimx
            do i=1,3
               integ(i,p+3,k)=fint(i,k)
            end do
         end do
      end do
c---------------------------------------------------------
c     -calculo FT-IR do estiramento atomo-atomo
c---------------------------------------------------------
      write(*,*)'-> FT-IR do estiramento sitio-sitio'
c
      call tcfstr(str,strmm,nt)
c
      do w=1,tmax
         if(nt(w).le.1)goto 3
      end do
c
      call tcfacf(str,nt,tcfmed,smedi)
      call fourier(op,tcfmed,fint)
c
      ni=0.
      do k=1,nimx
         do i=1,3
            integ(i,2,k)=fint(i,k)
            ni=ni+dni
         end do
      end do
c---------------------------------------------------------
c     -calculo FT-IR da flexao atomo-atomo-atomo
c---------------------------------------------------------
      write(*,*)'-> FT-IR angular sitio-sitio-sitio'
c
      call tcfflex(str,flexmm,nt)
c
      do w=1,tmax
         if(nt(w).le.1)goto 4
      end do
c
      call tcfacf(str,nt,flexmed,smedi)
      call fourier(op,flexmed,fint)
c
      ni=0.
      do k=1,nimx
         do i=1,3
            integ(i,3,k)=fint(i,k)
            ni=ni+dni
         end do
      end do
c----------------------------------------------------------------------
      write(iwrz,'(6x,a1,1x,a4,11x,a6,5x,a7,4x,a9,5x,22(a6,6x))')
     1     '#','Freq','dosvel','stretch','libration',(atx(p),p=1,nx)
c
      ni=0.
      do k=1,nimx
         write(iwrz,'((f12.4,5x,22(f12.5)))')ni,(integ(3,p,k),p=1,nx+3)
         ni=ni+dni
      end do
c
      write(iwtt,'(3x,a1,1x,a4,11x,a6,5x,a7,4x,a9,5x,22(a6,6x))')
     2     '#','Time','dosvel','stretch','libration',(atx(p),p=1,nx)
c
      do w=1,www
         dt=dtime*float(timestep(w)-timestep(1))
         write(iwtt,'(f9.4,5x,22(f12.5))')
     1        dt,dip(w),tcfmed(w),flexmed(w),(vel(p,w),p=1,nx)
      end do
c
      write(iwrt,*)
      write(iwrt,*)'Cálculo de autocorrelação (TCF)'
      write(iwrt,*)('-',i=1,79)
      write(iwrt,'(2x,a12,5x,a16)')'Periodo (ps)','autocorrelacao->'
      write(iwrt,*)('-',i=1,79)
      write(iwrt,'(3x,f10.7,6x,a6,6x,a11,6x,a6,22(3x,a3))')
     1 1/(2*light*nimax),'dosvel','estiramento','flexao',(atx(p),p=1,nx)
      write(iwrt,*)('-',i=1,79)
      write(iwrt,*)
c
      goto 5
c----------------------------------------------------------------------
 3    write(iwrz,'(5x,a4,6x,a1,3x,a6,22(2x,a1,3x,a6))')
     2     'Freq',',','dosvel',(',',atx(p),p=1,nx)
c
      ni=0.
      do k=1,nimx
         write(iwrz,'((f9.4,5x,22(a1,f11.5)))')
     1        ni,',',integ(3,1,k),(',',integ(3,p,k),p=4,nx+3)
         ni=ni+dni
      end do
c
      write(iwtt,'(7x,a4,3x,a1,3x,a6,22(2x,a1,3x,a6))')
     2     'Time',',','dosvel',(',',atx(p),p=1,nx)
c
      do w=1,www
         dt=dtime*float(timestep(w)-timestep(1))
         write(iwtt,'(f11.4,19(1x,a1,f10.4))')
     1           dt,',',dip(w),(',',vel(p,w),p=1,nx)
      end do
c
      write(iwrt,*)
      write(iwrt,*)'Cálculo de autocorrelação (TCF)'
      write(iwrt,*)('-',i=1,79)
      write(iwrt,'(2x,a12,5x,a16)')'Periodo (ps)','autocorrelacao->'
      write(iwrt,*)('-',i=1,79)
      write(iwrt,'(3x,f10.7,6x,a6,22(6x,a3))')
     1     1/(2*light*nimax),'dosvel',(atx(p),p=1,nx)
      write(iwrt,*)('-',i=1,79)
      write(iwrt,*)
c
      goto 5
c------------------------------------------------------------------------
 4    write(iwrz,
     1    '(5x,a4,6x,a1,3x,a6,2x,a1,2x,a7,22(2x,a1,3x,a6))')
     2    'Freq',',','dosvel',',','stretch',(',',atx(p),p=1,nx)
c
      ni=0.
      do k=1,nimx
         write(iwrz,'((f9.4,5x,22(a1,f11.5)))')
     1        ni,',',integ(3,1,k),',',integ(3,2,k),
     2        (',',integ(3,p,k),p=4,nx+3)
         ni=ni+dni
      end do
c
      write(iwtt,
     1    '(7x,a4,3x,a1,3x,a6,2x,a1,2x,a7,22(2x,a1,3x,a6))')
     2     'Time',',','dosvel',',','stretch',(',',atx(p),p=1,nx)
c
      do w=1,www
         dt=dtime*float(timestep(w)-timestep(1))
         write(iwtt,'(f11.4,21(1x,a1,f10.4))')dt,',',dip(w),
     1        ',',tcfmed(w),(',',vel(p,w),p=1,nx)
      end do
c
      write(iwrt,*)
      write(iwrt,*)'Cálculo de autocorrelação (TCF)'
      write(iwrt,*)('-',i=1,79)
      write(iwrt,'(2x,a12,5x,a16)')'Periodo (ps)','autocorrelacao->'
      write(iwrt,*)('-',i=1,79)
      write(iwrt,'(3x,f10.7,6x,a6,6x,a11,22(6x,a3))')
     1 1/(2*light*nimax),'dosvel','estiramento',(atx(p),p=1,nx)
      write(iwrt,*)('-',i=1,79)
      write(iwrt,*)
c-----------------------------------------------------------------------
 5    return
c
      end subroutine infrared
c
      subroutine fourier(op,smed,integ)
c
      implicit none
c
      integer pmax
c
      parameter (pmax=4000)
c
      integer k,w,i,op,nimx,www
      real(kind=4) light,ni,dni,dt,integ(3,pmax),escc(2),nimax
      real(kind=4) essc(2,2),esss(2),smed(fmax),pi
c
      data light/3.e-2/,dni/1.e-2/
c
      www=ifix((1+(ddw-d0w))/float(dtw))
c
      pi=acos(-1.0)
c
      nimax=1/((timestep(2)-timestep(1))*2*dtime*light)
      nimx=int(nimax/dni)
c
      if(nimx.gt.pmax)then
         nimx=pmax
         dni=nimax/nimx
      end if
c
      call smooth(www,dni,nimax,smed)
c---------------------------------------------------------
         select case(op)
c-calculo da transformada de Fourier discreta -> DFT
      case(1)
         ni=0.
         do k=1,nimx
            escc(1)=0.
            escc(2)=0.
            do w=1,www
               dt=(timestep(w)-timestep(1))*dtime!/(www+1)
               escc(1)=escc(1)+smed(w)*cos(2*pi*light*ni*dt)
               escc(2)=escc(2)+smed(w)*sin(2*pi*light*ni*dt)
            end do
            integ(1,k)=(escc(1)/www)
            integ(2,k)=(escc(2)/www)
            integ(3,k)=sqrt(integ(1,k)**2+integ(2,k)**2)
            ni=ni+dni
         end do
c-calculo da transformada de Fourier pelo metodo Simpson
      case(2)
         ni=0.
         do k=1,nimx
            dt=(timestep(www)-timestep(1))*dtime
            essc(1,1)=smed(1)
            essc(1,2)=smed(www)*cos(2*pi*light*ni*dt)
            essc(2,1)=smed(1)
            essc(2,2)=smed(www)*sin(2*pi*light*ni*dt)
c
            escc(1)=0.
            escc(2)=0.
            do w=2,www-1
               dt=(timestep(w)-timestep(1))*dtime
               esss(1)=smed(w)*cos(2*pi*light*ni*dt)
               esss(2)=smed(w)*sin(2*pi*light*ni*dt)
               do i=1,2
                  if(mod(w,2).eq.0)esss(i)=esss(i)*4
                  if(mod(w,2).ne.0)esss(i)=esss(i)*2
                  escc(i)=escc(i)+esss(i)
               end do
            end do
            dt=((timestep(www)-timestep(1))*dtime)/www
            integ(1,k)=dt*((essc(1,1)+essc(1,2))+escc(1))/(3*acos(-1.0))
            integ(2,k)=dt*((essc(2,1)+essc(2,2))+escc(2))/(3*acos(-1.0))
            integ(3,k)=sqrt(integ(1,k)**2+integ(2,k)**2)
            ni=ni+dni
         end do
      end select
c---------------------------------------------------------
      return
c
      end subroutine fourier
c
      end module ftir
