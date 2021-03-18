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
      module tcf
**********************************************************************
*     Cálculo TCF dos estiramento e deformação                       *
*     Flaviano Williams Fernandes, 20 de agosto de 2014              *
**********************************************************************
      use sizes
      use system
      use cmmolec
      use utilsrange
c
      contains
c
      subroutine tcfstr(sqtstr,strmm,n)
c
      implicit none
c
      integer strmax,dblemol,t0,nt,pp,ppp,ans,ans1(nmax),ans2(nmax)
      integer n(fmax)
c
      parameter (dblemol=2*nmmax,strmax=500)
c
      integer w,j,i,k,p,g,gg,t,s,ss,strmmax
      real(kind=4) sqtstr(fmax,dblemol),str(dblemol,iz)
      real(kind=4) strmm,strm,strmx,dsv,dstr0
      real(kind=4) nn(strmax),rc1(iz),rc2(iz),rx(2,nmax,iz)
c     
      common/tcfstrdata/ strmmax,strmx,dstr0,nn
c
      write(*,*)'-> Calculo do Estiramento'
c
      strmmax=strmax
c---------------------------------------------------------
c-calculo da coordenada relativa
      nt=0
      dstr0=0.
      strmx=0.
      do w=1,tmax
         n(w)=1
         select case(opz(1))
      case(1)
         do j=1,bnd
            s=0
            ss=0
            do i=1,qmol
               if(i.lt.oph(1,j,1))s=s+qatom(i)*nmolec(i)
               if(i.lt.oph(1,j,2))ss=ss+qatom(i)*nmolec(i)
            end do
c---------------------------------------------------------
c-contagem de atomos na regiao z1<z<z2
            pp=1
            do t=1,nmolec(oph(1,j,1))
               g=s+(t-1)*qatom(oph(1,j,1))
               do k=1,iz
                  rc1(k)=r(w,g+natk(1,j,1),k)
               end do
               call testrange(w,rc1,ans1(pp))
               if(ans1(pp).eq.1)then
                  do k=1,iz
                     rx(1,pp,k)=r(w,g+natk(1,j,1),k)
                  end do
                  pp=pp+1
               end if
            end do
            pp=pp-1
c
            ppp=1
            do t=1,nmolec(oph(1,j,2))
               gg=ss+(t-1)*qatom(oph(1,j,2))
               do k=1,iz
                  rc2(k)=r(w,gg+natk(1,j,2),k)
               end do
               call testrange(w,rc2,ans2(ppp))
               if(ans2(ppp).eq.1)then
                  do k=1,iz
                     rx(2,ppp,k)=r(w,gg+natk(1,j,2),k)
                  end do
                  ppp=ppp+1
               end if
            end do
            ppp=ppp-1
c---------------------------------------------------------
            do p=1,pp
               t0=p+1
               if(oph(1,j,1).ne.oph(1,j,2))t0=1
               do t=t0,ppp
                  strm=0.
                  do i=1,iz
                     strm=strm+(rx(2,t,i)-rx(1,p,i))**2
                  end do
                  strm=sqrt(strm)
                  if(strm.le.strcut.and.ans1(p)*ans2(t).eq.1)then
                     do i=1,iz
                        str(n(w),i)=rx(2,t,i)-rx(1,p,i)
                     end do
                     n(w)=n(w)+1
                  end if
               end do
            end do
            continue
         end do
      case(2)
         do j=1,bnd
            do p=1,nmolec(oph(1,j,1))
               g=(p-1)*qatom(oph(1,j,1))
               do k=1,iz
                  rc1(k)=rm(w,oph(1,j,1),p,k)
               end do
               call testrange(w,rc1,ans)
               if(ans.eq.1)then
                  do i=1,iz
                     str(n(w),i)=
     1                    r(w,g+natk(1,j,2),i)-r(w,g+natk(1,j,1),i)
                  end do
                  n(w)=n(w)+1
               end if
            end do
         end do
      end select
c     
      n(w)=n(w)-1
c
      if(n(w).le.0)n(w)=1
c
      nt=nt+n(w)
c---------------------------------------------------------
c-modulando vetores
      do j=1,n(w)
         sqtstr(w,j)=0.
         do i=1,iz
            sqtstr(w,j)=sqtstr(w,j)+str(j,i)**2
         end do
         sqtstr(w,j)=sqrt(sqtstr(w,j))
         strmx=max(strmx,sqtstr(w,j))
      end do
      end do
c
      write(*,*)nt,'arrays'
c---------------------------------------------------------
c-calculo da probabilidade
      call tcfprb(strmax,dstr0,strmx,sqtstr,n,nn,strmm,dsv)
c
      write(iwrt,*)'Estiramento ->',strmm,dsv
      write(*,*)'Estiramento ->',strmm
c---------------------------------------------------------
c-recalculando arrays
      do w=1,tmax
         strm=0.
         if(n(w).gt.0)then
            do j=1,n(w)
               strm=strm+sqtstr(w,j)
            end do
            strm=strm/float(n(w))
            do j=1,n(w)
               sqtstr(w,j)=sqtstr(w,j)-strm
            end do
         end if
      end do
c---------------------------------------------------------
      return
c     
      end subroutine tcfstr
c     
      subroutine tcfflex(sqtstr,strmm,n)
      
      implicit none
c
      integer w,j,i,k,p,g,gg,ggg,z,s,ss,sss,t,q,nt
      integer dblemol,strmax,strmmax
c
      parameter (dblemol=2*nmmax,strmax=500)
c
      integer n(fmax)
      integer ans1,ans2
      real(kind=4) sqtstr(fmax,dblemol),str(2,dblemol,iz)
      real(kind=4) strmm,strm,strx,dsv,strmx,dflx0
      real(kind=4) nn(strmax),rc1(iz),rc2(iz)
c     
      common/tcfflexdata/ strmmax,strmx,dflx0,nn
c
      write(*,*)'-> Calculo da deformacao!'
c
      strmmax=strmax
c---------------------------------------------------------
c-calculo dos vetores coplanares
      nt=0
      dflx0=-1.
      strmx=1.
      do w=1,tmax
         n(w)=1
         select case(opz(2))
      case(1)
         do j=1,flx
            s=0
            ss=0
            sss=0
            do i=1,qmol
               if(i.lt.oph(2,j,1))s=s+qatom(i)*nmolec(i)
               if(i.lt.oph(2,j,2))ss=ss+qatom(i)*nmolec(i)
               if(i.lt.oph(2,j,3))sss=sss+qatom(i)*nmolec(i)
            end do
            do p=1,nmolec(oph(2,j,2))
               gg=ss+(p-1)*qatom(oph(2,j,2))
               do k=1,nmolec(oph(2,j,1))
                  g=s+(k-1)*qatom(oph(2,j,1))
                  strx=0.
                  do i=1,iz
                     strx=strx+(r(w,g+natk(2,j,1),i)
     1                    -r(w,gg+natk(2,j,2),i))**2
                  end do
                  strx=sqrt(strx)
                  if(strx.gt.strcut)goto 2
                  do t=k+1,nmolec(oph(2,j,3))
                     ggg=sss+(t-1)*qatom(oph(2,j,3))
                     strx=0.
                     do i=1,iz
                        strx=strx+(r(w,ggg+natk(2,j,3),i)
     1                       -r(w,gg+natk(2,j,2),i))**2
                     end do
                     strx=sqrt(strx)
                     if(strx.gt.strcut)goto 1
                     do i=1,iz
                        rc1(i)=rm(w,oph(1,j,1),t,i)
                        rc2(i)=rm(w,oph(1,j,2),t,i)
                     end do
                     call testrange(w,rc1,ans1)
                     call testrange(w,rc2,ans2)
                     if(ans1*ans2.eq.1)then
                        do i=1,iz
                           str(1,n(w),i)=r(w,g+natk(2,j,1),i)
     1                          -r(w,gg+natk(2,j,2),i)
                           str(2,n(w),i)=r(w,ggg+natk(2,j,3),i)
     1                          -r(w,gg+natk(2,j,2),i)
                        end do
                        n(w)=n(w)+1
                     end if
 1                   continue
                  end do
 2                continue
               end do
            end do
         end do
      case(2)
         do j=1,flx
            do p=1,nmolec(oph(2,j,1))
               g=(p-1)*qatom(oph(2,j,1))
               do k=1,iz
                  rc1(k)=rm(w,oph(1,j,1),p,k)
               end do
               call testrange(w,rc1,ans1)
               if(ans1.eq.1)then
                  do z=1,2
                     do i=1,iz
                        str(z,n(w),i)=r(w,g+natk(2,j,2*z-1),i)
     1                       -r(w,g+natk(2,j,2),i)
                     end do
                  end do
                  n(w)=n(w)+1
               end if
            end do
         end do
      end select
c
      n(w)=n(w)-1
c
      nt=nt+n(w)
c---------------------------------------------------------
c-calculo do angulo entre vetores
         do j=1,n(w)
            sqtstr(w,j)=0.
            do i=1,iz
               sqtstr(w,j)=sqtstr(w,j)+str(1,j,i)*str(2,j,i)
            end do
            do q=1,2
               strx=0.
               do i=1,iz
                  strx=strx+str(q,j,i)**2
               end do
               strx=sqrt(strx)
               sqtstr(w,j)=sqtstr(w,j)/strx
            end do
         end do
c---------------------------------------------------------
      end do
c
      write(*,*)nt,'arrays'
c---------------------------------------------------------
c-calculo do desvio padrao e valor medio
c      call desvio(n(tmax)*tmax,valy,dsv,strmm)
c
c      strmm=strmm/(float(n*tmax))
c      strmm=acos(strmm)*180/acos(-1.d0)
c      dsv=acos(dsv)*180/acos(-1.d0)
c
c---------------------------------------------------------
c-calculo da probabilidade
      call tcfprb(strmax,dflx0,strmx,sqtstr,n,nn,strmm,dsv)
c
      strmm=acos(strmm)*180/acos(-1.)
      dsv=acos(dsv)*180/acos(-1.)
c
      write(iwrt,*)'Flexao ->',strmm,dsv
      write(*,*)'Flexao ->',strmm
c---------------------------------------------------------
c-recalculando arrays
      do w=1,tmax
         strm=0.
         if(n(w).gt.0)then
            do j=1,n(w)
               strm=strm+sqtstr(w,j)
            end do
            strm=strm/float(n(w))
            do j=1,n(w)
               sqtstr(w,j)=sqtstr(w,j)-strm
            end do
         end if
      end do
c---------------------------------------------------------
      return
      
      end subroutine tcfflex
c
      subroutine tcfprb(strmax,dstr0,strmx,sqtstr,n,nn,vm,dsv)
      
      implicit none

      integer strmax,i,j,w,nt,p,dblemol
c
      parameter (dblemol=2*nmmax)
c
      integer strprb(strmax),n(fmax)
      real(kind=4) prec,dstr,strmx,strv,sm,int,vm,dsv,dstr0
      real(kind=4) sqtstr(fmax,dblemol),f(2,strmax),integ(2)
      real(kind=4) nn(strmax)
c---------------------------------------------------------
c-calculo da probabilidade
      dstr=(strmx-dstr0)/strmax
      prec=dstr*5.e-1
      strv=0.+dstr0
      do i=1,strmax
         strprb(i)=0
         do w=1,tmax
            do j=1,n(w)
               if(sqtstr(w,j).ge.(strv-prec))then
                  if(sqtstr(w,j).lt.(strv+prec))then
                     strprb(i)=strprb(i)+1
                  end if
               end if
            end do
         end do
         strv=strv+dstr
      end do
c---------------------------------------------------------
c-contando as distancias relativas
      nt=0
      do i=1,strmax
         nt=nt+strprb(i)
      end do
c---------------------------------------------------------
c-normalizando a função probabilidade
      do i=1,strmax
         nn(i)=(float(strprb(i))/float(nt))
      end do
c
      sm=0.
      do i=2,strmax-1,2
         sm=sm+nn(i)*4+nn(i+1)*2
      end do
      int=dstr*(nn(1)+nn(strmax)+sm)/3
c
      do i=1,strmax
         nn(i)=nn(i)/int
      end do
c---------------------------------------------------------
c-calculo do valor médio de estiramento e desvio padrao
c      call desvio(n(tmax)*tmax,valy,dsv,strmm)
      strv=0.+dstr0
      do i=1,strmax
         f(1,i)=nn(i)*strv
         f(2,i)=nn(i)*strv**2
         strv=strv+dstr
      end do
c
      do p=1,2
         sm=0.
         do i=2,strmax-1,2
            sm=sm+f(p,i)*4+f(p,i+1)*2
         end do
         integ(p)=dstr*(f(p,1)+f(p,strmax)+sm)/3
      end do
c
      vm=integ(1)
      dsv=sqrt(integ(2)-integ(1)**2)
c
      return
c
      end subroutine tcfprb
c
      end module tcf
