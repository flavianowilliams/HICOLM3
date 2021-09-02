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
      module dipmol
**********************************************************************
*     Cálculo dos momentos de dipolo molecular e probabilidade       *
*     Flaviano Williams Fernandes, 20 de agosto de 2014              *
**********************************************************************
      use sizes
      use system
c
      contains
c
      subroutine dipole(dm)
c
      implicit none
c
      integer nwmax
c
      parameter (nwmax=int(nmax/3))
c
      integer i,p,j,w,g,s
      real(kind=4) dm(fmax,molmax,nmmax,iz),e,d
c
      data e/-1.602176487/,d/4.803/
c-----------------------------------------------------------------------
      s=0
      do i=1,qmol
         if(i.lt.opm)s=s+qatom(i)*nmolec(i)
      end do
c
      do w=1,tmax
         do p=1,nmolec(opm)
            g=s+(p-1)*qatom(opm)
            do j=1,qatom(opm)
               do i=1,iz
                  dm(w,opm,p,i)=0.d0
               end do
            end do
         end do
         do p=1,nmolec(opm)
            g=(p-1)*qatom(opm)
            do j=1,qatom(opm)
               do i=1,iz
                  dm(w,opm,p,i)=dm(w,opm,p,i)+qat(j)*r(w,g+j,i)
c                  dm(w,opm,p,i)=dm(w,opm,p,i)+qat(j)*v(w,g+j,i)
               end do
            end do
         end do
c-----------------------------------------------------------------------
c-conversao para Debyes
         do p=1,nmolec(opm)
            do i=1,iz
c               dm(w,opm,p,i)=e*dm(w,opm,p,i)!conversao para SI
               dm(w,opm,p,i)=d*dm(w,opm,p,i)!conversao para Debyes
            end do
         end do
      end do
c-----------------------------------------------------------------------
      return
c
      end subroutine dipole
c
      subroutine dipave(dm,integx)
c      
      implicit none
c
      integer mimax,w,p,i
c
      parameter (mimax=500)
c
      integer n(mimax),nt,nmimax
      real(kind=4) dm(fmax,molmax,nmmax,iz),vec(fmax,nmmax)
      real(kind=4) int,mi,dmi,ints,s,mimx,prec,dsv,integx
      real(kind=4) f(2,mimax),integ(2)
      real(kind=4) nn(mimax)
c
      common/dipavecalc/ nmimax,dmi,nn
c
c      data prec/1.e-3/
      nmimax=mimax
c------------------------------------------------------
c-modulando os vetores momento de dipolo
      mimx=0.
c     
      do w=1,tmax
         do p=1,nmolec(opm)
            vec(w,p)=0.
            do i=1,3
               vec(w,p)=vec(w,p)+dm(w,opm,p,i)**2
            end do
            vec(w,p)=sqrt(vec(w,p))
            mimx=max(mimx,vec(w,p))
         end do
      end do
c-----------------------------------------------------
c-calculando os dipolos em intervalos de momento
      dmi=mimx/mimax
      prec=dmi*1.e-1
c     
      mi=0.
      do i=1,mimax
         n(i)=0
         do w=1,tmax
            do p=1,nmolec(opm)
               if(vec(w,p).ge.(mi-prec).and.vec(w,p).lt.(mi+prec))
     1              n(i)=n(i)+1
            end do
         end do
         mi=mi+dmi
      end do
c-----------------------------------------------------
c-contando os momentos dipolo
      nt=0
      do i=1,mimax
         nt=nt+n(i)
      end do
c-----------------------------------------------------
c-normalizando a funcao probabilidade
      do i=1,mimax
         nn(i)=(float(n(i)))/float(nt)
      end do
c
      s=0.
      do i=2,mimax-1,2
         s=s+nn(i)*4+nn(i+1)*2
      end do
      int=dmi*(nn(1)+nn(mimax)+s)/3
c
      do i=1,mimax
         nn(i)=nn(i)/int
      end do
c-teste da normalizacao
      s=0.
      do i=1,mimax-1,2
         s=s+nn(i)*4+nn(i+1)*2
      end do
      ints=dmi*(nn(1)+nn(mimax)+s)/3
c-----------------------------------------------------
c-valor medio do momento dipolo e desvio padrao
      mi=0.
      do i=1,mimax
         f(1,i)=nn(i)*mi
         f(2,i)=nn(i)*mi**2
         mi=mi+dmi
      end do
c
      do p=1,2
         s=0
         do i=2,mimax-1,2
            s=s+f(p,i)*4+f(p,i+1)*2
         end do
         integ(p)=dmi*(f(p,1)+f(p,mimax)+s)/3
      end do
c
      integx=integ(1)
      dsv=sqrt(integ(2)-integ(1)**2)
c-----------------------------------------------------
      write(*,*)opm,'->',integ(1)
c     
      write(iwrt,*)
      write(iwrt,*)'Momento de dipolo'
      write(iwrt,*)('-',i=1,79)
      write(iwrt,*)'Valor medio e desvio padrao:',integx,dsv
      write(iwrt,'(a21,i10)')'Contagem:',nt
      write(iwrt,'(a21,i10)')'Moleculas:',tmax*nmolec(opm)
      write(iwrt,'(a21,f10.4)')'Teste da probabilidade:',ints
      write(iwrt,*)('-',i=1,79)
c     
      return
c     
      end subroutine dipave
c
      end module dipmol
