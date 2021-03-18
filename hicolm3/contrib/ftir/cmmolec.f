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
      module cmmolec
**********************************************************************
*     Remocao das projecoes imagem e calculo do centro de massa      *
*     Flaviano Williams Fernandes, 20 de agosto de 2014              *
**********************************************************************
      use error
      use sizes
      use system
      use utils
c
      implicit none
c
      integer natm(molmax,nmax)
      integer atmolec(molmax,atmmax)
      real(4) rm(fmax,molmax,nmmax,iz),vm(fmax,molmax,nmmax,iz)
c      character*8 atmolec(molmax,atmmax)
c
      save rm,vm,atmolec,natm
c
      contains
c
      subroutine molec(check,w,cmserr,cmsprec)
***************************************************************
*     Subrotina que calcula o centro de massa de cada molecula*
*     e desliga os efeitos da periodicidade devido as CC      *
***************************************************************
c
      implicit none
c
      integer i,j,h,p,w,g,t,s,c,cmserr,check
      real(kind=4) rsum(2,2,iz),rsq(iz),cmsprec
      real(kind=4) mtmol
      character*12 mark(nmax)
c
      cmsprec=5.e-3
c
      if(check.eq.1)goto 1
c--------------------------------------------------------------------
c-atribuindo valores iniciais
c
      if(w.eq.1)then
         do i=1,atmax
            mark(i)='nomodifi'
         end do
      end if
c--------------------------------------------------------------------
c-religando sítios da molecula
c
      if(w.eq.1)then
         s=0
         do h=1,qmol
            do p=1,nmolec(h)
               g=s+(p-1)*qatom(h)
               do j=1,qatom(h)
                  do i=1,iz
                     rsq(i)=r(w,g+j,i)-r(w,g+1,i)
                     if(abs(rsq(i)).gt.l(w,i))mark(g+j)='modified'
                  end do
                  call periodicidade(w,rsq)
                  do i=1,iz
                     r(w,g+j,i)=rsq(i)+r(w,g+1,i)
                  end do
               end do
            end do
            s=s+nmolec(h)*qatom(h)
         end do
      end if
c--------------------------------------------------------------------
c-desligando CC
c
      if(w.gt.1)then
         s=0
         do h=1,qmol
            do p=1,nmolec(h)
               g=s+(p-1)*qatom(h)
               do j=1,qatom(h)
                  do i=1,iz
                     rsq(i)=r(w,g+j,i)-r(1,g+j,i)
                     if(abs(rsq(i)).gt.l(w,i))mark(g+j)='modified'
                  end do
                  call periodicidade(w,rsq)
                  do i=1,iz
                     r(w,g+j,i)=rsq(i)+r(1,g+j,i)
                  end do
               end do
            end do
            s=s+nmolec(h)*qatom(h)
         end do
      end if
c--------------------------------------------------------------------
c-teste da eficiencia do calculo anterior
      s=0
      do h=1,qmol
         do p=1,nmolec(h)
            g=s+(p-1)*qatom(h)
            do j=1,qatom(h)
               do i=1,iz
                  if(abs(r(w,g+j,i)-r(w,g+1,i)).gt.l(w,i))
     1                 call erro(1,2,1)
                  if(w.gt.1)then
                     if(abs(r(w,g+j,i)-r(w-1,g+j,i)).gt.l(w,i))
     1                    write(*,*)r(w,g+j,i)
                     if(abs(r(w,g+j,i)-r(w-1,g+j,i)).gt.l(w,i))
     1                    call erro(1,2,1)
                  end if
               end do
            end do
         end do
         s=s+nmolec(h)*qatom(h)
      end do
c--------------------------------------------------------------------
c-imprimindo em REVCON
       call revcon(w,mark,2)
c--------------------------------------------------------------------
c-Calculo do centro de massa
c
      do i=1,iz
         do h=1,qmol
            do p=1,nmolec(h)
               rm(w,h,p,i)=0.
               vm(w,h,p,i)=0.
            end do
         end do
         do t=1,2
            do j=1,2
               rsum(t,j,i)=0.
            end do
         end do
      end do
c
      s=0
      do h=1,qmol
         c=0
         do p=1,nmolec(h)
            mtmol=0.
            g=s+(p-1)*qatom(h)
            do j=1,qatom(h)
               do i=1,iz
                  rm(w,h,p,i)=rm(w,h,p,i)+(r(w,g+j,i)*mmol(g+j))
                  vm(w,h,p,i)=vm(w,h,p,i)+(v(w,g+j,i)*mmol(g+j))
               end do
               atmolec(h,j)=at(g+j)
               mtmol=mtmol+mmol(g+j)
            end do
            do i=1,iz
               rsum(1,1,i)=rsum(1,1,i)+(rm(w,h,p,i))
               rsum(2,1,i)=rsum(2,1,i)+(vm(w,h,p,i))
               rm(w,h,p,i)=rm(w,h,p,i)/mtmol
               vm(w,h,p,i)=vm(w,h,p,i)/mtmol
            end do
            c=c+1
            natm(h,c)=c
         end do
         s=s+nmolec(h)*qatom(h)
      end do
c-teste do calculo do centro de massa
      do i=1,atmax
         do j=1,iz
            rsum(1,2,j)=rsum(1,2,j)+(r(w,i,j)*mmol(i))
            rsum(2,2,j)=rsum(2,2,j)+(v(w,i,j)*mmol(i))
         end do
      end do
      do t=1,1
         do j=1,iz
            if(abs(rsum(t,2,j)-rsum(t,1,j)).gt.cmsprec)then
               cmserr=cmserr+1
               exit
            end if
         end do
      end do
c--------------------------------------------------------------------
 1    return
c
      end subroutine molec
c
      subroutine revcon(w,mark,op)
c
      implicit none
c
      integer i,j,w,op
      character*12 mark(nmax)
c
      open(unit=irdd,file="REVCON.xyz",status="unknown")
c
      select case(op)
      case(1)
      write(irdd,'(a37)')'REVCON file generated by surfrag code'
      write(irdd,'(2i10)')keyi(1),keyi(2)
      write(irdd,'(3f12.4)')l(w,1),keyf(1),keyf(2)
      write(irdd,'(3f12.4)')keyf(3),l(w,2),keyf(4)
      write(irdd,'(3f12.4)')keyf(5),keyf(6),l(w,3)
      do i=1,atmax
         write(irdd,'(a8,i10,a31)')at(i),nat(i),mark(i)
         write(irdd,'(3e12.4)')(r(w,i,j),j=1,3)
         write(irdd,'(3e12.4)')(v(w,i,j),j=1,3)
      end do
      case(2)
         write(irdd,*)atmax
         write(irdd,'(a30,i5)')'REVCON file generated by frame',w
         do i=1,atmax
            write(irdd,'(a8,3f12.4,a31)')at(i),(r(w,i,j),j=1,3),mark(i)
         end do
      end select
c
      close(irdd)
c
      return
c
      end subroutine revcon
c
      end module cmmolec
