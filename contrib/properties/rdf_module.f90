!
!MIT License
!
!Copyright (c) 2020 flavianowilliams
!
!Permission is hereby granted, free of charge, to any person obtaining a copy
!of this software and associated documentation files (the "Software"), to deal
!in the Software without restriction, including without limitation the rights
!to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
!copies of the Software, and to permit persons to whom the Software is
!furnished to do so, subject to the following conditions:
!
!THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
!AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
!OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
!SOFTWARE.
!
module rdf_module

  implicit none

  integer nk,nkmx,nmolec,molectt,nat,nstp
  integer, allocatable :: idna(:)
  integer qmolec(10),qatom(10),namoltt(100000)
  character(2) spcat(2)
  character(2), allocatable :: atp(:)

  parameter(nkmx=1000000)

  real(8) rdfcut,drdfcut,a,b,c
  real(8), allocatable :: gr(:)
  real(8), allocatable :: v(:,:,:)
  real(8), allocatable :: xa(:),ya(:),za(:)

  contains

  subroutine rdf_prepare(nat,nstp)

    implicit none

    integer i,j,k,nat,nstp
    real(8) lxf
    character lxc

    read(2,*)
    read(2,*)
    read(2,*)
    read(2,*)lxc,lxc,lxc,nat,lxf,lxf,lxf,lxf,nstp
    read(2,*)
    read(2,*)

    allocate (xa(nat),ya(nat),za(nat))
    allocate (idna(nat),atp(nat))
    allocate (gr(nkmx))
    allocate (v(nstp,3,3))

    drdfcut=0.025d0

    write(*,*)'Types of molecules (qty):'
    read(*,*)nmolec
!    nmolec=1
    write(*,*)'Number of molecules per type:'
    read(*,*)(qmolec(i),i=1,nmolec)
!    qmolec(1)=309
    write(*,*)'Number of atoms per molecule per type:'
    read(*,*)(qatom(i),i=1,nmolec)
!    qatom(1)=3
    write(*,*)'Choose two species for the RDF calculus:'
    read(*,*)spcat(1),spcat(2)
!    spcat(1)='OW'
!    spcat(2)='HW'

    molectt=1
    do i=1,nmolec
       do j=1,qmolec(i)
          namoltt(molectt)=qatom(i)
          molectt=molectt+1
       end do
    end do

    molectt=molectt-1

    do k=1,nkmx
       gr(k)=0.d0
    end do

!    write(*,*)'R_cutoff','dR_cutoff'
!    read(*,*)rdfcut,drdfcut

!    nk=int(rdfcut/drdfcut)

    do i=1,6
       read(1,*)
    end do

    read(1,*)(lxf,k=1,2),((v(i,j,k),k=1,3),j=1,3),a,b,c

    rdfcut=0.5d0*min(a,min(b,c))

    nk=int(rdfcut/drdfcut)

    do i=2,nstp
       read(1,*)(lxf,k=1,2),((v(i,j,k),k=1,3),j=1,3),a,b,c
       rdfcut=min(rdfcut,0.5d0*min(a,min(b,c)))
       nk=max(nk,int(rdfcut/drdfcut))
    end do

    close(1)

    return

  end subroutine rdf_prepare

  subroutine rdf_calc(nstp,nat)

    implicit none

    integer nstp,nat,j,k,i,spct(3)
    real(8) lxf

    call rdf_prepare(nat,nstp)

    spct(1)=0
    spct(2)=0
    do i=1,nat
       read(2,*)(lxf,k=1,5),atp(i),lxf,lxf,xa(i),ya(i),za(i)
       if(atp(i).eq.spcat(1))spct(1)=spct(1)+1
       if(atp(i).eq.spcat(2))spct(2)=spct(2)+1
    end do

    call rdf(1)

    do i=2,nstp
       do j=1,nat
          read(2,*)(lxf,k=1,5),atp(j),lxf,lxf,xa(j),ya(j),za(j)
       end do
       call rdf(i)
    end do

    call rdf_final(nstp,spct)

    return

  end subroutine rdf_calc

  subroutine rdf(stp)

    implicit none

    integer i,j,l,ii,jj,nx,nxx,s,ss,stp
    real(8) xvz,yvz,zvz,dr,rr,drr

    drr=drdfcut

    rr=0.d0
    do l=1,nk
       s=0
       do i=1,molectt
          do j=1,namoltt(i)
             nx=s+j
             ss=0
             do ii=i+1,molectt
                do jj=1,namoltt(ii)
                   nxx=ss+jj+(s+namoltt(i))
                   if(atp(nx).eq.spcat(1).and.atp(nxx).eq.spcat(2))then
                      call mic(stp,nx,nxx,xvz,yvz,zvz)
                      dr=sqrt(xvz**2+yvz**2+zvz**2)
                      if(rr.gt.(dr-0.5d0*drr))then
                         if(rr.le.(dr+0.5d0*drr))gr(l)=gr(l)+2.d0
                      end if
                   end if
                end do
                ss=ss+namoltt(ii)
             end do
          end do
          s=s+namoltt(i)
       end do
       rr=rr+drr
    end do

    return

  end subroutine rdf

  subroutine rdf_final(nstp,spct)

    implicit none

    integer k,nstp,spct(3)
    real(8) rr,drr,vol

    write(3,'(a26,2a5)')'# RDF calculus of species:',spcat(1),spcat(2)

    drr=drdfcut

    rr=0.0d0
    do k=1,nk
       vol=4*3.14d0*((rr+0.5d0*drr)**3-(rr-0.5d0*drr)**3)/3.d0
       gr(k)=gr(k)*a*b*c/(vol*nstp*spct(1)*spct(2))
       rr=rr+drr
    end do

    rr=0.0d0
    do k=1,nk
       write(3,*)rr,gr(k)
       rr=rr+drr
    end do

    return

  end subroutine rdf_final

  subroutine mic(stp,i,j,xvz,yvz,zvz)
    !**************************************************************************
    !subrotina responsavel por aplicar a tecnica minimum image convention     *
    !**************************************************************************
    implicit none

    integer i,j,stp
    real(8) xx,yy,zz,xvz,yvz,zvz

    xx=v(stp,1,1)+v(stp,2,1)+v(stp,3,1)
    yy=v(stp,1,2)+v(stp,2,2)+v(stp,3,2)
    zz=v(stp,1,3)+v(stp,2,3)+v(stp,3,3)

    xvz=(xa(j)-xa(i))-xx*int(2.d0*(xa(j)-xa(i))/xx)
    yvz=(ya(j)-ya(i))-yy*int(2.d0*(ya(j)-ya(i))/yy)
    zvz=(za(j)-za(i))-zz*int(2.d0*(za(j)-za(i))/zz)

    return

  end subroutine mic

end module rdf_module

