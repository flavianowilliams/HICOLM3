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
module interaction_module
  !*******************************************************************************************
  !*******************************************************************************************

  use neighbourlist_module
  use bonds_module
  use angles_module
  use coulomb_module
  use dihedrals_module
  use vanderwaals_module

  implicit none

  private
  public :: interaction

  type, extends(neighbourlist) :: interaction
     type(bonds)                   :: bnd
     type(angles)                  :: tht
     type(coulomb)                 :: coul
     type(dihedrals)               :: dih
     type(vanderwaals)             :: vdw
     real(8), private              :: enpot
     real(8), private              :: virtot
     real(8), private              :: encorr
     real(8), private              :: vircorr
     real(8), private              :: fbi(3)
     real(8), private              :: fbj(3)
     real(8), private              :: fbk(3)
     real(8), private              :: fbl(3)
   contains
     procedure :: interaction_prepare
     procedure :: set_forcefield
     procedure :: set_force2
     procedure :: set_force3
     procedure :: set_force4
     generic   :: set_force => set_force2, set_force3, set_force4
     procedure :: set_enpot
     procedure :: get_enpot
     procedure :: set_virtot
     procedure :: get_virtot
     procedure :: set_vdwcorr
     procedure :: get_encorr
     procedure :: get_vircorr
  end type interaction

contains

  subroutine interaction_prepare(this)
    implicit none
    class(interaction), intent(inout) :: this
    allocate(this%fax(this%get_natom()),this%fay(this%get_natom()),this%faz(this%get_natom()))
  end subroutine interaction_prepare

  subroutine set_forcefield(this)
    implicit none
    class(interaction), intent(inout) :: this
    integer                           :: i,j,k,l,m,n,o,ni,nj,nk,nl,nx
    real(8)                           :: xvz,yvz,zvz,dr,enpot,virtot,theta,dr1,dr2
    real(8)                           :: prm(3),drij(3),drik(3),drjk(3),drkl(3),ri(3),rj(3)
    real(8)                           :: rk(3),rl(3)
    real(8)                           :: vc1x,vc1y,vc1z,vc2x,vc2y,vc2z,phi
    character(7)                      :: ptrm
    do i=1,this%get_natom()
       this%fax(i)=0.d0
       this%fay(i)=0.d0
       this%faz(i)=0.d0
    end do
    call this%coul%coulomb_prepare&
         (this%get_coulop(),this%get_kconv(),this%get_rcutoff(),this%get_pi())
    nx=0
    enpot=0.d0
    virtot=0.d0
    do i=1,this%get_nmol()
       do j=1,this%ntmol(i)
          do k=1,this%bondscnt(i)
             ni=nx+this%molbond(i,k,1)
             nj=nx+this%molbond(i,k,2)
             call this%mic(ni,nj,xvz,yvz,zvz)
             dr=sqrt(xvz**2+yvz**2+zvz**2)
             do l=1,2
                prm(l)=this%parbnd(i,k,l)
             end do
             ptrm=this%tbonds(i,k)
             call this%bnd%set_bonds(dr,prm,ptrm)
             call this%set_force(ni,nj,xvz,yvz,zvz,this%bnd%get_force())
             call this%bnd%set_virbond(this%bnd%get_force(),dr)
             enpot=enpot+this%bnd%get_enbond()
             virtot=virtot+this%bnd%get_virbond()
          end do
          do k=1,this%bendscnt(i)
             ni=nx+this%molbend(i,k,2)
             nj=nx+this%molbend(i,k,1)
             nk=nx+this%molbend(i,k,3)
             call this%mic(ni,nj,drij(1),drij(2),drij(3))
             call this%mic(ni,nk,drik(1),drik(2),drik(3))
             dr1=sqrt(drij(1)**2+drij(2)**2+drij(3)**2)
             dr2=sqrt(drik(1)**2+drik(2)**2+drik(3)**2)
             theta=acos((drij(1)*drik(1)+drij(2)*drik(2)+drij(3)*drik(3))/(dr1*dr2))
             do l=1,2
                prm(l)=this%parbend(i,k,l)
             end do
             ptrm=this%tbends(i,k)
             call this%tht%set_angles(theta,prm,ptrm)
             call this%set_force&
                  (ni,nj,nk,drij,drik,dr1,dr2,theta,this%tht%get_force())
             call this%tht%set_virbend(this%fbj,this%fbk,drij,drik)
             enpot=enpot+this%tht%get_enbend()
             virtot=virtot+this%tht%get_virbend()
          end do
          do k=1,this%torscnt(i)
             ni=nx+this%moltors(i,k,1)
             nj=nx+this%moltors(i,k,2)
             nk=nx+this%moltors(i,k,3)
             nl=nx+this%moltors(i,k,4)
             call this%mic(ni,nj,drij(1),drij(2),drij(3))
             call this%mic(nj,nk,drjk(1),drjk(2),drjk(3))
             call this%mic(nk,nl,drkl(1),drkl(2),drkl(3))
             vc1x=drij(2)*drjk(3)-drij(3)*drjk(2)
             vc1y=drij(3)*drjk(1)-drij(1)*drjk(3)
             vc1z=drij(1)*drjk(2)-drij(2)*drjk(1)
             vc2x=drjk(2)*drkl(3)-drjk(3)*drkl(2)
             vc2y=drjk(3)*drkl(1)-drjk(1)*drkl(3)
             vc2z=drjk(1)*drkl(2)-drjk(2)*drkl(1)
             dr1=sqrt(vc1x**2+vc1y**2+vc1z**2) !-|rij x rjk|
             dr2=sqrt(vc2x**2+vc2y**2+vc2z**2) !-|rjk x rkn|
             phi=acos((vc1x*vc2x+vc1y*vc2y+vc1z*vc2z)/(dr1*dr2))
             do l=1,3
                prm(l)=this%partors(i,k,l)
             end do
             ptrm=this%ttors(i,k)
             call this%dih%set_dihedrals(phi,prm,ptrm)
             call this%set_force&
                  (ni,nj,nk,nl,drij,drjk,drkl,dr1,dr2,phi,this%dih%get_force())
             ri(1)=this%xa(ni)
             ri(2)=this%ya(ni)
             ri(3)=this%za(ni)
             rj(1)=this%xa(nj)
             rj(2)=this%ya(nj)
             rj(3)=this%za(nj)
             rk(1)=this%xa(nk)
             rk(2)=this%ya(nk)
             rk(3)=this%za(nk)
             rl(1)=this%xa(nl)
             rl(2)=this%ya(nl)
             rl(3)=this%za(nl)
             call this%dih%set_virtors(this%fbi,this%fbj,this%fbk,this%fbl,ri,rj,rk,rl)
             virtot=virtot+this%dih%get_virtors()
          end do
          nx=nx+this%nxmol(i)
       end do
       nx=0
       if(this%nxmol(i).lt.4)goto 2
       do j=1,this%ntmol(i)
          do k=1,this%nxmol(i)
             ni=nx+k
             do l=k+1,this%nxmol(i)
                nj=nx+l
                do m=1,this%bendscnt(i)
                   do n=1,3
                      do o=n+1,3
                         if(ni.eq.this%molbend(i,m,n).and.nj.eq.this%molbend(i,m,o))goto 1
                         if(ni.eq.this%molbend(i,m,o).and.nj.eq.this%molbend(i,m,n))goto 1
                      end do
                   end do
                end do
                call this%mic(ni,nj,xvz,yvz,zvz)
                dr=sqrt(xvz**2+yvz**2+zvz**2)
                if(abs(this%qat(ni)*this%qat(nj)).gt.1.d-8)then
                   call this%coul%set_coulomb(dr,this%qat(ni),this%qat(nj))
                   call this%set_force&
                        (ni,nj,xvz,yvz,zvz,this%coul%get_force()*this%sf_coul(i))
                   call this%coul%set_vircoul(this%coul%get_force()*this%sf_coul(i),dr)
                   enpot=enpot+this%coul%get_encoul()*this%sf_coul(i)
                   virtot=virtot+this%coul%get_vircoul()
                end if
                do m=1,this%get_nvdw()
                   if(this%tpa(ni).eq.this%spcvdw(m,1).and.&
                        this%tpa(nj).eq.this%spcvdw(m,2).or.&
                        this%tpa(ni).eq.this%spcvdw(m,2).and.&
                        this%tpa(nj).eq.this%spcvdw(m,1))then
                      do n=1,2
                         prm(n)=this%parvdw(m,n)
                      end do
                      ptrm='charmm'
                      call this%vdw%set_vanderwaals(dr,prm,ptrm)
                      call this%set_force&
                           (ni,nj,xvz,yvz,zvz,this%vdw%get_force()*this%sf_vdw(i))
                      call this%vdw%set_virvdw(this%vdw%get_force()*this%sf_vdw(i),dr)
                      enpot=enpot+this%vdw%get_envdw()*this%sf_vdw(i)
                      virtot=virtot+this%vdw%get_virvdw()
                   end if
                end do
1               continue
             end do
          end do
          nx=nx+this%nxmol(i)
       end do
2      continue
    end do
    do i=1,this%get_natom()
       do j=1,this%nlist(i)
          ni=i
          nj=this%ilist(i,j)
          call this%mic(ni,nj,xvz,yvz,zvz)
          dr=sqrt(xvz**2+yvz**2+zvz**2)
          if(abs(this%qat(ni)*this%qat(nj)).gt.1.d-8)then
             call this%coul%set_coulomb(dr,this%qat(ni),this%qat(nj))
             call this%set_force(ni,nj,xvz,yvz,zvz,this%coul%get_force())
             call this%coul%set_vircoul(this%coul%get_force(),dr)
             enpot=enpot+this%coul%get_encoul()
             virtot=virtot+this%coul%get_vircoul()
          end if
          do k=1,this%get_nvdw()
             if(this%tpa(ni).eq.this%spcvdw(k,1).and.this%tpa(nj).eq.this%spcvdw(k,2).or.&
                  this%tpa(ni).eq.this%spcvdw(k,2).and.this%tpa(nj).eq.this%spcvdw(k,1))then
                do l=1,2
                   prm(l)=this%parvdw(k,l)
                end do
                ptrm=this%tvdw(k)
                call this%vdw%set_vanderwaals(dr,prm,ptrm)
                call this%set_force(ni,nj,xvz,yvz,zvz,this%vdw%get_force())
                call this%vdw%set_virvdw(this%vdw%get_force(),dr)
                enpot=enpot+this%vdw%get_envdw()
                virtot=virtot+this%vdw%get_virvdw()
             end if
          end do
       end do
    end do
    call this%set_enpot(enpot+this%get_encorr())
    call this%set_virtot(virtot+this%get_vircorr())
  end subroutine set_forcefield

  subroutine set_force2(this,ni,nj,xvz,yvz,zvz,fr)
    implicit none
    class(interaction), intent(inout) :: this
    integer, intent(in)               :: ni,nj
    real(8), intent(in)               :: fr,xvz,yvz,zvz
    this%fbi(1)=-fr*xvz
    this%fbi(2)=-fr*yvz
    this%fbi(3)=-fr*zvz
    this%fbj(1)=+fr*xvz
    this%fbj(2)=+fr*yvz
    this%fbj(3)=+fr*zvz
    this%fax(ni)=this%fax(ni)+this%fbi(1)
    this%fay(ni)=this%fay(ni)+this%fbi(2)
    this%faz(ni)=this%faz(ni)+this%fbi(3)
    this%fax(nj)=this%fax(nj)+this%fbj(1)
    this%fay(nj)=this%fay(nj)+this%fbj(2)
    this%faz(nj)=this%faz(nj)+this%fbj(3)
  end subroutine set_force2

  subroutine set_force3(this,i1,i2,i3,drij,drik,dr1,dr2,theta,fa)
    implicit none
    class(interaction), intent(inout) :: this
    integer, intent(in)               :: i1,i2,i3
    integer                           :: ix(3),i,j
    real(8), intent(in)               :: fa,dr1,dr2,theta
    real(8), intent(in)               :: drij(3),drik(3)
    real(8)                           :: derij(3,3)
    ix(1)=i1
    ix(2)=i2
    ix(3)=i3
    do j=1,3
       do i=1,3
          derij(i,j)=(kronij(ix(i),ix(2))-kronij(ix(i),ix(1)))*drik(j)/(dr1*dr2) &
               +(kronij(ix(i),ix(3))-kronij(ix(i),ix(1)))*drij(j)/(dr1*dr2) &
               -cos(theta)*((kronij(ix(i),ix(2))-kronij(ix(i),ix(1)))*drij(j)/dr1**2 &
               +(kronij(ix(i),ix(3))-kronij(ix(i),ix(1)))*drik(j)/dr2**2)
       end do
    end do
    this%fbi(1)=fa*derij(1,1)/sin(theta)
    this%fbi(2)=fa*derij(1,2)/sin(theta)
    this%fbi(3)=fa*derij(1,3)/sin(theta)
    this%fbj(1)=fa*derij(2,1)/sin(theta)
    this%fbj(2)=fa*derij(2,2)/sin(theta)
    this%fbj(3)=fa*derij(2,3)/sin(theta)
    this%fbk(1)=fa*derij(3,1)/sin(theta)
    this%fbk(2)=fa*derij(3,2)/sin(theta)
    this%fbk(3)=fa*derij(3,3)/sin(theta)
    this%fax(i1)=this%fax(i1)+this%fbi(1)
    this%fay(i1)=this%fay(i1)+this%fbi(2)
    this%faz(i1)=this%faz(i1)+this%fbi(3)
    this%fax(i2)=this%fax(i2)+this%fbj(1)
    this%fay(i2)=this%fay(i2)+this%fbj(2)
    this%faz(i2)=this%faz(i2)+this%fbj(3)
    this%fax(i3)=this%fax(i3)+this%fbk(1)
    this%fay(i3)=this%fay(i3)+this%fbk(2)
    this%faz(i3)=this%faz(i3)+this%fbk(3)
  end subroutine set_force3

  subroutine set_force4(this,i1,i2,i3,i4,drij,drjk,drkl,dr1,dr2,phi,fd)
    implicit none
    class(interaction), intent(inout) :: this
    integer, intent(in)               :: i1,i2,i3,i4
    integer                           :: ix(4),i,j
    real(8)                           :: dvc(4,3),phi2
    real(8), intent(in)               :: fd,dr1,dr2,phi
    real(8), intent(in)               :: drij(3),drjk(3),drkl(3)
    phi2=max(phi,1.d-8)
    ix(1)=i1 !-i
    ix(2)=i2 !-j
    ix(3)=i3 !-k
    ix(4)=i4 !-n
    do i=1,4
       do j=1,3
          dvc(i,j)=dfunc1(i1,i2,i3,i4,drij,drjk,drkl,ix(i),j)/(dr1*dr2) &
               -0.5d0*cos(phi2)*(dfunc2(i1,i2,i3,drij,drjk,ix(i),j)/dr1**2 &
               +dfunc2(i2,i3,i4,drjk,drkl,ix(i),j)/dr2**2)
       end do
    end do
    this%fbi(1)=fd*dvc(1,1)/sin(phi2)
    this%fbi(2)=fd*dvc(1,2)/sin(phi2)
    this%fbi(3)=fd*dvc(1,3)/sin(phi2)
    this%fbj(1)=fd*dvc(2,1)/sin(phi2)
    this%fbj(2)=fd*dvc(2,2)/sin(phi2)
    this%fbj(3)=fd*dvc(2,3)/sin(phi2)
    this%fbk(1)=fd*dvc(3,1)/sin(phi2)
    this%fbk(2)=fd*dvc(3,2)/sin(phi2)
    this%fbk(3)=fd*dvc(3,3)/sin(phi2)
    this%fbl(1)=fd*dvc(4,1)/sin(phi2)
    this%fbl(2)=fd*dvc(4,2)/sin(phi2)
    this%fbl(3)=fd*dvc(4,3)/sin(phi2)
    this%fax(i1)=this%fax(i1)+this%fbi(1)
    this%fay(i1)=this%fay(i1)+this%fbi(2)
    this%faz(i1)=this%faz(i1)+this%fbi(3)
    this%fax(i2)=this%fax(i2)+this%fbj(1)
    this%fay(i2)=this%fay(i2)+this%fbj(2)
    this%faz(i2)=this%faz(i2)+this%fbj(3)
    this%fax(i3)=this%fax(i3)+this%fbk(1)
    this%fay(i3)=this%fay(i3)+this%fbk(2)
    this%faz(i3)=this%faz(i3)+this%fbk(3)
    this%fax(i4)=this%fax(i4)+this%fbl(1)
    this%fay(i4)=this%fay(i4)+this%fbl(2)
    this%faz(i4)=this%faz(i4)+this%fbl(3)
  end subroutine set_force4

  subroutine set_enpot(this,enpot)
    implicit none
    class(interaction), intent(inout) :: this
    real(8), intent(in)               :: enpot
    this%enpot=enpot
  end subroutine set_enpot

  double precision function get_enpot(this)
    implicit none
    class(interaction), intent(inout) :: this
    get_enpot=this%enpot
  end function get_enpot

  subroutine set_virtot(this,virtot)
    implicit none
    class(interaction), intent(inout) :: this
    real(8), intent(in)               :: virtot
    this%virtot=virtot
  end subroutine set_virtot

  double precision function get_virtot(this)
    implicit none
    class(interaction), intent(inout) :: this
    get_virtot=this%virtot
  end function get_virtot

  subroutine set_vdwcorr(this)
    implicit none
    class(interaction), intent(inout) :: this
    integer                           :: natp1,natp2,i,j
    real(8)                           :: envdw_corr,virvdw_corr,es,vs
    real(8)                           :: prm(2)
    envdw_corr=0.d0
    virvdw_corr=0.d0
    do i=1,this%get_nvdw()
       select case(this%tvdw(i))
       case('charmm')
          do j=1,2
             prm(j)=this%parvdw(i,j)
          end do
          es=prm(1)*&
               (prm(2)**12-6.d0*(this%get_rcutoff()*prm(2))**6)/(9.d0*this%get_rcutoff()**9)
          vs=12.d0*prm(1)*&
               (prm(2)**12-3.d0*(this%get_rcutoff()*prm(2))**6)/(9.d0*this%get_rcutoff()**9)
       case('lj')
          do j=1,2
             prm(j)=this%parvdw(i,j)
          end do
          es=4.d0*prm(1)*&
               (prm(2)**12-3.d0*(this%get_rcutoff()*prm(2))**6)/(9.d0*this%get_rcutoff()**9)
          vs=24.d0*prm(1)*(2.d0*prm(2)**12-3.d0*(this%get_rcutoff()*prm(2))**6)/&
               (9.d0*this%get_rcutoff()**9)
       end select
       natp1=0
       natp2=0
       do j=1,this%get_natom()
          if(this%tpa(j).eq.this%spcvdw(i,1))natp1=natp1+1
          if(this%tpa(j).eq.this%spcvdw(i,2))natp2=natp2+1
       end do
       envdw_corr=envdw_corr+es*natp1*natp2
       virvdw_corr=virvdw_corr+vs*natp1*natp2
    end do
    this%encorr=2.0d0*this%get_pi()*envdw_corr/this%get_volume()
    this%vircorr=-2.0d0*this%get_pi()*virvdw_corr/this%get_volume()
  end subroutine set_vdwcorr

  double precision function get_encorr(this)
    implicit none
    class(interaction), intent(inout) :: this
    get_encorr=this%encorr
  end function get_encorr

  double precision function get_vircorr(this)
    implicit none
    class(interaction), intent(inout) :: this
    get_vircorr=this%vircorr
  end function get_vircorr

  integer function kronij(i,j)
    integer, intent(in)               :: i,j
    kronij=int((float(i+j)-abs(i-j))/(float(i+j)+abs(i-j)))
  end function kronij

  double precision function dfunc1(i1,i2,i3,i4,drij,drjk,drkn,i,j)
    implicit none
    integer i,j,i1,i2,i3,i4
    real(8) drij(3),drjk(3),drkn(3)
    dfunc1=drij(j)*acomt(drjk,drjk,j)*(kronij(i,i3)-kronij(i,i4)) &
         +drij(j)*acomt(drjk,drkn,j)*(kronij(i,i3)-kronij(i,i2)) &
         +drjk(j)*acomt(drij,drjk,j)*(kronij(i,i4)-kronij(i,i3)) &
         +drjk(j)*acomt(drjk,drkn,j)*(kronij(i,i2)-kronij(i,i1)) &
         +drkn(j)*acomt(drij,drjk,j)*(kronij(i,i3)-kronij(i,i2)) &
         +drkn(j)*acomt(drjk,drjk,j)*(kronij(i,i1)-kronij(i,i2)) &
         +2.d0*drjk(j)*acomt(drij,drkn,j)*(kronij(i,i2)-kronij(i,i3))
  end function dfunc1

  double precision function dfunc2(i1,i2,i3,dri1,dri2,i,j)
    implicit none
    integer i,j,i1,i2,i3
    real(8) dri1(3),dri2(3)
    dfunc2=dri1(j)*acomt(dri2,dri2,j)*(kronij(i,i2)-kronij(i,i1)) &
         +dri1(j)*acomt(dri1,dri2,j)*(kronij(i,i2)-kronij(i,i3)) &
         +dri2(j)*acomt(dri1,dri1,j)*(kronij(i,i3)-kronij(i,i2)) &
         +dri2(j)*acomt(dri1,dri2,j)*(kronij(i,i1)-kronij(i,i2))
    dfunc2=2.d0*dfunc2
  end function dfunc2

  double precision function acomt(x1,x2,k)
    implicit none
    integer i,k
    real(8) sum,x1(3),x2(3)
    sum=0.d0
    do i=1,3
       sum=sum+(1.d0-kronij(i,k))*x1(i)*x2(i)
    end do
    acomt=sum
  end function acomt

end module interaction_module
