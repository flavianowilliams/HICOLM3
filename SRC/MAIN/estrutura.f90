!
! This file is part of the HICOLM distribution (https://github.com/flavianowilliams/HICOLM).
!
! Copyright (c) 2019 Flaviano Williams Fernandes.
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, version 3.
!
! This program is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
! General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program. If not, see <http://www.gnu.org/licenses/>.
!
module estrutura
  !******************************************************************************************
  ! Propriedades relacionada ao espaco real:                                                *
  ! - Definicao do grupo de simetria;                                                       *
  ! - Definicao da celula unitaria;                                                         *
  ! - Condicoes de contorno periódicas;                                                     *
  ! - Condicoes de contorno inversa;                                                        *
  ! - Convencao de minima imagem.                                                           *
  ! - Impressao de variaveis canonicas em HICOLM.XSF e HICOLM.AXSF                              *
  !                                                                                         *
  ! Flaviano Williams Fernandes, 05 de dezembro de 2018                                     *
  !******************************************************************************************

  use sistema
  use input
  use alloc_arrays

  implicit none

  integer nmin(natmax,10),na(natmax),gsym
  real(8) rmin(2,natmax),alpha,beta,gamma

  save nmin,rmin,gsym

contains

  subroutine coordenadas(t1)
    !***************************************************************************************
    ! Propriedades:                                                                        *
    ! - Definicao do grupo de simetria;                                                    *
    ! - Definicao da celula unitaria;                                                      *
    !***************************************************************************************

    implicit none

    integer i,j
    real(8) vl(3),t1,tf,sum

    !-imprimindo coordenadas em arquivo AXSF

    open(7,file='HICOLM.AXSF',status='unknown')

    !-contagem de tempo

    call cpu_time(t1)

    !-volume da celula unitaria

    vl(1)=(v(1,2)*v(2,3)-v(1,3)*v(2,2))
    vl(2)=(v(1,3)*v(2,1)-v(1,1)*v(2,3))
    vl(3)=(v(1,1)*v(2,2)-v(1,2)*v(2,1))

    volume=0.d0
    do j=1,3
       volume=volume+v(3,j)*vl(j)
    end do

    volume=abs(volume)

    !-angulos da celula unitaria

    sum=0.d0
    do j=1,3
       sum=sum+v(2,j)*v(3,j)
    end do
    alpha=acos(sum/(b*c))
    sum=0.d0
    do j=1,3
       sum=sum+v(1,j)*v(3,j)
    end do
    beta=acos(sum/(a*c))
    sum=0.d0
    do j=1,3
       sum=sum+v(1,j)*v(2,j)
    end do
    gamma=acos(sum/(a*b))

    !-aplicando condicoes de contorno inversa

    call ccp_inv

    !-Definindo grupo de simetria

    call cell_symmetry

    !-imprimindo informacoes do espaço real

    write(6,*)('#',i=1,93)
    write(6,*)('STRUCTURE ',i=1,9)
    write(6,*)('#',i=1,93)
    write(6,*)
    write(6,'(a18,i5)')'Total of atoms:',natom
    write(6,*)
    write(6,*)'Real space:'
    write(6,*)
    write(6,'(a16,3f15.8)')'Lattice constts:',a*rconv,b*rconv,c*rconv
    write(6,*)
    write(6,'(a16,3f15.8)')'Lattice vectors:',(v(1,i)*rconv,i=1,3)
    write(6,'(16x,3f15.8)')(v(2,i)*rconv,i=1,3)
    write(6,'(16x,3f15.8)')(v(3,i)*rconv,i=1,3)

    write(6,*)
    write(6,'(a14,f14.4)')'       VOLUME:',volume*rconv**3
    write(6,*)
    write(6,'(a14,i10)')'     Symmetry:',gsym

    write(6,*)

    !-contagem de tempo

    call cpu_time(tf)

    t1=tf-t1

    return

 end subroutine coordenadas

  subroutine ccp(i,j,rx,ry,rz,xvz,yvz,zvz)
    !****************************************************************************************
    ! Condicoes de contorno periodicas                                                      *
    !****************************************************************************************
    implicit none

    integer i,j,rx,ry,rz
    real(8) xx,yy,zz,xvz,yvz,zvz

    xx=xa(j)+rx*v(1,1)+ry*v(2,1)+rz*v(3,1)
    yy=ya(j)+rx*v(1,2)+ry*v(2,2)+rz*v(3,2)
    zz=za(j)+rx*v(1,3)+ry*v(2,3)+rz*v(3,3)

    xvz=xx-xa(i)
    yvz=yy-ya(i)
    zvz=zz-za(i)

    return

  end subroutine ccp

  subroutine mic(i,j,xvz,yvz,zvz)
    !***************************************************************************************
    ! Convencao de minima imagem                                                           *
    !***************************************************************************************

    implicit none

    integer i,j
    real(8) xx,yy,zz,xvz,yvz,zvz

    xx=v(1,1)+v(2,1)+v(3,1)
    yy=v(1,2)+v(2,2)+v(3,2)
    zz=v(1,3)+v(2,3)+v(3,3)

    xvz=(xa(j)-xa(i))-xx*int(2.d0*(xa(j)-xa(i))/xx)
    yvz=(ya(j)-ya(i))-yy*int(2.d0*(ya(j)-ya(i))/yy)
    zvz=(za(j)-za(i))-zz*int(2.d0*(za(j)-za(i))/zz)

    return

  end subroutine mic

  subroutine ccp_inv()
    !***************************************************************************************
    ! Condicoes de contorno inversa                                                        *
    !***************************************************************************************
    implicit none

    integer i
    real(8) xvz,yvz,zvz

    do i=1,natom
       xvz=v(1,1)+v(2,1)+v(3,1)
       yvz=v(1,2)+v(2,2)+v(3,2)
       zvz=v(1,3)+v(2,3)+v(3,3)
       xa(i)=xa(i)-xvz*nint(xa(i)/xvz)
       ya(i)=ya(i)-yvz*nint(ya(i)/yvz)
       za(i)=za(i)-zvz*nint(za(i)/zvz)
    end do

    return

  end subroutine ccp_inv

  subroutine ccpmm(i,j,xvz,yvz,zvz)
    !***************************************************************************************
    !subrotina responsavel por aplicar condicoes de contorno                               *
    !***************************************************************************************

    implicit none

   integer i,j,rx,ry,rz
   real(8) xvz,yvz,zvz,xx,yy,zz,xvzz,yvzz,zvzz

   xvzz=xa(j)-xa(i)
   yvzz=ya(j)-ya(i)
   zvzz=za(j)-za(i)

   do rx=-1,1
      do ry=-1,1
         do rz=-1,1
            xx=rx*v(1,1)+ry*v(2,1)+rz*v(3,1)
            yy=rx*v(1,2)+ry*v(2,2)+rz*v(3,2)
            zz=rx*v(1,3)+ry*v(2,3)+rz*v(3,3)
            xvz=xvzz-xx*aint(2.d0*xvzz/xx)
            yvz=yvzz-yy*aint(2.d0*yvzz/yy)
            zvz=zvzz-zz*aint(2.d0*zvzz/zz)
         end do
      end do
   end do

   return

 end subroutine ccpmm

 subroutine cell_symmetry
   !***************************************************************************************
   ! Grupos de simetria                                                                   *
   !***************************************************************************************

   implicit none

   real(8) prec

   !-definindo precisão

   prec=1.d-2

   if(gsymopt.eq.0)goto 1

   gsym=gsymopt

   goto 2

   !-triclinico

1  gsym=1

   !-cubico

   if(abs(2.d0*alpha-pi).le.prec)then
      if(abs(2.d0*beta-pi).le.prec)then
         if(abs(2.d0*gamma-pi).le.prec)then
            if(abs(a-b).le.prec)then
               if(abs(a-c).le.prec)gsym=23
            end if
         end if
      end if
   end if

   !-Hexagonal

   if(abs(2.d0*alpha-acos(-1.d0)).le.prec)then
      if(abs(2.d0*beta-acos(-1.d0)).le.prec)then
         if(abs(a-b).le.prec)then
            if(abs(2.d0*gamma-pi/6.d0).le.prec.or.&
                 abs(2.d0*gamma-pi/3.d0).le.prec)gsym=6
         end if
      end if
   end if

2  return

 end subroutine cell_symmetry

 subroutine geometria()
   !***************************************************************************************
   ! - Impressao de variaveis canonicas em HICOLM.XSF                                       *
   !***************************************************************************************

   implicit none

   integer i,k,j

   !-abrindo ficheiro

   open(1,file='HICOLM.XSF',status='unknown')

   !-imprimindo arquivo XSF

   write(1,*)'BEGIN_INFO'
   write(1,*)'  #'
   write(1,*)'  # This is a XCRYSDEN-Structure-File'
   write(1,*)'  # aimed for Visualization of Fermi Surface'
   write(1,*)'  #'
   write(1,*)'  # Case: Tight-binding calculation'
   write(1,*)'  #'
   write(1,*)'  # Launch as: xcrysden --xsf HICOLM.XSF'
   write(1,*)'  #'
   write(1,*)'END_INFO'

   k=1

   write(1,*)'# estrutura final'
   write(1,'(a7)')'CRYSTAL'

   write(1,'(a7)')'PRIMVEC'
   do i=1,3
      write(1,'(3(3x,f14.8))')(v(i,j)*rconv,j=1,3)
   end do

   write(1,'(a9)')'PRIMCOORD'
   write(1,'(2i5)')natom,k

   do i=1,natom
      write(1,'(i5,3f14.8,2x,3f14.8,2x,3f14.8)') &
           idna(i),xa(i)*rconv,ya(i)*rconv,za(i)*rconv, &
           fax(i)*econv/rconv,fay(i)*econv/rconv,faz(i)*econv/rconv, &
           vax(i)*rconv/tconv,vay(i)*rconv/tconv,vaz(i)*rconv/tconv
   end do

   close(1)

    return

  end subroutine geometria

  subroutine history(ihist)
    !***************************************************************************************
    ! Impressao de variaveis canonicas por frame em HICOLM.AXSF                              *
    !***************************************************************************************

    implicit none

    integer i,j,k,ihist
    real(8) fnorm

    k=1

    if(ihist.eq.1)then
       write(7,*)'# Historico de coordenadas'
       write(7,'(a9,i5)')'ANIMSTEPS',int((ntrialmax-nrelax)/nhist)
       write(7,'(a7)')'CRYSTAL'
    end if

    write(7,'(a7,i5)')'PRIMVEC',ihist
    do i=1,3
       write(7,'(3(3x,f14.8))')(v(i,j)*rconv,j=1,3)
    end do
    write(7,'(a9,i5)')'PRIMCOORD',ihist
    write(7,'(2i5)')natom,k

    do i=1,natom
       fnorm=sqrt(fax(i)**2+fay(i)**2+faz(i)**2)*econv/rconv
       write(7,'(i5,7f14.8)')idna(i),xa(i)*rconv,ya(i)*rconv,za(i)*rconv,&
            fax(i)*econv/rconv,fay(i)*econv/rconv,faz(i)*econv/rconv,fnorm
    end do

    ihist=ihist+1

    return

  end subroutine history

end module estrutura
