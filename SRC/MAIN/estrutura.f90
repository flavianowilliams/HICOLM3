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

  subroutine structure_prepare(t1)
    !***************************************************************************************
    ! Propriedades:                                                                        *
    ! - Definicao do grupo de simetria;                                                    *
    ! - Definicao da celula unitaria;                                                      *
    !***************************************************************************************

    implicit none

    integer i,j
    real(8) vl(3),t1,tf,sum

    !-contagem de tempo

    call cpu_time(t1)

    !-redefinindo parametros de rede e posicoes atomicas

    if(reuse.gt.0)call frame

    !-imprimindo coordenadas em arquivo AXSF

    open(7,file='HICOLM.AXSF',status='unknown')

    !-convertendo unidades de medida

    call convert

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

    call ccp

    !-Definindo grupo de simetria

    call cell_symmetry

    !-checando viabilidade geometrica de cada molecula

    call structure_check

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

  end subroutine structure_prepare

 subroutine structure_check

   implicit none

   integer i,j,k,np,ni,nj,nk,nn
   real(8) shift,dx1,dy1,dz1,dx2,dy2,dz2

   shift=1.d-2

   !-aplicando translacao nas moleculas

   np=0
   do i=1,nmolec
      do j=1,ntmolec(i)
         do k=1,torscnt(i)
            ni=np+moltors(i,k,1)
            nj=np+moltors(i,k,2)
            nk=np+moltors(i,k,3)
            nn=np+moltors(i,k,4)
            call mic(ni,nj,dx1,dy1,dz1)
            call mic(nk,nn,dx2,dy2,dz2)
            if(abs(dx1).lt.1.e-8.and.abs(dx2).lt.1.e-8)then
               za(ni)=za(ni)-shift
               za(nn)=za(nn)+shift
            end if
            if(abs(dy1).lt.1.e-8.and.abs(dy2).lt.1.e-8)then
               za(ni)=za(ni)-shift
               za(nn)=za(nn)+shift
            end if
            if(abs(dz1).lt.1.e-8.and.abs(dz2).lt.1.e-8)then
               za(ni)=za(ni)-shift
               za(nn)=za(nn)+shift
            end if
         end do
         np=np+nxmolec(i)
      end do
   end do

   !-aplicando condicoes de contorno periodica

   call ccp

   return

 end subroutine structure_check

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

  subroutine ccp()
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

  end subroutine ccp

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

  subroutine frame
    !***************************************************************************************
    ! Leitura do ficheiro de entrada:                                                      *
    ! - Coordenadas espaciais;                                                             *
    ! - Velocidades;                                                                       *
    ! - Forca.                                                                             *
    !***************************************************************************************
    implicit none

    integer i,j

    open(1,file='HICOLM.XSF',status='old')

    do i=1,13
       read(1,*)
    end do

    do i=1,3
       read(1,'(3(3x,f14.8))')(v(i,j),j=1,3)
    end do

    read(1,*)
    read(1,*)natom

    select case(reuse)
    case(1)
       do i=1,natom
          read(1,'(i5,3f14.8,2(2x,3f14.8))') &
               idna(i),xa(i),ya(i),za(i)
       end do
    case(2)
       do i=1,natom
          read(1,'(i5,3f14.8,2(2x,3f14.8))') &
               idna(i),xa(i),ya(i),za(i),fax(i),fay(i),faz(i)
       end do
    case(3)
       do i=1,natom
          read(1,'(i5,3f14.8,2(2x,3f14.8))') &
               idna(i),xa(i),ya(i),za(i),fax(i),fay(i),faz(i),vax(i),vay(i),vaz(i)
       end do
    end select

    !-fechando arquivo XSF com as coordenadas atomicas

    close(1)

    return

  end subroutine frame

end module estrutura
