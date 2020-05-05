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
module tight_binding
  !*****************************************************************************************
  ! Tight Binding                                                                          *
  !                                                                                        *
  ! Flaviano Williams Fernandes, 05 de dezembro de 2018                                    *
  !*****************************************************************************************

  use input
  use fitting
  use optimizing
  use propriedades

contains

  subroutine tb(t3)

    implicit none

    real(8) t0,t3

    !-contagem do tempo absoluto

    call cpu_time(t0)

    !-preparando tight-binding

    call tb_prepare

    !-imprimindo dados de entrada

    call tb_print

    !-convertendo unidades de medida

    call tb_convert

    !-obtendo bandas de energia do calculo DFT

    call DFT

    !-otimizando parametros SK

    call optimize

    !-calculo das propriedades eletronicas

    call solving

    !-contagem do tempo absoluto

    call cpu_time(t3)

    t3=t3-t0

    return

  end subroutine tb

  subroutine tb_print
    !*************************************************************************
    ! Impressão dos dados de entrada                                         *
    !*************************************************************************
    implicit none

    integer i,j,k,l,m1,nx
    real(8) t1(3)

    !-parametros do calculo tight-bindind

    nx=0
    do i=1,nlo
       do j=1,i
          do m1=1,lo(j)
             call labels(lo(j),lo(i),m1,k)
             do l=1,nprllm(k)
                nx=nx+2
             end do
          end do
       end do
    end do
    do i=1,nlo
       do j=1,mlo(i)
          do k=1,npril(lo(i),mo(i,j))
             nx=nx+1
          end do
       end do
    end do

    do i=1,nlo
       do j=1,mlo(i)
          if(npril(lo(i),mo(i,j)).gt.1)then
             nx=nx+1
             goto 1
          end if
       end do
    end do

    nx=nx+1

1   write(6,'(5x,a36)')'Parametros do calculo Tight-binding:'
    write(6,'(5x,36a1)')('-',j=1,36)
    write(6,'(5x,a7,i5)')'diagop:',diagop
    write(6,'(5x,a7,f5.2)')' rcutt:',rcutt*rconv
    write(6,'(5x,a7,f5.2)')' lcutt:',lcutt*econv
    write(6,'(5x,36a1)')('-',j=1,36)
    write(6,*)
    write(6,*)'     Orbitais atomicos:'
    write(6,'(5x,36a1)')('-',j=1,36)
    write(6,'(5x,a3,3x,a2)')'l','m'
    write(6,'(5x,36a1)')('-',j=1,36)
    do j=1,nlo
       write(6,'(6x,i2,5(3x,i2))')(lo(j)-1),((mo(j,k)-1),k=1,mlo(j))
    end do
    write(6,'(5x,36a1)')('-',j=1,36)
    write(6,*)
    write(6,*)'     Numero total parametros:',nx
    write(6,*)
    do i=1,nlo
       do j=1,mlo(i)
          if(npril(lo(i),mo(i,j)).ge.1)then
             write(6,41)'Parametro de densidade:',dens0*rconv**(-0.5d0)
             goto 2
          end if
       end do
    end do
2   write(6,*)
    write(6,*)'     Termos on-site:'
    write(6,43)('-',j=1,21)
    write(6,44)'l','m','on-site'
    write(6,43)('-',j=1,21)
    do i=1,nlo
       do j=1,mlo(i)
          write(6,47)lo(i)-1,mo(i,j)-lo(i),&
               (onst0(lo(i),mo(i,j),k)*econv,k=1,npril(lo(i),mo(i,j)))
       end do
    end do
    write(6,43)('-',j=1,21)
    write(6,*)
    write(6,*)'     Parametros de hoppings:'
    write(6,43)('-',j=1,41)
    write(6,45)'l','l','m','offsite'
    write(6,43)('-',j=1,41)
    do i=1,nlo
       do j=1,i
          do m1=1,lo(j)
             call labels(lo(j),lo(i),m1,k)
             t1(1)=tpr0(k,1)*econv
             t1(2)=tpr0(k,2)
             t1(3)=tpr0(k,3)*(econv/rconv)
             write(6,48)lo(j)-1,lo(i)-1,m1-1,(t1(l),l=1,nprllm(k))
          end do
       end do
    end do
    write(6,43)('-',j=1,41)
    write(6,*)
    write(6,*)'     Parametros de overlaps:'
    write(6,43)('-',j=1,41)
    write(6,45)'l','l','m','overlap'
    write(6,43)('-',j=1,41)
    do i=1,nlo
       do j=1,i
          do m1=1,lo(j)
             call labels(lo(j),lo(i),m1,k)
             t1(1)=spr0(k,1)
             t1(2)=spr0(k,2)
             t1(3)=spr0(k,3)*kconv
             write(6,48)lo(j)-1,lo(i)-1,m1-1,(t1(l),l=1,nprllm(k))
          end do
       end do
    end do
    write(6,43)('-',j=1,41)
    write(6,*)

    return

41  format(5x,a23,f12.4)
43  format(5x,43a1)
44  format(5x,2(a2,2x),a8,2x)
45  format(5x,3(a1,2x),3x,a8,2x)
47  format(5x,2(i2,2x),4(f8.4,2x))
48  format(5x,3(i1,2x),10(f8.4,2x))

  end subroutine tb_print

  subroutine tb_convert
    !**************************************************************************************
    ! Conversão para o sistema padrão u.a. de unidades                                    *
    !**************************************************************************************

    implicit none

    integer i,j,k

    err=err*0.5d0/rconv

    !-Parametros SK

    do i=1,lmax
       do j=1,2*i+1
          do k=1,npril(i,j)
             onst0(i,j,k)=onst0(i,j,k)/econv
          end do
       end do
    end do

    do i=1,nparamt
       tpr0(mskt(i),1)=tpr0(mskt(i),1)/econv
       tpr0(mskt(i),2)=tpr0(mskt(i),2)
       tpr0(mskt(i),3)=tpr0(mskt(i),3)/(econv/rconv)
    end do

    do i=1,nparamt
       spr0(mskt(i),1)=spr0(mskt(i),1)
       spr0(mskt(i),2)=spr0(mskt(i),2)
       spr0(mskt(i),3)=spr0(mskt(i),3)/(1.0d0/rconv)
    end do

    dens0=dens0/rconv**(-0.5d0)

    lcutt=lcutt/rconv
    rcutt=rcutt/rconv

    !-energia de Fermi

    fermi=fermi/econv

    !-parametros no espaço K

    do i=1,3
       kptt0(i)=kptt0(i)/kconv
    end do

    return

  end subroutine tb_convert

end module tight_binding
