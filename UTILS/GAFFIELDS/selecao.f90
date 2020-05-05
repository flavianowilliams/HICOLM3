module selecao

  use input
  use populacao
  use objetivo

  integer iwst,ibst,nx

  real(8) sumfit,fitnessmn,fitnessmx,wndwn,fitmx

contains

  subroutine avaliacao

    implicit none

    integer i,j,aux
    real(8) fitnessmn,fitnessmx,xn(nparam)

    sumfit=0.d0
    do i=1,npop
       do j=1,nparam
          xn(j)=ind(i,j)
       end do
!       if(opcod.eq.1)call aptidao_bin(xn,fitness(i))
       call aptidao(xn,fitness(i))
       sumfit=sumfit+fitness(i)
    end do

    do i=1,npop
       prob(i)=fitness(i)/sumfit
    end do

    !-rankeando o fitness pelo metodo bubble sort

    do i=1,npop
       ik0(i)=i
    end do

    do i=1,npop
       do j=1,npop-1
          if(fitness(ik0(j)).lt.fitness(ik0(j+1)))then
             aux=ik0(j)
             ik0(j)=ik0(j+1)
             ik0(j+1)=aux
          end if
       end do
    end do

    iwst=ik0(npop)
    ibst=ik0(1)
    fitnessmn=fitness(iwst)
    fitnessmx=fitness(ibst)

    !-aplicando metodo windowing

!    wndwn=fitnessmn-mod(fitnessmn,2.d0)
!
!    do i=1,npop
!       fitness(i)=fitness(i)-wndwn
!    end do
!    fitnessmn=fitness(iwst)
!    fitnessmx=fitness(ibst)
!    sumfit=sumfit-wndwn*npop
    wndwn=0.d0

    return

  end subroutine avaliacao

  subroutine roleta()

    implicit none

    integer i,j,ix
    real(8) parcial,dummy,rand

    do i=1,npop
       ix0(i)=i
    end do

    do i=1,npop

       call RANDOM_NUMBER(dummy)

       rand=dummy*sumfit

       ix=1
       parcial=0.d0
       do j=1,npop
          parcial=parcial+fitness(ik0(j))
          if(parcial.ge.rand)exit
          ix=ix+1
       end do

       ix0(i)=ix

    end do

    return

  end subroutine roleta

end module selecao
