module utils

  use input
  use selecao

contains

  subroutine dsv(dvrst)

    implicit none

    integer i
    real(8) med,dvrst

    med=0.d0
    do i=1,npop
       med=med+fitness(i)
    end do

    med=med/npop

    dvrst=0.d0
    do i=1,npop
       dvrst=dvrst+(fitness(i)-med)**2
    end do

    dvrst=dvrst/npop

    return

  end subroutine dsv

end module utils
