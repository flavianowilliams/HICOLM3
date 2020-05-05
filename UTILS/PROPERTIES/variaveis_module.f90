module variaveis_module

  use input_module

contains

  subroutine variaveis(nstp,nat)

    implicit none

    integer i,j,k,nstp,nat,varmax

    parameter (varmax=13)

    real(8) var(varmax),sum(varmax),med(varmax),dsv(varmax)

    do i=1,varmax
       sum(i)=0.d0
    end do

    write(2,'(7x,13a12)')'#TIME',('STRESS',i=1,6),'VOLUME','TEMPERA.','PRESSURE','EKINET',&
         'ETOTAL','DENSITY'

    do i=1,nstp
       read(1,*)
       read(1,*)(var(j),j=1,6)
       read(1,*)(var(j),j=7,12)
       read(1,*)(var(j),j=13,varmax)
       do k=1,5
          read(1,*)
       end do
       do k=1,nat
          read(1,*)
          read(1,*)
          read(1,*)
       end do
       do j=1,varmax
          sum(j)=sum(j)+var(j)
       end do
       write(2,'(7x,13f12.4)')var(7),(var(j),j=1,6),(var(j),j=8,13)
    end do

    do j=1,varmax
       med(j)=sum(j)/nstp
    end do

    rewind(1)

    read(1,*)

    do j=1,varmax
       sum(j)=0.d0
    end do

    do i=1,nstp
       read(1,*)
       read(1,*)(var(j),j=1,6)
       read(1,*)(var(j),j=7,12)
       read(1,*)(var(j),j=13,varmax)
       do k=1,5
          read(1,*)
       end do
       do k=1,nat
          read(1,*)
          read(1,*)
          read(1,*)
       end do
       do j=1,varmax
          sum(j)=sum(j)+(var(j)-med(j))**2
       end do
    end do

    do j=1,varmax
       dsv(j)=sqrt(sum(j)/(nstp-1))
    end do

    write(2,*)
    write(2,'(a7,6f12.4)')'## MED1',(med(i),i=1,6)
    write(2,'(a7,6f12.4)')'## DSV1',(dsv(i),i=1,6)
    write(2,*)
    write(2,'(a7,6f12.4)')'## MED2',(med(i),i=7,12)
    write(2,'(a7,6f12.4)')'## DSV2',(dsv(i),i=7,12)
    write(2,*)
    write(2,'(a7,6f12.4)')'## MED3',(med(i),i=13,varmax)
    write(2,'(a7,6f12.4)')'## DSV3',(dsv(i),i=13,varmax)

    return

  end subroutine variaveis

end module variaveis_module
