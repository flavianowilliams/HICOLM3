module populacao

  use input

  integer(8), allocatable :: dec(:),decl(:)

  integer, allocatable :: imm(:),imn(:)
  integer, allocatable :: ix0(:),ik0(:)

  real(8), allocatable :: ind(:,:),indl(:,:)
  real(8), allocatable :: fitness(:),prob(:)
  real(8), allocatable :: bin(:,:),binl(:,:)

contains

  subroutine alloc_populacao

    implicit none

    !-alocando arrays

    allocate(ind(npop,nparam),indl(npop,nparam))
    allocate(bin(npop,ncron),binl(npop,ncron),dec(npop),decl(npop))
    allocate(imm(npop*ncron),imn(npop*ncron))
    allocate(fitness(npop),prob(npop))
    allocate(ix0(npop),ik0(npop))

    return

  end subroutine alloc_populacao

  subroutine individuos

    implicit none

    integer i,j
    real(8) dummy

    write(*,*)'Populacao:'
    write(*,*)
    write(*,'(4x,a21,i5)')'Tamanho da populacao:',npop
    write(*,'(4x,a31,f7.4)')'Intervalo de valores dos genes:',dind
    write(*,'(4x,a26,f7.4)')'Porcentagem de cruzamento:',txgr
    write(*,'(4x,a25,f7.4)')'Probabilidade de mutacao:',pmut
    write(*,*)

    !-gerando aleatoriamente populacao inicial na base dez

    do i=1,npop
       do j=1,nparam
          call RANDOM_NUMBER(dummy)
          !          ind(i,j)=ind0(j)+(2.d0*dummy-1.0d0)*dind*ind0(j)
          ind(i,j)=nimin(j)+(nimax(j)-nimin(j))*dummy
       end do
    end do

    return

  end subroutine individuos

  subroutine codificacao

    implicit none

    integer i,j
!    integer(8) xx
!    real(8) dni,num(ncron)

!    do i=1,npop
!       do j=1,ncron
!          bin(i,j)=0.d0
!       end do
!    end do

!    if(opcod.eq.1)goto 1 !-binario
!    if(opcod.eq.2)goto 2 !-real

!1   dni=nimax(1)-nimin(1)

    !-mapeando a populacao inicial

!    do i=1,npop
!       dec(i)=nint((ind(i,1)-nimin(1))*(2**ncron-1)/dni)
!    end do

!    do i=1,npop
!       xx=dec(i)
!       do j=1,ncron
!          num(j)=mod(xx,2)
!          xx=int(float(xx)/2.d0)
!       end do
!       do j=1,ncron
!          bin(i,ncron+1-j)=num(j)
!       end do
!    end do
!
!    goto 10

    do i=1,npop
       do j=1,ncron
          bin(i,j)=ind(i,j)
       end do
    end do

    return

  end subroutine codificacao

  subroutine frames_PES

    implicit none

    integer i,j

    open(4,file='GAFFIELDS.AXSF',status='unknown')

    write(4,'(a9,i5)')'ANIMSTEPS',mgeo

    do i=1,mgeo
       write(4,'(a5,i5)')'ATOMS',i
       do j=1,natom
          write(4,'(i5,6f14.8)')idax(j),xa(i,j),ya(i,j),za(i,j),fx(i,j),fy(i,j),fz(i,j)
       end do
    end do

    close(4)

    return

  end subroutine frames_PES

  subroutine frame_ending

    implicit none

    integer j

    open(6,file='GAFFIELDS.XSF',status='unknown')
    open(7,file='GAFFIELDS_chk.XSF',status='unknown')

    write(6,'(a5)')'ATOMS'

    do j=1,natom
       write(6,'(i5,6f14.8)')idax(j),xa(mgeo,j),ya(mgeo,j),za(mgeo,j),&
            fx(mgeo,j),fy(mgeo,j),fz(mgeo,j)
    end do

    write(7,'(a5)')'ATOMS'

    do j=1,natom
       write(7,'(i5,6f14.8)')idax(j),xa(mgeo,j),ya(mgeo,j),za(mgeo,j),&
            dx(mgeo,j),dy(mgeo,j),dz(mgeo,j)
    end do

    close(6)
    close(7)

    return

  end subroutine frame_ending

end module populacao
