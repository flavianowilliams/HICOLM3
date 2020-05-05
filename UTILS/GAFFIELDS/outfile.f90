module outfile

  use populacao
  use selecao
  use operadores

contains

  subroutine output

    implicit none

    integer i,j,nx
    real(8) xn(nparam)

    open(3,file='population.dat',status='unknown')
    
!    select case(opcod)
!    case(1)
!       do i=1,npop
!          write(*,'(i3,2x,f14.8,2x,i10,2x,f14.8,2x,f14.8,2x,100i1)')&
!               i,ind(i,1),dec(i),fitness(i),prob(i),(nint(bin(i,j)),j=1,ncron)
!       end do
!    case(2)
       do i=1,npop
          write(3,'(i3,2x,f14.8,2x,f14.8,2x,100f10.4)')&
               i,fitness(i),prob(i),(ind(i,j),j=1,ncron)
          do j=1,nparam
             xn(j)=ind(ibst,j)
          end do
       end do
       nx=1
       do i=1,nmolec
          do j=1,nprmstr(i)
             if(chkvar(1,i,j).eqv..true.)then
                prmstr(i,j)=xn(nx)
                nx=nx+1
             end if
          end do
          do j=1,nprmscr(i)
             if(chkvar(2,i,j).eqv..true.)then
                prmscr(i,j)=xn(nx)
                nx=nx+1
             end if
          end do
       end do
       write(*,*)
       write(*,'(1x,a33)')'Calculo finalizado com sucesso!!!'
       write(*,*)
       write(*,'(1x,a5,f7.2)')'RMSD:',1.d0/fitness(ibst)
       write(*,*)
       write(*,'(a22,2x,15f12.6)')'Individuo selecionado:',(ind(ibst,i),i=1,nparam)
       write(*,*)
       write(*,*)'Campo de forca:'
       write(*,*)
       do i=1,nmolec
          write(*,*)'Molecula:',idmolec(i)
          write(*,*)
          write(*,*)'Bonds:'
          write(*,*)'Tipo:',ntstr(i)
          write(*,*)'Parametros:',(prmstr(i,j),j=1,nprmstr(i))
          write(*,*)
          write(*,*)'Bends:'
          write(*,*)'Tipo:',ntscr(i)
          write(*,*)'Parametros:',(prmscr(i,j),j=1,nprmscr(i))
          write(*,*)
       end do
!    end select

    write(*,*)
    write(*,'(a21,f12.4,i4)')'Fitness minima:',fitness(iwst),iwst
    write(*,'(a21,f12.4,i4)')'Fitness maxima:',fitness(ibst),ibst
    write(*,'(a21,f12.4)')'Fitness acumu.:',sumfit

    write(*,*)

    return

  end subroutine output

end module outfile
