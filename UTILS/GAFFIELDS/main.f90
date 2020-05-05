program ag

  use input
  use populacao
  use selecao
  use operadores
  use utils
  use outfile

  implicit none

  integer i,j
  real(8) dvrst,t0,t1

  call cpu_time(t0)

  open(5,file='evolucao.dat',status='unknown')

  !-lendo dados de entrada

  write(*,*)

  call parametros

  write(*,*)

  !-alocando arrays

  call alloc_inicial
  call alloc_populacao

  !-PES

  call coords

  !-convertendo unidades de medida

  call convert

  !-gerando populacao inicial

  write(*,*)'-> Gerando populacao inicial'
  write(*,*)

  call individuos

  write(*,*)

  !-codificando populacao inicial

  write(*,*)'-> Codificando populacao inicial'
  write(*,*)

  call codificacao

  write(*,*)

  !-avaliando individuos da geracao

  write(*,*)'-> Avaliando individuos da geracao inicial'
  write(*,*)

  call funcaomax
  call avaliacao

  write(*,*)

  write(*,*)'-> Imprimindo dados da geracao inicial'
  write(*,*)

  i=1
  do while (1.d0/abs(fitness(ibst)).gt.tol)
     call roleta
     call geracao
     call dsv(dvrst)
     write(5,*)i,1.d0/fitness(ibst),fitness(ibst),murat,dvrst
     if(i.eq.nstp)then
        write(*,*)'Encerrado por limite maximo de geracoes:',nstp
        write(*,'(4x,a14,f7.4)')'Erro estimado:',1.d0/fitness(ibst)
        write(*,*)
        write(*,*)'Individuo selecionado:',(ind(ibst,j),j=1,nparam)
        write(*,*)
        goto 1
     end if
     i=i+1
  end do

  call output

  !-imprimindo PES em arquivo AXSF

1 call frames_PES

  call frame_ending

  call cpu_time(t1)

  t0=t1-t0

  write(*,'(a17,f5.0,a7)')'CPU elapsed time:',t0,'seconds'

  write(*,*)

end program ag
