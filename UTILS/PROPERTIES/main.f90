program tbmc

  use input_module
  use rdf_module
  use variaveis_module

  implicit none

  open(1,file='HICOLM.md',status='old')
  open(2,file='HICOLM.dat',status='unknown')

  call entrada

  select case(nopt)
  case(1)
     call rdf_calc(nstp,nat)
  case(2)
     call variaveis(nstp,nat)
  end select

end program tbmc
