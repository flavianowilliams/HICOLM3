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
program HICOLM
  !*******************************************************************************************
  !programa principal responsavel por carregar os modulos essenciais para efetuar os cálculos*
  !*******************************************************************************************
  use input                ! parametros de entrada
  use alloc_arrays         ! alocando arrays essenciais
  use estrutura            ! definindo espaço real
  use brillouin            ! definindo zona de brillouin
  use molecular_dynamics   ! metodo Dinâmica Molecular
  use optimize             ! ferramenta de otimizacao do sistema
  use outfile              ! imprimindo informacoes em arquivo out

  implicit none

  real(8) :: t0,t1,t2,t3,t4
  character(10) host,time
  character(8) date

  call cpu_time(t0)

  open(3,file='HICOLM.df',status='unknown')  ! print dataframe
  open(6,file='HICOLM.out',status='unknown') ! imprimindo arquivo out

  t1=0.d0
  t2=0.d0
  t3=0.d0
  t4=0.d0

  !-elapsed time information

  call date_and_time(Date=date)
  call date_and_time(TIME=time)
  !  call hostnm(host)
  host='undefined'

  !-cabecalho
  !
  write(6,*)
  write(6,'(18x,a57)')'HICOLM: Multi-Methods for Molecules and Condensed Systems'
  write(6,*)
  write(6,'(39x,a14)')'Version: 2.2.1'
  write(6,*)
  write(6,'(''Host: '',2x,a10)')host
  write(6,'(''Date: '',2x,a8)')date
  write(6,'(''Time: '',2x,a10)')time
  write(6,*)

  !========================================================
  !-obtendo valores iniciais
  call cte
  !========================================================
  !-definindo parametros de entrada
  call entrada(t1)
  !========================================================
  !-alocando arrays
  call system_arrays
  !========================================================
  !-convertendo unidades de medida
!  call convert
  !========================================================
  !-definindo coordenadas atomicas
  call structure_prepare(t2)
  !========================================================
  !-definindo vetores da rede reciproca
  call reciprocal
  !========================================================
  select case(method)
  case('@MDPREPARE')
     !-calculando Dinâmica molecular
     call opt
  case('@MDRUNNING')
     !-otimizando sistema
     call md(t3)
  end select
  !========================================================
  !-calculando propriedades
  !if(ebndtot.eq.1)call solving(t4)
  !========================================================
  !-imprimindo informacoes em ficheiro de saida
  call printfile(t0,t1,t2,t3,t4)
  !========================================================

end program HICOLM
