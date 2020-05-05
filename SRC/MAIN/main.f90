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
program HICOLM
  !*******************************************************************************************
  !programa principal responsavel por carregar os modulos essenciais para efetuar os cálculos*
  !*******************************************************************************************
  use input                ! parametros de entrada
  use alloc_arrays         ! alocando arrays essenciais
  use estrutura            ! definindo espaço real
  use brillouin            ! definindo zona de brillouin
  use molecular_dynamics   ! metodo Dinâmica Molecular
  use outfile              ! imprimindo informacoes em arquivo out

  implicit none

  real(8) :: t0,t1,t2,t3,t4
  character(10) host,time
  character(8) date

  call cpu_time(t0)

  open(6,file='HICOLM.out',status='unknown') ! imprimindo arquivo out

  t1=0.d0
  t2=0.d0
  t3=0.d0
  t4=0.d0

  !-elapsed time information

  call date_and_time(Date=date)
  call date_and_time(TIME=time)
  call hostnm(host)

  !-cabecalho
  !
  write(6,'(5x,a57)')'HICOLM: Multi-Methods for Molecules and Condensed Systems'
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
  !-definindo coordenadas atomicas
  call coordenadas(t2)
  !========================================================
  !-definindo vetores da rede reciproca
  call reciprocal
  !========================================================
  !-calculando Dinâmica molecular
  call md(t3)
  !========================================================
  !-calculando propriedades
  !if(ebndtot.eq.1)call solving(t4)
  !========================================================
  !-imprimindo informacoes em ficheiro de saida
  call printfile(t0,t1,t2,t3,t4)
  !========================================================

end program HICOLM
