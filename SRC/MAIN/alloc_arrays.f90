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
module alloc_arrays
  !******************************************************************************************
  ! Alocação estática de arrays                                                             *
  !                                                                                         *
  ! Flaviano Williams Fernandes, 05 de dezembro de 2018                                     *
  !******************************************************************************************

  use input

  integer, allocatable :: bondsmlc(:),bendsmlc(:),torsmlc(:)

  real(8) strmax,str(6),cxx(6,6),exx(6),sxx(6,6)

contains

  subroutine system_arrays
    !****************************************************************************************
    ! Alocacao de variaveis globais                                                         *
    !****************************************************************************************
    implicit none

    integer ierr

    !-alocando arrays referentes ao Campo de Força
    !
    allocate(bondsmlc(nmolec),bendsmlc(nmolec),torsmlc(nmolec),stat=ierr)
    !
    !-alocando variaveis canonicas
    !
    if(ierr.ne.0)stop 'system_arrays: allocation failed'

    return

  end subroutine system_arrays

end module alloc_arrays
