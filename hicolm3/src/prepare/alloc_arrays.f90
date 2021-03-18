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
