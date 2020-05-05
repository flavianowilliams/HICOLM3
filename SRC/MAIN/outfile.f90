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
module outfile

  use input
  use estrutura

contains

  subroutine printfile(t0,t1,t2,t3,t4)

    implicit none

    integer i
    real(8) t0,t1,t2,t3,t4,tf

    call cpu_time(tf)

    tf=tf-t0

    write(6,*)('-',i=1,70)
    write(6,'(''         input = '',f8.4)') t1
    write(6,'(''   coordinates = '',f8.4)') t2
    write(6,'(''  m. mechanics = '',f8.4)') t3
    write(6,'(''       solving = '',f8.4)') t4
    write(6,*)('-',i=1,70)
    write(6,'(''  Elapsed time = '',f8.2,'' seconds'')') tf
    write(6,*)
    write(6,'(93a1)')('#',i=1,93)
    write(6,*)('END!',i=1,23)
    write(6,'(93a1)')('#',i=1,93)
    write(6,*)

  !========================================================
  !-aviso de direitos autorais

    write(6,'(5x,a46)')'Copyright (C) 2019 Flaviano Williams Fernandes'
    write(6,*)
    write(6,'(5x,a68)')'This program is free software: you can redistribute it and/or modify'
    write(6,'(5x,a68)')'it under the terms of the GNU General Public License as published by'
    write(6,'(5x,a40)')'the Free Software Foundation, version 3.'
    write(6,*)
    write(6,'(5x,a63)')'This program is distributed in the hope that it will be useful,'
    write(6,'(5x,a62)')'but WITHOUT ANY WARRANTY; without even the implied warranty of'
    write(6,'(5x,a61)')'MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the'
    write(6,'(5x,a44)')'GNU General Public License for more details.'
    write(6,*)
    write(6,'(5x,a65)')'You should have received a copy of the GNU General Public License'
    write(6,'(5x,a66)')'along with this program.  If not, see <https://www.gnu.org/licenses/>.'
    write(6,*)

    return

  end subroutine printfile

end module outfile
