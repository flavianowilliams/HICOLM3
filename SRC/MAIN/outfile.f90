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

    write(6,'(5x,a12)')'MIT License'
    write(6,*)
    write(6,'(5x,a36)')'Copyright (c) 2020 flavianowilliams'
    write(6,*)
    write(6,'(5x,a77)')'Permission is hereby granted, free of charge, to any person obtaining a copy'
    write(6,'(5x,a78)')'of this software and associated documentation files (the "Software"), to deal'
    write(6,'(5x,a77)')'in the Software without restriction, including without limitation the rights'
    write(6,'(5x,a74)')'to use, copy, modify, merge, publish, distribute, sublicense, and/or sell'
    write(6,'(5x,a70)')'copies of the Software, and to permit persons to whom the Software is'
    write(6,'(5x,a57)')'furnished to do so, subject to the following conditions:'
    write(6,*)
    write(6,'(5x,a75)')'THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR'
    write(6,'(5x,a73)')'IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,'
    write(6,'(5x,a76)')'FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE'
    write(6,'(5x,a71)')'AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER'
    write(6,'(5x,a77)')'LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,'
    write(6,*)

    return

  end subroutine printfile

end module outfile
