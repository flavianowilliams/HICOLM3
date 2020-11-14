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
module scale_factor_module
  !******************************************************************************************
  ! Fator de escalonamento intramolecular de Van der Waals e coulombiano:                   *
  ! - Energia potencial;                                                                    *
  ! - Forças atômicas;                                                                      *
  ! - Virial;                                                                               *
  ! - Stress.                                                                               *
  !                                                                                         *
  ! Flaviano Williams Fernandes, 05 de dezembro de 2018                                     *
  !******************************************************************************************

  use input
  use estrutura

  contains

  subroutine sf_calc(envdw,encoul,virvdw,vircoul)

    implicit none

    integer i,j,k,ni,nj
    real(8) xvz,yvz,zvz
    real(8) envdw,encoul,virvdw,vircoul

    do i=1,nmolec
       do j=1,ntmolec(i)
          do ni=1,nxmolec(i)
             do nj=ni+1,nxmolec(i)
                do k=1,bondscnt(i)
                   if(ni.eq.molbond(i,k,1).and.nj.ne.molbond(i,k,2))then
                      !call mic(ni,nj,xvz,yvz,zvz)
                      !                      call vdw_sf(envdw,virvdw)
                      !                      call coul_sf(encoul,vircoul)
                      print*,ni,nj
                   elseif(ni.eq.molbond(i,k,2).and.nj.ne.molbond(i,k,1))then
                      !call mic(ni,nj,xvz,yvz,zvz)
                      !                      call vdw_sf(envdw,virvdw)
                      !                      call coul_sf(encoul,vircoul)
                   end if
                end do
             end do
          end do
       end do
    end do

    return

  end subroutine sf_calc

end module scale_factor_module
