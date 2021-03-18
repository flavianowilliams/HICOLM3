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
  use vdw_module
  use coulomb_module

  contains

  subroutine sf_calc(envdw,encoul,virvdw,vircoul)

    implicit none

    integer i,j,k,ni,nj,ix,ixx,np,i1,i2
    real(8) xvz,yvz,zvz
    real(8) envdw,encoul,virvdw,vircoul

    np=0
    do i=1,nmolec
       do j=1,ntmolec(i)
          do i1=1,nxmolec(i)
             ni=np+i1
             do i2=i1+1,nxmolec(i)
                nj=np+i2
                do k=1,bendscnt(i)
                   do ix=1,3
                      do ixx=ix+1,3
                         if(ni.eq.molbend(i,k,ix).and.nj.eq.molbend(i,k,ixx))goto 1
                         if(ni.eq.molbend(i,k,ixx).and.nj.eq.molbend(i,k,ix))goto 1
                      end do
                   end do
                end do
                do k=1,(torscnt(i)+itorscnt(i))
                   do ix=1,4
                      do ixx=ix+1,4
                         if(i1.eq.moltors(i,k,ix).and.i2.eq.moltors(i,k,ixx))goto 1
                         if(i1.eq.moltors(i,k,ixx).and.i2.eq.moltors(i,k,ix))goto 1
                      end do
                   end do
                end do
                call mic(ni,nj,xvz,yvz,zvz)
                call vdw_sf(ni,nj,i,xvz,yvz,zvz,envdw,virvdw)
                call coulomb_sf(ni,nj,i,xvz,yvz,zvz,encoul,vircoul)
1               continue
             end do
          end do
          np=np+nxmolec(i)
       end do
    end do

    return

  end subroutine sf_calc

end module scale_factor_module
