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
module amber

contains

  subroutine amber_vdw(p1,p2,prms)

    implicit none

    integer i
    real(8) epsi(2),ri(2),prms(2)
    character(2) p1,p2,pa

    !-atribuindo valores iniciais

    do i=1,2
       epsi(i)=0.d0
       ri(i)=0.d0
    end do

    !-calculando parametros de LJ de acordo com a regra de Lorentz-Berttelot

    pa=p1

    do i=1,2
       select case(pa)
       case('H ')
          epsi(i)=0.0157d0
          ri(i)=0.6000d0
       case('HO')
          epsi(i)=0.0000d0
          ri(i)=0.0000d0
       case('HS')
          epsi(i)=0.0157d0
          ri(i)=0.6000d0
       case('HC')
          epsi(i)=0.0157d0
          ri(i)=1.4870d0
       case('H1')
          epsi(i)=0.0157d0
          ri(i)=1.3870d0
       case('H2')
          epsi(i)=0.0157d0
          ri(i)=1.2870d0
       case('H3')
          epsi(i)=0.0157d0
          ri(i)=1.1870d0
       case('HP')
          epsi(i)=0.0157d0
          ri(i)=1.1000d0
       case('HA')
          epsi(i)=0.0150d0
          ri(i)=1.4590d0
       case('H4')
          epsi(i)=0.0150d0
          ri(i)=1.4090d0
       case('H5')
          epsi(i)=0.0150d0
          ri(i)=1.3590d0
       case('HW')
          epsi(i)=0.0000d0
          ri(i)=0.0000d0
       case('HZ')
          epsi(i)=0.0150d0
          ri(i)=1.4590d0
       case('O ')
          epsi(i)=0.2100d0
          ri(i)=1.6612d0
       case('O2')
          epsi(i)=0.2100d0
          ri(i)=1.6612d0
       case('OW')
          epsi(i)=0.1520d0
          ri(i)=1.7683d0
       case('OH')
          epsi(i)=0.2104d0
          ri(i)=1.7210d0
       case('OS')
          epsi(i)=0.1700d0
          ri(i)=1.6837d0
       case('C*')
          epsi(i)=0.0860d0
          ri(i)=1.9080d0
       case('CT')
          epsi(i)=0.1094d0
          ri(i)=1.9080d0
       case('C ')
          epsi(i)=0.0860d0
          ri(i)=1.9080d0
       case('N ')
          epsi(i)=0.1700d0
          ri(i)=1.8240d0
       case('N3')
          epsi(i)=0.1700d0
          ri(i)=1.8240d0
       case('NY')
          epsi(i)=0.1700d0
          ri(i)=1.8240d0
       case('S ')
          epsi(i)=0.2500d0
          ri(i)=2.0000d0
       case('SH')
          epsi(i)=0.2500d0
          ri(i)=2.0000d0
       case('P ')
          epsi(i)=0.2000d0
          ri(i)=2.1000d0
       case('IM')
          epsi(i)=0.1000d0
          ri(i)=2.4700d0
       case('Li')
          epsi(i)=0.0183d0
          ri(i)=1.1370d0
       case('IP')
          epsi(i)=0.00277d0
          ri(i)=1.8680d0
       case('Na')
          epsi(i)=0.00277d0
          ri(i)=1.8680d0
       case('K ')
          epsi(i)=0.000328d0
          ri(i)=2.6580d0
       case('Rb')
          epsi(i)=0.00017d0
          ri(i)=2.9560d0
       case('Cs')
          epsi(i)=0.0000806d0
          ri(i)=3.3950d0
       case('MG')
          epsi(i)=0.8947d0
          ri(i)=0.7926d0
       case('C0')
          epsi(i)=0.459789d0
          ri(i)=1.7131d0
       case('Zn')
          epsi(i)=0.0125d0
          ri(i)=1.1000d0
       case('F ')
          epsi(i)=0.0610d0
          ri(i)=1.7500d0
       case('Cl')
          epsi(i)=0.2650d0
          ri(i)=1.9480d0
       case('Br')
          epsi(i)=0.3200d0
          ri(i)=2.2200d0
       case('I ')
          epsi(i)=0.4000d0
          ri(i)=2.3500d0
       case('IB')
          epsi(i)=0.1000d0
          ri(i)=5.0000d0
       case('LP')
          epsi(i)=0.0000d0
          ri(i)=0.0000d0
       end select
       pa=p2
    end do

    prms(1)=sqrt(epsi(1)*epsi(2))
    prms(2)=ri(1)+ri(2)

!    prms(2)=prms(2)/2**(1.d0/6.d0) !convertendo parametro p/ outra versao de LJ

    return

  end subroutine amber_vdw

  subroutine amber_bonds(p1,p2,prms)

    implicit none

    integer i
    real(8) prms(2)
    character(2) p1,p2,pa,pb

    !-atribuindo valores iniciais

    prms(1)=0.d0
    prms(2)=0.d0

    !-calculando parametros de ligação

    pa=p1
    pb=p2

    do i=1,2
       select case(pa)
       case('OW')
          select case(pb)
          case('HW')
             prms(1)=553.0d0
             prms(2)=0.9572d0
          end select
       case('C ')
          select case(pb)
          case('C ')
             prms(1)=310.0d0
             prms(2)=1.525d0
          case('CA')
             prms(1)=469.0d0
             prms(2)=1.409d0
          case('CB')
             prms(1)=447.0d0
             prms(2)=1.419d0
          case('CM')
             prms(1)=410.0d0
             prms(2)=1.444d0
          case('CT')
             prms(1)=317.0d0
             prms(2)=1.522d0
          case('N ')
             prms(1)=490.0d0
             prms(2)=1.335d0
          case('N*')
             prms(1)=424.0d0
             prms(2)=1.383d0
          case('NA')
             prms(1)=418.0d0
             prms(2)=1.388d0
          case('NC')
             prms(1)=457.0d0
             prms(2)=1.358d0
          case('O ')
             prms(1)=570.0d0
             prms(2)=1.229d0
          case('O2')
             prms(1)=656.0d0
             prms(2)=1.250d0
          case('OH')
             prms(1)=450.0d0
             prms(2)=1.364d0
          case('OS')
             prms(1)=450.0d0
             prms(2)=1.323d0
          case('H4')
             prms(1)=367.0d0
             prms(2)=1.080d0
          case('H5')
             prms(1)=367.0d0
             prms(2)=1.080d0
          end select
       case('CA')
          select case(pb)
          case('CA')
             prms(1)=469.0d0
             prms(2)=1.400d0
          case('CB')
             prms(1)=469.0d0
             prms(2)=1.404d0
          case('CM')
             prms(1)=427.0d0
             prms(2)=1.433d0
          case('CN')
             prms(1)=469.0d0
             prms(2)=1.400d0
          case('CT')
             prms(1)=317.0d0
             prms(2)=1.510d0
          case('HA')
             prms(1)=367.0d0
             prms(2)=1.080d0
          case('H4')
             prms(1)=367.0d0
             prms(2)=1.080d0
          case('N2')
             prms(1)=481.0d0
             prms(2)=1.340d0
          case('NA')
             prms(1)=427.0d0
             prms(2)=1.381d0
          case('NC')
             prms(1)=483.0d0
             prms(2)=1.339d0
          case('OH')
             prms(1)=450.0d0
             prms(2)=1.364d0
          end select
       case('CB')
          select case(pb)
          case('CB')
             prms(1)=520.0d0
             prms(2)=1.370d0
          case('N*')
             prms(1)=436.0d0
             prms(2)=1.374d0
          case('NB')
             prms(1)=414.0d0
             prms(2)=1.391d0
          case('NC')
             prms(1)=461.0d0
             prms(2)=1.354d0
          case('CN')
             prms(1)=447.0d0
             prms(2)=1.419d0
          end select
       case('CD')
          select case(pb)
          case('HA')
             prms(1)=367.0d0
             prms(2)=1.080d0
          case('CD')
             prms(1)=469.0d0
             prms(2)=1.400d0
          case('CM')
             prms(1)=549.0d0
             prms(2)=1.350d0
          case('CT')
             prms(1)=317.0d0
             prms(2)=1.510d0
          end select
       case('CK')
          select case(pb)
          case('H5')
             prms(1)=367.0d0
             prms(2)=1.080d0
          case('N*')
             prms(1)=440.0d0
             prms(2)=1.371d0
          case('NB')
             prms(1)=529.0d0
             prms(2)=1.304d0
          end select
       case('CM')
          select case(pb)
          case('CM')
             prms(1)=549.0d0
             prms(2)=1.350d0
          case('CT')
             prms(1)=317.0d0
             prms(2)=1.510d0
          case('HA')
             prms(1)=367.0d0
             prms(2)=1.080d0
          case('H4')
             prms(1)=367.0d0
             prms(2)=1.080d0
          case('H5')
             prms(1)=367.0d0
             prms(2)=1.080d0
          case('N*')
             prms(1)=448.0d0
             prms(2)=1.365d0
          case('OS')
             prms(1)=480.0d0
             prms(2)=1.240d0
          end select
       case('CQ')
          select case(pb)
          case('H5')
             prms(1)=367.0d0
             prms(2)=1.080d0
          case('NC')
             prms(1)=502.0d0
             prms(2)=1.324d0
          end select
       case('CT')
          select case(pb)
          case('CT')
             prms(1)=310.0d0
             prms(2)=1.526d0
          case('HC')
             prms(1)=340.0d0
             prms(2)=1.090d0
          case('H1')
             prms(1)=340.0d0
             prms(2)=1.090d0
          case('H2')
             prms(1)=340.0d0
             prms(2)=1.090d0
          case('H3')
             prms(1)=340.0d0
             prms(2)=1.090d0
          case('HP')
             prms(1)=340.0d0
             prms(2)=1.090d0
          case('N*')
             prms(1)=337.0d0
             prms(2)=1.475d0
          case('N2')
             prms(1)=337.0d0
             prms(2)=1.463d0
          case('OH')
             prms(1)=320.0d0
             prms(2)=1.410d0
          case('OS')
             prms(1)=320.0d0
             prms(2)=1.410d0
          case('N ')
             prms(1)=337.0d0
             prms(2)=1.449d0
          case('N3')
             prms(1)=367.0d0
             prms(2)=1.471d0
          case('NT')
             prms(1)=367.0d0
             prms(2)=1.471d0
          case('S ')
             prms(1)=227.0d0
             prms(2)=1.810d0
          case('SH')
             prms(1)=237.0d0
             prms(2)=1.810d0
          case('CY')
             prms(1)=400.0d0
             prms(2)=1.458d0
          case('CZ')
             prms(1)=400.0d0
             prms(2)=1.459d0
          end select
       case('C*')
          select case(pb)
          case('HC')
             prms(1)=367.0d0
             prms(2)=1.080d0
          case('CB')
             prms(1)=388.0d0
             prms(2)=1.459d0
          case('CT')
             prms(1)=317.0d0
             prms(2)=1.495d0
          case('CW')
             prms(1)=546.0d0
             prms(2)=1.352d0
          end select
       case('CC')
          select case(pb)
          case('CT')
             prms(1)=317.0d0
             prms(2)=1.504d0
          case('CV')
             prms(1)=512.0d0
             prms(2)=1.375d0
          case('CW')
             prms(1)=518.0d0
             prms(2)=1.371d0
          case('NA')
             prms(1)=422.0d0
             prms(2)=1.385d0
          case('NB')
             prms(1)=410.0d0
             prms(2)=1.394d0
          end select
       case('CN')
          select case(pb)
          case('NA')
             prms(1)=428.0d0
             prms(2)=1.380d0
          end select
       case('CR')
          select case(pb)
          case('H5')
             prms(1)=367.0d0
             prms(2)=1.080d0
          case('NA')
             prms(1)=477.0d0
             prms(2)=1.343d0
          case('NB')
             prms(1)=488.0d0
             prms(2)=1.335d0
          end select
       case('CV')
          select case(pb)
          case('H4')
             prms(1)=367.0d0
             prms(2)=1.080d0
          case('NB')
             prms(1)=410.0d0
             prms(2)=1.394d0
          end select
       case('CW')
          select case(pb)
          case('H4')
             prms(1)=367.0d0
             prms(2)=1.080d0
          case('NA')
             prms(1)=427.0d0
             prms(2)=1.381d0
          end select
       case('CY')
          select case(pb)
          case('NY')
             prms(1)=600.0d0
             prms(2)=1.150d0
          end select
       case('CZ')
          select case(pb)
          case('CZ')
             prms(1)=600.0d0
             prms(2)=1.206d0
          case('HZ')
             prms(1)=400.0d0
             prms(2)=1.056d0
          end select
       case('O2')
          select case(pb)
          case('P ')
             prms(1)=525.0d0
             prms(2)=1.480d0
          end select
       case('OH')
          select case(pb)
          case('P ')
             prms(1)=230.0d0
             prms(2)=1.610d0
          end select
       case('OS')
          select case(pb)
          case('P ')
             prms(1)=230.0d0
             prms(2)=1.610d0
          end select
       case('H ')
          select case(pb)
          case('N2')
             prms(1)=434.0d0
             prms(2)=1.010d0
          case('N*')
             prms(1)=434.0d0
             prms(2)=1.010d0
          case('NA')
             prms(1)=434.0d0
             prms(2)=1.010d0
          case('N ')
             prms(1)=434.0d0
             prms(2)=1.010d0
          case('N3')
             prms(1)=434.0d0
             prms(2)=1.010d0
          case('NT')
             prms(1)=434.0d0
             prms(2)=1.010d0
          end select
       case('HO')
          select case(pb)
          case('OH')
             prms(1)=553.0d0
             prms(2)=0.960d0
          case('OS')
             prms(1)=553.0d0
             prms(2)=0.960d0
          end select
       case('HS')
          select case(pb)
          case('SH')
             prms(1)=274.0d0
             prms(2)=1.336d0
          end select
       case('S ')
          select case(pb)
          case('S ')
             prms(1)=166.0d0
             prms(2)=2.038d0
          end select
       case('F ')
          select case(pb)
          case('CT')
             prms(1)=367.0d0
             prms(2)=1.380d0
          case('CA')
             prms(1)=386.0d0
             prms(2)=1.359d0
          end select
       case('Cl')
          select case(pb)
          case('CT')
             prms(1)=232.0d0
             prms(2)=1.766d0
          case('CA')
             prms(1)=193.0d0
             prms(2)=1.727d0
          end select
       case('Br')
          select case(pb)
          case('CT')
             prms(1)=159.0d0
             prms(2)=1.944d0
          case('CA')
             prms(1)=172.0d0
             prms(2)=1.890d0
          end select
       case('I ')
          select case(pb)
          case('CT')
             prms(1)=148.0d0
             prms(2)=2.166d0
          case('CA')
             prms(1)=171.0d0
             prms(2)=2.075d0
          end select
       case('LP')
          select case(pb)
          case('O ')
             prms(1)=600.0d0
             prms(2)=0.200d0
          case('OH')
             prms(1)=600.0d0
             prms(2)=0.200d0
          case('OS')
             prms(1)=600.0d0
             prms(2)=0.200d0
          case('N3')
             prms(1)=600.0d0
             prms(2)=0.200d0
          case('NT')
             prms(1)=600.0d0
             prms(2)=0.200d0
          case('NB')
             prms(1)=600.0d0
             prms(2)=0.200d0
          case('NC')
             prms(1)=600.0d0
             prms(2)=0.200d0
          case('NS')
             prms(1)=600.0d0
             prms(2)=0.700d0
          case('SH')
             prms(1)=600.0d0
             prms(2)=0.700d0
          end select
       end select
       pa=p2
       pb=p1
    end do

    return

  end subroutine amber_bonds

  subroutine amber_bends(p1,p2,p3,prms)

    implicit none

    integer i
    real(8) prms(2)
    character(2) p1,p2,p3,pa,pb,pc

    prms(1)=0.d0
    prms(2)=0.d0

    pa=p1
    pb=p2
    pc=p3

    do i=1,2
       select case(pa)
       case('HW')
          select case(pb)
          case('OW')
             select case(pc)
             case('HW')
                prms(1)=100.0d0
                prms(2)=104.24d0
             end select
          case('HW')
             select case(pc)
             case('OW')
                prms(1)=0.0d0
                prms(2)=127.74d0
             end select
          end select
       case('HC')
          select case(pb)
          case('CT')
             select case(pc)
             case('HC')
                prms(1)=35.0d0
                prms(2)=109.50d0
             case('CT')
                prms(1)=50.0d0
                prms(2)=109.50d0
             end select
          end select
       case('C ')
          select case(pb)
          case('C ')
             select case(pc)
             case('O ')
                prms(1)=80.0d0
                prms(2)=120.00d0
             case('OH')
                prms(1)=80.0d0
                prms(2)=120.00d0
             case('H4')
                prms(1)=50.0d0
                prms(2)=120.00d0
             end select
          case('CB')
             select case(pc)
             case('CB')
                prms(1)=63.0d0
                prms(2)=119.20d0
             case('NB')
                prms(1)=70.0d0
                prms(2)=130.00d0
             end select
          case('CT')
             select case(pc)
             case('H1')
                prms(1)=50.0d0
                prms(2)=109.50d0
             case('HP')
                prms(1)=50.0d0
                prms(2)=109.50d0
             case('HC')
                prms(1)=50.0d0
                prms(2)=109.50d0
             case('N ')
                prms(1)=63.0d0
                prms(2)=110.10d0
             case('N3')
                prms(1)=80.0d0
                prms(2)=111.20d0
             case('CT')
                prms(1)=63.0d0
                prms(2)=111.10d0
             case('OS')
                prms(1)=60.0d0
                prms(2)=109.50d0
             end select
          case('CA')
             select case(pc)
             case('CA')
                prms(1)=63.0d0
                prms(2)=120.00d0
             case('HA')
                prms(1)=50.0d0
                prms(2)=120.00d0
             end select
          case('CM')
             select case(pc)
             case('CM')
                prms(1)=63.0d0
                prms(2)=120.70d0
             case('CT')
                prms(1)=70.0d0
                prms(2)=119.70d0
             case('HA')
                prms(1)=50.0d0
                prms(2)=119.70d0
             case('H4')
                prms(1)=50.0d0
                prms(2)=119.70d0
             end select
          case('N ')
             select case(pc)
             case('CT')
                prms(1)=50.0d0
                prms(2)=121.90d0
             case('H ')
                prms(1)=50.0d0
                prms(2)=120.00d0
             end select
          case('N*')
             select case(pc)
             case('CM')
                prms(1)=70.0d0
                prms(2)=121.60d0
             case('CT')
                prms(1)=70.0d0
                prms(2)=117.60d0
             case('H ')
                prms(1)=50.0d0
                prms(2)=119.20d0
             end select
          case('NA')
             select case(pc)
             case('C ')
                prms(1)=70.0d0
                prms(2)=126.40d0
             case('CA')
                prms(1)=70.0d0
                prms(2)=125.20d0
             case('H ')
                prms(1)=50.0d0
                prms(2)=116.80d0
             end select
          case('NC')
             select case(pc)
             case('CA')
                prms(1)=70.0d0
                prms(2)=120.50d0
             case('LP')
                prms(1)=150.0d0
                prms(2)=120.00d0
             end select
          case('OH')
             select case(pc)
             case('HO')
                prms(1)=50.0d0
                prms(2)=113.00d0
             case('LP')
                prms(1)=150.0d0
                prms(2)=120.00d0
             end select
          case('OS')
             select case(pc)
             case('CT')
                prms(1)=60.0d0
                prms(2)=117.00d0
             case('LP')
                prms(1)=150.0d0
                prms(2)=109.50d0
             end select
          end select
       case('CA')
          select case(pb)
          case('C ')
             select case(pc)
             case('CA')
                prms(1)=63.0d0
                prms(2)=120.00d0
             case('OH')
                prms(1)=70.0d0
                prms(2)=120.00d0
             end select
          case('HW')
             select case(pc)
             case('OW')
                prms(1)=0.0d0
                prms(2)=127.74d0
             end select
          case('CA')
             select case(pc)
             case('HA')
                prms(1)=50.0d0
                prms(2)=120.00d0
             end select
          end select
       case('O ')
          select case(pb)
          case('C ')
             select case(pc)
             case('OH')
                prms(1)=80.0d0
                prms(2)=120.00d0
             end select
          end select
       end select
       pa=p3
       pc=p1
    end do

    return

  end subroutine amber_bends

  subroutine amber_dihedrals(p1,p2,p3,p4,prms)

    implicit none

    real(8) prms(4)
    character(2) p1,p2,p3,p4

    call amber_dihedrals_general(p2,p3,prms)
    call amber_dihedrals_proper(p1,p2,p3,p4,prms)

    return

  end subroutine amber_dihedrals

  subroutine amber_dihedrals_proper(p1,p2,p3,p4,prms)

    implicit none

      integer i
      real(8) prms(4)
      character(2) p1,p2,p3,p4,pa,pb,pc,pd

    pa=p1
    pb=p2
    pc=p3
    pd=p4

    do i=1,2
       select case(pa)
       case('N ')
          select case(pb)
          case('CT')
             select case(pc)
             case('C ')
                select case(pd)
                case('N ')
                   prms(1)=1.0d0
                   prms(2)=2.00d0
                   prms(3)=180.0d0
                   prms(4)=2.0d0
                end select
             end select
          end select
       case('HC')
          select case(pb)
          case('CT')
             select case(pc)
             case('CT')
                select case(pd)
                case('HC')
                   prms(1)=1.0d0
                   prms(2)=0.15d0
                   prms(3)=0.0d0
                   prms(4)=3.0d0
                end select
             end select
          end select
       case('HO')
          select case(pb)
          case('OH')
             select case(pc)
             case('C ')
                select case(pd)
                case('O ')
                   prms(1)=1.0d0
                   prms(2)=1.9d0
                   prms(3)=0.0d0
                   prms(4)=1.0d0
                end select
             end select
          end select
       case('C ')
          select case(pb)
          case('N ')
             select case(pc)
             case('CT')
                select case(pd)
                case('C ')
                   prms(1)=1.0d0
                   prms(2)=0.8d0
                   prms(3)=0.0d0
                   prms(4)=1.0d0
                end select
             end select
          end select
       case('CT')
          select case(pb)
          case('CT')
             select case(pc)
             case('N ')
                select case(pd)
                case('C ')
                   prms(1)=1.0d0
                   prms(2)=0.53d0
                   prms(3)=0.0d0
                   prms(4)=1.0d0
                end select
             end select
          end select
       case('CT')
          select case(pb)
          case('CT')
             select case(pc)
             case('C ')
                select case(pd)
                case('N ')
                   prms(1)=1.0d0
                   prms(2)=0.07d0
                   prms(3)=0.0d0
                   prms(4)=2.0d0
                end select
             end select
          end select
       case('H ')
          select case(pb)
          case('N ')
             select case(pc)
             case('C ')
                select case(pd)
                case('O ')
                   prms(1)=1.0d0
                   prms(2)=2.0d0
                   prms(3)=0.0d0
                   prms(4)=1.0d0
                end select
             end select
          end select
       end select
       pa=p4
       pb=p3
       pc=p2
       pd=p1
    end do

    return

  end subroutine amber_dihedrals_proper

  subroutine amber_dihedrals_general(p2,p3,prms)

    implicit none

    integer i
    real(8) prms(4)
    character(2) p2,p3,pb,pc

    prms(1)=0.d0
    prms(2)=0.d0
    prms(3)=0.d0
    prms(4)=0.d0

    pb=p2
    pc=p3

    do i=1,2
       select case(pb)
       case('C ')
          select case(pc)
          case('C ')
             prms(1)=4.0d0
             prms(2)=14.50d0
             prms(3)=180.0d0
             prms(4)=2.0d0
          case('CA')
             prms(1)=4.0d0
             prms(2)=14.50d0
             prms(3)=180.0d0
             prms(4)=2.0d0
          case('CB')
             prms(1)=4.0d0
             prms(2)=12.00d0
             prms(3)=180.0d0
             prms(4)=2.0d0
          case('CM')
             prms(1)=4.0d0
             prms(2)=8.70d0
             prms(3)=180.0d0
             prms(4)=2.0d0
          case('CT')
             prms(1)=6.0d0
             prms(2)=10.00d0
             prms(3)=180.0d0
             prms(4)=2.0d0
          case('N ')
             prms(1)=4.0d0
             prms(2)=10.00d0
             prms(3)=180.0d0
             prms(4)=2.0d0
          case('N*')
             prms(1)=4.0d0
             prms(2)=5.80d0
             prms(3)=180.0d0
             prms(4)=2.0d0
          case('NA')
             prms(1)=4.0d0
             prms(2)=5.20d0
             prms(3)=180.0d0
             prms(4)=2.0d0
          case('NC')
             prms(1)=2.0d0
             prms(2)=8.00d0
             prms(3)=180.0d0
             prms(4)=2.0d0
          case('O ')
             prms(1)=4.0d0
             prms(2)=11.20d0
             prms(3)=180.0d0
             prms(4)=2.0d0
          case('OS')
             prms(1)=2.0d0
             prms(2)=5.40d0
             prms(3)=180.0d0
             prms(4)=2.0d0
          case('OH')
             prms(1)=2.0d0
             prms(2)=4.60d0
             prms(3)=180.0d0
             prms(4)=2.0d0
          end select
       case('CA')
          select case(pc)
          case('CA')
             prms(1)=4.0d0
             prms(2)=14.50d0
             prms(3)=180.0d0
             prms(4)=2.0d0
          case('CB')
             prms(1)=4.0d0
             prms(2)=14.00d0
             prms(3)=180.0d0
             prms(4)=2.0d0
          case('CM')
             prms(1)=4.0d0
             prms(2)=10.20d0
             prms(3)=180.0d0
             prms(4)=2.0d0
          case('CN')
             prms(1)=4.0d0
             prms(2)=14.50d0
             prms(3)=180.0d0
             prms(4)=2.0d0
          case('CT')
             prms(1)=6.0d0
             prms(2)=0.00d0
             prms(3)=0.0d0
             prms(4)=2.0d0
          case('N2')
             prms(1)=4.0d0
             prms(2)=9.60d0
             prms(3)=180.0d0
             prms(4)=2.0d0
          case('NA')
             prms(1)=4.0d0
             prms(2)=6.00d0
             prms(3)=180.0d0
             prms(4)=2.0d0
          case('NC')
             prms(1)=2.0d0
             prms(2)=9.60d0
             prms(3)=180.0d0
             prms(4)=2.0d0
          case('OH')
             prms(1)=2.0d0
             prms(2)=1.80d0
             prms(3)=180.0d0
             prms(4)=2.0d0
          end select
       case('CB')
          select case(pc)
          case('CB')
             prms(1)=4.0d0
             prms(2)=21.80d0
             prms(3)=180.0d0
             prms(4)=2.0d0
          case('CN')
             prms(1)=4.0d0
             prms(2)=12.00d0
             prms(3)=180.0d0
             prms(4)=2.0d0
          case('N*')
             prms(1)=4.0d0
             prms(2)=6.60d0
             prms(3)=180.0d0
             prms(4)=2.0d0
          case('NB')
             prms(1)=2.0d0
             prms(2)=5.10d0
             prms(3)=180.0d0
             prms(4)=2.0d0
          case('NC')
             prms(1)=2.0d0
             prms(2)=8.30d0
             prms(3)=180.0d0
             prms(4)=2.0d0
          end select
       case('CC')
          select case(pc)
          case('CT')
             prms(1)=6.0d0
             prms(2)=0.00d0
             prms(3)=0.0d0
             prms(4)=2.0d0
          case('CV')
             prms(1)=4.0d0
             prms(2)=20.60d0
             prms(3)=180.0d0
             prms(4)=2.0d0
          case('CW')
             prms(1)=4.0d0
             prms(2)=21.50d0
             prms(3)=180.0d0
             prms(4)=2.0d0
          case('NA')
             prms(1)=4.0d0
             prms(2)=5.60d0
             prms(3)=180.0d0
             prms(4)=2.0d0
          case('NB')
             prms(1)=2.0d0
             prms(2)=4.80d0
             prms(3)=180.0d0
             prms(4)=2.0d0
          end select
       case('CD')
          select case(pc)
          case('CD')
             prms(1)=4.0d0
             prms(2)=4.00d0
             prms(3)=180.0d0
             prms(4)=2.0d0
          case('CT')
             prms(1)=6.0d0
             prms(2)=0.00d0
             prms(3)=0.0d0
             prms(4)=2.0d0
          case('CM')
             prms(1)=4.0d0
             prms(2)=26.00d0
             prms(3)=180.0d0
             prms(4)=2.0d0
          end select
       case('CK')
          select case(pc)
          case('N*')
             prms(1)=4.0d0
             prms(2)=6.80d0
             prms(3)=180.0d0
             prms(4)=2.0d0
          case('NB')
             prms(1)=2.0d0
             prms(2)=20.00d0
             prms(3)=180.0d0
             prms(4)=2.0d0
          end select
       case('CM')
          select case(pc)
          case('CM')
             prms(1)=4.0d0
             prms(2)=26.60d0
             prms(3)=180.0d0
             prms(4)=2.0d0
          case('CT')
             prms(1)=6.0d0
             prms(2)=0.00d0
             prms(3)=0.0d0
             prms(4)=3.0d0
          case('N*')
             prms(1)=4.0d0
             prms(2)=7.40d0
             prms(3)=180.0d0
             prms(4)=2.0d0
          case('OS')
             prms(1)=2.0d0
             prms(2)=2.10d0
             prms(3)=180.0d0
             prms(4)=2.0d0
          end select
       case('CN')
          select case(pc)
          case('NA')
             prms(1)=4.0d0
             prms(2)=6.10d0
             prms(3)=180.0d0
             prms(4)=2.0d0
          end select
       case('CQ')
          select case(pc)
          case('NC')
             prms(1)=2.0d0
             prms(2)=13.60d0
             prms(3)=180.0d0
             prms(4)=2.0d0
          end select
       case('CT')
          select case(pc)
          case('CT')
             prms(1)=9.0d0
             prms(2)=1.40d0
             prms(3)=0.0d0
             prms(4)=3.0d0
          case('CY')
             prms(1)=3.0d0
             prms(2)=0.00d0
             prms(3)=0.0d0
             prms(4)=1.0d0
          case('CZ')
             prms(1)=3.0d0
             prms(2)=0.00d0
             prms(3)=0.0d0
             prms(4)=1.0d0
          case('N ')
             prms(1)=6.0d0
             prms(2)=0.00d0
             prms(3)=0.0d0
             prms(4)=2.0d0
          case('N*')
             prms(1)=6.0d0
             prms(2)=0.00d0
             prms(3)=0.0d0
             prms(4)=2.0d0
          case('N2')
             prms(1)=6.0d0
             prms(2)=0.00d0
             prms(3)=0.0d0
             prms(4)=2.0d0
          case('NT')
             prms(1)=6.0d0
             prms(2)=1.80d0
             prms(3)=0.0d0
             prms(4)=3.0d0
          case('N3')
             prms(1)=9.0d0
             prms(2)=1.40d0
             prms(3)=0.0d0
             prms(4)=3.0d0
          case('OH')
             prms(1)=3.0d0
             prms(2)=0.50d0
             prms(3)=0.0d0
             prms(4)=3.0d0
          case('OS')
             prms(1)=3.0d0
             prms(2)=1.15d0
             prms(3)=0.0d0
             prms(4)=3.0d0
          case('S ')
             prms(1)=3.0d0
             prms(2)=1.00d0
             prms(3)=0.0d0
             prms(4)=3.0d0
          case('SH')
             prms(1)=3.0d0
             prms(2)=0.75d0
             prms(3)=0.0d0
             prms(4)=3.0d0
          end select
       case('C*')
          select case(pc)
          case('CB')
             prms(1)=4.0d0
             prms(2)=6.70d0
             prms(3)=180.0d0
             prms(4)=2.0d0
          case('CT')
             prms(1)=6.0d0
             prms(2)=0.00d0
             prms(3)=0.0d0
             prms(4)=2.0d0
          case('CW')
             prms(1)=4.0d0
             prms(2)=26.10d0
             prms(3)=180.0d0
             prms(4)=2.0d0
          end select
       case('CR')
          select case(pc)
          case('NA')
             prms(1)=4.0d0
             prms(2)=9.30d0
             prms(3)=180.0d0
             prms(4)=2.0d0
          case('NB')
             prms(1)=2.0d0
             prms(2)=10.00d0
             prms(3)=180.0d0
             prms(4)=2.0d0
          end select
       case('CV')
          select case(pc)
          case('NB')
             prms(1)=2.0d0
             prms(2)=4.80d0
             prms(3)=180.0d0
             prms(4)=2.0d0
          end select
       case('CW')
          select case(pc)
          case('NA')
             prms(1)=4.0d0
             prms(2)=6.00d0
             prms(3)=180.0d0
             prms(4)=2.0d0
          end select
       case('OH')
          select case(pc)
          case('P ')
             prms(1)=3.0d0
             prms(2)=0.75d0
             prms(3)=0.0d0
             prms(4)=3.0d0
          end select
       case('OS')
          select case(pc)
          case('P ')
             prms(1)=3.0d0
             prms(2)=0.75d0
             prms(3)=0.0d0
             prms(4)=3.0d0
          end select
       end select
       pb=p3
       pc=p2
    end do

    return

  end subroutine amber_dihedrals_general

end module amber
