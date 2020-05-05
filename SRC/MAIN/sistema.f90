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
module sistema
  !********************************************************************************************
  ! modulo de dimensionamento estático de arrays                                              *
  !                                                                                           *
  ! Flaviano Williams Fernandes, 05 de dezembro de 2018                                       *
  !********************************************************************************************
  public

  integer ntpmax,natmax,nkmax,nbandmax,nparammax,nprllmmax,nkrlxmax,dkrlxmax
  integer lmax,mmax,molecmax,bendmax,bondmax,torsmax
  !
  !-----------------------------------------------------------------
  !-nparammax: Limite de parametros SK
  !- nbandmax: Limite de bandas de calculo DFT
  !-    nkmax: Limite de pontos K
  !-   natmax: Limite de atomos
  !-   ntpmax: Limite de tipos atomicos
  !-   lpmmax: Limite de subniveis magneticos
  !-     mmax: Limite do subnivel magnetico
  !- molecmax: Limite de moleculas
  !-  bendmax: Limite de parametros de deformacao angular
  !-  bondmax: Limite de parametros de estiramento
  !-  torsmax: Limite de parametros de torsão
  !- nkrlxmax: Limite de pontos K para ajuste dos parametros SK
  !- dkrlxmax: Limite de pontos K para ajuste dos parametros SK
  !----------------------------------------------------------------
  !
  parameter (natmax=3000)
  parameter (ntpmax=25)
  parameter (nkmax=15)
  parameter (nprllmmax=4)
  parameter (nparammax=10)

  parameter (nbandmax=5)

  parameter (lmax=4)
  parameter (mmax=9)

  parameter (nkrlxmax=10)
  parameter (dkrlxmax=10)

  parameter (molecmax=1000)

  parameter (bendmax=50)
  parameter (bondmax=50)
  parameter (torsmax=50)

end module sistema
