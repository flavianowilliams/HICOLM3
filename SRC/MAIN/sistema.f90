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
  parameter (natmax=2000)
  parameter (ntpmax=10)
  parameter (nkmax=150000)
  parameter (nprllmmax=4)
  parameter (nparammax=10)

  parameter (nbandmax=5000)

  parameter (lmax=4)
  parameter (mmax=9)

  parameter (nkrlxmax=100)
  parameter (dkrlxmax=100)

  parameter (molecmax=500)

  parameter (bendmax=20)
  parameter (bondmax=30)
  parameter (torsmax=10)

end module sistema
