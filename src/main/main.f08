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
program HICOLM
  !*******************************************************************************************
  !programa principal responsavel pelo instanciamento dos objetos                            *
  !*******************************************************************************************

  use prepare_module           ! prepare class to prepare the physical environment
  use moleculardynamics_module ! prepare class to molecular dynamics procedure

  implicit none

  integer       :: i,j,i0
  real(8)       :: t0,t1,t2,t3,t4
  real(8)       :: sf_coul,sf_vdw
  real(8)       :: drx,dry,drz,drmax
  character(10) :: host,time
  character(8)  :: date,in
  logical       :: lval

  type(prepare)           :: prp ! instanciating prepare object
  type(moleculardynamics) :: md  ! instanciating molecular dynamics object

  call cpu_time(t0)

  !-abrindo ficheiro de dados

  open(5,file='HICOLM.in',status='old')               ! reading input datas
  open(6,file='HICOLM.out',status='unknown')          ! printing output informations

  t1=0.d0
  t2=0.d0
  t3=0.d0
  t4=0.d0

  !-elapsed time information

  call date_and_time(Date=date)
  call date_and_time(TIME=time)
  call hostnm(host)
!  host='undefined'

  !-cabecalho
  !
  write(6,*)
  write(6,'(18x,a57)')'HICOLM: Multi-Methods for Molecules and Condensed Systems'
  write(6,*)
  write(6,'(39x,a14)')'Version: x.x.x'
  write(6,*)
  write(6,'(''Host: '',2x,a10)')host
  write(6,'(''Date: '',2x,a8)')date
  write(6,'(''Time: '',2x,a10)')time
  write(6,*)

  ! assign 1-4 scale factor according AMBER force field
  !
  sf_coul=1.d0/1.2d0
  sf_vdw=1.d0/2.d0

  !========================================================
  !
  !-read MD option
  !
  lval=.false.
  do while (lval.eqv..false.)
     read(5,*,end=1)in
     if(in.eq.'@PREPARE')then
        open(9,file='HICOLM.xyz',status='old')   ! read atomic coordinates
!
        prp=prepare()                            ! definindo valores default
!
        call prp%constants_prepare()             ! definindo constantes
        call prp%molecules()                     ! atribuindo qde moléculas e sitios atomicos
        call prp%set_natom()                     ! calculando qde de sitios atomicos
        call prp%set_latticevectors()            ! lendo coordenadas da celula unitaria
        call prp%set_lattice_constants()         ! calculando constantes de rede
        call prp%set_lattice_angles()            ! calculando angulos de rede
        call prp%set_volume()                    ! calculando volume da supercelula
        call prp%set_symmetry()                  ! calculando grupo de simetria
        call prp%sites()                         ! atribuindo coordenadas atomicas e Z
        call prp%translate()                     ! aplicando translação do sistema coordenadas
        call prp%molecule_prepare()              ! atribuindo informacoes moleculares
        call prp%set_massmol()                   ! calculando massa molecular
        call prp%set_mmolar()                    ! calculando massa molecular
        call prp%set_scale_factor(sf_coul,sf_vdw)! atribuindo fatores escalonamento 1-4
        call prp%set_internal_coordinates()      ! atribuindo coordenadas internas
        call prp%set_parbnd()                    ! atribuindo potenciais de ligacao
        call prp%set_parbend()                   ! atribuindo potenciais angulares
        call prp%set_partors()                   ! atribuindo potenciais de torsão
        call prp%set_paritors()                  ! atribuindo diedros improprios
        call prp%set_parvdw()                    ! atribuindo potenciais de Van der Waals
        call prp%set_extra_parbnd()              ! alterando potenciais de ligacao
        call prp%set_extra_parbend()             ! alterando potenciais angulares
        call prp%set_extra_partors()             ! alterando diedros
        call prp%set_extra_parvdw()              ! alterando potenciais de Van der Waals
        call prp%check()                         ! checando parametros de entrada
        call prp%set_global()                    ! imprimindo propriedades globais
        call prp%print_sys()                     ! imprimindo estrutura em HICOLM.sys
        call prp%print_top()                     ! imprimindo topologia em HICOLM.top
        call prp%print_out()                     ! imprimindo valores em HICOLM.out
        lval=.true.
     elseif(in.eq.'@MD')then
        open(3,file='HICOLM.AXSF',status='unknown')        ! printing coordinates per frame
        open(4,file='atoms.csv',status='unknown')          ! imprimindo informacoes atomicas
        open(7,file='thermodynamics.csv',status='unknown') ! imprimindo informacoes termodin.
        open(8,file='lattice.csv',status='unknown')        ! imprimindo informacoes da rede
!
        md=moleculardynamics()                   ! set default values
!
        call md%constants_prepare()              ! definindo constantes
        call md%set_input()                      ! lendo parametros de entrada em HICOLM.in
        call md%set_molecules()                  ! lendo tipos e qde de moleculas
        call md%set_latticevectors()             ! lendo coordenadas da celula unitaria
        call md%set_natom()                      ! calculando qde de sitios atomicos
        call md%set_atoms()                      ! lendo coordenadas atomicas
        call md%set_topology()                   ! lendo parametros do campo de forca
        if(md%get_restart().gt.0)call md%read_geometry() ! reiniciando simulação
        call md%ccp()                            ! aplicando condicoes de contorno periodicas
        call md%convert_units()                  ! convertendo unidades de medida
        call md%set_mmolar()                     ! calculando massa molecular
        call md%set_global()                     ! calculando a carga total do sistema
        call md%set_zat()                        ! atribuindo numero atomico
        call md%set_mass()                       ! atribuindo massa atomica
        call md%set_qat()                        ! atribuindo cargas atomicas
        call md%set_tpa()                        ! atribuindo tipos atomicos
        call md%set_lattice_constants()          ! calculando constantes de rede
        call md%set_lattice_angles()             ! calculando angulos de rede
        call md%set_symmetry()                   ! calculando grupo de simetria
        call md%set_volume()                     ! calculando volume da supercelula
        call md%print_out()                      ! imprimindo valores em HICOLM.out
        call md%print()                          ! imprimindo parametros da MD
        call md%set_velocity()                   ! atribuindo velocidades iniciais
        call md%neighbour_prepare()              ! preparando lista de vizinhos de Verlet
        call md%verlet_list()                    ! atribuindo lista de vizinhos de Verlet
        call md%interaction_prepare()            ! preparando campo de forca
        call md%set_forcefield()                 ! calculo das interacoes moleculares
        call md%set_time(0.d0)                   ! setando instante inicial
        write(6,'(4x,111a1)')('-',i=1,84)
        write(6,10)'##','STEP','TIME','VOLUME','TEMPERATURE' ,'PRESSURE','E(TOTAL)'
        write(6,'(4x,111a1)')('-',i=1,84)
        drmax=md%get_drcutoff()
        i0=0
        do i=1,md%get_nframes()
           if((i-i0).ge.nint(md%get_drcutoff()/drmax))then
              call md%verlet_list()
              drmax=0.d0
              do j=1,md%get_natom()
                 drx=md%vax(j)*md%get_timestep()+&
                      0.5d0*md%fax(j)*md%get_timestep()**2/md%mass(j)
                 dry=md%vay(j)*md%get_timestep()+&
                      0.5d0*md%fay(j)*md%get_timestep()**2/md%mass(j)
                 drz=md%vaz(j)*md%get_timestep()+&
                      0.5d0*md%faz(j)*md%get_timestep()**2/md%mass(j)
                 drmax=max(drmax,sqrt(drx**2+dry**2+drz**2))
              end do
              i0=i
           end if
           call md%set_forcefield()
           call md%print_geometry(i)
!           if(mod(i,25).eq.0)write(6,20)'MD',i,time*tconv,volume*rconv**3,&
!                temp*teconv,press*pconv,(ekinet+enpot+envdw_corr)*econv
           if(mod(i,25).eq.0)write(6,20)&
                'MD',i,md%get_time()*md%get_tconv(),md%get_volume()*md%get_rconv()**3,&
                (md%get_enpot())*md%get_econv()
           call md%set_time(i*md%get_timestep())
        end do
        write(6,'(4x,111a1)')('-',i=1,84)
        write(6,*)
        lval=.true.
     end if
  end do
  !===========================================================================================
  !
  !-print copyright
  !
    write(6,'(93a1)')('#',i=1,93)
    write(6,*)('END!',i=1,23)
    write(6,'(93a1)')('#',i=1,93)
    write(6,*)
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
    write(6,'(5x,a78)')'LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,'
    write(6,'(5x,a78)')'OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE'
    write(6,'(5x,a10)')'SOFTWARE.'
    write(6,*)
    !=========================================================================================

  stop

1 write(6,*)'ERROR: Method does not found!'
  write(6,*)'Hint: You must choose one of the following methods: @PREPARE or @MD'
  stop

10 format(5x,a2,6x,a4,6x,a5,9x,a6,6x,a10,5x,a8,6x,a8)
20 format(5x,a2,2x,i8,2x,es12.4,2x,es12.4,3(2x,es12.4))

end program HICOLM
