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
  use optimize_module          ! prepare class to molecular dynamics procedure

  implicit none

  integer       :: i,j,k,i0,i1,nx
  real(8)       :: t0,t1,t2,t3
  real(8)       :: sf_coul,sf_vdw
  real(8)       :: drmax
  real(8)       :: dgg0,dgg,q,gg,alpha
  character(10) :: host,time
  character(8)  :: date
  character(9)  :: in
  logical       :: lval

  type(prepare)           :: prp ! instanciating prepare object
  type(moleculardynamics) :: md  ! instanciating molecular dynamics object
  type(optimize)          :: opt ! instanciating optimize object

  call cpu_time(t0)

  !-abrindo ficheiro de dados

  open(5,file='INPUT',status='old')          ! reading input datas
  open(6,file='hicolm.out',status='unknown') ! printing output informations

  !-elapsed time information

  call date_and_time(Date=date)
  call date_and_time(TIME=time)
  call hostnm(host)
!  host='undefined'
  !
  !-cabecalho
  !
  write(6,*)
  write(6,'(18x,a57)')'HICOLM: Multi-Methods for Molecules and Condensed Systems'
  write(6,*)
  write(6,'(39x,a14)')'Version: x.x.x'
  write(6,*)
  !
  !===========================================================================================
  !
  !-print copyright
  !
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
  !
  ! general information
  !
  write(6,'(''Host: '',2x,a10)')host
  write(6,'(''Date: '',2x,a8)')date
  write(6,*)
  !
  !=========================================================================================
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
        call prp%set_system()                    ! definindo parametros de estrut.
        call prp%set_natom()                     ! calculando qde de sitios atomico
        call prp%set_freeze()                    ! definindo atomos congelados
        call prp%set_lattice_constants()         ! calculando constantes de rede
        call prp%set_lattice_angles()            ! calculando angulos de rede
        call prp%set_volume()                    ! calculando volume da supercelula
        call prp%set_symmetry()                  ! calculando grupo de simetria
        call prp%sites()                         ! atribuindo coordenadas atomicas e Z
        call prp%translate()                     ! aplicando translação do sistema coordenadas
        call prp%set_internal_coordinates()      ! atribuindo coordenadas internas
        call prp%set_topology()
        call prp%check_vdw()
        call prp%check_parbnd()
        call prp%check_parbend()
        call prp%check_partors()
        call prp%check_paritors()
        call prp%set_massmol()                   ! calculando massa molecular
        call prp%set_mmolar()                    ! calculando massa molecular
        call prp%set_scale_factor(sf_coul,sf_vdw)! atribuindo fatores escalonamento 1-4
        call prp%check()                         ! checando parametros de entrada
        call prp%set_global()                    ! imprimindo propriedades globais
        call prp%print_sys()                     ! imprimindo estrutura em SYSTEM
        call prp%print_top()                     ! imprimindo topologia em TOPOLOGY
        call prp%print_out()                     ! imprimindo valores em hicolm.out
        lval=.true.
     elseif(in.eq.'@MD')then
        call cpu_time(t1)
        open(4,file='atoms.csv',status='unknown') ! imprimindo informacoes atomicas
        open(7,file='thermodynamics.csv',status='unknown') ! imprimindo informacoes termodin.
        open(8,file='lattice.csv',status='unknown')        ! imprimindo informacoes da rede
        !
        md=moleculardynamics()                   ! set default values
        !
        call md%constants_prepare()              ! definindo constantes
        call md%set_input()                      ! lendo parametros de entrada em INPUT
        call md%set_latticevectors()             ! lendo coordenadas da celula unitaria
        call md%set_molecules()                  ! lendo tipos e qde de moleculas
        call md%set_natom()                      ! calculando qde de sitios atomicos
        call md%set_nfree()                      ! atribuindo graus de liberdade
        call md%set_topology()                   ! lendo parametros do campo de forca
        call md%set_rcutoff2()                   ! atribuindo raio de corte
        call md%set_canonicalvariables()         ! atribuindo valores iniciais para var. can.
        call md%check()                          ! checando parametros de entrada
        call md%convert_units()                  ! convertendo unidades de medida
        call md%interaction_init()               ! preparando campo de forca
        call md%ccp()                            ! aplicando condicoes de contorno periodicas
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
        call md%set_vdwcorr()                    ! calculando correção de Van der Waals
        call md%print()                          ! imprimindo parametros da MD
        call md%neighbour_prepare()              ! preparando lista de vizinhos de Verlet
        call md%verlet_list()                    ! atribuindo lista de vizinhos de Verlet
        call md%set_verlchk()                    ! atribuindo passos p/ atualizacao da lista
        call md%print_out()                      ! imprimindo valores em hicolm.out
        call md%set_forcefield()                 ! calculo das interacoes moleculares
        call md%set_ekinetic()                   ! atribuindo energia cinecita
        call md%set_pressure()                   ! atribuindo pressao
        call md%set_temperature()                ! atribuindo temperatura
        call md%set_sigma()                      ! atribuindo sigma para o termostato Berend.
        call md%set_qmass()                      ! atribuindo massa do calorimetro
        call md%set_pmass()                      ! atribuindo massa do barostato
        call md%set_nhist()                      ! definindo quantidade de frames
        call md%set_time(0.d0)                   ! setando instante inicial
        call cpu_time(t2)

        write(6,*)('#',i=1,93)
        write(6,*)('MD (CYCLE) ',i=1,8)
        write(6,*)('#',i=1,93)
        write(6,*)
        write(6,'(2x,a27,f12.6,1x,a8)')'  Correction of VdW energy:',&
             md%get_encorr()*md%get_econv(),'kcal/mol'
        write(6,'(2x,a27,f12.6,1x,a3)')'Correction of VdW pressure:',&
             (-md%get_vircorr()/(3.d0*md%get_volume()))*md%get_pconv(),'atm'
        write(6,*)
        write(6,'(4x,111a1)')('-',i=1,84)
        write(6,10)'##','STEP','TIME','VOLUME','TEMPERATURE' ,'PRESSURE','E(TOTAL)'
        write(6,'(4x,111a1)')('-',i=1,84)
        drmax=md%get_drcutoff()
        i0=0
        i1=md%get_verlchk()
        do i=1,md%get_nstep()
           if(md%get_ensble().eq.'nve')then
              call md%set_nve()
           elseif(md%get_ensble().eq.'nvt')then
              if(md%get_ensble_mt().eq.'berendsen')then
                 call md%set_nvt_berendsen()
              elseif(md%get_ensble_mt().eq.'hoover')then
                 call md%set_nvt_nosehoover()
              end if
           elseif(md%get_ensble().eq.'npt')then
              if(md%get_ensble_mt().eq.'berendsen')then
                 call md%set_npt_berendsen()
              elseif(md%get_ensble_mt().eq.'hoover')then
                 call md%set_npt_nosehoover()
              end if
           end if
           if(i1.ge.md%get_verlchk())then
              call md%verlet_list()
              call md%set_verlchk()
              i0=i
           end if
           call md%set_verlchk()
           call md%print_geometry(i)
           call md%print_frames(i)
           call md%print_dataframes(i)
           if(mod(i,25).eq.0)write(6,20)&
                'MD',i,md%get_time()*md%get_tconv(),md%get_volume()*md%get_rconv()**3,&
                md%get_temperature()*md%get_teconv(),md%get_pressure()*md%get_pconv(),&
                md%get_etotal()*md%get_econv()
           call md%set_time(i*md%get_timestep())
           i1=i-i0
        end do
        write(6,'(4x,111a1)')('-',i=1,84)
        write(6,*)
        call cpu_time(t3)
        write(6,'(5x,a19)')'CPU time evaluation'
        write(6,'(4x,84a1)')('-',i=1,84)
        write(6,'(5x,a25,1x,f8.4,1x,a7)')'To the preparing system =',(t2-t1),'seconds'
        write(6,'(5x,14x,a11,1x,f8.4,1x,a7)')'MD(cycle) =',(t3-t2)/md%get_nstep(),'seconds'
        write(6,'(4x,84a1)')('-',i=1,84)
        write(6,'(5x,11x,a14,1x,i8,1x,a7)')'Elapsed time =',nint(t3-t1),'seconds'
        write(6,*)
        lval=.true.
     elseif(in.eq.'@OPTIMIZE')then
        call cpu_time(t1)
        open(4,file='optimize.csv',status='unknown')       ! imprimindo data frame da otimz.
        !
        opt=optimize()                            ! set default values
        !
        call opt%constants_prepare()              ! definindo constantes
        call opt%set_latticevectors()             ! lendo coordenadas da celula unitaria
        call opt%set_molecules()                  ! lendo tipos e qde de moleculas
        call opt%set_natom()                      ! calculando qde de sitios atomicos
        call opt%set_nfree()                      ! atribuindo graus de liberdade
        call opt%set_inopt()                      ! lendo parametros de entrada em INPUT
        call opt%set_topology()                   ! lendo parametros do campo de forca
        call opt%set_rcutoff2()                   ! atribuindo raio de corte
        call opt%check()                          ! checando parametros de entrada
        call opt%set_canonicalvariables()         ! atribuindo valores iniciais para var. can.
        call opt%convert_units()                  ! convertendo unidades de medida
        call opt%set_lattice_constants()          ! calculando constantes de rede
        call opt%set_lattice_angles()             ! calculando angulos de rede
        call opt%set_symmetry()                   ! calculando grupo de simetria
        call opt%set_volume()                     ! calculando volume da supercelula
        call opt%ccp()                            ! aplicando condicoes de contorno periodicas
        call opt%set_mmolar()                     ! calculando massa molecular
        call opt%set_global()                     ! calculando a carga total do sistema
        call opt%set_zat()                        ! atribuindo numero atomico
        call opt%set_qat()                        ! atribuindo cargas atomicas
        call opt%set_tpa()                        ! atribuindo tipos atomicos
        call opt%gd_init()                        ! preparando otimizacao
        call opt%print()                          ! imprimindo parametros da otimizacao
        call opt%neighbour_prepare()              ! preparando lista de vizinhos de Verlet
        call opt%verlet_list()                    ! atribuindo lista de vizinhos de Verlet
        call opt%set_loop()                       ! calculando residuo e hessiana
        call cpu_time(t2)
        !
        write(6,*)('#',i=1,93)
        write(6,*)('OPTIMIZING ',i=1,8)
        write(6,*)('#',i=1,93)
        write(6,*)
        write(6,'(4x,111a1)')('-',i=1,84)
        write(6,30)'##','STEP','LISEARCH','E(TOTAL)','MAXFORCE'
        write(6,'(4x,111a1)')('-',i=1,84)
        nx=1
        dgg0=0.0d0
        gg=0.d0
        do i=1,opt%get_nmatrix()
           q=0.d0
           do j=1,opt%get_nmatrix()
              q=q+opt%hess(i,j)*opt%res(j)
           end do
           gg=gg+opt%res(i)*q
           dgg0=dgg0+opt%res(i)**2
        end do
        if(gg.eq.0.d0)stop 'residue reached null value!'
        do i=1,opt%get_nstep()
           gg=0.d0
           dgg=0.d0
           do j=1,opt%get_nmatrix()
              q=0.d0
              do k=1,opt%get_nmatrix()
                 q=q+opt%hess(j,k)*opt%res(k)
              end do
              gg=gg+opt%res(j)*q
              dgg=dgg+opt%res(j)**2
           end do
           if(gg.lt.0.d0)then
              print*,'Hessian did not positive definite. Randomrizing positions...'
              call opt%random_coordinates()
              call opt%ccp()
              call opt%verlet_list()
              call opt%set_loop()
              gg=0.d0
              dgg=0.d0
              do j=1,opt%get_nmatrix()
                 q=0.d0
                 do k=1,opt%get_nmatrix()
                    q=q+opt%hess(j,k)*opt%res(k)
                 end do
                 gg=gg+opt%res(j)*q
                 dgg=dgg+opt%res(j)**2
              end do
           end if
           if(gg.eq.0.d0)stop 'residue reached null value!'
           alpha=dgg/gg
           do j=1,opt%get_natom()
              opt%xa(j)=opt%xa(j)+alpha*opt%res(3*j-2)*opt%freeze(j)
              opt%ya(j)=opt%ya(j)+alpha*opt%res(3*j-1)*opt%freeze(j)
              opt%za(j)=opt%za(j)+alpha*opt%res(3*j)*opt%freeze(j)
           end do
           call opt%ccp()
           call opt%verlet_list()
           call opt%set_loop()
           call opt%set_lsearch(alpha)
           call opt%set_maxforce()
           call opt%print_dataframes(nx)
           write(6,40)'SD',i,alpha*opt%get_rconv()**2/opt%get_econv()&
                ,opt%get_enpot()*opt%get_econv()&
                ,opt%get_maxforce()*(opt%get_econv()/opt%get_rconv())
           if(opt%get_maxforce().le.opt%get_tolerance())then
              call opt%print_geometry()
              write(6,*)
              write(6,*)&
                   'SUCCESS! The optimization procedure reached the convergence criteria.'
              write(6,*)
              goto 2
           end if
           dgg0=dgg
           nx=nx+1
        end do
        call opt%print_geometry()
        write(6,*)
        write(6,*)'Warning: The optimization did not converge to the convergence criteria.'
        write(6,*)
        lval=.true.
     end if
  end do

  goto 2

1 write(6,*)'ERROR: Method does not found!'
  write(6,*)'Hint: You must choose one of the following methods: @PREPARE, @OPTIMIZE or @MD'
  stop

10 format(5x,a2,6x,a4,6x,a5,9x,a6,6x,a10,5x,a8,6x,a8)
20 format(5x,a2,2x,i8,2x,es12.4,2x,es12.4,3(2x,es12.4))
30 format(5x,a2,6x,a4,3x,a10,4x,a10,6x,a8,6x,a8)
40 format(5x,a2,2x,i8,2x,es12.4,2x,es12.4,2x,es12.4)

2 write(6,'(93a1)')('#',i=1,93)
  write(6,*)('END!',i=1,23)
  write(6,'(93a1)')('#',i=1,93)
  write(6,*)

end program HICOLM
