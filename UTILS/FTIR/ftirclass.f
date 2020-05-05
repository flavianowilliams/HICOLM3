      program ftirclass
**********************************************************************
*     Programa FTIR-Class para o calculo do espectro vibracional     *
*     Flaviano Williams Fernandes, 20 de agosto de 2014              *
**********************************************************************
      use sizes
      use system
      use cmmolec
      use utilsrange
      use dipmol
      use corr
      use tcf
      use ftir
      use printfile
c
      implicit none
c
      integer i,w,check,cmserr
      real(kind=4) t0,dt,cmsprec
      real(kind=4) strmm,flexmm,int,dm(fmax,molmax,nmmax,iz)
c
      open(unit=ird,file="TBMD.md",status="old")
      open(unit=iwrz,file="infrared.dat",status="unknown")
      open(unit=iwtt,file="vacf.dat",status="unknown")
      open(unit=iwrx,file="probability.dat",status="unknown")
      open(unit=iwrt,file="ftirclass.out",status="unknown")
c
      call cpu_time(t0)
c
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      call input
c
      write(*,*)'-> calculando CMS e desligando CC'
c
      check=0
      cmserr=0
      do w=1,fmax
         call coord(check,w)
         call molec(check,w,cmserr,cmsprec)
         if(check.eq.1)exit
      end do
c
      if(ddw.gt.tmax)ddw=tmax   !acertando o limite máximo p/ ddw
c
      write(*,*)'-> Calculo CMS e desligando CC: ok'
c--------------------------------------------------------------------
      dt=dtime*float(timestep(tmax)-timestep(1))
c
      write(iwrt,*)
      write(iwrt,*)'Desligando CC?:',' SIM'
      write(iwrt,*)
      write(iwrt,*)'Parametros do arquivo TBMD.md:'
      write(iwrt,*)'-------------------------------------'
      write(iwrt,'(a7,2x,i5)')'Frames:',tmax
      write(iwrt,'(a14,2x,f9.6)')'Timestep (ps):',dtime
      write(iwrt,'(a24,2x,f7.3)')'Intervalo de tempo (ps):',dt
      write(iwrt,*)'-------------------------------------'
      write(iwrt,*)
      write(iwrt,*)'Parametros do calculo de autocorrelacao:'
      write(iwrt,*)'----------------------------------------'
      write(iwrt,'(a10,i7)')'  Inicial:',d0w
      write(iwrt,'(a10,i7)')' Amostras:',dtw
      write(iwrt,'(a10,i7)')'Intervalo:',ddw
      write(iwrt,*)'----------------------------------------'
      write(iwrt,*)
c
      write(iwrt,*)'Checando calculo do centro de massa!'
      write(iwrt,*)'Precisao:',cmsprec
      write(iwrt,*)
c
      if(cmserr.eq.0)then
         write(iwrt,*)'CMS ok!'
      else
         write(iwrt,*)'CMS fora de precisao (frames):',cmserr
      end if
c
      write(iwrt,*)
c--------------------------------------------------------------------
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      write(iwrt,*)('%',i=1,79)
      write(*,*)'-> calculando espectro FT-IR'
      call infrared(1,strmm,flexmm)
      write(*,*)'calculo do espectro FT-IR: ok!'
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      write(*,*)'-> calculando momento dipolo'
      call dipole(dm)
      call dipave(dm,int)
      write(*,*)'-> calculo do momento dipolo: ok!'
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      write(*,*)'-> imprimindo probabilidade'
      call prbwrt
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c-imprimindo valores em output
      call output(strmm,flexmm,int)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c-calculando o tempo de CPU
      write(iwrt,*)('%',i=1,79)
      call info(t0)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      end program ftirclass
