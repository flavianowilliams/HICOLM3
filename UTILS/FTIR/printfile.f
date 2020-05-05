      module printfile
**********************************************************************
*     Imprimindo probabilidades                                      *
*     Flaviano Williams Fernandes, 20 de agosto de 2014              *
**********************************************************************
      use sizes
      use system
c
      contains
c
      subroutine prbwrt
c
      implicit none
c
      integer maxn,mimax,i,strmax,nmax,flxmax
c
      parameter (maxn=500)
c
      real(kind=4) dmi,mi,strmx,dstr,strv,dstr0,dflx0,dflx,flxmx,flxv
      real(kind=4) nn(maxn),str(maxn),flx(maxn)
c
      common/dipavecalc/ mimax,dmi,nn
      common/tcfstrdata/ strmax,strmx,dstr0,str
      common/tcfflexdata/ flxmax,flxmx,dflx0,flx
c-----------------------------------------------------
      nmax=min(maxn,min(mimax,min(strmax,flxmax)))
c
      dstr=(strmx-dstr0)/nmax
      dflx=(flxmx-dflx0)/nmax
c
      mi=0.
      strv=0.+dstr0
      flxv=0.+dflx0
      do i=1,nmax
         write(iwrx,'(3(a4,4x,f7.4,f12.4))')
     1        '1',mi,nn(i),'2',strv,str(i),'3',flxv,flx(i)
         mi=mi+dmi
         strv=strv+dstr
         flxv=flxv+dflx
      end do
c-----------------------------------------------------
      return
c
      end subroutine prbwrt
c
      subroutine output(strmm,flexmm,int)

      implicit none
c
      integer i
      real(kind=4) strmm,flexmm,int
c
      write(iwrt,*)
      write(iwrt,*)
      write(iwrt,*)'Valores de saida'
      write(iwrt,*)'Molecula',opm
      write(iwrt,*)('-',i=1,41)
      write(iwrt,'(4x,a8,7x,a7,7x,a5)')'Grandeza','Unidade','Valor'
      write(iwrt,*)('-',i=1,41)
      write(iwrt,'(a13,4x,a8,7x,f5.2)')'Dipolo','Debye',int
      write(iwrt,'(a13,4x,a9,6x,f5.2)')'Estiramento','Angstrom',strmm
      write(iwrt,'(a13,4x,a8,6x,f6.2)')'Angulo flexao','Graus',flexmm
      write(iwrt,*)('-',i=1,41)
c
      return
c
      end subroutine output
c
      end module printfile
