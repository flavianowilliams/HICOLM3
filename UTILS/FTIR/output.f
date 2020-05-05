      subroutine output(vacf,strmm,flexmm,int)

      implicit none

      integer ird,irdd,iwrx,iwrt,iwtt,iwrz,iwrh
      integer opm,ddw,bnd,flx,i
      integer oph(2,10,3),natk(2,10,3),opz(2)
      real(kind=4) dni,strcut,vacf,strmm,flexmm,int

      common/acfdata/ opm,opz,oph,ddw,bnd,flx,dni,strcut,natk
      common/units/ ird,irdd,iwrh,iwrx,iwtt,iwrz,iwrt

      write(iwrt,*)
      write(iwrt,*)
      write(iwrt,*)'Valores de saida'
      write(iwrt,*)'Molecula',opm
      write(iwrt,*)('-',i=1,41)
      write(iwrt,'(4x,a8,7x,a7,7x,a5)')'Grandeza','Unidade','Valor'
      write(iwrt,*)('-',i=1,41)
      write(iwrt,'(a13,4x,a11,4x,f5.2)')'Difusao','x10-5 cm2/s',vacf
      write(iwrt,'(a13,4x,a8,7x,f5.2)')'Dipolo','Debye',int
      write(iwrt,'(a13,4x,a9,6x,f5.2)')'Estiramento','Angstrom',strmm
      write(iwrt,'(a13,4x,a8,6x,f6.2)')'Angulo flexao','Graus',flexmm
      write(iwrt,*)('-',i=1,41)

      return

      end subroutine output
