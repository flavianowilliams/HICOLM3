      module sizes
**********************************************************************
*     Dimensionando arrays                                           *
*     Flaviano Williams Fernandes, 20 de agosto de 2014              *
**********************************************************************
      public
c********************************************
c-fmax:   limite de frames                  * 
c-nmax:   limite de átomos                  *  
c-molmax: limite da quantidade de espécies  *
c-nmmax:  limite de moléculas               *
c-iz:     dimensão do espaço vetorial       *
c-atmmax: limite de átomos por molécula     *
c********************************************
      integer fmax,nmax,iz,molmax,nmmax,atmmax
c
      parameter (fmax=1100)
      parameter (nmax=6000)
      parameter (molmax=4)
      parameter (nmmax=2000)
      parameter (atmmax=200)
      parameter (iz=3)
c
      integer ird,irdd,iwrx,iwrt,iwtt,iwrz,iwrh
      parameter (ird=1,irdd=2,iwrx=3,iwtt=4,iwrz=8,iwrh=6,iwrt=7)
c     
      end module sizes
