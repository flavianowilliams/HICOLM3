      subroutine prbwrt
c
      implicit none
c
      integer maxn,mimax,i,strmax,nmax,flxmax
c
      parameter (maxn=500)
c
      integer ird,irdd,iwrx,iwrt,iwtt,iwrz,iwrh
      real(kind=4) dmi,mi,strmx,dstr,strv,dstr0,dflx0,dflx,flxmx,flxv
      real(kind=4) nn(maxn),str(maxn),flx(maxn)
c
      common/units/ ird,irdd,iwrh,iwrx,iwtt,iwrz,iwrt
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
