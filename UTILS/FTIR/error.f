      module error
**********************************************************************
*     Ceck por erros nos parametros iniciais                         *
*     Flaviano Williams Fernandes, 20 de agosto de 2014              *
**********************************************************************
      use sizes

      contains

      subroutine erro(er,mo,rotn)
      
      implicit none
      
      integer er,mo,rotn
      
      if(mo.eq.1)goto 1
      if(mo.eq.2)goto 2
      
 1    select case(rotn)
      case(1)
         write(iwrt,*)'molec:'
      case(2)
         write(iwrt,*)'input:'
      case(3)
         write(iwrt,*)'dipave:'
      case(4)
         write(iwrt,*)'tcfstr'
      end select

      select case(er)
      case(1)
      write(iwrt,'(a6,a19)')'Aviso:','Falha em desligar CC'
      case(2)
      write(iwrt,'(a6,a27)')'Aviso:','Erro na contagem dos atomos'
      end select

      goto 3
      
 2    select case(rotn)
      case(1)
         write(*,*)'molec:'
      case(2)
         write(*,*)'input:'
      case(3)
         write(*,*)'dipave:'
      case(4)
         write(*,*)'tcfstr'
      end select

      select case(er)
      case(1)
         write(*,'(a5,a27)')'Erro:','Falha em desligar CC'
      case(2)
         write(*,'(a5,a16)')'Erro:','Divisao por zero'
      case(3)
         write(*,'(a5,a18)')'Erro:','Atomos n/ coincide'
      case(4)
         write(*,'(a5,a34)')'Erro:','Parametro nao coincide com HISTORY'
      case(5)
         write(*,'(a5,a27)')'Erro:','Erro na contagem dos atomos'
      end select
      
      stop
      
 3    continue
      
      return
      
      end subroutine erro

      end module error
