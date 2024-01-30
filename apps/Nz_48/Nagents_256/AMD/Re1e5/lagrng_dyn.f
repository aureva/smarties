c....subroutine to compute lagrangian averaging for dynamic model

      subroutine lagrng_dyn(a,b,a_old,b_old,
     +     un,vn,wn,i,j,k,me,t)
      
      include 'dimen.h'

      real*8, dimension(nx,ny,nz2):: a_old,b_old

      real*8 a,b,un,vn,wn,eps,Tn,TLM,TMM,
     +     xi,yi,zi,aint,bint,iw,jw,kw
      
      integer*4 i,j,k,ii,jj,kk,t,iplus,jplus

      TLM=a_old(i,j,k)
      TMM=b_old(i,j,k)

      if(TLM.le.0.0.or.TMM.lt.(1e-10))then
         eps=0.0
      else
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@     
         Tn=1.5*delta*(TLM*TMM)**(-1./8.)
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         eps=(dtl/Tn)/(1.+dtl/Tn)
      endif

      xi=-un*dtl
      yi=-vn*dtl
      zi=-wn*dtl
      
c@@@@@@@@@@@@@@@@@@@@@ New Code @@@@@@@@@@@@@@@@@@@@@@@@@
c     *********** x-comp *************  
      if(xi.lt.0.0)then 
         iw=(dx-abs(xi))*idx
         if(i.eq.1)then
            ii=nx
         else
            ii=i-1
         endif
      else
         iw=(xi)*idx
         ii=i
      endif  
      
c     *********** y-comp *************
      if(yi.lt.0.0)then 
         jw=(dy-abs(yi))*idy
         if(j.eq.1)then
            jj=ny
         else
            jj=j-1
         endif
      else
         jw=(yi)*idy
         jj=j
      endif  
      
c     *********** z-comp *************  
      if(me.eq.0.and.k.le.3)then
         if(zi.lt.0.0)then
            kk=2
            if(k.eq.2)then
               kw=0.0
            else
               if(mom_nodes.eq.0)then
                  kw=(dz/2.-abs(zi))/(dz/2.)
               else
                  kw=(dz-abs(zi))*idz
               endif
            endif
         else
            if(k.eq.2)then
               if(mom_nodes.eq.0)then
                  kw=(zi)/(dz/2.)
               else
                  kw=(zi)*idz
               endif
               kk=2
            else
               kw=(zi)*idz
               kk=3
            endif
         endif
      elseif(me.eq.nprocs-1.and.k.eq.nzb+1)then
         if(zi.lt.0.0)then
            kw=(dz-abs(zi))*idz
            kk=nzb
         else
            kw=0
            kk=nzb+1
         endif
      else
         if(zi.lt.0.0)then 
            kw=(dz-abs(zi))*idz
            kk=k-1
         else
            kw=(zi)*idz
            kk=k
         endif  
      endif
c     ********************************
      
      if(t.le.1.and.INITU.eq.0)then
         kw=0.0
         kk=k
      endif
      
c      if(t.le.1.and.INITS.eq.0)then
c         kw=0.0
c         kk=k
c      endif
          
      if(ii.eq.nx)then
         iplus=-(nx-1)
      else
         iplus=1
      endif
      
      if(jj.eq.ny)then
         jplus=-(ny-1)
      else
         jplus=1
      endif
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      call interp3D(ii,jj,kk,iw,jw,kw,iplus,jplus,a_old,aint)
      call interp3D(ii,jj,kk,iw,jw,kw,iplus,jplus,b_old,bint)
       
      a=eps*a+(1.-eps)*aint
      b=eps*b+(1.-eps)*bint

      return
      
      end
      










