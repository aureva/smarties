      SUBROUTINE tridag(a,b,c,r,u,n)
      INTEGER*4 n,NMAX
      REAL*8 a(n),b(n),c(n),r(n),u(n)
      PARAMETER (NMAX=500)
      INTEGER*4 j
      REAL*8 bet,gam(NMAX)

      if(b(1).eq.0.)pause 'tridag: rewrite equations'
      bet=b(1)
      u(1)=r(1)/bet
      do 11 j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(bet.eq.0.) then
           print *, 'tridag failed at k=',j  
           print *,'a, b, c, gam, and bet=',a(j),b(j),c(j),gam(j),bet
           pause                        
        end if   
        u(j)=(r(j)-a(j)*u(j-1))/bet
11    continue
      do 12 j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
12    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ]2#"0>Ya%.
