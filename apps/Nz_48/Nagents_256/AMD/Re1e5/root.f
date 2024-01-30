      subroutine root(rtnewt,a1,a2,b1,b2,c1,c2,d1,d2,e1,e2)
      INTEGER JMAX
      REAL*8 rtnewt,xinit,xacc,b,f1,f2,a1,a2,b1,b2,c1,c2,d1,d2,e1,e2
      PARAMETER (JMAX=100)
      PARAMETER (xinit=100.0,xacc=0.0001)
      INTEGER j
      REAL*8 df,dx,f,sign
      REAL*8 AA,BB,CC,DD,EE,FF,AAp,BBp,CCp,DDp
 
      rtnewt=xinit

      AA = a1*d2
      BB = -b1*d2-a2*d1
      CC = -a1*e2+a2*e1
      DD = b1*e2-a2*c1+b2*d1
      EE = a1*c2-b2*e1
      FF = -b1*c2+b2*c1

      AAp = AA*5.
      BBp = BB*4.
      CCp = CC*3.
      DDp = DD*2.
      
      do 11 j=1,JMAX
c.....................................
	b=rtnewt

        f = AA*b**5. + BB*b**4. + CC*b**3. +
     +      DD*b**2. + EE*b + FF
c        f = AA*b*b*b*b*b + BB*b*b*b*b + CC*b*b*b +
c     +      DD*b*b + EE*b + FF

        df = AAp*b**4. + BBp*b**3. + CCp*b**2. +
     +       DDp*b + EE
c        df = AAp*b*b*b*b + BBp*b*b*b + CCp*b*b +
c     +       DDp*b + EE
c.....................................
        dx=f/df

	rtnewt=rtnewt-dx

        if(abs(dx).lt.xacc) return
11    continue

222	rtnewt=0.

      end
