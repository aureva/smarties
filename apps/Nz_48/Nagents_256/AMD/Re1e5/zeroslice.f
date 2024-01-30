      subroutine  zeroslice(au,av,aw,ap,u2,v2,w2,w3,atxx,atxz,atyy,
     +     atyz,atzz,atxy,p2,auw,avw,auv,adudz,adudx,advdz,adwdz,adwdx,
     +     e,aCs2,aCs,abeta1,atxz_s,aESGS,
     +     aFSu1,aFSu2,aFSu3,aFSu4,aFSv1,aFSv2,aFSv3,aFSv4)
      
      implicit none
      include 'dimen.h'
      real*8,dimension(anx,nz2) :: au,av,aw,u2,v2,w2,w3,atxx,atxz,atyy,
     +     atyz,atzz,atxy,p2,auw,avw,auv,ap,adudz,adudx,e,
     +     aCs2,aCs,abeta1,adwdz,adwdx,advdz,aESGS,
     +     aFSu1,aFSu2,aFSu3,aFSu4,aFSv1,aFSv2,aFSv3,aFSv4	 
      real*8,dimension(nx,ny) :: atxz_s

      au = 0.0
      av = 0.0
      aw = 0.0
      ap = 0.0
      u2 = 0.0
      v2 = 0.0
      w2 = 0.0
      w3 = 0.0
      atxx = 0.0
      atxz = 0.0
      atyy = 0.0
      atyz = 0.0
      atzz = 0.0
      atxy = 0.0
      p2 = 0.0
      auw = 0.0
      avw = 0.0
      auv = 0.0
      adudz = 0.0
      adudx = 0.0
      advdz = 0.0
      adwdz = 0.0
      adwdx = 0.0
      e = 0.0
      aCs2 = 0.0
      aCs = 0.0
      abeta1 = 0.0
      atxz_s   = 0.0
      aESGS = 0.0
	  
      aFSu1=0.0
      aFSu2=0.0	  
      aFSu3=0.0	
      aFSu4=0.0
      aFSv1=0.0
      aFSv2=0.0	  
      aFSv3=0.0	
      aFSv4=0.0	    

      return
      end
