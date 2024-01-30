      Subroutine openfiles ()
      include 'dimen.h'
      character*11 ff
      character*6 pp
      character*7 ss
      
      ff='unformatted'
      pp='append'
      ss='unknown'
      
c...  Open output files  

c...  Forcing files
      Open (unit=61,file='output/Ugeo.bin',status=ss,form=ff,
     +     position=pp)
      Open (unit=62,file='output/Vgeo.bin',status=ss,form=ff,
     +     position=pp)
      Open (unit=63,file='output/Uadv.bin',status=ss,form=ff,
     +     position=pp)
      Open (unit=64,file='output/Vadv.bin',status=ss,form=ff,
     +     position=pp)
      Open (unit=65,file='output/Tadv.bin',status=ss,form=ff,
     +     position=pp)
      Open (unit=66,file='output/Qadv.bin',status=ss,form=ff,
     +position=pp)
    
      Open (unit=11,file='output/ke.out',status=ss,position=pp)
      Open (unit=39,file='output/tke.bin',status=ss,form=ff,position=pp)
      Open (unit=71,file='output/au.bin',status=ss,form=ff,position=pp)
      Open (unit=72,file='output/av.bin',status=ss,form=ff,position=pp)
      Open (unit=73,file='output/aw.bin',status=ss,form=ff,position=pp)
      Open (unit=74,file='output/u2.bin',status=ss,form=ff,position=pp)
      Open (unit=75,file='output/v2.bin',status=ss,form=ff,position=pp)
      Open (unit=76,file='output/w2.bin',status=ss,form=ff,position=pp)
      Open (unit=40,file='output/w3.bin',status=ss,form=ff,position=pp)
      Open (unit=57,file='output/atxx.bin',status=ss,form=ff,
     +     position=pp)
      Open (unit=78,file='output/atxz.bin',status=ss,form=ff,
     +     position=pp)
      Open (unit=79,file='output/atyy.bin',status=ss,form=ff,
     +     position=pp)
      Open (unit=80,file='output/atyz.bin',status=ss,form=ff,
     +     position=pp)
      Open (unit=87,file='output/atzz.bin',status=ss,form=ff,
     +     position=pp)
      Open (unit=90,file='output/atxy.bin',status=ss,form=ff,
     +     position=pp)
      Open (unit=81,file='output/p2.bin',status=ss,form=ff,
     +     position=pp)
      Open (unit=82,file='output/auw.bin',status=ss,form=ff,
     +     position=pp)
      Open (unit=83,file='output/avw.bin',status=ss,form=ff,
     +     position=pp)
      Open (unit=77,file='output/auv.bin',status=ss,form=ff,
     +     position=pp)
      Open (unit=85,file='output/ap.bin',status=ss,form=ff,
     +     position=pp)
      Open (unit=88,file='output/dudz.bin',status=ss,form=ff,
     +     position=pp)
      Open (unit=91,file='output/dudx.bin',status=ss,form=ff,
     +     position=pp)
      Open (unit=95,file='output/dvdz.bin',status=ss,form=ff,
     +     position=pp)
      Open (unit=93,file='output/dwdz.bin',status=ss,form=ff,
     +     position=pp)
      Open (unit=92,file='output/dwdx.bin',status=ss,form=ff,
     +     position=pp)
      Open (unit=299,file='output/spectru.bin',status=ss,form=ff,
     +     position=pp)
      Open (unit=298,file='output/spectrv.bin',status=ss,form=ff,
     +     position=pp)
      Open (unit=297,file='output/spectrw.bin',status=ss,form=ff,
     +     position=pp)
      Open(unit=398,file='output/Cs2_ALL.bin',status=ss,form=ff,
     +     position=pp)
      Open(unit=399,file='output/Cs_ALL.bin',status=ss,form=ff,
     +     position=pp)
      Open (unit=978,file='output/beta1.bin',status=ss,form=ff,
     +     position=pp)
      Open (unit=888,file='output/ESGS.bin',status=ss,form=ff,
     +     position=pp)

	 
      Open (unit=889,file='output/FSu1.bin',status=ss,form=ff,
     +     position=pp)
      Open (unit=890,file='output/FSu2.bin',status=ss,form=ff,
     +     position=pp)
      Open (unit=891,file='output/FSu3.bin',status=ss,form=ff,
     +     position=pp)
      Open (unit=892,file='output/FSu4.bin',status=ss,form=ff,
     +     position=pp)
      Open (unit=893,file='output/FSv1.bin',status=ss,form=ff,
     +     position=pp)
      Open (unit=894,file='output/FSv2.bin',status=ss,form=ff,
     +     position=pp)
      Open (unit=895,file='output/FSv3.bin',status=ss,form=ff,
     +     position=pp)
      Open (unit=896,file='output/FSv4.bin',status=ss,form=ff,
     +     position=pp)
      Open (unit=897,file='output/t_ground.out',status=ss,position=pp)
 

      open(unit=132,file='output/ustar.bin',status=ss,form=ff,
     +     position=pp)

      Open (unit=134,file='output/meanQ.out',status=ss,position=pp)	 
	 
      if(nframe.gt.0)then
         Open(unit=994,file='output/u_frame.bin',status=ss,
     +        form=ff,position=pp)
         Open(unit=995,file='output/v_frame.bin',status=ss,
     +        form=ff,position=pp)
         Open(unit=998,file='output/w_frame.bin',status=ss,
     +        form=ff,position=pp)
         if(S_Flag.eq.1)then
            Open(unit=997,file='output/t_frame.bin',status=ss,
     +           form=ff,position=pp)
         endif
         if(Q_Flag.eq.1)then
            Open(unit=901,file='output/q_frame.bin',status=ss,
     +           form=ff,position=pp)
         endif
      endif

      If (S_Flag.eq.1) then
         Open (unit=89,file='output/dtdz.bin',status=ss,form=ff,
     +        position=pp)
         Open (unit=114,file='output/dtdx.bin',status=ss,form=ff,
     +        position=pp)
         Open (unit=115,file='output/dtdy.bin',status=ss,form=ff,
     +        position=pp)
         Open (unit=100,file='output/at.bin',status=ss,form=ff,
     +        position=pp)
         Open (unit=102,file='output/t2.bin',status=ss,form=ff,
     +        position=pp)
         Open (unit=41,file='output/t3.bin',status=ss,form=ff,
     +        position=pp)
         Open (unit=112,file='output/flux_t1.bin',status=ss,form=ff,
     +        position=pp)
         Open (unit=116,file='output/flux_t2.bin',status=ss,form=ff,
     +        position=pp)
         Open (unit=104,file='output/flux_t3.bin',status=ss,form=ff,
     +        position=pp)
         Open(unit=106,file='output/awt.bin',status=ss,form=ff,
     +        position=pp)
         Open(unit=113,file='output/aut.bin',status=ss,form=ff,
     +        position=pp)
         Open(unit=117,file='output/avt.bin',status=ss,form=ff,
     +        position=pp)
         open (unit=300,file='output/spectrt.bin',status=ss,form=ff,
     +        position=pp)
         open(unit=1526,file='output/prandtl.bin',status=ss,form=ff,
     +        position=pp)
         open(unit=1527,file='output/cs2pr.bin',status=ss,form=ff,
     +        position=pp)
         Open (unit=979,file='output/beta2.bin',status=ss,form=ff,
     +        position=pp)
         Open (unit=999,file='output/ET.bin',status=ss,form=ff,
     +        position=pp)
         open(unit=133,file='output/qz_surf.bin',status=ss,form=ff,
     +        position=pp)
         open(unit=135,file='output/ts_surf.bin',status=ss,form=ff,
     +        position=pp)	 
      end if

      If (Q_Flag.eq.1) then
         Open (unit=903,file='output/dqdz.bin',status=ss,form=ff,
     +        position=pp)
         Open (unit=904,file='output/dqdx.bin',status=ss,form=ff,
     +        position=pp)
         Open (unit=905,file='output/dqdy.bin',status=ss,form=ff,
     +        position=pp)
         Open (unit=906,file='output/aq.bin',status=ss,form=ff,
     +        position=pp)
         Open (unit=907,file='output/q2.bin',status=ss,form=ff,
     +        position=pp)
         Open (unit=908,file='output/q3.bin',status=ss,form=ff,
     +        position=pp)
         Open (unit=909,file='output/flux_q1.bin',status=ss,form=ff,
     +        position=pp)
         Open (unit=910,file='output/flux_q2.bin',status=ss,form=ff,
     +        position=pp)
         Open (unit=911,file='output/flux_q3.bin',status=ss,form=ff,
     +        position=pp)
         Open(unit=912,file='output/awq.bin',status=ss,form=ff,
     +        position=pp)
         Open(unit=913,file='output/auq.bin',status=ss,form=ff,
     +        position=pp)
         Open(unit=914,file='output/avq.bin',status=ss,form=ff,
     +        position=pp)
         open(unit=915,file='output/spectrq.bin',status=ss,form=ff,
     +        position=pp)
         open(unit=916,file='output/Schmidt.bin',status=ss,form=ff,
     +        position=pp)
         open(unit=917,file='output/cs2sc.bin',status=ss,form=ff,
     +        position=pp)
         Open (unit=918,file='output/beta3.bin',status=ss,form=ff,
     +        position=pp)
         Open (unit=919,file='output/EQ.bin',status=ss,form=ff,
     +        position=pp)
         open(unit=920,file='output/qm_surf.bin',status=ss,form=ff,
     +        position=pp)
      end if

!      if(initu.eq.0.and.inits.eq.0)then
         
         Close(unit=61,status='Delete')
         Close(unit=62,status='Delete')
         Close(unit=63,status='Delete')
         Close(unit=64,status='Delete')
         Close(unit=65,status='Delete')
         Close(unit=66,status='Delete')

         Close (unit=11,status='Delete')
         Close(unit=39,status='Delete')
         Close(unit=71,status='Delete')
         Close(unit=72,status='Delete')
         Close (unit=73,status='Delete')
         Close (unit=74,status='Delete')
         Close (unit=77,status='Delete')
         Close (unit=75,status='Delete')
         Close (unit=76,status='Delete')
         Close (unit=40,status='Delete')
         Close (unit=57,status='Delete')
         Close (unit=78,status='Delete')
         Close (unit=79,status='Delete')
         Close (unit=80,status='Delete')
         Close (unit=87,status='Delete')
         Close (unit=81,status='Delete')
         Close (unit=82,status='Delete')
         Close (unit=83,status='Delete')
         Close (unit=85,status='Delete')
         Close (unit=88,status='Delete')
         Close (unit=90,status='Delete')
         Close (unit=91,status='Delete')
         Close (unit=95,status='Delete')
         Close (unit=92,status='Delete')
         Close (unit=93,status='Delete')
         Close (unit=297,status='Delete')
         Close (unit=298,status='Delete')
         Close (unit=299,status='Delete')
         Close (unit=398,status='Delete')
         Close (unit=399,status='Delete')
         Close (unit=978,status='Delete')
         Close (unit=888,status='Delete')
         Close (unit=132,status='Delete')

         If (S_Flag.eq.1) then
            Close (unit=89,status='Delete')
            Close (unit=100,status='Delete')
            Close (unit=102,status='Delete')
            Close (unit=41,status='Delete')
            Close (unit=104,status='Delete')
            Close (unit=106,status='Delete')
            Close (unit=112,status='Delete')
            Close (unit=113,status='Delete')
            Close (unit=114,status='Delete')
            Close (unit=115,status='Delete')
            Close (unit=116,status='Delete')
            Close (unit=117,status='Delete')
            Close (unit=1526,status='Delete')
            Close (unit=1527,status='Delete')
            Close (unit=300,status='Delete')
            Close (unit=979,status='Delete')
            Close (unit=999,status='Delete')
            Close (unit=133,status='Delete')
         end if

         If (Q_Flag.eq.1) then
            Close (unit=903,status='Delete')
            Close (unit=904,status='Delete')
            Close (unit=905,status='Delete')
            Close (unit=906,status='Delete')
            Close (unit=907,status='Delete')
            Close (unit=908,status='Delete')
            Close (unit=909,status='Delete')
            Close (unit=910,status='Delete')
            Close (unit=911,status='Delete')
            Close (unit=912,status='Delete')
            Close (unit=913,status='Delete')
            Close (unit=914,status='Delete')
            Close (unit=915,status='Delete')
            Close (unit=916,status='Delete')
            Close (unit=917,status='Delete')
            Close (unit=918,status='Delete')
            Close (unit=919,status='Delete')
            Close (unit=920,status='Delete')
         end if
         
         if(nframe.gt.0)then
            Close(unit=994,status='Delete')
            Close(unit=995,status='Delete')
            Close(unit=998,status='Delete')
            if(S_Flag.eq.1)then
               Close(unit=997,status='Delete')
            endif
            if(Q_Flag.eq.1)then
               Close(unit=901,status='Delete')
            endif
            
         endif

         Open (unit=61,file='output/Ugeo.bin',status=ss,
     +        form=ff,position=pp)
         Open (unit=62,file='output/Vgeo.bin',status=ss,
     +        form=ff,position=pp)
         Open (unit=63,file='output/Uadv.bin',status=ss,
     +        form=ff,position=pp)
         Open (unit=64,file='output/Vadv.bin',status=ss,
     +        form=ff,position=pp)
         Open (unit=65,file='output/Tadv.bin',status=ss,
     +        form=ff,position=pp)
         Open (unit=66,file='output/Qadv.bin',status=ss,
     +        form=ff,position=pp)

         Open (unit=11,file='output/ke.out',status=ss,position=pp)
         Open (unit=39,file='output/tke.bin',status=ss,form=ff,
     +        position=pp)
         Open (unit=71,file='output/au.bin',status=ss,form=ff,
     +        position=pp)
         Open (unit=72,file='output/av.bin',status=ss,form=ff,
     +        position=pp)
         Open (unit=73,file='output/aw.bin',status=ss,form=ff,
     +        position=pp)
         Open (unit=74,file='output/u2.bin',status=ss,form=ff,
     +        position=pp)
         Open (unit=75,file='output/v2.bin',status=ss,form=ff,
     +        position=pp)
         Open (unit=76,file='output/w2.bin',status=ss,form=ff,
     +        position=pp)
         Open (unit=40,file='output/w3.bin',status=ss,form=ff,
     +        position=pp)
         Open (unit=57,file='output/atxx.bin',status=ss,form=ff,
     +        position=pp)
         Open (unit=78,file='output/atxz.bin',status=ss,form=ff,
     +        position=pp)
         Open (unit=79,file='output/atyy.bin',status=ss,form=ff,
     +        position=pp)
         Open (unit=80,file='output/atyz.bin',status=ss,form=ff,
     +        position=pp)
         Open (unit=87,file='output/atzz.bin',status=ss,form=ff,
     +        position=pp)
         Open (unit=90,file='output/atxy.bin',status=ss,form=ff,
     +        position=pp)
         Open (unit=81,file='output/p2.bin',status=ss,form=ff,
     +        position=pp)
         Open (unit=82,file='output/auw.bin',status=ss,form=ff,
     +        position=pp)
         Open (unit=83,file='output/avw.bin',status=ss,form=ff,
     +        position=pp)
         Open (unit=77,file='output/auv.bin',status=ss,form=ff,
     +        position=pp)
         Open (unit=85,file='output/ap.bin',status=ss,form=ff,
     +        position=pp)
         Open (unit=88,file='output/dudz.bin',status=ss,form=ff,
     +        position=pp)
         Open (unit=91,file='output/dudx.bin',status=ss,form=ff,
     +        position=pp)
         Open (unit=95,file='output/dvdz.bin',status=ss,form=ff,
     +        position=pp)
         Open (unit=93,file='output/dwdz.bin',status=ss,form=ff,
     +        position=pp)
         Open (unit=92,file='output/dwdx.bin',status=ss,form=ff,
     +        position=pp)
         Open (unit=299,file='output/spectru.bin',status=ss,form=ff,
     +        position=pp)
         Open (unit=298,file='output/spectrv.bin',status=ss,form=ff,
     +        position=pp)
         Open (unit=297,file='output/spectrw.bin',status=ss,form=ff,
     +        position=pp)
         Open(unit=398,file='output/Cs2_ALL.bin',status=ss,form=ff,
     +        position=pp)
         Open(unit=399,file='output/Cs_ALL.bin',status=ss,form=ff,
     +        position=pp)
         Open (unit=978,file='output/beta1.bin',status=ss,form=ff,
     +        position=pp)
        Open (unit=888,file='output/ESGS.bin',status=ss,form=ff,
     +        position=pp)

         open(unit=132,file='output/ustar.bin',status=ss,form=ff,
     +        position=pp)

        open(1525,file='output/surface_flux.out',status=ss,position=pp)
        open(1524,file='output/power.out',status=ss,position=pp)
        open(1523,file='output/omega.out',status=ss,position=pp)
		
        open(1528,file='output/power_vawt.out',status=ss,position=pp)
        open(1529,file='output/omega_vawt.out',status=ss,position=pp)		
		
         if(nframe.gt.0)then
            Open(unit=994,file='output/u_frame.bin',status=ss,
     +           form=ff,position=pp)
            Open(unit=995,file='output/v_frame.bin',status=ss,
     +           form=ff,position=pp)
            Open(unit=998,file='output/w_frame.bin',status=ss,
     +           form=ff,position=pp)
            if(S_Flag.eq.1)then
               Open(unit=997,file='output/t_frame.bin',status=ss,
     +              form=ff,position=pp)
            endif
            if(Q_Flag.eq.1)then
               Open(unit=901,file='output/q_frame.bin',status=ss,
     +              form=ff,position=pp)
            endif
         endif

         If (S_Flag.eq.1) then
            Open (unit=89,file='output/dtdz.bin',status=ss,form=ff,
     +           position=pp)
            Open (unit=114,file='output/dtdx.bin',status=ss,form=ff,
     +           position=pp)
            Open (unit=115,file='output/dtdy.bin',status=ss,form=ff,
     +           position=pp)
            Open (unit=100,file='output/at.bin',status=ss,form=ff,
     +           position=pp)
            Open (unit=102,file='output/t2.bin',status=ss,form=ff,
     +           position=pp)
            Open (unit=41,file='output/t3.bin',status=ss,form=ff,
     +        position=pp)
            Open (unit=112,file='output/flux_t1.bin',status=ss,form=ff,
     +           position=pp)
            Open (unit=116,file='output/flux_t2.bin',status=ss,form=ff,
     +           position=pp)
            Open (unit=104,file='output/flux_t3.bin',status=ss,form=ff,
     +           position=pp)
            Open(unit=106,file='output/awt.bin',status=ss,form=ff,
     +           position=pp)
            Open(unit=113,file='output/aut.bin',status=ss,form=ff,
     +           position=pp)
            Open(unit=117,file='output/avt.bin',status=ss,form=ff,
     +           position=pp)
            open (unit=300,file='output/spectrt.bin',status=ss,form=ff,
     +           position=pp)
            open(unit=1526,file='output/prandtl.bin',status=ss,form=ff,
     +           position=pp)
            open(unit=1527,file='output/cs2pr.bin',status=ss,form=ff,
     +           position=pp)
            Open (unit=979,file='output/beta2.bin',status=ss,form=ff,
     +           position=pp)
            Open (unit=999,file='output/ET.bin',status=ss,form=ff,
     +           position=pp)
            open(unit=133,file='output/qz_surf.bin',status=ss,form=ff,
     +           position=pp)
	 
c            open(1525,file='output/obukov.out',status=ss,position=pp)
            open(1999,file='output/t_star.out',status=ss,position=pp)

         end if

         If (Q_Flag.eq.1) then
            Open (unit=903,file='output/dqdz.bin',status=ss,form=ff,
     +           position=pp)
            Open (unit=904,file='output/dqdx.bin',status=ss,form=ff,
     +           position=pp)
            Open (unit=905,file='output/dqdy.bin',status=ss,form=ff,
     +           position=pp)
            Open (unit=906,file='output/aq.bin',status=ss,form=ff,
     +           position=pp)
            Open (unit=907,file='output/q2.bin',status=ss,form=ff,
     +           position=pp)
            Open (unit=908,file='output/q3.bin',status=ss,form=ff,
     +           position=pp)
            Open (unit=909,file='output/flux_q1.bin',status=ss,form=ff,
     +           position=pp)
            Open (unit=910,file='output/flux_q2.bin',status=ss,form=ff,
     +           position=pp)
            Open (unit=911,file='output/flux_q3.bin',status=ss,form=ff,
     +           position=pp)
            Open(unit=912,file='output/awq.bin',status=ss,form=ff,
     +           position=pp)
            Open(unit=913,file='output/auq.bin',status=ss,form=ff,
     +           position=pp)
            Open(unit=914,file='output/avq.bin',status=ss,form=ff,
     +           position=pp)
            open(unit=915,file='output/spectrq.bin',status=ss,form=ff,
     +           position=pp)
            open(unit=916,file='output/Schmidt.bin',status=ss,form=ff,
     +           position=pp)
            open(unit=917,file='output/cs2sc.bin',status=ss,form=ff,
     +           position=pp)
            Open (unit=918,file='output/beta3.bin',status=ss,form=ff,
     +           position=pp)
            Open (unit=919,file='output/EQ.bin',status=ss,form=ff,
     +           position=pp)
            open(unit=920,file='output/qm_surf.bin',status=ss,form=ff,
     +           position=pp)
         end if
         
!      endif
      
      end

      Subroutine closefiles ()
      include 'dimen.h'

         Close(unit=61)
         Close(unit=62)
         Close(unit=63)
         Close(unit=64)
         Close(unit=65)
         Close(unit=66)

         Close(unit=11)
         Close(unit=39)
         Close(unit=71)
         Close(unit=72)
         Close(unit=73)
         Close(unit=74)
         Close(unit=77)
         Close(unit=75)
         Close(unit=76)
         Close(unit=40)
         Close(unit=57)
         Close(unit=78)
         Close(unit=79)
         Close(unit=80)
         Close(unit=87)
         Close(unit=81)
         Close(unit=82)
         Close(unit=83)
         Close(unit=85)
         Close(unit=88)
         Close(unit=90)
         Close(unit=91)
         Close(unit=95)
         Close(unit=92)
         Close(unit=93)
         Close(unit=297)
         Close(unit=298)
         Close(unit=299)
         Close(unit=398)
         Close(unit=399)
         Close(unit=978)
         Close(unit=888)
         Close(unit=132)

         If (S_Flag.eq.1) then
            Close (unit=89)
            Close (unit=100)
            Close (unit=102)
            Close (unit=41)
            Close (unit=104)
            Close (unit=106)
            Close (unit=112)
            Close (unit=113)
            Close (unit=114)
            Close (unit=115)
            Close (unit=116)
            Close (unit=117)
            Close (unit=1526)
            Close (unit=1527)
            Close (unit=300)
            Close (unit=979)
            Close (unit=999)
            Close (unit=133)
         end if

         If (Q_Flag.eq.1) then
            Close (unit=903)
            Close (unit=904)
            Close (unit=905)
            Close (unit=906)
            Close (unit=907)
            Close (unit=908)
            Close (unit=909)
            Close (unit=910)
            Close (unit=911)
            Close (unit=912)
            Close (unit=913)
            Close (unit=914)
            Close (unit=915)
            Close (unit=916)
            Close (unit=917)
            Close (unit=918)
            Close (unit=919)
            Close (unit=920)
         end if

         if(nframe.gt.0)then
            Close(unit=994)
            Close(unit=995)
            Close(unit=998)
            if(S_Flag.eq.1)then
               Close(unit=997)
            endif
            if(Q_Flag.eq.1)then
               Close(unit=901)
            endif
         endif

       end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine openfiles1

      include 'dimen.h'
      character*11 ff,pp,ss,aa,fm1,fm2,fm3
      character*33 path
      integer irecl,fnum,i,fid
      integer*4 fsize
      real*8 recx
      inquire(iolength=irecl) recx

      Open (500,file='flist.txt',form='formatted')

      ff  ='unformatted'
      pp  ='rewind'
      ss  ='replace'
      aa  ='sequential'
      fm1 ='output/'
      fm3 ='.bin'

!      if(initu.eq.1.and.inits.eq.1)then
!      pp  ='append'
!      ss  ='unknown'
!      endif

      recx = (Nxe-Nxs+1)*(Nye-Nys+1)*(Nze-Nzs+1)
      recx = recx*Nsteps/p_count
      fsize= nint(recx)


c...  Open output files

      read(500,*)fnum

      do i=1,fnum
         read(500,*)fid,fm2
         path=trim(fm1)//trim(fm2)//trim(fm3)

         Open (fid,file=path,access='sequential',
     +   status=ss,form=ff,position=pp)
      end do
      
      close(500)

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccc	  
      subroutine post(fid,ac,me,nall)
      implicit none
      include 'dimen.h'

      integer*4 fid,i,k
      real*8, dimension(nx,ny,nz2):: ac
      real*8, dimension(nx,ny,nzb)::  at1	  
         
      IF (me>0) then
         call MPI_SEND(ac(1,1,2),nzb*nx*ny,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )
      ELSE
         do k=1,nzb
            at1(:,:,k)=ac(:,:,k+1)
         enddo
         write(fid) at1		 
         do i=1,nprocs-1
               call MPI_RECV(at1(1,1,0*nzb+1),nzb*nx*ny,
     +              MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )
         write(fid) at1
         at1=0.0d0		 
         end do
         call flush(fid)
      END IF


      ac = 0.0d0
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      subroutine post_dimen(me)
      implicit none
      include 'dimen.h'
      integer*4 i
      if(me.eq.0)then
      Open(unit=455,file='output/dimen.dat',form='formatted')
      write(455,'(I4)')Nx
      write(455,'(I4)')Ny
      write(455,'(I4)')Nz
      write(455,'(E13.5)')l_r
      write(455,'(I4)')model
      write(455,'(I4)')averaging
c      write(455,'(I4)')Nxbs
c      write(455,'(I4)')Nxbe
      write(455,'(I8)')Nsteps
      write(455,'(I8)')p_count

      write(455,'(E13.5)')u_star
c      write(455,'(E13.5)')zo1
c      write(455,'(E13.5)')zoH1
      write(455,'(E13.5)')z_i
      write(455,'(E13.5)')l_z
      write(455,'(E13.5)')z_d
      write(455,'(E13.5)')dt
      write(455,'(E13.5)')fgr
      write(455,'(E13.5)')tfr

      write(455,'(E13.5)')T_Scale
      write(455,'(E13.5)')theta_0
      write(455,'(E13.5)')Ugal
      write(455,'(E13.5)')Ugeo_inf
      write(455,'(E13.5)')Vgeo_inf

      write(455,'(I4)')num_turbine

      do i = 1,num_turbine
      write(455,'(I4)')wtx(i)
      write(455,'(I4)')wty(i)
      write(455,'(E13.5)')Zhub(i)
      end do
      close(455)
      end if

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      subroutine post_instant3(fid,ac,me,nall)
      implicit none
      include 'dimen.h'

      integer*4 fid,i,k
      real*8, dimension(nx,ny,nz2):: ac
      real*8, dimension(nx,ny,nzb)::  at1	  
         
      IF (me>0) then
         call MPI_SEND(ac(1,1,2),nzb*nx*ny,
     +        MPI_DOUBLE_PRECISION,0,me,nall, ierr )
      ELSE
         do k=1,nzb
            at1(:,:,k)=ac(:,:,k+1)
         enddo
         write(fid) at1		 
         do i=1,nprocs-1
               call MPI_RECV(at1(1,1,0*nzb+1),nzb*nx*ny,
     +              MPI_DOUBLE_PRECISION,i,i,nall,status2,ierr )
         write(fid) at1
         at1=0.0d0		 
         end do
         call flush(fid)
      END IF

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      






