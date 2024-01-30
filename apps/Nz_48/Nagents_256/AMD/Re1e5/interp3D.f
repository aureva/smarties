c  recieves an array and a set of three global points and returns the value
c  interpolated using linear interpolation for the new point
c  note that ot=0 for momentum and 1 for scalars...

	subroutine interp3D(ii,jj,kk,iw,jw,kw,iplus,jplus,U,ui)
	
	implicit none
	include 'dimen.h'
	real*8 ui,U(Nx,Ny,Nz2),iw,jw,kw,
     +	     ui_hgh,ui_low   
	integer*4 ii,jj,kk,iplus,jplus
	
c       use periodic boundary condition in x and y
	
	if(kw.eq.0.0)then
	    ui =    U(ii,jj,kk)*(1.-iw)*(1.-jw)+
     +	            U(ii+iplus,jj,kk)*iw*(1.-jw)+
     +              U(ii+iplus,jj+jplus,kk)*iw*jw+
     +	            U(ii,jj+jplus,kk)*(1.-iw)*jw
	else
	   ui_low = U(ii,jj,kk)*(1.-iw)*(1.-jw)+
     +	            U(ii+iplus,jj,kk)*iw*(1.-jw)+
     +              U(ii+iplus,jj+jplus,kk)*iw*jw+
     +	            U(ii,jj+jplus,kk)*(1.-iw)*jw
     +          
	   ui_hgh = U(ii,jj,kk+1)*(1.-iw)*(1.-jw)+
     +	            U(ii+iplus,jj,kk+1)*iw*(1.-jw)+
     +              U(ii+iplus,jj+jplus,kk+1)*iw*jw+
     +	            U(ii,jj+jplus,kk+1)*(1.-iw)*jw

	   ui = ui_low*(1.-kw)+ui_hgh*kw
	 endif

        return
	
        end



