!--------- 2D RDF Subroutine-----------------------------------------------
! Reference: 
! Mathias & Tavan "Angular resolution & range of dipole-dipole correlations in water" 
! J. Chem. Phys. (2004) 120,9 4393 
!
! Dan Elton 4/13, updated and incorporated into RDF code 7/14
!
! This program is supposed to reproduce plots as shown in fig. 2 of Mathias & Tavan's paper, 
! given an input file with OHH coordinates in the order OHHOHHOHHOHH... 
!   
! The program outputs 2D arrays which are histograms
!
! The output format is chosen for easy plotting in OCTAVE or MATHEMATICA
!-----------------------------------------------------------------------------
module RDF2D
 use global_vars
Implicit None

!2DRDF variables
real,dimension(:,:),allocatable,save ::  g_hist, dip_hist, D_hist, num_hist
real,save      :: div_size
integer,save      :: numxbins, numzbins 
 

contains 
!-----------------------------------------------------------------------
!--------------------2D RDFS ------------------------------------------
!-----------------------------------------------------------------------
 subroutine m_2DRDF(Oxy, Hydro, TTM3Fdips)
    Implicit None
    real,dimension(numzbins,numzbins) :: g_hist, dip_hist, D_hist 
    Real,dimension(3,Nmol),intent(in)   :: Oxy, TTM3Fdips
    real,dimension(3,Nmol*2),intent(in) :: Hydro
    real :: rho, vol, distance, tmp, xd, zd
    integer :: ia, ja, ix, i, countx, countz 
    real, dimension(3) :: v1, v2, v3, d1, d2, r, x, z, summ, rhat,  zhat
    real :: dip_fac, D_fac, cos_fac

    ia = 0
    ja = 0
    rho = Nmol / (  length(1)*length(2)*length(3) )

    do ia=1,Nmol-1 ! do for every molecule
        ! get the dipole vector for this molecule
        if (TTM3F) then
	   	d1 = TTM3Fdips(1:3,ia)
	else if (TIP4P2005) then
		v1=hydro(:,2*ia-0)-Oxy(:,ia)
   		if (PBC) v1=v1 - length*anint(v1/length)!PBC
		v2=hydro(:,2*ia-1)-Oxy(:,ia)
       		if (PBC) v2=v2 - length*anint(v2/length)!PBC
		summ = v1+v2!find vector to M-site
   	 	v3=(summ/sqrt(dot_product(summ,summ)))*.1546
      		d1=(v1*qH+v2*qH+v3*qO)*ang2m*e2coul/3.33564e-30!conv. to Debye
	  	!write(*,*) sqrt(dot_product(d1,d1))
	else 
       		v1=hydro(:,2*ia-0)-Oxy(:,ia)
		if (PBC) v1=v1 - length*anint(v1/length)!PBC
      		v2=hydro(:,2*ia-1)-Oxy(:,ia)
       		if (PBC) v2=v2 - length*anint(v2/length)!PBC
       		d1=v1*qH + v2*qH
	endif

       do ja=ia+1,Nmol ! do for every other molecule
	  if (TTM3F) then
	   	d2 = TTM3Fdips(1:3,ja)
          else if (TIP4P2005) then
		v1=hydro(:,2*ja-0)-Oxy(:,ja)
	  	if (PBC) v1=v1 - length*anint(v1/length)!PBC
          	v2=hydro(:,2*ja-1)-Oxy(:,ja)
  	 	if (PBC) v2=v2 - length*anint(v2/length)!PBC
		summ = v1+v2!find vector to M-site
   	 	v3=(summ/sqrt(dot_product(summ,summ)))*.1546
  		d2=(v1*qH+v2*qH+v3*qO)*ang2m*e2coul/3.33564e-30!conv. to Debye
	  else 
          	v1=hydro(:,2*ja-0)-Oxy(:,ja)
	  	if (PBC) v1=v1 - length*anint(v1/length)!PBC
          	v2=hydro(:,2*ja-1)-Oxy(:,ja)
  	 	if (PBC) v2=v2 - length*anint(v2/length)!PBC
          	d2=v1*qH + v2*qH
	  endif 

	  !Find the 'x' and 'z' components (in the frame of the molecule!) 

	  r = Oxy(:,ja) - Oxy(:,ia) !order here is important
	  if (PBC) r = r - length*anint(r/length)!PBC
 	  
    	  zhat = d1/sqrt(dot_product(d1,d1))
	
 	  x = R - dot_product(R,zhat)*zhat
	    
	  z = r - x 
		 	  
            
	  distance = sqrt(dot_product(r,r))
          rhat = r/distance

	  xd = sqrt( sum(x**2) )  
          zd = sqrt( sum(z**2) ) * dot_product(z,d1)/abs(dot_product(z,d1))!preserve sign
          
	  dip_fac = dot_product(d1,d2)
	  D_fac = 3*dot_product(d1,rhat)*dot_product(d2,rhat) - dip_fac

	  if (NORMFLAG) then
	        dip_fac = dip_fac/(  sqrt(dot_product(d1,d1))*sqrt(dot_product(d2,d2)) )
		D_fac = D_fac/(  sqrt(dot_product(d1,d1))*sqrt(dot_product(d2,d2)) )
	  endif 
          if( (xd .lt. MINVAL(halflength)).and.(abs(zd) .lt. MINVAL(halflength)) ) then
                countx = ceiling(xd/div_size)
		countz = ceiling(zd/div_size) + floor(real(numzbins)/2)
	
		g_hist(countz,countx)   =  g_hist(countz,countx)   + 1
		dip_hist(countz,countx) =  dip_hist(countz,countx) + dip_fac
		D_hist(countz,countx)   =  D_hist(countz,countx)   + D_fac
		
          endif
       enddo
    enddo

end subroutine m_2DRDF



!------------------------------------------------------------------
!------------ initial allocations --------------------------------
!------------------------------------------------------------------
subroutine allocate_2DRDF
	numxbins = ceiling(MINVAL(halflength)/div_size)
	numzbins = 2*numxbins + 1 

	allocate(g_hist(numzbins,numxbins))
	allocate(dip_hist(numzbins,numxbins))
	allocate(D_hist(numzbins,numxbins))
	allocate(num_hist(numzbins,numxbins))
	g_hist = 0 
end subroutine allocate_2DRDF


!------------------------------------------------------------------
!------------  compute the appropriate averages and normalize ----
!------------------------------------------------------------------
subroutine normalize_2DRDF 
 Implicit none 
 real :: rho, vol, histgas
 integer :: i, j

 do i = 1,numzbins 
	do j = 1,numxbins
		if (g_hist(i,j) .ne. 0) then
			dip_hist(i,j) = dip_hist(i,j)/g_hist(i,j)
			D_hist(i,j)   = D_hist(i,j)/g_hist(i,j) 
		endif
	enddo
 enddo 

 num_hist = g_hist
 rho = Nmol / (   length(1)*length(2)*length(3)  ) !"gas" number density

 do i = 1,numxbins
        vol=pi*((i*div_size)**2-((i-1)*div_size)**2)*div_size 
        histgas=rho*vol
	g_hist(:,i) = 2*g_hist(:,i)/(histgas*real(nsteps)*real(Nmol))
 enddo

end subroutine normalize_2DRDF

!------------------------------------------------------------------
!---------------- open output files & output data ----------------
!------------------------------------------------------------------
subroutine write_2DRDF
 Implicit none 
 integer :: i, j

 Open(41,file=trim(fileheader)//"g_hist.dat",status="unknown")
 Open(42,file =trim(fileheader)//"dip_hist.dat",status="unknown")
 Open(43,file=trim(fileheader)//"D_hist.dat",status="unknown")
 Open(44,file=trim(fileheader)//"num_hist.dat",status="unknown")

 do i = 1, numzbins
	WRITE(41,*) (g_hist(i,j),   j=1,numxbins)
	WRITE(42,*) (dip_hist(i,j), j=1,numxbins)
	WRITE(43,*) (D_hist(i,j),   j=1,numxbins)
	WRITE(44,*) (num_hist(i,j), j=1,numxbins)
 enddo

 close(41)
 close(42)
 close(43)
 close(44)

 deallocate(g_hist)
 deallocate(dip_hist)
 deallocate(D_hist)
 deallocate(num_hist)
end subroutine write_2DRDF

end module 
