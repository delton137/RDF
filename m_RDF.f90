!-------------------------  m_RDF module --------------------------------------
!--- contains subroutines for OO RDF, OH RDF, dip-dip/cosine function and 
!--- the distance dependent Kirkwood function G_K(r)
!-----------------------------------------------------------------------------

module m_RDF
use global_vars
Implicit None

contains 
!-----------------------------------------------------------------------
!------------------------ O-O RDF -------------------------------------
!-----------------------------------------------------------------------
  subroutine RDF_OO(Nmol, xa, hist)
    implicit none
    ! DUMMY ARGUMENTS
    integer,intent(in) :: Nmol
    real,dimension(3,Nmol),intent(in) :: xa
    real,dimension(Ndiv),intent(out) :: hist
    
    ! INTERNAL VARS
    real :: rho, vol, distance, tmp, histgas
    integer :: ia, ja, ix, i, Ntot ! loop flags
    
    ! ARRAYS
    real,dimension(Ndiv) :: histliquid ! pre-normalization histogram

    histliquid = 0.0d0

    rho = Nmol / (length(1)*length(2)*length(3))

! ADJUST COORDINATES FOR PERIODIC BOUNDING
    do ia=1,Nmol-1                 ! do for every molecule
       do ja=ia+1,Nmol             ! to every other molecule
          distance=0.0d0            
          do ix=1,3
             if(abs(xa(ix,ia)-xa(ix,ja)) >= minval(halflength)) then
                if(xa(ix,ia)-xa(ix,ja) > 0) then        
                   tmp=xa(ix,ja)+minval(halflength)*2.0d0       
                else                                    
                   tmp=xa(ix,ja)-minval(halflength)*2.0d0       
                end if
                distance=distance+(xa(ix,ia)-tmp)**2  
             else                                       
                tmp=xa(ix,ja)                           
                distance=distance+(xa(ix,ia)-tmp)**2    
             endif
          enddo
          distance=sqrt(distance)                      
 ! ACTUALIZE HISTOGRAM
          if(distance .lt. minval(halflength)) then
             call updatehist(distance,histliquid)
          endif

       enddo
    enddo
! NORMALIZE HISTOGRAM
    do i=1,Ndiv    
       vol=(4.0d0/3.0d0)*pi*((i*delta)**3-((i-1)*delta)**3)            
       histgas=rho*vol 
       hist(i)= hist(i) + histliquid(i)/(histgas*Nmol)
    enddo
  end subroutine RDF_OO


!-------------------------------------------------------------------------------
!-------------------- OH RDF --------------------------------------------------
!-------------------------------------------------------------------------------
  subroutine RDF_OH(Nmol, Oxy, Hydro, hist)
    implicit none

    ! DUMMY ARGUMENTS
    integer,intent(in) :: Nmol
    real,dimension(3,Nmol),intent(in) :: Oxy
    real,dimension(3,2*Nmol),intent(in) :: Hydro
    real,dimension(Ndiv),intent(inout) :: hist
    
    ! INTERNAL VARS
    real :: rho, vol, distance, tmp, histgas
    integer :: ia, ja, ix, i ! loop flags
    ! ARRAYS
    real,dimension(Ndiv) :: histliquid ! pre-normalization histogram

    ! INITIALIZE VARIABLES
    histliquid = 0.0d0
    rho = 2*Nmol / (length(1)*length(2)*length(3))

    ! ADJUST FOR PERIODIC BOUNDING AND CALCULATE DISTANCE
    do ia=1,Nmol
       do ja=1,2*Nmol 
          distance=0.0d0
          do ix=1,3
             if(abs(Oxy(ix,ia)-Hydro(ix,ja)) >= minval(halflength)) then 
                if(Oxy(ix,ia)-Hydro(ix,ja) > 0) then  
                   tmp=Hydro(ix,ja)+minval(halflength)*2.0d0  
                else                                  
                   tmp=Hydro(ix,ja)-minval(halflength)*2.0d0       
                end if
                distance=distance+(Oxy(ix,ia)-tmp)**2 
             else                                       
                tmp=Hydro(ix,ja)                            
                distance=distance+(Oxy(ix,ia)-tmp)**2   
             endif
          enddo
          distance=sqrt(distance)                     

 ! ACTUALIZE HISTOGRAM
          if(distance .lt. minval(halflength)) then
             call updatehist(distance,histliquid)     
          endif
       enddo
    enddo

! NORMALIZE HISTOGRAM
    do i=1,Ndiv   
       vol=(4.0d0/3.0d0)*pi*((i*delta)**3-((i-1)*delta)**3) 
       histgas=rho*vol
       hist(i)= hist(i) + histliquid(i)/(histgas*2*Nmol)
    enddo
    
end subroutine RDF_OH


!-----------------------------------------------------------------------
!---------------------- Cosine function -------------------------------
!-----------------------------------------------------------------------
subroutine RDF_OO_dip(Nmol, xao, xah, histPOS, histNEG, histNumMolsPOS, histNumMolsNEG, TTM3Fdips,TTM3F,coshist)
    implicit none
    ! DUMMY ARGUMENTS
    integer,intent(in) :: Nmol
    logical,intent(in) :: TTM3F
    real,dimension(3,Nmol),intent(in) :: xao, TTM3Fdips
    real,dimension(3,Nmol*2),intent(in) :: xah
    real,dimension(Ndiv),intent(inout) :: histPOS, histNEG, coshist
    Integer,dimension(Ndiv),intent(inout) :: histNumMolsPOS, histNumMolsNEG

    ! INTERNAL VARS
    real :: rho, vol, distance, tmp, lever, avg_dip
    integer :: ia, ja, ix, i, j, Ntot
    real, dimension(3) :: v1, v2, v3, d1, d2, tmp3, summ
    real :: dip_fac, cos_fac
    real,parameter :: e2coul=1.60217646e-19
    real,parameter :: ang2m=1e-10

    ! ARRAYS
    integer :: count
	
    rho = Nmol / (length(1)*length(2)*length(3))

! ADJUST COORDINATES FOR PERIODIC BOUNDING
    do ia=1,Nmol-1                
       ! compute the dipole vector for this molecule
       ! fixed by D. Elton
  if (TTM3F) then
	   	d1 = TTM3Fdips(1:3,ia)
	else if (TIP4P) then
	   	v1=xah(:,2*ia-0)-xao(:,ia)
   		v1=v1 - length*anint(v1/length)!PBC
	        v2=xah(:,2*ia-1)-xao(:,ia)
   		v2=v2 - length*anint(v2/length)!PBC
   		summ = v1+v2!find vector to M-site
   	 	v3=(summ/sqrt(dot_product(summ,summ)))*rOM
      		d1=(v1*qH+v2*qH+v3*qO)*ang2m*e2coul/3.33564e-30!conv. to Debye
	  	!write(*,*) sqrt(dot_product(d1,d1))
	else
   		v1=v1 - length*anint(v1/length)!PBC
  		v2=xah(:,2*ia-1)-xao(:,ia)
   		v2=v2 - length*anint(v2/length)!PBC
   		d1=v1*qH+v2*qH
	endif


  do ja=ia+1,Nmol             ! to for every other molecule
	  if (TTM3F) then
	   	d2 = TTM3Fdips(1:3,ja)
    else if (TIP4P) then
	  	v1=xah(:,2*ja-0)-xao(:,ja)
	  	v1=v1 - length*anint(v1/length)!PBC
     	v2=xah(:,2*ja-1)-xao(:,ja)
  	 	v2=v2 - length*anint(v2/length)!PBC
   		summ = v1+v2!find vector to M-site
   	 	v3=(summ/sqrt(dot_product(summ,summ)))*rOM
  		d2=(v1*qH+v2*qH+v3*qO)*ang2m*e2coul/3.33564e-30!conv. to Debye
	  else 
     	v1=xah(:,2*ja-0)-xao(:,ja)
	  	v1=v1 - length*anint(v1/length)!PBC
     	v2=xah(:,2*ja-1)-xao(:,ja)
  	 	v2=v2 - length*anint(v2/length)!PBC
    	d2=v1*qH+v2*qH
	  endif 

    tmp3 = xao(:,ia)-xao(:,ja)
    tmp3  = tmp3  - length*anint(tmp3/length)!PBC
    distance = sum(tmp3**2)
    distance= sqrt(distance)
    dip_fac = dot_product(d1,d2)
    cos_fac = dip_fac/(  sqrt(dot_product(d1,d1))*sqrt(dot_product(d2,d2))  )

 !--------------- Put stuff in histogram bins  ----------------------------------
 ! if(distance <= minval(halflength)) then
                count=int(distance/delta)
		!lever=mod(distance,delta)/delta
		if (dip_fac .gt. 0) then
			histPOS(count)=histPOS(count) + dip_fac!*(1.0d0-lever)
			!histPOS(count+1)=histPOS(count) + dip_fac*lever
			histNumMolsPOS(count)=histNumMolsPOS(count) + 1
		endif
		if (dip_fac .lt. 0) then
			histNEG(count)=histNEG(count) + dip_fac!*(1.0d0-lever)
			!histNEG(count+1)=histNEG(count) + dip_fac*lever
			histNumMolsNEG(count)=histNumMolsNEG(count) + 1
		endif	
		coshist(count) = coshist(count) + cos_fac	
 ! endif
    enddo
enddo
end subroutine RDF_OO_dip



!-------------------------------------------------------------------------------
!------------------------Kirkwood G factor calculation ------------------------
!-------------------------------------------------------------------------------
  subroutine GKRcalc(Nmol, xao, xah, gKr, gka, gke, TTM3Fdips,TTM3F)
  implicit none
    ! DUMMY ARGUMENTS
    integer,intent(in) :: Nmol
    logical,intent(in) :: TTM3F
    real,dimension(3,Nmol),intent(in) :: xao, TTM3Fdips
    real,dimension(3,Nmol*2),intent(in) :: xah
    real,dimension(Ndiv) :: gKr, gka, gke 
    ! INTERNAL VARS
    real :: distance
    integer :: ia, ja, ix, i, j, Ntot
    real, dimension(3) :: v1, v2, v3, d1, d2, tmp3, summ
    real :: dip_fac, avgdip2, costheta, totgka, totgke 
    real,parameter :: e2coul=1.60217646e-19
    real,parameter :: ang2m=1e-10

    ! ARRAYS
    integer :: count
	
do ia=1,Nmol-1               
       ! compute the dipole vector for this molecule
       ! fixed by D. Elton
        if (TTM3F) then
	   	d1 = TTM3Fdips(1:3,ia)
	else if (TIP4P) then
		v1=xah(:,2*ia-0)-xao(:,ia)
   		v1=v1 - length*anint(v1/length)!PBC
		v2=xah(:,2*ia-1)-xao(:,ia)
       		v2=v2 - length*anint(v2/length)!PBC
		summ = v1+v2!find vector to M-site
   	 	v3=(summ/sqrt(dot_product(summ,summ)))*.1546
      		d1=(v1*qH+v2*qH+v3*qO)*ang2m*e2coul/3.33564e-30!conv. to Debye
	  	!write(*,*) sqrt(dot_product(d1,d1))
	else 
       		v1=v1 - length*anint(v1/length)!PBC
      		v2=xah(:,2*ia-1)-xao(:,ia)
       		v2=v2 - length*anint(v2/length)!PBC
       		d1=v1*qH+v2*qH
	endif
       avgdip2 = 0 
       tgka = 0
       tgke = 0 
	
       do ja = 1,Nmol         
	  if (TTM3F) then
	   	d2 = TTM3Fdips(1:3,ja)
          else if (TIP4P) then
		v1=xah(:,2*ja-0)-xao(:,ja)
	  	v1=v1 - length*anint(v1/length)!PBC
          	v2=xah(:,2*ja-1)-xao(:,ja)
  	 	v2=v2 - length*anint(v2/length)!PBC
		summ = v1+v2!find vector to M-site
   	 	v3=(summ/sqrt(dot_product(summ,summ)))*.1546
  		d2=(v1*qH+v2*qH+v3*qO)*ang2m*e2coul/3.33564e-30!conv. to Debye
	  else 
          	v1=xah(:,2*ja-0)-xao(:,ja)
	  	v1=v1 - length*anint(v1/length)!PBC
          	v2=xah(:,2*ja-1)-xao(:,ja)
  	 	v2=v2 - length*anint(v2/length)!PBC
          	d2=v1*qH+v2*qH
	  endif 
	  avgdip2 = avgdip2 + sum(d2**2)
          distance=0.0d0
          tmp3 = xao(:,ia)-xao(:,ja)
	  tmp3  = tmp3  - length*anint(tmp3/length)!PBC
	  distance = sum(tmp3**2)
          distance= sqrt(distance)
          dip_fac = dot_product(d1,d2)
	  
	  costheta = dot_product(d1,tmp3)/sqrt(dot_product(tmp3,tmp3)*dot_product(d1,d1))

	  count=int(distance/delta)
		if (costheta .lt. .577) then 
			tgka(count) = tgka(count) + dip_fac
		endif
		if (costheta .gt. .577) then 
			tgke(count) = tgke(count) + dip_fac
		endif             
	 	
  enddo

	avgdip2 = avgdip2/Nmol


enddo

end subroutine GKRcalc

!-----------------------------------------------------------------------
!----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine updatehist(x,histliquid)

    implicit none

    ! INPUT VARIABLES
    real :: x,lever
    real, dimension(Ndiv) :: histliquid

    ! INTERNAL VARIABLES
    integer :: i,j,count

    ! PERFORM LEVER SMOOTHING CALCULATIONS
    count=floor(x/delta)     
   ! lever=mod(x,delta)/delta
    !there is a factor of two because each distance is for 2 atoms
    if (count .ne. 0) then
    histliquid(count) = histliquid(count) + 2
    !histliquid(count)=histliquid(count)+2.0d0*(1.0d0-lever)
   ! histliquid(count+1)=histliquid(count+1)+2.0d0*lever
	endif

  end subroutine updatehist

end module 
