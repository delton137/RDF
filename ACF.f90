!--------- 3D Angular Correlation Function (ACF) ----------------------
!-- This calculates the ACF introduced by Liu & Wu 
!-- REF: JCP 139 041103 (5113) 
!-- Dan Elton, 2014        
!----------------------------------------------------------------------
module ACF
use global_vars
Implicit None
integer, parameter :: nxbins = 10
integer, parameter :: nybins = 10
integer, parameter :: nzbins = 10
integer, parameter :: nbins = 10
real,dimension(:,:,:), allocatable  :: rho_avg
real,dimension(:,:,:), allocatable  :: chi
real,dimension(:,:,:,:,:,:), allocatable  :: rho
real :: d_theta, d_alpha, d_gamma, d_x, d_y, d_z
 
contains 

!-----------------------------------------------------------------------
 subroutine make_ACF(Oxy, Hydro, TTM3Fdips)
    Implicit None
    real,dimension(3,Nmol),intent(in)   :: Oxy, TTM3Fdips
    real,dimension(3,Nmol*2),intent(in) :: Hydro
    real :: vol, distance, tmp, xd, zd
    integer :: ia, ja, ix, countx, county, countz, countt, counta, countg
    real, dimension(3) :: v1, v2, v3, d1, d2, r, x, y, z, summ, rhat
    real, dimension(3) :: Xprime, Yprime, Zprime, Rprime, x2, y2, z2
    real :: theta, alpha, gamma
    real, dimension(3,2) :: A

    write(*,*) "sum rho_avg = ", sum(rho_avg)
    write(*,*) "sum rho = ", sum(rho)
	

    do ia=1, Nmol-1 
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

	  !set up initial x-y-z frame unit vectors for molecule one
	  !these are called the local coordinates
	  x = hydro(:,2*ia-1) - hydro(:,2*ia-0)
   	  if (PBC) x = x - length*anint(x/length)!PBC
	  x = x/sqrt(dot_product(x,x))

 	  z = d1/sqrt(dot_product(d1,d1))

	  y = cross_product(x,z)
	  y = y/sqrt(dot_product(y,y))


	  !find 3 unit vectors describing the second molecule's frame
	  !they are denoted x2, y2, z2
	  x2 = hydro(:,2*ja-1) - hydro(:,2*ja-0)
   	  if (PBC) x2 = x2 - length*anint(x2/length)!PBC
	  x2 = x2/sqrt(dot_product(x2,x2))

	  z2 = d2/sqrt(dot_product(d2,d2))

	  y2 = cross_product(x2, z2)
	  y2 = y2/sqrt(dot_product(y2,y2))

	  !transform 2nd molecules's unit vectors into local coordinates, they are denoted Xprime, Yprime, Zprime
	  call proj_local_coords(x, y, z, x2, Xprime)
	  call proj_local_coords(x, y, z, y2, Yprime)
	  call proj_local_coords(x, y, z, z2, Zprime)

	  !Find the three 'Euler angles' (see fig. 1 of JCP 139 041103)
 	  !Ref: https://en.wikipedia.org/wiki/Euler_angles#Geometric_derivation
	
	  theta = acos(Zprime(3))
	  alpha = atan2(Zprime(1),-Zprime(2))
	  gamma = atan2(Xprime(3), Yprime(3))  

	  if ( ISNAN(theta) )then
		theta = acos(dot_product(d1,d2)/sqrt(dot_product(d1,d1)*dot_product(d2,d2)))
		write(*,*) "NaN error with theta"
 	  endif
	  if ( ISNAN(alpha) )then
		alpha = 0 
		write(*,*) "NaN error with alpha"
 	  endif
	  if ( ISNAN(gamma) )then
		gamma = 0 
		write(*,*) "NaN error with gamma"
 	  endif

	  !project the distance vector r into local coordinates
	  r = Oxy(:,ja) - Oxy(:,ia) 
   	  if (PBC) r = r - length*anint(r/length)!PBC

	  call proj_local_coords(x, y, z, r, Rprime)

	 !write(*,*) "r = ", sqrt(dot_product(r,r))
	 !write(*,*) "R = ", sqrt(dot_product(Rprime,Rprime))

	  !add results to histogram bins 
          countx = ceiling(abs(Rprime(1))/d_x)
          county = ceiling(abs(Rprime(2))/d_y)
	  countz = ceiling(Rprime(3)/d_z) + nzbins 
          countt = ceiling(theta/d_theta) + ceiling(real(nbins)/2)+1
          counta = ceiling(alpha/d_alpha) + ceiling(real(nbins)/2)+1
	  countg = ceiling(gamma/d_gamma) + ceiling(real(nbins)/2)+1

	 if( sqrt(dot_product(Rprime,Rprime)) .lt. minval(length)) then
	  
	  !write(*,*) countx, county, countz, countt, counta, countg

          if( countx .lt. 1 ) then 
		write(*,*) "error countx= ", countx, "Rprime(1) = ", Rprime(1)
	  else if( county .lt. 1 ) then
		write(*,*) "error county= ", county, "Rprime(2) = ", Rprime(2)
	  else if( countz .lt. 1 ) then
		write(*,*) "error countz= ", county, "Rprime(3) = ", Rprime(3)
	  else if( countt .lt. 1 ) then
	  	write(*,*) "error countt= ", countt, "theta = ", theta
          else if( counta .lt. 1 ) then
		write(*,*) "error counta= ", counta, "alpha = ", alpha
          else if( countg .lt. 1 ) then 
		write(*,*) "error countg= ", countg, "gamma = ", gamma
	  else

		rho(countx,county,countz,countt,counta,countg) = rho(countx,county,countz,countt,counta,countg) + 1.0
		rho_avg(countx,county,countz) = rho_avg(countx,county,countz) + 1.0

 	  endif
	 endif		
       enddo
    enddo

end subroutine make_ACF


!----------------------------------------------------------------------------------
!- Projects a vector R from global coordinates into the local (x, y, z) coordinates 
!- it does this by first projecting into the x-y plane (called rproj)
!- x, y, and z are presumed to be unit vectors
!----------------------------------------------------------------------------------
subroutine proj_local_coords(x, y, z, R, Rprime)
Implicit none
    real, dimension(3), intent(in) :: x, y, z, R
    real, dimension(3), intent(out) :: Rprime
    real, dimension(3) :: rproj 

    Rprime(3) = dot_product(R,z)

    rproj = R - Rprime(3)*z
	    
    Rprime(1) = dot_product(x, rproj) 
    Rprime(2) = dot_product(y, rproj) 
   
end subroutine proj_local_coords 	



!----------------------------------------------------------------------------------
!------- Do all the averaging at the end to construct the ACF --------------------
!----------------------------------------------------------------------------------
subroutine averaging_ACF
Implicit none
 real, dimension(nbins,nbins) :: density_these_bins
 real num_angular_bins, avg_density, surface_area, average_these_bins
 Integer i, j, k, l, m 

 num_angular_bins = real(nbins)**3
 write(*,*) "sum rho = ", sum(rho)
 write(*,*) "sum rho_avg = ", sum(rho_avg)

 chi(i,j,k)  = 0 


do l = 0, nbins-1

   surface_area = (cos(d_theta*l) - cos(d_theta*(l+1)) )*d_alpha*d_gamma 

   do k = 1, nzbins
	do j = 1, nybins
 		do i = 1, nxbins
			if (rho_avg(i,j,k) .eq. 0) then 
				chi(i,j,k) = 0
			else
				avg_density = rho_avg(i,j,k)/(4.0*pi*2.0*pi)

				density_these_bins = rho(i,j,k,l+1,:,:)/surface_area

				chi(i,j,k) = chi(i,j,k) + sqrt( sum(  (density_these_bins - avg_density)**2 ) )/avg_density
			endif
		enddo 
	enddo
    enddo
enddo


end subroutine averaging_ACF




!----------------------------------------------------------------------------------
!------- Write out data in a number of different formats -------------------------
!----------------------------------------------------------------------------------
subroutine write_ACF
Implicit none
 Integer i, j, k 


 Open(51,file=trim(fileheader)//"z",status="unknown")
 Open(52,file=trim(fileheader)//"x_y",status="unknown")
 Open(53,file=trim(fileheader)//"xyplane.dat",status="unknown")
 
 write(51,'(a)') '# This .xvg is formated for xmgrace "'
 write(51,'(a)') '@ xaxis label "z(\cE\C)" '
 write(51,'(a)') '@ yaxis label "\f{Symbol}c\f{Times-Roman}(0,0,z)" '
 write(51,'(a)') '@ TYPE nxy '
 write(51,'(a)') '@ legend on '
 write(51,'(a)') '@ legend box off '
 write(51,'(a)') '@ legend loctype view '
 write(51,'(a)') '@ legend 0.78, 0.8'
 write(51,'(a)') '@ legend length 2'
 write(51,'(a)') '@ s0 legend \" ", "\"" '
 do i = 1, 2*nzbins
	write(51,*) i*d_x-length(3), chi(0,0,i) 
 enddo

 write(52,'(a)') '# This .xvg is formated for xmgrace "'
 write(52,'(a)') '@ xaxis label "x,y(\cE\C)" '
 write(52,'(a)') '@ yaxis label "\f{Symbol}c\f{Times-Roman}" '
 write(52,'(a)') '@ TYPE nxy '
 write(52,'(a)') '@ legend on '
 write(52,'(a)') '@ legend box off '
 write(52,'(a)') '@ legend loctype view '
 write(52,'(a)') '@ legend 0.78, 0.8'
 write(52,'(a)') '@ legend length 2'
 write(52,'(a)') '@ s0 legend \" ", "\"" '
 do i = 1, nxbins
	write(52,*) i*d_x, chi(i,0,0) , chi(0,j,0) 
 enddo

 do k = 1, nzbins
 do i = 1, nxbins
	write(53,*) (chi(i,j,k), j=1,nybins) 
 enddo
 enddo

end subroutine write_ACF

!----------------------------------------------------------------------------------
!------------------- allocate and initialize -------------------------------------
!----------------------------------------------------------------------------------
subroutine init_ACF 
 allocate(rho_avg(nxbins, nybins, 2*nzbins))
 allocate(chi(nxbins+1, nybins+1, nzbins+1))
 allocate(rho(nxbins+1, nybins+1, 2*nzbins+1, 2*nbins+1, 2*nbins+1, 2*nbins+1))
 d_x = length(1)/nxbins 
 d_y = length(2)/nybins
 d_z = length(3)/nzbins
 d_theta = 2.0*pi/nbins
 d_alpha = 2.0*pi/nbins
 d_gamma = 2.0*pi/nbins
end subroutine init_ACF 


!----------------------------------------------------------------------------------
!------------------- function to compute cross product ---------------------------
!----------------------------------------------------------------------------------
function cross_product(x,z)
 Implicit None
 real, dimension(3), intent(in)  :: x,z
 real, dimension(3) :: cross_product

  cross_product(1) = x(2)*z(3) - z(2)*x(3) 
  cross_product(2) = z(1)*x(3) - x(1)*z(3)
  cross_product(3) = x(1)*z(2) - z(1)*x(2)	

end function cross_product


end module 
