!--------- 3D Angular Correlation Function (ACF) ----------------------
!-- This calculates the ACF introduced by Liu & Wu 
!-- REF: JCP 139 041103 (2013) 
!-- Dan Elton, 2014        
!----------------------------------------------------------------------
module ACF
use global_vars
Implicit None
integer, parameter :: nxbins = 50
integer, parameter :: nybins = 50
integer, parameter :: nzbins = 50
integer, parameter :: nbins = 10
real,dimension(:,:,:), allocatable, save :: rho_avg
real,dimension(:,:,:), allocatable, save :: chi
real,dimension(:,:,:,:,:,:), allocatable, save :: rho
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
    real, dimension(3) :: Xprime, Yprime, Zprime, Rprime, x2
    real :: theta, alpha, gamma
    real, dimension(3,2) :: A

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

 	  r = Oxy(:,ja) - Oxy(:,ia) 

	  !set up initial x-y-z frame unit vectors for molecule one
	  !these are called the local coordinates
	  x = hydro(:,2*ia-1) - hydro(:,2*ia-0)
   	  if (PBC) x = x - length*anint(x/length)!PBC

	  x = x/sqrt(dot_product(x,x))
 	  z = d1/sqrt(dot_product(d1,d1))

	  !y = x cross z 
	  y(1) = x(2)*z(3) - z(2)*x(3) 
	  y(2) = z(1)*x(3) - x(1)*z(3)
	  y(3) = x(1)*z(2) - z(1)*x(1)		

	  write(*,*) "check mag y =" , sqrt(dot_product(y,y))
	  y = y/sqrt(dot_product(y,y))
		

	  !find 2 unit vectors describing the second molecule's frame
	  !in the local coordintaes
	  call proj_local_coords(x, y, z, r, Rprime)

	  Zprime = Rprime/dot_product(Rprime, Rprime)

	  x2 = hydro(:,2*ja-1) - hydro(:,2*ja-0)
   	  if (PBC) x2 = x2 - length*anint(x2/length)!PBC

	  call proj_local_coords(x, y, z, x2, Xprime)

	  Xprime = Xprime/dot_product(Xprime, Xprime)

	  !Yprime = Xprime cross Zprime
	  Yprime(3) = Xprime(1)*Zprime(2) - Zprime(1)*Xprime(1)

	  !Find the three 'Euler angles' (see fig. 1 of JCP 139 041103)
 	  !Ref: https://en.wikipedia.org/wiki/Euler_angles#Geometric_derivation
	  theta = acos(Zprime(3))

	  write(*,*) "theta 1 = ", theta
	  write(*,*) "theta 2 = ", dot_product(d1,d2) /(  sqrt(dot_product(d1,d1))*sqrt(dot_product(d2,d2))  )

	  alpha = atan2(Zprime(1),-Zprime(2))
	  gamma = atan2(Xprime(3), Yprime(3))  
	
	  !add to histogram bins 
          countx = floor(Rprime(1)/d_x)
          county = floor(Rprime(2)/d_y)
	  countz = floor(Rprime(3)/d_z) + floor(real(nzbins)/2)
          countt = floor(theta/d_theta)
          counta = floor(alpha/d_alpha)
	  countg = floor(gamma/d_gamma)  

	  rho(countx,county,countz,countt,counta,countg) = rho(countx,county,countz,countt,counta,countg) + 1
	  rho_avg(countx,county,countz) = rho_avg(countx,county,countz) + 1
		
       enddo
    enddo

end subroutine make_ACF


!----------------------------------------------------------------------------------
!------- Do all the averaging at the end to construct the ACF --------------------
!----------------------------------------------------------------------------------
subroutine averaging_ACF
Implicit none
 Integer i, j, k 

 do i = 1, nxbins
	do j = 1, nybins
		do k = 1, nzbins
			chi(i,j,k) = sqrt(sum(  rho(i,j,k,:,:,:) - rho_avg(i,j,k) ) )/ rho_avg(i,j,k)
		enddo 
	enddo
 enddo

end subroutine averaging_ACF


!----------------------------------------------------------------------------------
!------- Write out data in a number of different formats -------------------------
!----------------------------------------------------------------------------------
subroutine write_ACF
Implicit none
 Integer i, j 


 Open(51,file=trim(fileheader)//"z",status="unknown")
 Open(52,file=trim(fileheader)//"x_y",status="unknown")
 Open(53,file=trim(fileheader)//"xyplane.dat",status="unknown")
 

 do i = 1, nzbins
	write(51,*) i*d_x-length(3)/2, chi(0,0,i) 
 enddo

 do i = 1, nxbins
	write(52,*) i*d_x, chi(i,0,0) , chi(0,j,0) 
 enddo

 do i = 1, nxbins
	write(52,*) (chi(i,j,0), j=1,nybins) 
 enddo

end subroutine write_ACF

!----------------------------------------------------------------------------------
!------------------- allocate and initialize -------------------------------------
!----------------------------------------------------------------------------------
subroutine init_ACF 
 allocate(rho_avg(nxbins, nybins, nzbins))
 allocate(chi(nxbins, nybins, nzbins))
 allocate(rho(nxbins, nybins, nzbins, nbins, nbins, nbins))
 d_x = length(1)/nxbins
 d_y = length(2)/nybins
 d_z = length(3)/nzbins
 d_theta = pi/nbins
 d_alpha = 2.0*pi/nbins
 d_gamma = 2.0*pi/nbins
end subroutine init_ACF 

!----------------------------------------------------------------------------------
!- Projects a vector R from global coordinates into the local (x, y, z) coordinates 
!- it does this by first projecting into the x-y plane
!- this may not be most efficient method, but it works. 
!----------------------------------------------------------------------------------
subroutine proj_local_coords(x, y, z, R, Rprime)
Implicit none
    real, dimension(3), intent(in) :: x, y, z, R
    real, dimension(3), intent(out) :: Rprime
    real, dimension(3,2) :: A
    real, dimension(2,3) :: At
    real, dimension(2,2) :: C
    real, dimension(3,2) :: C1
    real, dimension(3,3) :: C2
    real, dimension(3) :: rproj, zvec
	
    A(:,1) = (/  x(1), x(2), x(3) /)
    A(:,2) = (/  y(1), y(2), y(3) /)

    At = transpose(A) 

    Call inverse(matmul(At, A), C, 2)

    C1 = matmul(A, C)
    C2 = matmul(C1,At)

    rproj = matmul(C2,R)

    zvec = R - rproj

    Rprime(1) = dot_product(x, rproj)/dot_product(rproj, rproj) 
    Rprime(2) = dot_product(y, rproj)/dot_product(rproj, rproj) 
    Rprime(3) = sqrt(dot_product(zvec, zvec))*dot_product(zvec, z)/abs(dot_product(zvec, z))
   
end subroutine proj_local_coords 	


subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
implicit none 
integer n
real a(n,n), c(n,n)
real L(n,n), U(n,n), b(n), d(n), x(n)
real coeff
integer i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do
end subroutine inverse


end module 
