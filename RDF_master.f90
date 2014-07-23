!----------------------------------------------------------------------------- 
! Master RDF program, reads from an input file called called "RDF_master.inp"
! calls subroutines from m_RDF to calculate RDFs and dipole corr. functions
!-----------------------------------------------------------------------------
! Program rewritten by Dan Elton 4/13
!----------------------------------------------------------------------------- 
! .xtc support added by Dan Elton 3/14
!-----------------------------------------------------------------------------

Program RDF_Master
  use global_vars
  use m_RDF    
  use m_timer  
  use ACF
  use RDF2D
Implicit None

!------------------------ declare all the variables! -----------------------------
 character(len=3) :: sym     
 character(120) :: fileinp, TTM3F_dip_input,model
real,dimension(:,:),allocatable :: Oxy, Hydro,TTM3Fdips,runninggKr
real,dimension(:),allocatable :: gKr, gKr2, histOO, histOH, histHH,histOOdipPOS,histOOdipNEG, histOOdip, coshist, coshistn, gKrerr
real,dimension(:),allocatable :: gka, gke, diff
integer,dimension(:),allocatable :: histNumMolsPOS, histNumMolsNEG
real :: density, vol, histgas, junk, totdip, N, errpercent, sum_
real(8) ::  seconds
integer, parameter :: Lun1 = 21, Lun3 = 23, LunTTM3F = 24 
integer :: Nmol_HH, npts
integer :: i,j, ierror, istep, ia, ix, Ntot, errpts
logical :: first,DIPDIPFLAG,OOFLAG,OHFLAG,HHFLAG,COSFLAG,GKFLAG,GK2FLAG,ALLFLAG, ACFFLAG, RDF2DFLAG
 character(20), parameter  :: constants = "RDF_master.inp"
!Variables for reading XTC
 character(len=3) :: filetype    
Integer               :: ID1,STEP,RET, NAT, indx, maxsteps
Integer, parameter   :: MAGIC=1995, MAXAT=1000000001
real(4)                :: BOX(9),DT,TIME,PREC
real(4), dimension(:),allocatable :: X


call set_timer()

avgdip2 = 3
!----------------------- open input file and read -----  ---------------------------
Open(Lun3,file="RDF_master.inp",status="old",action="read",iostat=ierror)
read(Lun3,*) fileinp
read(Lun3,*) filetype 
read(Lun3,*) TTM3F_dip_input
read(Lun3,*) fileheader
read(Lun3,*) length
read(Lun3,*) Ndiv
read(Lun3,*) model
read(Lun3,*) qO
read(Lun3,*) qH
read(Lun3,*) nsteps
read(Lun3,*) ALLFLAG
read(Lun3,*) DIPDIPFLAG
read(Lun3,*) POSNEG
read(Lun3,*) COSFLAG
read(Lun3,*) GKFLAG
read(Lun3,*) GK2FLAG
read(Lun3,*) errpercent
read(Lun3,*) OOFLAG
read(Lun3,*) OHFLAG
read(Lun3,*) HHFLAG
read(Lun3,*) PBC
read(Lun3,*) Na
read(Lun3,*) ACFFLAG
read(Lun3,*) RDF2DFLAG
Close(Lun3)

!----------------------- Debug ---------------------------------------
if ( (GK2FLAG) .and. (.not.DIPDIPFLAG) ) then
	write(*,*) " G_K from dip-dip was selected but dip-dip was not enabled. exiting.."
endif
if ( (COSFLAG) .and. (.not.DIPDIPFLAG) ) then
	write(*,*) " cos(r) was selected but dip-dip was not enabled. exiting.."
endif

!----------------------------------------------------------------------------
!--------------------------- set up model ----------------------------------
!----------------------------------------------------------------------------

if (model == 'spce') then
	qO = -.8476
	qH = .4238
	write(*,*) "Model is SPC/E"
else if (model == 'tip3p') then
	qO = -0.834
	qH =  0.417 	
	write(*,*) "Model is TIP3P"
else if ((model == 'tip4eps') .or. (model == 'TIP4eps')) then
	qO = -1.054
	qH =  0.527
	rOM = .105
	TIP4P = .true. 
	write(*,*) "Model is TIP4eps"
else if (model == 'tip4p') then
	qO = -1.04
	qH = .52
	rOM = .15
	TIP4P = .true.
	write(*,*) "Model is TIP4P"
else if (model == 'tip4p2005') then
	qO = -1.1128
	qH = .5564
	rOM = .1546
	TIP4P = .true.
	write(*,*) "Model is TIP4P/2005"
else if (model == 'ttm3') then
	TTM3F = .true. 
	write(*,*) "Model is TTM3F - TTM3F SUPPORT NOT FULLY IMPLEMENTED YET!!"
else 
	write(*,*) "Model is generic 3 site"
endif 



!----------------------------------------------------------------------------
!-------------- read number of atoms and count number of timesteps ---
!----------------------------------------------------------------------------
if (filetype .eq. "xyz" ) then   
	open(Lun1,file=fileinp,status="old",action="read",iostat=ierror)
	if (TTM3F) Open(LunTTM3F,file=TTM3F_dip_input,status="unknown",action="read",iostat=ierror)
	if (ALLFLAG) then
		first = .true.
		nsteps = 0
		npts = 0
		do
  	 		if (first) then       
  	 			read(Lun1,*,iostat=ierror) Na
   	   			first = .false.
   	    	        else
   	    	 		Read(Lun1,*,iostat=ierror) 
   			endif  
   			if (ierror /= 0) then
   	    			exit
   			endif
  		npts = npts + 1 !Count number of lines
		enddo
	nsteps = floor(real(npts) / real((Na + 2)))! Number of timesteps
	endif
	if (.not. ALLFLAG) then
	   	read(Lun1,*,iostat=ierror) Na
	endif

	Nmol = Na / 3
endif

if (filetype .eq. "xtc" ) then   
        write(*,*) "Reading XTC. If this is a  4 site model we assume all 4 sites are in the XTC."
        write(*,*) "we assume the units are nm in the .xtc. They will be converted to Ang"
	if (ALLFLAG .and. GKFLAG) then 
		stop 
		write(*,*) "WARNING: program cannot count the number of steps with .xtc." 
		write(*,*) "To ensure gk(R) error calculation works, please set 'read all' to false" 
		write(*,*) "and make sure you set the number of steps correctly"
	endif

	if (.not. ALLFLAG) then 
		CALL XTCOPEN(ID1,fileinp,"R",MAXAT)
     		if (TIP4P) then
			Nmol = Na/4
    		else
	 		Nmol = Na/3
   		endif
		NAT = Na
	endif
endif
	
write(*,*) "number of steps: ", nsteps
write(*,*) "numer of errpts: ", errpts
write(*,*) "numer of molecules: ", Nmol
write(*,*) "numer of atoms: ", Na


Nmol_HH = Nmol * 2
errpts = floor((errpercent/100)*nsteps)
halflength = length / 2.0d0
delta = sqrt(halflength(1)**2 + halflength(2)**2 + halflength(3)**2 ) / real(Ndiv)

!----------------------------------------------------------------------------
!---------------------------- allocate all the arrays! --------------------- 
!----------------------------------------------------------------------------
rewind(Lun1)
allocate(Oxy(3,Nmol))
allocate(Hydro(3,Nmol*2))
allocate(TTM3Fdips(3,Nmol))
allocate(histOO(Ndiv))
allocate(histOOdipPOS(Ndiv))
allocate(histOOdipNEG(Ndiv))
allocate(histHH(Ndiv))
allocate(histOH(Ndiv))
allocate(histNumMolsPOS(Ndiv))
allocate(histNumMolsNEG(Ndiv))
allocate(histOOdip(Ndiv))
allocate(coshist(Ndiv))
allocate(coshistn(Ndiv))
if (GKFLAG) then
	allocate(gKrerr(Ndiv))
	allocate(diff(Ndiv))
	allocate(gKr(Ndiv))
	allocate(gka(Ndiv))
	allocate(gke(Ndiv))
	allocate(tgke(Ndiv))
	allocate(tgka(Ndiv))	
	allocate(runninggKr(Ndiv,errpts+2))
	gKr = 0
endif 
if (GK2FLAG) then 
	allocate(gkr2(Ndiv))
endif
if (filetype .eq. "xtc") allocate(X(3*Na))
histOOdipPOS = 0 
histOOdipNEG = 0 
histNumMolsPOS = 0 
histNumMolsNEG = 0 
histOOdip = 0
TTM3Fdips = 0
coshist = 0 
if (ACFFLAG) call init_ACF
if (RDF2DFLAG) call  allocate_2DRDF

!----------------------------------------------------------------------------
!---------------------------------- main loop ------------------------------ 
!----------------------------------------------------------------------------
Do istep = 1,nsteps
  !------------------------------------------------
  !----------- Reading xyz------------------------
  !------------------------------------------------
  if (filetype .eq. "xyz") then
        read(Lun1,*)
	read(Lun1,*)
	Do ia = 1, Nmol
		Read (Lun1,*) sym, (Oxy(ix,ia),ix=1,3)
		Read (Lun1,*) sym, (Hydro(ix,2*ia-0), ix=1,3)
		Read (Lun1,*) sym, (Hydro(ix,2*ia-1), ix=1,3)
		if (TTM3F) Read (LunTTM3F,*) (TTM3Fdips(ix,ia), ix=1,3)
	enddo
  endif
  !------------------------------------------------
  !----------- Reading xtc ------------------------
  !------------------------------------------------
  if (filetype .eq. "xtc") then

        CALL READXTC(ID1,NAT,STEP,TIME,BOX,X,PREC,RET)

	if (RET.EQ.0) then 
		write(*,*) "End of file at frame no. ", istep, " = ", nsteps
		nsteps = istep
		goto 1000
	endif 
        if ( .not. (NAT .eq. Na)) then
            write (*,*) 'ERROR: Number of Atoms in file is different than specified'
	    write(*,*)  NAT, " neq " , Na
            stop
        endif
	if (istep .eq. 1) write(*,*) "Box size from file is", BOX(1),BOX(5),BOX(9), " (nm)"
    

	!move coords from the X(:) array to the Oxy & Hydro arrays
	indx = 1
	!four site model case - the m sites are in the .xtc, these are considered the "oxy"
	if (TIP4P .eq. .true.) then
		do ia = 1, Nmol
			Oxy(:,ia)       = X(indx+0:indx+2)  
			Hydro(:,2*ia-0) = X(indx+3:indx+5)   
			Hydro(:,2*ia-1) = X(indx+6:indx+8)  
			indx = indx+12	
		enddo
	else
	!three site model case
		do ia = 1, Nmol
			Oxy(:,ia)       = X(indx+0:indx+2)   
			Hydro(:,2*ia-0) = X(indx+3:indx+5)   
			Hydro(:,2*ia-1) = X(indx+6:indx+8)  
			indx = ia*9  	
		enddo
	endif
	!Convert to Ang!!
	Hydro = 10*Hydro
	Oxy = 10*Oxy
  endif

  	if (OOFLAG) call RDF_OO(Nmol, Oxy, histOO)
        if (DIPDIPFLAG)call RDF_OO_dip(Nmol,Oxy,Hydro,histOOdipPOS,histOOdipNEG,histNumMolsPOS,histNumMolsNEG,TTM3Fdips,TTM3F,coshist)
  	if (HHFLAG) call RDF_OO(Nmol_HH, Hydro, histHH)
  	if (OHFLAG) call RDF_OH(Nmol, Oxy, Hydro, histOH)
	if (GKFLAG) then 
		call GKRcalc(Nmol, Oxy, Hydro, gKr, gka, gke, TTM3Fdips,TTM3F)
		if (istep .gt. (nsteps - errpts -  1)) then 
			do j = 1, Ndiv		       
				runninggKr(j,(istep - nsteps + errpts + 1 )) = gKr(j)/(Nmol*istep)
			enddo
		endif
	endif 
	if (ACFFLAG) call make_ACF(Oxy, Hydro, TTM3Fdips)
 	if (RDF2DFLAG) call m_2DRDF(Oxy, Hydro, TTM3Fdips)
   write(6,*) istep, nsteps
enddo! istep = 1,nsteps-1

1000 continue 


!----------------------------------------------------------------------------
!------------------------ Do all the normalizations! ----------------------- 
!----------------------------------------------------------------------------
!(there are two seperate functions for the pos. and neg. correlations)

if (DIPDIPFLAG) then
    density = Nmol / (length(1)*length(2)*length(3))
    do i = 1,Ndiv
	if ((histNumMolsNEG(i) + histNumMolsPOS(i)) .eq. 0) then 
		coshist(i) = 0 
	else 
		coshistn(i) = coshist(i)/ (histNumMolsNEG(i) + histNumMolsPOS(i) )
	endif
        
	vol=(4.0d0/3.0d0)*pi*((i*delta)**3-((i-1)*delta)**3) 
	histgas=density*vol
 	histOOdipPOS(i) = histOOdipPOS(i)/(histgas*real(nsteps)*Nmol)
	histOOdipNEG(i) = histOOdipNEG(i)/(histgas*real(nsteps)*Nmol)
	histOOdip(i) = ( histOOdipPOS(i)+histOOdipNEG(i) ) 
     enddo 
endif
! gkr stuff
if (GKFLAG ) then
      do i = 1, Ndiv
 	totgka = 0 
 	totgke = 0
		do j = 1,i 
			totgka = totgka + tgka(j)
			totgke = totgke + tgke(j) 
		enddo
	gKr(i) = gKr(i) + (totgka+totgke)/avgdip2 
	gka(i) = gka(i) + totgka/avgdip2
	gke(i) = gke(i) + totgke/avgdip2
	enddo

	gKr = gKr/(Nmol*nsteps)
	gka = gka/(Nmol*nsteps)
	gke = gke/(Nmol*nsteps)
endif

if (ACFFLAG) call averaging_ACF
if (RDF2DFLAG) call normalize_2DRDF

close(Lun1) 

!--------------------------- Get all the gKr errors! -----------------------------
if (GKFLAG) then
	gKrerr = 0 
	do i = 1, errpts
		diff = gKr - runninggKr(:,i)
		gKrerr = gkrerr + diff**2
	enddo
	gKrerr = sqrt(gKrerr/real(errpts))
endif 

!----------------------- Calculate G_K(r) by integration method! ------------------
if (GK2FLAG) then
	sum_ = 0
	do i = 1,Ndiv
		sum_ = sum_ + ( (i*delta )**3 - ((i-1)*delta)**3  )*histOOdip(i)
		gkr2(i) = sum_
	enddo
	gkr2 = (4.0d0/3.0d0)*pi*density*gkr2/avgdip2 + 1
endif


!--------------------------  Write out all the things! -------------------------
if (ACFFLAG) call write_ACF
if (RDF2DFLAG) call  write_2DRDF
if (DIPDIPFLAG) Open(17,file=trim(fileheader)//"dip.dat",status="unknown")
if (COSFLAG)    Open(43,file=trim(fileheader)//"cos.dat",status="unknown")
if (GKFLAG)     Open(16,file=trim(fileheader)//"gKr.dat",status="unknown") 
if (GK2FLAG)    Open(18,file=trim(fileheader)//"gKr_int.dat",status="unknown") 
if (OOFLAG)     Open(13,file=trim(fileheader)//"_OO.dat",status="unknown")
if (HHFLAG)     Open(14,file=trim(fileheader)//"_HH.dat",status="unknown")
if (OHFLAG)     Open(15,file=trim(fileheader)//"_OH.dat",status="unknown")

do i = 1, Ndiv
   if (DIPDIPFLAG)  then	
	if (POSNEG) then
	  write(17,110) real(i-1)*delta+delta/2, histOOdip(i), histOOdipPOS(i), -histOOdipNEG(i) 
	else 
	  write(17,100) real(i-1)*delta+delta/2, histOOdip(i)
	endif 
   endif
	if (COSFLAG) write(43,*) real(i-1)*delta+delta/2, coshistn(i)
	if (GKFLAG)  write(16,120) real(i-1)*delta+delta/2, gKr(i), gKrerr(i), gka(i), gke(i) 
	if (GK2FLAG) write(18,100) real(i-1)*delta+delta/2, gKr2(i) 
	if (OOFLAG)  write(13,*) real(i-1)*delta+delta/2, histOO(i)/real(nsteps)
	if (HHFLAG)  write(14,*) real(i-1)*delta+delta/2, histHH(i)/real(nsteps)
	if (OHFLAG)  write(15,*) real(i-1)*delta+delta/2, histOH(i)/real(nsteps)
enddo

100 Format (f12.6,f16.13) 
110 Format (f12.6,3(1x,f16.13))
120 Format (f12.4,4(1x,f10.4))


if (DIPDIPFLAG) Close(17)
if (COSFLAG)    Close(43)
if (GKFLAG)     Close(16)
if (OOFLAG)     Close(13)
if (HHFLAG)     Close(14)
if (OHFLAG)     Close(15)

call elapsed_time(seconds) ! finish timing and output elapsed seconds
End Program RDF_master
