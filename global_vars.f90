module global_vars
real,parameter :: pi = acos(-1.0d0) 
real,dimension(3), save  :: length, halflength
Integer,save  :: Ndiv
real,save  :: delta
real,save  :: qO, qH
real,save  :: rOM
integer :: nsteps
integer,save  :: DipType
integer,save  :: Na, Nmol
real,parameter :: e2coul=1.60217646e-19
real,parameter :: ang2m=1e-10
 character(120)   :: fileheader
logical   :: POSNEG,TIP4P,TIP4P2005,TIP4P2005F,SPCE,TTM3F,GEN,PBC,NORMFLAG
real :: totgka, totgke , avgdip2
real,dimension(:),allocatable :: tgke, tgka

end module global_vars
