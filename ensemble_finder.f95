	program inputreader
	IMPLICIT NONE
	integer :: nfrs,n,nmol
	CHARACTER(32) :: anly
	open(unit=5,file="protname.txt",status='old')
	read(5,*)
	read(5,*)n
	read(5,*)nfrs
	close(5)
	!nfrs=5000 !Testing
	WRITE(*,*)n,nfrs
	OPEN(unit=10,file='anly',status='old')
	READ(10,*)anly
	CLOSE(10)
	anly=adjustl(anly)
	IF(TRIM(anly) .EQ. "le4pd")call lproject(n,nfrs)
	IF(TRIM(anly) .EQ. "pca")call project(n,nfrs)
	IF(TRIM(anly) .EQ. "DNA")call dnaproject(n,nfrs)
	STOP
	End program inputreader

	subroutine lproject(n,nfrs)
	IMPLICIT NONE
	INTEGER :: n,nfrs,nmol,mode1,mode2,ncoords,nbins,ibin,jbin,c1bin,c2bin
	INTEGER :: i,j,k,a,counter,countery,counterz,maxframe(1),minframe(1)
	DOUBLE PRECISION :: xix(3*n,nfrs),xiy(3*n,nfrs),xiz(3*n,nfrs),xixavg(3*n),xiyavg(3*n),xizavg(3*n)
        DOUBLE PRECISION :: xim(3*n,nfrs),xfluct(3*n),yfluct(3*n),zfluct(3*n),ximavg(3*n),lx(3*n,nfrs)
	DOUBLE PRECISION :: ly(3*n,nfrs),lz(3*n,nfrs),xiavg(3*n),xifluct(3*n),lml(3*n),dummy,dummy2
	DOUBLE PRECISION :: QINV(3*n,3*n),rx(n,nfrs),ry(n,nfrs),rz(n,nfrs),avbl,xi(3*n,nfrs),Q(3*n,3*n)
	DOUBLE PRECISION :: xi1(nfrs),xi2(nfrs),xmin,xmax,ymin,ymax,proj1,proj2
	DOUBLE PRECISION, ALLOCATABLE :: coords(:,:),coords2(:,:)
        CHARACTER(120) :: aa,protname,bb,cc

	!5' end
	xi1=0.0
	xi2=0.0
	proj1=0.0
	proj2=0.0
	!nfrs=INT(nfrs/2)
	OPEN(unit=24,file='ncoords',status='old')
	READ(24,*)ncoords
	CLOSE(24)
	WRITE(*,*)'Number of structures: ',ncoords
	WRITE(*,*)'Number of frames: ',nfrs

	OPEN(unit=11,file='mode',status='old')
	READ(11,*)bb
	CLOSE(11)
	bb=adjustl(bb)

	ALLOCATE(coords(ncoords,2))
	OPEN(unit=24,file='anly_'//TRIM(bb)//'.dat',status='old') !Read in the trajectory data
	counter=1
	DO k=1,nfrs
	  READ(24,*)xi1(counter),xi2(counter) !Theta, phi
	  counter=counter+1
	END DO
	CLOSE(24)
	WRITE(*,*)'nfrs :',counter-1
	!nfrs=counter-1

	!proj1: projection along the first mode
	!proj2: projection along the second mode

	!proj1=proj1-xmin
	!proj2=proj2-ymin
	
	OPEN(unit=24,file='coords_'//TRIM(bb)//'.dat',status='old') 
	DO i=1,ncoords
	  READ(24,*)coords(i,1),coords(i,2) !Read in the selected distance and twist values
	END DO
	CLOSE(24)
	WRITE(*,*)coords

	!Now that I have the projections along each mode, I need to find the corresponding 
	!structures from the projections of the trajectory.
	DO i=1,ncoords
	  WRITE(aa,*)i
	  aa=adjustl(aa)
	  WRITE(*,*)TRIM(aa)
	  OPEN(unit=24,file='frames_'//TRIM(aa)//'.ndx',status='unknown')
	  WRITE(24,'(A)')'[ State '//TRIM(aa)//']'
	  DO k=1,nfrs !counter-1
	    nbins=100
	    ibin=1+INT((xi1(k)/180)*nbins) !Theta bin
	    jbin=1+INT((xi2(k)/180)*nbins) !Phi bin
	    c1bin=1+INT((coords(i,1)/180)*nbins)
	    c2bin=1+INT((coords(i,2)/180)*nbins)
	    IF((ibin .EQ. c1bin) .AND. (jbin .EQ. c2bin) )THEN !Granularity of 0.5% (start here for now)
	      WRITE(24,*)k !Frame number
	      WRITE(*,*)k,xi1(k),xi2(k)
	      !WRITE(*,*)k
	    END IF
	  END DO
	  CLOSE(24)
	END DO


	end subroutine lproject

	subroutine project(n,nfrs)
	IMPLICIT NONE
	INTEGER :: n,nfrs,nmol,mode1,mode2
	INTEGER :: i,j,k,a,counter,countery,counterz,maxframe(1),minframe(1)
	DOUBLE PRECISION :: xix(3*n,nfrs),xiy(n-1,nfrs),xiz(n-1,nfrs),xixavg(n-1),xiyavg(n-1),xizavg(n-1)
        DOUBLE PRECISION :: xim(n-1,nfrs),xfluct(n-1),yfluct(n-1),zfluct(n-1),ximavg(n-1),lx(n-1,nfrs)
	DOUBLE PRECISION :: ly(n-1,nfrs),lz(n-1,nfrs),xiavg(n-1),xifluct(n-1),lml(n-1),dummy
	DOUBLE PRECISION :: QINV(n-1,n-1),rx(n,nfrs),ry(n,nfrs),rz(n,nfrs),avbl,xi(n-1,nfrs),Q(n-1,n-1)
	DOUBLE PRECISION :: xi1(nfrs),xi2(nfrs),xmin,xmax,ymin,ymax,proj1,proj2
        CHARACTER(120) :: aa,protname,bb,cc

	xi1=0.0
	xi2=0.0
	proj1=0.0
	proj2=0.0

	OPEN(unit=24,file='mode1',status='old')
	READ(24,*)mode1
	WRITE(*,*)'First mode :',mode1
	CLOSE(24)
	WRITE(aa,*)mode1
	aa=adjustl(aa)

	OPEN(unit=24,file='mode2',status='old')
	READ(24,*)mode2
	WRITE(*,*)'First mode :',mode2
	CLOSE(24)
	WRITE(bb,*)mode2
	bb=adjustl(bb)

	OPEN(unit=24,file='../pc_traj_'//TRIM(aa)//'.dat',status='old') !Read in the first mode
	DO k=1,nfrs
	  READ(24,*)dummy,xi1(k)
	END DO
	CLOSE(24)

	OPEN(unit=24,file='../pc_traj_'//TRIM(bb)//'.dat',status='old') !Read in the second mode
	DO k=1,nfrs
	  READ(24,*)dummy,xi2(k)
	END DO
	CLOSE(24)

	open(unit=10,file="xmax")	
	read(10,*)xmax
	close(10)
	open(unit=10,file="ymax")	
	read(10,*)ymax
	close(10)
	open(unit=10,file="xmin")	
	read(10,*)xmin
	close(10)
	open(unit=10,file="ymin")	
	read(10,*)ymin
	close(10)
	open(unit=10,file="x")	
	read(10,*)proj1
	close(10)
	open(unit=10,file="y")	
	read(10,*)proj2
	close(10)

	!proj1: projection along the first mode
	!proj2: projection along the second mode

	!proj1=proj1-xmin
	!proj2=proj2-ymin

	!Now that I have the projections along each mode, I need to find the corresponding 
	!structures from the projections of the trajectory.
	OPEN(unit=24,file='frame',status='old')
	READ(24,'(A)')cc
	cc=adjustl(cc)
	CLOSE(24)
	OPEN(unit=24,file='pca_frames_'//TRIM(aa)//'_'//TRIM(bb)//'_frame_'//TRIM(cc),status='unknown')
	DO k=1,nfrs
	  !IF(MOD(k,1000) .EQ. 0)WRITE(*,*)ABS((xi1(k)-proj1)/proj1),ABS((xi2(k)-proj2)/proj2)
	  !WRITE(*,*)k
	  IF((ABS((xi1(k)-proj1)/proj1) .LE. 0.005) .AND. (ABS((xi2(k)-proj2)/proj2) .LE. 0.005) )THEN !Granularity of 1.0% (start here for now)
	    WRITE(24,*)k
	    !WRITE(*,*)k
	  END IF
	END DO
	CLOSE(24)


	end subroutine project

	subroutine dnaproject(n,nfrs)
	IMPLICIT NONE
	INTEGER :: n,nfrs,nmol,mode1,mode2,ncoords
	INTEGER :: i,j,k,a,counter,countery,counterz,maxframe(1),minframe(1)
	DOUBLE PRECISION :: xix(3*n,nfrs),xiy(3*n,nfrs),xiz(3*n,nfrs),xixavg(3*n),xiyavg(3*n),xizavg(3*n)
        DOUBLE PRECISION :: xim(3*n,nfrs),xfluct(3*n),yfluct(3*n),zfluct(3*n),ximavg(3*n),lx(3*n,nfrs)
	DOUBLE PRECISION :: ly(3*n,nfrs),lz(3*n,nfrs),xiavg(3*n),xifluct(3*n),lml(3*n),dummy,dummy2
	DOUBLE PRECISION :: QINV(3*n,3*n),rx(n,nfrs),ry(n,nfrs),rz(n,nfrs),avbl,xi(3*n,nfrs),Q(3*n,3*n)
	DOUBLE PRECISION :: xi1(nfrs),xi2(nfrs),xmin,xmax,ymin,ymax,proj1,proj2
	DOUBLE PRECISION, ALLOCATABLE :: coords(:,:),coords2(:,:)
        CHARACTER(120) :: aa,protname,bb,cc

	!5' end
	xi1=0.0
	xi2=0.0
	proj1=0.0
	proj2=0.0
	!nfrs=INT(nfrs/2)
	OPEN(unit=24,file='ncoords',status='old')
	READ(24,*)ncoords
	CLOSE(24)
	WRITE(*,*)'Number of structures: ',ncoords
	WRITE(*,*)'Number of frames: ',nfrs
	!READ(*,*)ncoords
	!ncoords=3

	ALLOCATE(coords(ncoords,2))
	OPEN(unit=24,file='anly.dat',status='old') !Read in the trajectory data
	counter=1
	DO k=1,nfrs
	  !IF(MOD(k,2) .EQ. 0)THEN
	    READ(24,*)dummy,dummy2
	  !ELSE
	    READ(24,*)xi1(counter),xi2(counter) !Distance, twist
	    counter=counter+1
	  !END IF
	END DO
	CLOSE(24)
	WRITE(*,*)'nfrs :',counter-1
	!nfrs=counter-1

	!proj1: projection along the first mode
	!proj2: projection along the second mode

	!proj1=proj1-xmin
	!proj2=proj2-ymin
	
	OPEN(unit=24,file='coords.dat',status='old') 
	DO i=1,ncoords
	  READ(24,*)coords(i,1),coords(i,2) !Read in the selected distance and twist values
	END DO
	CLOSE(24)
	WRITE(*,*)coords

	!Now that I have the projections along each mode, I need to find the corresponding 
	!structures from the projections of the trajectory.
	DO i=1,ncoords
	  WRITE(aa,*)i
	  aa=adjustl(aa)
	  WRITE(*,*)TRIM(aa)
	  OPEN(unit=24,file='frames_'//TRIM(aa)//'.ndx',status='unknown')
	  WRITE(24,'(A)')'[ State '//TRIM(aa)//']'
	  DO k=1,nfrs !counter-1
	    IF((ABS((xi1(k)-coords(i,1))/coords(i,1)) .LE. 0.005) .AND. (ABS((xi2(k)-coords(i,2))/coords(i,2)) .LE. 0.005) )THEN !Granularity of 0.5% (start here for now)
	      WRITE(24,*)k !Frame number
	      WRITE(*,*)k,xi1(k),xi2(k)
	      !WRITE(*,*)k
	    END IF
	  END DO
	  CLOSE(24)
	END DO

	!For Trinucleotides:
	!##########################################################################################
	!3' end
	!xi1=0.0
	!xi2=0.0
	!proj1=0.0
	!proj2=0.0
	!nfrs=INT(nfrs/2)
	!OPEN(unit=24,file='ncoords2',status='old')
	!READ(24,*)ncoords
	!CLOSE(24)
	!WRITE(*,*)'Number of structures: ',ncoords
	!READ(*,*)ncoords
	!ncoords=3

	!ALLOCATE(coords2(ncoords,2))
	!OPEN(unit=24,file='anly.dat',status='old') !Read in the trajectory data
	!counter=1
	!DO k=1,nfrs
	!  IF(MOD(k,2) .EQ. 0)THEN
	!    READ(24,*)xi1(counter),xi2(counter) !Distance, twist
	!    counter=counter+1
	!  ELSE
	!    READ(24,*)dummy,dummy2
	!  END IF
	!END DO
	!CLOSE(24)
	!WRITE(*,*)'nfrs :',counter-1

	!proj1: projection along the first mode
	!proj2: projection along the second mode

	!proj1=proj1-xmin
	!proj2=proj2-ymin
	
	!OPEN(unit=24,file='coords2.dat',status='old') 
	!DO i=1,ncoords
	!  READ(24,*)coords2(i,1),coords2(i,2) !Read in the selected distance and twist values
	!END DO
	!CLOSE(24)
	!WRITE(*,*)coords2

	!Now that I have the projections along each mode, I need to find the corresponding 
	!structures from the projections of the trajectory.
	!DO i=1,ncoords
	!  WRITE(aa,*)i
	!  aa=adjustl(aa)
	!  WRITE(*,*)TRIM(aa)
	!  OPEN(unit=24,file='frames2_'//TRIM(aa)//'.ndx',status='unknown')
	!  WRITE(24,'(A)')'[ State '//TRIM(aa)//']'
	!  DO k=1,counter-1
	!    IF((ABS((xi1(k)-coords2(i,1))/coords2(i,1)) .LE. 0.005) .AND. (ABS((xi2(k)-coords2(i,2))/coords2(i,2)) .LE. 0.005) )THEN !Granularity of 0.1% (start here for now)
	!      WRITE(24,*)k !Frame number
	!      WRITE(*,*)k,xi1(k),xi2(k)
	      !WRITE(*,*)k
	!    END IF
	!  END DO
	!  CLOSE(24)
	!END DO

	DEALLOCATE(coords)
	!DEALLOCATE(coords2)
	end subroutine dnaproject

