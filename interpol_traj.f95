	program inputreader
	IMPLICIT NONE
	integer :: nfrs,n,nmol,natoms
	CHARACTER(32) :: anly
	open(unit=5,file="protname.txt",status='old')
	read(5,*)
	read(5,*)n
	read(5,*)nfrs
	READ(5,*)natoms
	close(5)
	nfrs=10 !Number of interpolated frames
	WRITE(*,*)n,nfrs
	call traj(n,nfrs,natoms)
	STOP
	End program inputreader

	subroutine traj(n,nfrs,natoms)
	IMPLICIT NONE
	INTEGER :: n,nfrs,i,j,k,a,iatom,natoms,ncoords,resnum(n),atomnum(n),mode
	DOUBLE PRECISION :: rx(n,nfrs),ry(n,nfrs),rz(n,nfrs),rxcom(nfrs),rycom(nfrs),rzcom(nfrs),x(natoms,nfrs),y(natoms,nfrs)
	DOUBLE PRECISION :: rxavg(n),ryavg(n),rzavg(n),xim(3*n,nfrs),rdeg,degr,pi,z(natoms,nfrs)
	DOUBLE PRECISION :: QINV(3*n,3*n),xix(3*n,nfrs),xiy(3*n,nfrs),xiz(3*n,nfrs),theta(3*n,nfrs),phi(3*n,nfrs)
	CHARACTER(32) :: atom,res,protname,aa,bb
	CHARACTER(164) :: line
	CHARACTER(3),Dimension(n) :: resname

	OPEN(unit=12,file='mode',status='old')
	READ(12,*)mode
	CLOSE(12)
	WRITE(*,*)mode

	rxavg=0.0
	rxcom=0.0
	ryavg=0.0
	rycom=0.0
	rzavg=0.0
	rzcom=0.0
	QINV=0.0
	xix=0.0
	xiy=0.0
	xiz=0.0
	xim=0.0
	theta=0.0
	phi=0.0
	x=0.0
	y=0.0
	z=0.0

	pi=3.141592653589
	degr=((2.0*pi)/360.0) !deg to rad
	rdeg=1.0/degr !rad to deg

	OPEN(unit=12,file='3N_QHAINVmatrix.dat',status='old')
	DO a=1,3*n
	  DO i=1,3*n
	    READ(12,*)QINV(a,i)
	  END DO
	END DO
	CLOSE(12)

	OPEN(unit=12,file='com_avg_coor',status='old')
	DO i=1,n
	  READ(12,*)rxavg(i),ryavg(i),rzavg(i)
	END DO
	CLOSE(12)

	!Read in interpolated frames from PDB file
	WRITE(aa,*)mode
	aa=adjustl(aa)
        OPEN(unit=13,file='min'//TRIM(aa)//'.pdb')

        !OPEN(unit=13,file='test.pdb')
	!WRITE(*,*)'interpol1.pdb'
	!DO k=1,nfrs !Loop over conformations
          k=1 !Read minimum
	  DO i=1,5 !Read through TRJCONV header
            READ(13,*)
	  END DO
	  j=0
	  !WRITE(*,*)k
	  DO a=1,natoms !Loop over atoms
	    READ(13,'(A)')line
	    !WRITE(*,*)line
	    READ(line(14:16),'(A)')res
            read(line(33:39),'(F8.3)')x(a,k)
	    read(line(41:47),'(F8.3)')y(a,k)
	    read(line(49:55),'(F8.3)')z(a,k) 
	    IF (res .EQ. "CA")THEN !Select out the coordinates of the alpha-carbons
	      !WRITE(*,*)line
	      j=j+1
	      read(line(33:39),'(F8.3)')rx(j,k) !Adjust due to chain information
	      read(line(41:47),'(F8.3)')ry(j,k) !Adjust due to chain information
	      read(line(49:55),'(F8.3)')rz(j,k) !Adjust due to chain information
	      rx(j,k)=rx(j,k)/10.0 !Convert to nm
	      ry(j,k)=ry(j,k)/10.0 !Convert to nm
	      rz(j,k)=rz(j,k)/10.0 !Convert to nm
              !WRITE(*,*)rx(j,k),ry(j,k),rz(j,k)
	      rxcom(k)=rxcom(k)+rx(j,k)
	      rycom(k)=rycom(k)+ry(j,k)
	      rzcom(k)=rzcom(k)+rz(j,k)
	      !rxavg(j)=rxavg(j)+rx(j,k)
	      !ryavg(j)=ryavg(j)+ry(j,k)
	      !rzavg(j)=rzavg(j)+rz(j,k)
	      read(line(24:26),'(i3)')resnum(j)
	      read(line(8:11),'(i4)')atomnum(j)
	      read(line(18:20),'(A)')resname(j)
	    END IF
	  END DO
	  !READ(13,*)
	  rxcom(k)=rxcom(k)/REAL(n)
	  rycom(k)=rycom(k)/REAL(n)
	  rzcom(k)=rzcom(k)/REAL(n)	
	  DO j=1,n !Centre-of-mass coordinates
	    rx(j,k)=rx(j,k)-rxcom(k)-rxavg(j)
	    ry(j,k)=ry(j,k)-rycom(k)-ryavg(j)
	    rz(j,k)=rz(j,k)-rzcom(k)-rzavg(j)
	    !rxavg(j)=rxavg(j)+rx(j,k)
	    !ryavg(j)=ryavg(j)+ry(j,k)
	    !rzavg(j)=rzavg(j)+rz(j,k)
	    !WRITE(*,*)rx(j,k),ry(j,k),rz(j,k)
	  END DO
	  CLOSE(13)

          OPEN(unit=13,file='max'//TRIM(aa)//'.pdb')
	  k=nfrs  !Read maximum
	  DO i=1,5 !Read through TRJCONV header
            READ(13,*)
	  END DO
	  j=0
	  !WRITE(*,*)k
	  DO a=1,natoms !Loop over atoms
	    READ(13,'(A)')line
	    !WRITE(*,*)line
	    READ(line(14:16),'(A)')res
            read(line(33:39),'(F8.3)')x(a,k)
	    read(line(41:47),'(F8.3)')y(a,k)
	    read(line(49:55),'(F8.3)')z(a,k) 
	    IF (res .EQ. "CA")THEN !Select out the coordinates of the alpha-carbons
	      !WRITE(*,*)line
	      j=j+1
	      read(line(33:39),'(F8.3)')rx(j,k) !Adjust due to chain information
	      read(line(41:47),'(F8.3)')ry(j,k) !Adjust due to chain information
	      read(line(49:55),'(F8.3)')rz(j,k) !Adjust due to chain information
	      rx(j,k)=rx(j,k)/10.0 !Convert to nm
	      ry(j,k)=ry(j,k)/10.0 !Convert to nm
	      rz(j,k)=rz(j,k)/10.0 !Convert to nm
              !WRITE(*,*)rx(j,k),ry(j,k),rz(j,k)
	      rxcom(k)=rxcom(k)+rx(j,k)
	      rycom(k)=rycom(k)+ry(j,k)
	      rzcom(k)=rzcom(k)+rz(j,k)
	      !rxavg(j)=rxavg(j)+rx(j,k)
	      !ryavg(j)=ryavg(j)+ry(j,k)
	      !rzavg(j)=rzavg(j)+rz(j,k)
	      read(line(24:26),'(i3)')resnum(j)
	      read(line(8:11),'(i4)')atomnum(j)
	      read(line(18:20),'(A)')resname(j)
	    END IF
	  END DO
	  !READ(13,*)
	  rxcom(k)=rxcom(k)/REAL(n)
	  rycom(k)=rycom(k)/REAL(n)
	  rzcom(k)=rzcom(k)/REAL(n)	
	  DO j=1,n !Centre-of-mass coordinates
	    rx(j,k)=rx(j,k)-rxcom(k)-rxavg(j)
	    ry(j,k)=ry(j,k)-rycom(k)-ryavg(j)
	    rz(j,k)=rz(j,k)-rzcom(k)-rzavg(j)
	    !rxavg(j)=rxavg(j)+rx(j,k)
	    !ryavg(j)=ryavg(j)+ry(j,k)
	    !rzavg(j)=rzavg(j)+rz(j,k)
	    !WRITE(*,*)rx(j,k),ry(j,k),rz(j,k)
	  END DO
	  CLOSE(13)

        !Now interpolate between these frames
        DO k=2,9
	  DO j=1,natoms
	    x(j,k)=x(j,1)+(k-1)*((x(j,nfrs)-x(j,1))/(nfrs-1))
	    y(j,k)=y(j,1)+(k-1)*((y(j,nfrs)-y(j,1))/(nfrs-1))
	    z(j,k)=z(j,1)+(k-1)*((z(j,nfrs)-z(j,1))/(nfrs-1))
	  END DO
	END DO

	!Write interpolated coordinates to file and find Ca
	DO k=2,nfrs-1
	  WRITE(bb,*)k
	  bb=adjustl(bb)
	  OPEN(unit=13,file='min'//TRIM(aa)//'.pdb')
          OPEN(unit=14,file='interpol_'//TRIM(aa)//'.pdb',position='append',status='unknown')
	  WRITE(14,'(A)')'MODEL '//TRIM(bb)
	  DO i=1,5
            READ(13,*)
	  END DO
	  j=0
	  !WRITE(*,*)k
	  DO a=1,natoms !Loop over atoms
	    READ(13,'(A)')line
	    !WRITE(*,*)line
	    READ(line(14:16),'(A)')res
	    IF (res .EQ. "CA")THEN !Select out the coordinates of the alpha-carbons
	      !WRITE(*,*)line
	      j=j+1
	      read(line(33:39),'(F8.3)')rx(j,k) 
	      read(line(41:47),'(F8.3)')ry(j,k) 
	      read(line(49:55),'(F8.3)')rz(j,k) 
	      rx(j,k)=rx(j,k)/10.0 !Convert to nm
	      ry(j,k)=ry(j,k)/10.0 !Convert to nm
	      rz(j,k)=rz(j,k)/10.0 !Convert to nm
	      rx(j,k)=x(a,k)/10.0
	      ry(j,k)=y(a,k)/10.0 
	      rz(j,k)=z(a,k)/10.0 
              !WRITE(*,*)rx(j,k),ry(j,k),rz(j,k)
	      rxcom(k)=rxcom(k)+rx(j,k)
	      rycom(k)=rycom(k)+ry(j,k)
	      rzcom(k)=rzcom(k)+rz(j,k)
	      !rxavg(j)=rxavg(j)+rx(j,k)
	      !ryavg(j)=ryavg(j)+ry(j,k)
	      !rzavg(j)=rzavg(j)+rz(j,k)
	      read(line(24:26),'(i3)')resnum(j)
	      read(line(8:11),'(i4)')atomnum(j)
	      read(line(18:20),'(A)')resname(j)
	      WRITE(*,*)rx(j,k),ry(j,k),rz(j,k)
	    END IF
	    !Write the interpolated coordinates to file
	    write(line(32:38),'(F7.3)')x(a,k)
	    write(line(40:46),'(F7.3)')y(a,k)
	    write(line(48:54),'(F7.3)')z(a,k)
	    WRITE(14,'(A)')line
	  END DO
	  CLOSE(13)
	  rxcom(k)=rxcom(k)/REAL(n)
	  rycom(k)=rycom(k)/REAL(n)
	  rzcom(k)=rzcom(k)/REAL(n)	
	  DO j=1,n !Centre-of-mass coordinates
	    rx(j,k)=rx(j,k)-rxcom(k)-rxavg(j)
	    ry(j,k)=ry(j,k)-rycom(k)-ryavg(j)
	    rz(j,k)=rz(j,k)-rzcom(k)-rzavg(j)
	    !rxavg(j)=rxavg(j)+rx(j,k)
	    !ryavg(j)=ryavg(j)+ry(j,k)
	    !rzavg(j)=rzavg(j)+rz(j,k)
	    !WRITE(*,*)rx(j,k),ry(j,k),rz(j,k)
	  END DO
	  WRITE(14,'(A)')'TER'
	  WRITE(14,'(A)')'ENDMDL'
	END DO !End loop over interpolated frames
	CLOSE(14)
	
	DO k=1,nfrs	  
	  !Now calculate mode trajectories
	  a=mode
	  do j=1,3*n !residue loop
	    IF(MOD(j,3) .EQ. 1)xix(mode,k)=qinv(mode,j)*rx(j,k)+xix(mode,k)
	    IF(MOD(j,3) .EQ. 2)xiy(mode,k)=qinv(mode,j)*ry(j,k)+xiy(mode,k)
	    IF(MOD(j,3) .EQ. 0)xiz(mode,k)=qinv(mode,j)*rz(j,k)+xiz(mode,k)
	  end do
	  xim(a,k)=(xix(a,k)**2+xiy(a,k)**2+xiz(a,k)**2)**.5
!	  calculate theta, phi
	  theta(mode,k)=acos(xiz(mode,k)/xim(mode,k))
	  phi(mode,k)=atan(xiy(mode,k)/xix(mode,k))
	  if(xix(a,k).lt.0.0)phi(a,k)=phi(a,k)+pi
	  theta(a,k)=theta(a,k)*rdeg
	  if(phi(a,k).lt.0.0)phi(a,k)=phi(a,k)+2.0*pi
	  phi(a,k)=phi(a,k)*rdeg
	  WRITE(*,*)theta(a,k),phi(a,k)
	  WRITE(aa,*)mode
	  aa=adjustl(aa)
	END DO
	OPEN(unit=25,file='interpol'//TRIM(aa)//'_traj.dat',status='unknown')
	DO k=1,nfrs
	  !WRITE(*,*)'Writing the trajectory for mode ',a
	  WRITE(25,*)theta(mode,k),phi(mode,k)
	END DO
        CLOSE(25)
	
	end subroutine
