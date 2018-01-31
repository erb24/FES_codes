	program inputreader
	integer nfrs,n,natoms
	CHARACTER(32) :: protname
	open(unit=5,file="protname.txt",status='old')
	read(5,*)protname
	read(5,*)n
	read(5,*)nfrs
	read(5,*)natoms
	close(5)
	nbins=60
	write(*,*)protname,n,nfrs,natoms
	call average(protname,n,nfrs,natoms)
	End program inputreader

	subroutine average(protname,n,nfrs,natoms)
	IMPLICIT NONE
	INTEGER :: n,nfrs,i,j,k,a,iatom,natoms,ncoords,resnum(n),atomnum(n)
	DOUBLE PRECISION :: rx,ry,rz,rxavg(n),ryavg(n),rzavg(n)
	CHARACTER(32) :: atom,res,protname,aa,mode
	CHARACTER(164) :: line
	CHARACTER(3),Dimension(n) :: resname

	OPEN(unit=11,file="ncoords",status='old')
	READ(11,*)ncoords
	CLOSE(11)

	OPEN(unit=11,file="mode",status='old')
	READ(11,*)mode
	CLOSE(11)
	mode=adjustl(mode)

	WRITE(*,*)ncoords
	
	protname=adjustl(protname)
	DO i=1,ncoords !Loop over clicks
	  rxavg=0.0
	  ryavg=0.0
	  rzavg=0.0
	  WRITE(aa,*)i
	  aa=adjustl(aa)
	  OPEN(unit=13,file=TRIM(protname)//'_'//TRIM(aa)//'.pdb')
	  WRITE(*,*)TRIM(protname)//'_'//TRIM(aa)//'.pdb'
	  !OPEN(unit=13,file='test',status='old')
	  OPEN(unit=12,file="nframes_"//TRIM(aa),status='old')
	  READ(12,*)nfrs
	  CLOSE(12)
	  DO k=1,nfrs-1 !Loop over conformations
	    j=0
	    !WRITE(*,*)k
	    DO a=1,5
	      READ(13,*)
	    END DO
	    DO a=1,natoms !Loop over atoms
	      READ(13,'(A)')line
	      !WRITE(*,*)line
	      READ(line(14:16),'(A)')res
	      IF (res .EQ. "CA")THEN !Select out the coordinates of the alpha-carbons
	        !WRITE(*,*)line
	        j=j+1
	        read(line(33:39),'(F8.3)')rx !Adjust due to chain information
	        read(line(41:47),'(F8.3)')ry !Adjust due to chain information
	        read(line(49:55),'(F8.3)')rz !Adjust due to chain information
	        rxavg(j)=rxavg(j)+rx
	        ryavg(j)=ryavg(j)+ry
	        rzavg(j)=rzavg(j)+rz
	        read(line(24:26),'(i3)')resnum(j)
	        read(line(8:11),'(i4)')atomnum(j)
	        read(line(18:20),'(A)')resname(j)
	      END IF
	    END DO
	    DO a=1,2
	      READ(13,*)	
	    END DO
	  END DO !End loop over frames
	  CLOSE(13)
	  DO j=1,n
	    rxavg(j)=rxavg(j)/REAL(nfrs-1)
	    ryavg(j)=ryavg(j)/REAL(nfrs-1)
	    rzavg(j)=rzavg(j)/REAL(nfrs-1)
	  END DO
	  !Write averaged structure to file
	  OPEN(unit=12,file='avg_struct_mode_'//TRIM(mode)//'_conf_'//TRIM(aa)//'.pdb',status='unknown')
	  DO j=1,n
	    WRITE(line(8:11),'(i4)')atomnum(j)
	    WRITE(line(13:16),'(A)')"CA"
	    WRITE(line(18:20),'(A)')resname(j)
	    WRITE(line(24:26),'(i3)')resnum(j)
	    write(line(33:38),'(F6.3)')rxavg(j)
	    write(line(41:46),'(F6.3)')ryavg(j) 
	    write(line(49:54),'(F6.3)')rzavg(j) 
	    !WRITE(line(78:78),'(A)')"C"
	    WRITE(12,'(A)')line
	    WRITE(*,*)line
	  END DO
	  CLOSE(12)
	END DO !End loop over conformations

	end subroutine
	
