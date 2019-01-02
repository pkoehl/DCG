c===================================================================================
c===================================================================================
c
c	EnsembleEigendrv.f		Version 1 7/20/2016		Patrice Koehl
c									Jiahui Guan
c
c	This is the driver routine for EnsembleEigen; it is called from the
c	C++ code for DCG
c
c===================================================================================
c===================================================================================
c
	subroutine eigendrv(N, Neigen, Mat, eigen, eigen0, eigenvect, 
     1			WorkD, WorkL)
c
c===================================================================================
c===================================================================================
c	Input:
c		N: 				size of the matrix
c		Neigen:				number of eigenvalues to compute
c		mat:				Ensemble matrix
c		WorkD:				work array (double)
c		WorkL:				work array (bool)
c		
c	Output:
c		eigen:				the Neigen smallest eigenvalues
c						of the matrix
c		eigenvect:			the corresponding eigenvectors
c		eigen0:				the largest eigenvalue of the matrix
c===================================================================================
c===================================================================================
c
	integer		N,Neigen,i,Neig
	integer		m_krylov
	integer		k_V, k_d,k_workl,k_workd,k_work1
c
	real*8		eigen0,tol,X
	real*8		mat(*)
	real*8		WorkD(*)
	real*8		eigen(*),eigenvect(*)
c
	character	which*2
c
	logical		WorkL(*)
c
c===================================================================================
c	m_krylov needs to be at least 2*Neigen
c===================================================================================
c
	m_krylov = 2*Neigen
c 
c===================================================================================
c	Allocate space to the different arrays
c===================================================================================
c
	k_V     = 1
	k_d     = k_V + N*m_krylov
	k_workl = k_d + 2*m_krylov
	k_workd = k_workl + m_krylov*(m_krylov+8)
	k_work1 = k_workd + 3*N
c
c===================================================================================
c	First compute an estimate of the largest eigenvalue using Power Method
c===================================================================================
c
	tol = 1.d-8
	call PowerMethod(N,mat,WorkD(k_V),WorkD(k_workd),tol,eigen0)
c
c===================================================================================
c	Now compute the Neigen smallest eigenvalue (in fact, computes 2x Neigen)
c===================================================================================
c
	which = 'SM'
	call EnsembleEigen(N,mat,m_krylov,Neigen,WorkD(k_d),WorkD(k_V),
     1	eigen,eigenvect,WorkD(k_work1),WorkD(k_workl),WorkD(k_workd),
     2  WorkL,which)
c
	return
	end
c
c===================================================================================
c===================================================================================
c
c	EnsembleEigen.f		Version 1 7/20/2016		Patrice Koehl
c								Jiahui Guan
c
c	This code uses ARPACK to find a few eigenvalues 
c	(lambda) and corresponding eigenvectors (x) for the standard 
c	eigenvalue problem:
c          
c                        A*x = lambda*x
c 
c     	where A is an n by n real symmetric matrix.
c
c     The only thing that must be supplied in order to use this routine
c     is a function that performs a matrix-vector product
c
c                         w <-  Av
c
c	This specific implementation is designed to compute a few eigenvalues and the
c	corresponding eigenvectors of an Ensemble matrix to compute the DCG of a data
c	set, with the Ensemble matrix being stored with a sparse row condensed format.
c
c===================================================================================
c	Subroutine called
c===================================================================================
c
c	From ARPACK:
c
c	dsaupd  ARPACK reverse communication interface routine.
c	dseupd  ARPACK routine that returns Ritz values and (optionally)
c             	Ritz vectors.
c	dnrm2   Level 1 BLAS that computes the norm of a vector.
c	daxpy   Level 1 BLAS that computes y <- alpha*x+y.
c
c	For matrix-vector multiplication:
c
c	mat_mult : 	performs A*x, where A is a symmetric matrix
c
c===================================================================================
c===================================================================================
c
	subroutine EnsembleEigen(ndim,mat,m,nval,d,V,eigenval,
     1		eigenvect,work1,workl,workd,select,which)
c
c===================================================================================
c	Local Arrays
c===================================================================================
c
	integer		ndim,m,nval,i
c
	real*8		V(ndim,m)
	real*8		workl(*),workd(*),work1(*)
	real*8		d(m,2)
c
	logical		select(*)
c
	integer		iparam(11),ipntr(11)
c
c===================================================================================
c	Local Scalars
c===================================================================================
c
	character	bmat*1,which*2
c
	integer		ido,n,nev,ncv,lworkl,info,ierr
	integer		k,j,nx,ishfts,maxitr,mode1,nconv
c
	logical		rvec
c
	real*8		tol,sigma
c
c===================================================================================
c	Parameters
c===================================================================================
c
	real*8		zero
	parameter	(zero=0.0d0)
c  
c===================================================================================
c	BLAS & LAPACK routines used
c===================================================================================
c
	real*8		dnrm2
	external	dnrm2,daxpy
c
c===================================================================================
c	Intrinsic function
c===================================================================================
c
	intrinsic 	abs
c
c===================================================================================
c	Input for Hessian matrix
c===================================================================================
c
	real*8		mat(*)
c
c===================================================================================
c	Results of the calculation: eigenvalues, eigenvectors
c===================================================================================
c
	real*8		eigenval(*),eigenvect(ndim,*)
c
c ============================================================================
c 	Set dimensions for this problem
c 	(see header in dsaupd for the meaning of each variable)
c ============================================================================
c
	nx = ndim
	n = nx
c
	nev = m/2
	ncv = m
c
	bmat  = 'I'
c	which = 'LM'
c
	lworkl = ncv*(ncv+8)
	tol = zero 
	info = 0
	ido = 0
c
	ishfts = 1
	maxitr = 300 
	mode1 = 1
c
	iparam(1) = ishfts
	iparam(3) = maxitr
	iparam(7) = mode1
c
c ============================================================================
c 	Main loop to get Lanczos decomposition
c ============================================================================
c
100	continue
c
c               Call dsaupd repeatedly to build Krylov subspace
c               ===================================================
c
         	call dsaupd (ido,bmat,n,which,nev,tol,work1,ncv,V,
     1			ndim,iparam,ipntr,workd,workl,lworkl,info) 
c
c               Take action indicated by parameter IDO
c               ===================================================
c
		if(ido.eq.-1.or.ido.eq.1) then
c
c                       Perform matrix-vector multiplication
c                       ====================================
c
			call mat_mult(ndim,mat,workd(ipntr(1)),
     1			workd(ipntr(2)))
c
			goto 100
c
		endif
c
c       We have either reached convergence, or there is an error
c       =========================================================
c
	if(info.lt.0) then
c
		print *, ' '
		print *, ' Error with dsaupd, info = ', info
		print *, ' Check documentation in dsaupd '
		print *, ' '
c
	else 
c
c       ===============================================================
c	No fatal error: proceed and compute eigenvalues (Ritz values)
c	and eigenvectors (Ritz vectors) of Hessian
c       ===============================================================
c
		rvec = .true.
c
		call dseupd(rvec,'All',select,d,V,ndim,sigma,
     1			bmat,n,which,nev,tol,work1,ncv,V,ndim,
     2			iparam,ipntr,workd,workl,lworkl,ierr)
c
		if(ierr.ne.0) then
c
c		Error condition: check the documentation of dseupd
c       	======================================================
c
			print *, ' '
			print *, ' Error with dseupd, info = ', ierr
             		print *, ' Check the documentation of _deupd. '
             		print *, ' '
c
		else
c
c       	======================================================================
c		Only nconv (=iparam(5)) values have converged; corresponding
c		eigenvalues are returned in the first column of the two dimensional
c		array d and the corresponding eigenvectors are returned in the first
c		nconv columns of the matrix V, if requested
c       	======================================================================
c
			nconv =  iparam(5)
			nval = nconv
c
			do 300 j = 1,nconv
c
				eigenval(j) = d(j,1)
				do 200 k = 1,nx
					eigenvect(k,j) = v(k,j)
200				continue
c
300			continue
c
		endif
c
	endif
c
c ============================================================================
c	Done with program
c ============================================================================
c
	return
	end

c===================================================================================
c===================================================================================
c
c	mat_mult.f
c
c	This subroutine computes the product of a real, symmetric
c	matrix M, with a vector A.
c		(M A = B)
c
c===================================================================================
c===================================================================================
c
	subroutine mat_mult(N,Mat,A,B)
c
c===================================================================================
c===================================================================================
c	Input:
c		
c		- N:		Mat is NxN, and A and B are of size N
c		- Mat		the matrix 
c		- A:		left vector
c
c	Output:
c		- B:		result vector
c===================================================================================
c===================================================================================
c
	integer	i,j,k,N
c
	real*8	Mat(N,N),A(N),B(N)
c
c===================================================================================
c	Now loop over all non zero elements of M:
c===================================================================================
c
	do 200 i = 1,N
		B(i) = 0.d0
		do 100 j = 1,N
			B(i) = B(i) + Mat(i,j)*A(j)
100		continue
200	continue
c
	return
	end
c===================================================================================
c===================================================================================
c
c	PowerMethod.f
c
c	This subroutine finds the largest eigenvalue of a symmetric, real matrix
c
c===================================================================================
c===================================================================================
c
	subroutine PowerMethod(N,Mat,work1,work2,tol,eig)
c
c===================================================================================
c	Input:
c		
c		- N:		Mat is NxN
c		- Mat:		the matrix (full storage)
c		- work1:	work array (double) of size N
c		- work2:	work array (double) of size N
c		- tol:		tolerance on converged eigenvalue
c
c	Output:
c		- eig:		largest eigenvalue of H
c
c===================================================================================
c
	integer	N
	integer	i,iter,iter_max,iseed
c
	real*8	tol
	real*8	norm,fact,dnrm2
	real*8	eig0,eig
	real*8	Mat(*),work1(*),work2(*)
c
c===================================================================================
c	Initialize
c===================================================================================
c
	iter_max = 100
c
	iseed = -1
	call random_vect(N,iseed,work1)
c
c===================================================================================
c	Starts power iteration
c===================================================================================
c
	eig0 = 0.d0
c
	do 500 iter = 1,iter_max
c
c		====================================================================
c		Compute Hb, where b is the current estimate of the largest eigenvector
c		====================================================================
c
		call mat_mult(N,Mat,work1,work2)
c
c		====================================================================
c		Normalize new vector
c		====================================================================
c
		norm = 0.d0
		do 100 i = 1,N
			norm = norm + work2(i)*work2(i)
100		continue
		norm = sqrt(norm)
		fact = 1.d0/norm
		do 200 i = 1,N
			work2(i) = work2(i)*fact
200		continue
c
c		====================================================================
c		Check for convergence
c		====================================================================
c
		eig = norm
		if(dabs(eig-eig0) .le. tol) goto 500
c
c		====================================================================
c		Prepare for next iteration
c		====================================================================
c
		do 300 i = 1,N
			work1(i) = work2(i)
300		continue
		eig0 = eig
c
400	continue
c
500	continue
c
c===================================================================================
c	Warns if the procedure did not converge
c===================================================================================
c
c	if(iter.ge.iter_max) then
c		write(6,*) ' '
c		write(6,*) 'WARNING: PowerMethod may not have converged'
c		write(6,*) ' '
c	endif
c
	return
	end
c===================================================================================
c===================================================================================
c
c	Random_vect.f
c
c	This subroutine fills an array with random values in [0,1[ 
c
c===================================================================================
c===================================================================================
c
	subroutine random_vect(N,iseed,vect)
c
c===================================================================================
c	Input:
c		
c		- N:		size of the vector
c		- iseed:	seed for random number generator
c
c	Output:
c		
c		- vect:		random vector
c
c===================================================================================
c
	integer	i
	integer	N,iseed
c
	real*4	ran2
c
	real*8	x,norm
	real*8	vect(*)
c
c===================================================================================
c	Fill in vector: calls random number generator ran2
c===================================================================================
c
	norm = 0.d0
	do 100 i = 1,N
		x = dble(ran2(iseed))
		vect(i) = x
		norm = norm + x*x
100	continue
c
	norm = 1.d0/norm
c
	do 200 i = 1,N
		vect(i) = vect(i)*norm
200	continue
c
	return
	end
c==========================================================================================
c==========================================================================================
c	Portable Random number generator.
c
c	From Numerical Recipes, W.H.Press, B.P.Flannery
c				S.A.Teukolsky, W.T.Vetterling
c		Cambridge Univ. Press
c
c	Function that returns a uniform random deviate between 0.0 and 1.0.
c	Set idum to any negative value to initialize the sequence.
c==========================================================================================
c==========================================================================================
c
	function ran2(idum)
C	-------------------
C
	integer	m,ia,ic
	real	rm,ran2
c
	parameter (m=714025,ia=1366,ic=150889,rm=1./m)
C
	integer	ir(97), iy,iff,idum,j
	data iff /0/
c
	save iy,ir
C
C	BEGIN.
C
	if(idum.lt.0.or.iff.eq.0) then
		iff=1
		idum=mod(ic-idum,m)
		do 11 j=1,97		! init the shuffle table
			idum=mod(ia*idum+ic,m)
			ir(j)=idum
11		continue
		idum=mod(ia*idum+ic,m)
		iy=idum
	endif
	j=1+(97*iy)/m
	if(j.gt.97.or.j.lt.1) then
		write(6,*) 'j=',j
		write(6,*) 'iy,m :',iy,m
		stop
	endif
	iy=ir(j)
C
C	RETURNED VALUE
C
	ran2=iy*rm
C
	idum=mod(ia*idum+ic,m)
	ir(j)=idum
C
	return
	end
