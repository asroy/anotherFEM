module GMRES_mod
      use KindDefinition_mod, only : DP
      use GlobalVariableForDebug_mod
      implicit none

contains

!=================================== GMRES ===================================
!
! GMRES with right preconditioning
! Solve Ax + b = zero.
! To do this solve A*M(inv)*u + b = 0 (or Abar*u + b = 0)
! and then get x by solving Mx = u
!
!=============================================================================
subroutine gmres_ (dq,A,rhs,ia,iau,ja,ALU,iaFill,iauFill,jaFill,nnodes,nnz,nnzFill,blockSize,tolerance, nrestart,nsearch)
      implicit none
      integer,  intent(in)   :: nnodes, blockSize, nnz, nnzFill, nrestart, nsearch
      real(dp), intent(in)   :: rhs(blockSize,nnodes)
      real(dp), intent(out)  ::  dq(blockSize,nnodes)
      real(dp), intent(in)   :: A(blockSize,blockSize,nnz)
      real(dp), intent(out)  :: ALU(blockSize,blockSize,nnzFill)
      integer,  intent(in)   :: ia(nnodes+1), iau(nnodes), ja(nnz)
      integer,  intent(in)   :: iaFill(nnodes+1), jaFill(nnzFill)
      integer,  intent(inout):: iauFill(nnodes)
      real(dp), intent(in)   :: tolerance

      integer :: NEQ

      real(dp), dimension(:),   allocatable :: AP, X
      real(dp), dimension(:,:), allocatable :: phi, tres
      real(dp), dimension(:,:), allocatable :: F
      real(dp), dimension(:),   allocatable :: G           ! Right-hand side for least squares problem
      real(dp), dimension(:),   allocatable :: C           ! Cosine for Givens rotations
      real(dp), dimension(:),   allocatable :: S           ! Sin for Givens rotations
      real(dp), dimension(:,:), allocatable :: H           ! Hessenberg matrix

      real(dp), external :: SNRM2, SDOT

      integer :: i,j,m,l,ir
      integer :: nx,np
      integer :: lrec
      integer :: idolu,iprecon,ind
      integer :: jit,lit

      real(dp) :: bnorm,resnrm
      real(dp) :: res,tempc,apnorm
      real(dp) :: relres

      allocate(G(nsearch))
      allocate(c(nsearch))
      allocate(s(nsearch))
      allocate(H(nsearch,nsearch))
      allocate(AP(blockSize*nnodes))
      allocate(x(blockSize*nnodes))
      allocate(phi(blockSize,nnodes))
      allocate(tres(blockSize,nnodes))
      allocate(F(blockSize*nnodes,nsearch+1))


      nx = 1
      np = 1

      NEQ = blockSize*nnodes
!
! Initialize u and store in dq
!
      do i = 1,nnodes
        do m = 1,blockSize
          dq(m,i) = 0.
        end do
      end do
!
! First copy dq(blockSize,nnodes) into x(blockSize*nnodes)
!
      call DQX(nnodes,dq,X,blockSize)  ! Set initial guess u0 to dq u0 <- dq (which is zero)
!
      LREC=NEQ
      CALL SCOPY(LREC,X,1,F(1,NX),1)  ! Copy u0 to F(.,1)
!
! Compute r0 = b + Abar*u0
! Result is in AP
!
      idolu = 1
      iprecon = 0 ! Dont apply preconditioner; just fill and do the LU decompostion
      ind = 1
      call FCNEVAL(nnodes,idolu,dq,A,rhs,ia,iau,ja,ALU,iaFill,iauFill,jaFill,AP,ind,phi,iprecon,tres,nnz,nnzFill,blockSize)
!
! bnorm = ||r0||
      BNORM=SNRM2(NEQ,AP,1)
      LIT = 0
      JIT = 0
      RELRES = 1.
      IF (BNORM .le. tolerance) GO TO 2000
      RESNRM=BNORM
      write(6,*)' nrestart=',nrestart,' BNORM=',bnorm
      write(fpd1,*)' nrestart=',nrestart,' BNORM=',bnorm
!
! Now solve Abar*u + r0 = 0
!
      DO l=1,nrestart  ! Loop over cycles (nrestart 1 means no restart)
      RES=RESNRM
      APNORM=RESNRM
      DO 300 J=1,nsearch  ! Loop over dimension of Krylov subspace
      G(J)=RES            ! This is the right-hand side of the least squares problem
                          ! It will start with beta*e1 = ||r0|| but new entries will
                          ! get added and it will continually change as Givens rotations are applied
!
! Normalize v
!
      CALL SSCAL(NEQ,1./APNORM,AP,1)
!
! Copy v into F(*,NP+J) Note that NP = 1 so we will copy into 2nd column and higher
!
      CALL SCOPY(LREC,AP,1,F(1,NP+J),1)
!
! Copy v into X
!
      CALL SCOPY(NEQ,AP,1,X,1)
!
! copy x into dq
!
      call XDQ(nnodes,dq,X,blockSize)
!
! Calculate Abar*v, v is stored in dq, Abar*v will be stored in AP
!
      idolu = 0
      iprecon = 1
      ind = 0
      call FCNEVAL(nnodes,idolu,dq,A,rhs,ia,iau,ja,ALU,iaFill,iauFill,jaFill,AP,ind,phi,iprecon,tres,nnz,nnzFill,blockSize)
!
! Orthoganolize against previous vectors using modified Gram-Schmidt
! Here, AP keeps getting modified to subtract off components of AP
! that point in the direction of all previous vectors in subspace
!
      DO I=1,J
        CALL SCOPY(LREC,F(1,NP+I),1,X,1)
        H(I,J)=SDOT(NEQ,X,1,AP,1)
        CALL SAXPY(NEQ,-H(I,J),X,1,AP,1)
      end do
!
! Apply Givens rotations to H for all the previous cos and sin
!
      DO 200 I=1,J-1
        TEMPC=C(I)*H(I,J)+S(I)*H(I+1,J)
        H(I+1,J)=-S(I)*H(I,J)+C(I)*H(I+1,J)
        H(I,J)=TEMPC
 200  CONTINUE
!
! Get current cos and sin and apply Givens to this row including RHS
!
      APNORM=SNRM2(NEQ,AP,1)
      TEMPC=APNORM
      CALL SROTG(H(J,J),TEMPC,C(J),S(J))
      RES=-S(J)*G(J)
      G(J)=C(J)*G(J)
      RESNRM=ABS(RES)
      RELRES=RESNRM/BNORM
        write(6,'("j = ",i5," resnrm = ",e19.10," resres = ",e19.10)')j,resnrm,relres
        write(fpd1,'("j = ",i5," resnrm = ",e19.10," resres = ",e19.10)')j,resnrm,relres
      JIT = J
      IF (RELRES .LE. tolerance .OR. J .EQ. nsearch) GO TO 550
 300  CONTINUE ! End of loop over Krylov subspace
 550  CONTINUE
      IF (ABS(H(J,J)) .EQ. 0.) THEN
      GO TO 2000
      END IF
!
! Solve Hy=g by backsubstitution
! After this, G holds y (recall that x = x0 + Z and Z = V*y)
!
      G(J)=G(J)/H(J,J)
      DO IR=J-1,1,-1
        G(IR)=(G(IR)-SDOT(J-IR,G(IR+1),1,H(IR,IR+1),nsearch))/H(IR,IR)
      end do
!
! We now know y so compute z=Vy and add to x
! u <- u + Vy (this gets put back into F(,.1)
! Recall that the rest of the V's are stored in F(.,2)
!
      CALL SCOPY(LREC,F(1,NX),1,X,1)  ! Copy u0 into X (remember u0 is stored in F(.,1)
      DO I=1,J
        CALL SCOPY(LREC,F(1,NP+I),1,AP,1)
        CALL SAXPY(NEQ,-G(I),AP,1,X,1)
      end do
      CALL SCOPY(LREC,X,1,F(1,NX),1) ! Copy u into F(,.1)
      call XDQ(nnodes,dq,X,blockSize) ! Copy X(.) into dq(.,.) (this puts u into dq)
!
! At this point, X, dq and F(,.1)contain u
! We need to compute x = M(inv)*u
! Recall that we solved A*M(inv)*u + r0 = 0 for u
! We now know u so we need to compute x from Mx = u
! We do this by applying the preconditioner
!
      call QTOQ(nnodes,phi,dq,blockSize) ! copy dq into phi
                                ! (this puts u into phi to use as the right-hand side for the preconditioner)
!
! After qtoq, x, dq, and phi all contain u
!

!
!  calculate M(inv)*phi
!  phi is the solution to preconditioned linear system, M(inv)*phi the real solution,
!
      idolu = 0
      iprecon = 1
      call ILUPRE(nnodes,nnz,idolu,dq,A,ALU,ia,ja,iau,iaFill,iauFill,jaFill,phi,iprecon,blockSize,nnzFill)

!
! After ilupre, dq contains the actual x whereas the phi and X arrays contain the current value of u.
! Also note that F(.,1) also contains the current value of u
!
! Now dq contains x so we can compute Ax + r0 to see if it converged
! To do this don't apply preconditioner. Result comes back in AP
!
      idolu = 0
      iprecon = 0
      ind = 1
      call FCNEVAL(nnodes,idolu,dq,A,rhs,ia,iau,ja,ALU,iaFill,iauFill,jaFill,AP,ind,phi,iprecon,tres,nnz,nnzFill,blockSize)
!

      RESNRM=SNRM2(NEQ,AP,1)
      write(6,*)"Residual at end of restart cycle",L,RESNRM
      write(fpd1,*)"Residual at end of restart cycle",L,RESNRM
      RELRES=RESNRM/BNORM
      LIT = L
      IF (RELRES .LE. tolerance) GO TO 2000
      end do ! End of restart loop
 2000 CONTINUE ! If we reached tolerance kick out to here
!
! At this point, F(,.1) has u in it. It does not have x in it.
! It doesn't do any harm to copy into x but we don't need it
! so comment it out
!     CALL SCOPY(LREC,F(1,NX),1,X,1) ! Copy F to X
      write(6,105)lit,jit,relres
      write(fpd1,105)lit,jit,relres
  105 format(1h 'gmres iters=',i6,' srch directions=',i6,' rms=',e14.7)

      deallocate( G, H, S, C )
      deallocate( AP, X, phi, tres, F )

      return
end subroutine gmres_
!=================================== FCNEVAL =================================
!
! Compute Ax or Ax + b
! Actually computing A*M(inv)*dq + B
! Here also, M is stored in ALU and we solve M(inv)*dq by M*t = dq
!
!=============================================================================
subroutine FCNEVAL(nnodes,idolu,dq,A,B,ia,iau,ja,ALU,iaFill,iauFill,jaFill,AP,ind,phi,iprecon,tres,nnz,nnzFill,blockSize)
      implicit none
      integer,  intent(in)    :: nnodes,idolu,ind,iprecon,nnz,nnzFill,blockSize
      real(dp), intent(inout) :: dq(blockSize,nnodes)
      real(dp), intent(in)    :: A(blockSize,blockSize,nnz)
      real(dp), intent(inout) :: ALU(blockSize,blockSize,nnzFill)
      real(dp), intent(out)   :: phi(blockSize,nnodes)
      real(dp), intent(in)    :: B(blockSize,nnodes)
      real(dp), intent(out)   :: tres(blockSize,nnodes),AP(blockSize*nnodes)
      integer,  intent(in)    :: ia(nnodes+1),iau(nnodes),ja(nnz)
      integer,  intent(in)    :: iaFill(nnodes+1),jaFill(nnzFill)
      integer,  intent(inout) :: iauFill(nnodes)

      integer :: i,j,k,m
      integer :: istart,iend,icol,ii,index
      real(dp) :: rnd
!
!
! First form Ax
! Put preconditioner here.
! Note, in the preconditioner, we solve [A]{dq} = {phi}
! so first, lets go ahead and copy dq into phi just because
! it might make it a little clearer
!
      call QTOQ(nnodes,phi,dq,blockSize) ! phi <- dq
!
! After ilupre, dq holds M(inv)*dq where dq is the original
! version that came into the routine.
! We will then need to multiply A*dq(current) = A*M(inv)*dq(original) and return
!
      call ILUPRE(nnodes,nnz,idolu,dq,A,ALU,ia,ja,iau,iaFill,iauFill,jaFill,phi,iprecon,blockSize,nnzFill)
!
! Compute A*dq(current) = A*M(inv)*dq(original)
!
      do i = 1,nnodes
        do m = 1,blockSize
          tres(m,i) = 0.
        end do
      end do

      do i = 1,nnodes
        istart = ia(i) ! Row of matrix starts here
        iend   = ia(i+1) - 1 ! Row of matrix ends here
        do ii = istart,iend
          icol = ja(ii)
            do j=1,blockSize
              do k=1,blockSize
                tres(j,i) = tres(j,i) + A(j,k,ii)*dq(k,icol)
              end do
            end do
        end do
      end do
!
! If ind = 0, we are calculating Ax so don't add the "residual
! If ind = 1, we are calculating Ax + b so do add the "residual"
!
      rnd = float(ind)
      do i = 1,nnodes
        do m = 1,blockSize
          tres(m,i) = tres(m,i) + rnd*B(m,i)
        end do
      end do
!
! Copy tres into AP
!
      do i=1,nnodes
        index = blockSize*i - (blockSize - 1)
        do m = 1,blockSize
          AP(index) = tres(m,i)
          index = index + 1
        end do
      end do
!
      return

end subroutine FCNEVAL
!=================================== QTOQ =====================================
!
! Simply copy X into dq where both arrays are blockSize x nnodes
!
!=============================================================================
subroutine QTOQ(nnodes,dq,X,blockSize)
      implicit none
      integer,  intent(in)  :: nnodes,blockSize
      real(dp), intent(out) :: dq(blockSize,nnodes)
      real(dp), intent(in)  :: X(blockSize,nnodes)
!
      integer i,m

      do i = 1,nnodes
!
        do m = 1,blockSize
          dq(m,i) = X(m,i)
        end do
!
      end do
      return
end subroutine QTOQ
!================================ ILUPRE ====================================
!
! Preconditioning using ILU(k)
! we solve [A]{dq} = {phi} to give us a "new" value stored in dq
!
!============================================================================
subroutine ILUPRE(nnodes,nnz,idolu,dq,A,ALU,ia,ja,iau,iaFill,iauFill,jaFill,phi,iprecon,blockSize,nnzFill)
      implicit none
      integer,  intent(in)    :: nnodes,nnz,idolu,iprecon,blockSize,nnzFill
      integer,  intent(in)    :: ia(nnodes+1), iau(nnodes), ja(nnz)
      integer,  intent(in)    :: iaFill(nnodes+1),jaFill(nnzFill)
      integer,  intent(inout) :: iauFill(nnodes)
      real(dp), intent(inout) :: dq(blockSize,nnodes)
      real(dp), intent(inout) :: alu(blockSize,blockSize,nnzFill)
      real(dp), intent(in)    :: a(blockSize,blockSize,nnz)
      real(dp), intent(in)    :: phi(blockSize,nnodes)

      integer :: i,j,k
      integer :: icode
      integer :: m
      integer :: jstart,jend,kstart,kend,jcol,jcolnew,kcol,ii,jj


!
! First we must do the ILU decompositions so lets assemble the A
! matrix
!
      if(idolu.eq.1)then
!
! Now copy the A matrix into the new matrix (ALU)
!
         do i = 1,nnzFill
           do j = 1,blockSize
             do k = 1,blockSize
               ALU(j,k,i) = 0.
             end do
           end do
         end do

!
! Loop over rows and fill the new matrix values with the original ones
!
         do i = 1,nnodes
           jstart = ia(i)        ! Start of row in original matrix
           jend   = ia(i+1) - 1  ! End of row in original matrix
           kstart = iaFill(i)           ! Start of row in new matrix
           kend   = iaFill(i+1) - 1     ! End of row in new matrix
!
! Now loop over the columns in the old row
! For each element of A, find where in the new row it belongs
!
           do j = jstart,jend ! Original row
             jcol = ja(j)  ! Original column
             jcolnew = jcol
             do k = kstart,kend ! New row
               kcol = jaFill(k) ! New column
! If new column matches the old column, put this value of A into Afill
               if(kcol.eq.jcolnew)then
                 do ii = 1,blockSize
                   do jj = 1,blockSize
                     ALU(ii,jj,k) = A(ii,jj,j)
                   end do
                 end do
                 goto 5464 ! If we've found it, kick out
               end if
             end do ! k loop
 5464   continue
           end do ! jloop
         end do ! loop over rows to fill ALU
!
! For no preconditioning set ALU to identity
!       do 999 i = 1,nnodes
!         ii = iau(i)
!         do m = 1,blockSize
!           ALU(m,m,ii) = 1.
!         end do
! 999  continue
!
! Now do LU decomposition
!
        call GBLKILU(nnodes,nnzFill,jaFill,iaFill,ALU,iauFill,icode,blockSize)
        if(icode.ne.0)then
         write(6,333)
  333    format(1h ,'Error in BLKILU')
         stop
        end if
      end if
!
! If we want to apply the preconditioner, do the back substitution
!
      if(iprecon.eq.1)then

      dq = 0.
!
! Now do backsubstitution to get dq
!
      call GBLKSOL(nnodes,nnzFill, phi, dq, alu, jaFill, iaFill, iauFill, blockSize)
      end if
!
      return
end subroutine ilupre
!====================== GBLKILU ====================================
!
! ILU(0) for general sized matrices
! Based on SAAD's book
! neq - number of equations (nnodes)
! blockSize - block size
!
!===================================================================
subroutine GBLKILU(neq,nnz,ja,ia,alu,iau,icode,blockSize)
        integer,  intent(in)    :: neq,nnz,blockSize
        integer,  intent(out)   :: icode
        integer,  intent(in)    :: ja(nnz),ia(neq+1)
        integer,  intent(inout) :: iau(neq)
        real(dp), intent(out)   :: alu(blockSize,blockSize,nnz)

        integer, dimension(:), allocatable :: iw
        real(dp) :: tmat(blockSize,blockSize)

        real(dp) :: sum
        integer :: segment,row,column
        integer :: i,j,k,m,n,kn
        integer :: j1,j2
        integer :: jrow
        integer :: L1,L2
        integer :: jw,jj

        real(dp) :: tmp


        allocate( iw(neq) )
!-----------------------------------------------------------
!
! initialize work vector to zero's
!
        do i=1, neq
           iw(i) = 0
        enddo
!-----------------------------------------------------------
!-------------- MAIN LOOP ----------------------------------
!-----------------------------------------------------------
        do 500 k = 1, neq

!
!------------------------ k = row number -------------------
!
           j1 = ia(k)      ! Index in alu marking the beginning of row k
           j2 = ia(k+1)-1  ! index in alu marking the End of row k
           do j=j1, j2
              iw(ja(j)) = j
           enddo

           j=j1            ! Set j to beginning index of row k
 150       jrow = ja(j)    ! Column
!
!----------------------- Exit if diagonal element is reached.
!
           if (jrow .ge. k) goto 200
!
!----------------------- Compute the multiplier for jrow.
!
! Initialize tmat
!
           do L1 = 1,blockSize
             do L2 = 1,blockSize
               tmat(L1,L2) = 0.
             end do
           end do
!
! Multiply alu(j) by alu(iau(jrow))
!
           do L1 = 1,blockSize
              do L2 = 1,blockSize
                do m = 1,blockSize
                  tmat(L1,L2) =  tmat(L1,L2) + alu(L1,m,j)*alu(m,L2,iau(jrow))
                end do
              enddo
           enddo
!
! Replace alu(j) with tmat
!
           do L1 = 1,blockSize
             do m = 1,blockSize
              alu(L1,m,j) = tmat(L1,m)
           end do
           enddo

!
!-----------------------perform linear combination
!
           do jj = iau(jrow)+1, ia(jrow+1)-1  ! Loop over columns past the diagonal
              jw = iw(ja(jj))
              if (jw .ne. 0)then              ! If there is an element here
                 do L1=1,blockSize
                    do L2=1,blockSize
                      do m = 1,blockSize
                        alu(L1,L2,jw) = alu(L1,L2,jw) - tmat(L1,m)*alu(m,L2,jj)
                      end do
                    enddo
                 enddo
              endif
           enddo


           j=j+1
           if (j .le. j2) goto 150
!
!----------------------- define and store lower part
!
 200       iau(k) = j


           if (jrow .ne. k) goto 600

           iau(k) = j
!
! Compute inverse of alu(j)
! First do lu decomposition
!
        do m = 1,blockSize
          do n = 1,blockSize
            tmat(m,n) = alu(m,n,j)
          end do
        end do

        do n = 2,blockSize
          tmat(1,n) = tmat(1,n)/tmat(1,1)
        end do
!
! Now loop over each "segment" starting with row 2
!
          do segment = 2,blockSize
!
! First do the column corresponding to the segment
!
            do row = segment,blockSize
              sum = 0.
              do n = 1,segment-1
                sum = sum + tmat(row,n)*tmat(n,segment)
              end do
              tmat(row,segment) = tmat(row,segment) - sum
            end do
!
! Now do the upper triangular portion (the row)
!
            do column = segment+1,blockSize
              sum = 0.
              do n = 1,segment-1
                sum = sum + tmat(segment,n)*tmat(n,column)
              end do
              tmat(segment,column) = (tmat(segment,column) - sum)/tmat(segment,segment)
            end do

          end do
!
! The lu decomposition is now in tmat
! Get the inverse of alu(j) by repeated backsubstitution (once for each column of identity matrix)
! First set alu(j) to the identity matrix
!
           do m = 1,blockSize
             do n = 1,blockSize
               alu(m,n,j) = 0.
             end do
             alu(m,m,j) = 1.
           end do
!
! Now solve for the entries in the inverse of alu(j)
! Do this one column at a time
!
        do m = 1,blockSize ! This will be the column of alu
        alu(1,m,j) = alu(1,m,j)/tmat(1,1)
        do n = 2,blockSize
          sum = alu(n,m,j)
          do kn = 1,n-1
            sum = sum - tmat(n,kn)*alu(kn,m,j)
          end do
          alu(n,m,j) = sum/tmat(n,n)
        end do
!
! Backward (note that B(neqn,m) is already obtained so we start at the next row up)
!
        do n = blockSize-1,1,-1
          sum = alu(n,m,j)
          do kn = blockSize,n+1,-1
            alu(n,m,j) = alu(n,m,j) - tmat(n,kn)*alu(kn,m,j)
          end do
        end do
      end do
!
!
!----------------------- zero all entries of iw.
!
          do i = j1, j2
            iw(ja(i)) = 0
          enddo
500   continue

      deallocate( iw )

      icode = 0
      return

600   icode = k
      write(*,*)' zero pvt ',k
      return
end subroutine GBLKILU
!=========================== GBLKSOL ===================================
!
! L-U solve for ILU preconditioning processed by ILU.
!
!=======================================================================
! parameters
! n     = dimension of problem
! y     = (input) right hand side (can be the same as x)
! x     = (output) solution to L U x = y.
! alu, jalu, ial
!       = matrices L and U stored in sparse format
! iau   = integer array pointing to diagonal entries in L-U matrices
!-----------------------------------------------------------s
subroutine GBLKSOL(n,nnz, y, x, alu, jalu, ial, iau, blockSize)
      implicit none
      integer :: blockSize
      integer :: n,nnz
      integer :: jalu(nnz), ial(n+1), iau(n)
      real(dp) ::  alu(blockSize,blockSize,nnz), x(blockSize,n), y(blockSize,n)

      real(dp) :: t(blockSize)
      integer :: i,k,m
      integer :: m1,m2
      integer :: l1
!
      do i = 1, n
          m1 = ial(i)
          m2 = iau(i) -1
          do l1=1,blockSize
            x(l1,i) = y(l1,i)
            do k=m1, m2
              do m = 1,blockSize
                  x(l1,i) = x(l1,i) - alu(l1,m,k)*x(m,jalu(k))
              end do
            enddo
          enddo
      enddo
!
      do i = n, 1, -1
          m1 = iau(i) + 1
          m2 = ial(i+1) - 1
          do l1=1,blockSize
            do k=m1, m2
              do m = 1,blockSize
                x(l1,i) = x(l1,i) - alu(l1,m,k)*x(m,jalu(k))
              end do
            enddo
          enddo

          do l1=1,blockSize
            t(l1) = 0.
            do m = 1,blockSize
            t(l1) = t(l1) + alu(l1,m,iau(i))*x(m,i)
            end do
          enddo

          do m = 1,blockSize
            x(m,i) = t(m)
          end do
      enddo
      return
end subroutine GBLKSOL
!============================ getIAJAFILL ==============================
!
! For ILU(k) get the new ia, ja, and iau arrays to accomodate the fill
!
!=======================================================================
subroutine getIAJAFill(nnodes, nnz, ia, iau, ja, nnzFill, iaFill, iauFill, jaFill, fillLevel)

      integer, intent(in)                             :: fillLevel
      integer, intent(in)                             :: nnodes, nnz
      integer, intent(in)                             :: ia(nnodes+1), iau(nnodes), ja(nnz)
      integer, intent(out)                            :: nnzFill
      integer, dimension(:), allocatable, intent(out) :: iaFill, iauFill, jaFill

      type :: compressedBigA
        integer :: nnzrow
        integer, dimension(:), allocatable :: jarow
        integer, dimension(:), allocatable :: entry
      end type compressedBigA
!
      type(compressedBigA), dimension(:), allocatable :: bigAC

      integer :: i,j,k
      integer :: bwidth
      integer :: jstart,jend
      integer :: kstart,kend
      integer :: jcol
      integer :: nonZeroOnRow
      integer :: icount,jcount
      integer, dimension(:), allocatable :: bigArow

      integer :: nonzeros,nnzk,iTotalMemoryAllocated
      integer :: searchStart ! Starting location for search
      real(dp) :: ratio
!     external integer bigAvalue

!
! If ILU(0) just copy ia, ja, and iau into iaFill, jaFill, and iauFill and return
! Also set nnzFill to nnz
!
      if(fillLevel.eq.0)then
        nnzFill = nnz
        allocate(iaFill (nnodes+1))      ! The new ia
        allocate(iauFill(nnodes))        ! The new iau
        allocate(jaFill(nnzFill)) ! Allocate enough space for adding a fill for every node
        do i = 1,nnodes
          iaFill(i) = ia(i)
          iauFill(i) = iau(i)
        end do
        iaFill(nnodes+1) = ia(nnodes+1)
        do i = 1,nnz
          jaFill(i) = ja(i)
        end do
        return
      end if
!
! Get the bandwidth of the original matrix
!
      bwidth = bandwidth( ia, ja, nnodes, nnz )
!
! Now lets compute the fill pattern
!
! Allocate memory
!
      allocate(bigArow(nnodes)) ! Used for working 1 row at a time
      allocate(bigAC(nnodes))   ! Stores only the non-infinite entries
!
! Put zeros in the locations where there are elements present in the original matrix
! First fill bigAC for the first row. We will fill the other rows as we go
!
      nonzeros = ia(2) - ia(1)
      bigAC(1)%nnzrow = nonzeros
      allocate(bigAC(1)%entry(nonzeros))
      allocate(bigAC(1)%jarow(nonzeros))
      jstart = ia(1)
      jend   = ia(2) - 1
      jcount = 0
      do j = jstart,jend
        jcount = jcount + 1
        bigAC(1)%entry(jcount) = 0
        bigAC(1)%jarow(jcount) = ja(j)
      end do

      write(6,*)"Finished setting non-zero elements to zero"

      write(6,*)" Starting bigA "
!
! First set bigArow to a big number
!
      do i = 1,nnodes
        bigArow(i) = 500000000
      end do
!
! Fill the remainder of rows
!
      nodeLoop: do i = 2,nnodes
!
        if((i/100)*100.eq.i)write(6,*)"Working row ",i
!
! For row, put a zero where there are non-zero entries
!
        jstart = ia(i)
        jend = ia(i+1) - 1
        do j = jstart,jend
          jcol = ja(j)
          bigArow(jcol) = 0
        end do
!
! Determine starting and ending of row based on band structure
!
        kstart = max(1,i-bwidth)      ! First non-zero column in this row of bigB
        kend   = min(nnodes,i+bwidth)
!       do k = 1,i-1 ! Loops over columns 1,i-1 of row i
        do k = kstart,i-1 ! Loops over columns 1,i-1 of row i
          if(bigArow(k).le.fillLevel)then
            nnzk = bigAC(k)%nnzrow
            searchStart = 1
!           do j = k+1,nnodes
            do j = k+1,kend
              bigArow(j) = min(bigArow(j), bigArow(k) + bigAvalue(bigAC, k,j,nnzk,nnodes,bwidth,searchStart) + 1)
            end do
          end if
        end do
!
! Now find entries on current row that are less than infinity
!
        jcount = 0
!       do j = 1,nnodes
        do j = kstart,kend
          if(bigArow(j).le.fillLevel)jcount = jcount + 1
        end do
        bigAC(i)%nnzrow = jcount
        allocate(bigAC(i)%entry(jcount))
        allocate(bigAC(i)%jarow(jcount))
        jcount = 0
!       do j = 1,nnodes
        do j = kstart,kend
          if(bigArow(j).le.fillLevel)then
            jcount = jcount + 1
            bigAC(i)%jarow(jcount) = j
            bigAC(i)%entry(jcount) = bigArow(j)
          end if
        end do
!
! Reset bigArow to infinity before going to next row
!
        do j = 1,bigAC(i)%nnzrow
          jcol = bigAC(i)%jarow(j)
          bigArow(jcol) = 500000000
        end do

      end do nodeLoop

      write(6,*)"Finished computing bigA"
!
! Allocate iew arrays for compressed row storage
!
        allocate(iaFill(nnodes+1))       ! The new ia
        allocate(iauFill(nnodes))        ! The new iau
!
! Figure out how many non-zeros there are
! Remember that at this point, bigAC contains all the entries less than infinity
! so we need to only include the ones that are less than or equal to our fill level
!
        nnzFill = 0
        iaFill(1) = 1
        do i = 1,nnodes
          nonZeroOnRow = bigAC(i)%nnzrow
          jcount = 0 ! Use to determine how many on this row have an entry < fillLevel
          do j = 1,nonZeroOnRow
            if(bigAC(i)%entry(j).le.fillLevel)then
              jcount = jcount + 1
              nnzFill = nnzFill + 1
            end if
          end do
          iaFill(i+1) = iaFill(i) + jcount
        end do
!
! Allocate jaFill
!
        allocate(jaFill(nnzFill)) ! Allocate enough space for adding a fill for every node
        write(6,*)"nnzFill = ",nnzFill
!
! Now fill jaFill
!
        icount = 0
        iTotalMemoryAllocated = 0
        do i = 1,nnodes
          nonZeroOnRow = bigAC(i)%nnzrow
          iTotalMemoryAllocated = iTotalMemoryAllocated + nonZeroOnRow
          do j = 1,nonZeroOnRow
            if(bigAC(i)%entry(j).le.fillLevel)then
              icount = icount + 1
              jaFill(icount) = bigAC(i)%jarow(j)
              if(bigAC(i)%jarow(j).eq.i)iauFill(i) = icount
            end if
          end do
        end do

        write(6,'("Finished data structures for compressed row storage")')
        write(6,'("Memory allocated for bigAC = ",i10)')iTotalMemoryAllocated
        write(6,'("Memory the old way = ",i10)')nnodes*(2*bwidth + 1)
        ratio = real(nnodes*(2*bwidth + 1))/real(iTotalMemoryAllocated)
        write(6,'("Ratio = ",e19.10)')ratio
        write(6,'("nnzFill = ",i7)')nnzFill

!
! Free allocated memory
!
      do i = 1,nnodes
        deallocate(bigAC(i)%entry)
        deallocate(bigAC(i)%jarow)
      end do
      deallocate(bigArow) ! Used for working 1 row at a time
      deallocate(bigAC)   ! Stores only the non-infinite entries

      return


      contains

      !============================ BIGAVALUE ============================
      !
      ! Given a row and column in the big A matrix, determine the entry
      ! Assumes the entry is infinity unless a smaller one is actually found
      !
      !===================================================================
      integer function bigAvalue(bigAC, k,j,nnzk,nnodes,bwidth,searchStart)
            implicit none
            integer j,k,nnzk
            integer nnodes
            integer bwidth
            integer :: searchStart
            type(compressedBigA) :: bigAC(nnodes)

            integer :: i, L
            integer kstart,kend
            integer :: middle,half
            integer :: Lstart,Lend,lcount

            bigAvalue = 500000000
            kstart = max(1,k-bwidth)      ! First non-zero column in this row of bigB
            kend   = min(nnodes,k+bwidth)
            if(j.lt.kstart.or.j.gt.kend)return

      !     if(searchStart.ne.-1)then
            if(searchStart.eq.1)then
              middle = nnzk/2
              if(bigAC(k)%jarow(middle).ge.j)then
                half = middle/2
                if(bigAC(k)%jarow(half).ge.j)then
                  Lstart = 1
                  Lend = half
                else
                  Lstart = half
                  Lend   = middle
                end if
              else
                half = (middle + nnzk)/2
                if(bigAC(k)%jarow(half).ge.j)then
                  Lstart = middle
                  Lend   = half
                else
                  Lstart = half
                  Lend   = nnzk
                end if
              end if
            else
              Lstart = searchStart
              Lend   = nnzk
            end if

            lcount = 0
            do L = Lstart,Lend
              lcount = lcount + 1
              if(bigAC(k)%jarow(L).gt.j)then
                searchStart = L
                exit
              end if
              if(bigAC(k)%jarow(L).eq.j)then
                bigAvalue = bigAC(k)%entry(L)
                searchStart = L
      !     write(6,*)"function found at lcount = ",lcount," nnzk = ",nnzk
      !     write(6,*)" middle = ",middle," Lstart Lend = ",Lstart,Lend
                exit
              end if
            end do
            return
      end function bigAvalue

end subroutine getIAJAFill

!===================================================================
!
! Finds the maximum bandwidth
!
!===================================================================
      integer function bandwidth(ia, ja, nnodes, nnz)
      integer, intent(in) :: nnodes, nnz
      integer, intent(in) :: ia(nnodes+1), ja(nnz)

      integer i,jstart,jend,lowerBand,upperBand

      bandwidth = 0
      do i = 1,nnodes
        jstart = ia(i)
        jend   = ia(i+1) - 1
        lowerBand  = i - ja(jstart)
        upperBand = ja(jend) - i
        bandwidth = max(lowerBand,upperBand,bandwidth)
      end do

      end function bandwidth

!=================================== DQX =====================================
!
! Copy dq into the 'X' array for GMRES
! The reason we have to go to this trouble is because we have to renumber
!
!=============================================================================
      subroutine DQX(nnodes,dq,X,blockSize)
!
      integer nnodes,blockSize
      integer :: i,m
      integer :: index
      real(dp) :: dq(blockSize,nnodes),X(nnodes*blockSize)
      do i = 1,nnodes
        index = blockSize*i - (blockSize - 1)
        do m = 1,blockSize
          X(index) = dq(m,i)
          index = index + 1
        end do
!
      end do
      return
      end subroutine dqx
!
!=================================== XDQ =====================================
!
! Copy 'X' back into dq
! The reason we have to go to this trouble is because we have to renumber
!
!=============================================================================
      subroutine XDQ(nnodes,dq,X,blockSize)
      integer nnodes,blockSize
      real(dp) :: dq(blockSize,nnodes),X(nnodes*blockSize)
!
      integer :: i,m
      integer :: index
      do i = 1,nnodes
        index = blockSize*i - (blockSize - 1)
!
        do m = 1,blockSize
          dq(m,i) = X(index)
          index = index + 1
        end do
!
      end do
      return
      end subroutine XDQ
end module
