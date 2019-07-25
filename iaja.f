! Copyright (c) 2017, 
! Eyal Shalev (eyal@gsi.gov.il)
! Vladimir Lyakhovsky
! All rights reserved to Geological Survey of Israel (GSI) 
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!    * Redistributions of source code must retain the above copyright
!      notice, this list of conditions and the following disclaimer.
!    * Redistributions in binary form must reproduce the above copyright
!      notice, this list of conditions and the following disclaimer in the
!      documentation and/or other materials provided with the distribution.
!    * Neither the name of Harel Levin or Gal Oren, nor the
!      names of its contributors may be used to endorse or promote products
!      derived from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
! ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
! WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL Harel Levin & Gal Oren BE LIABLE FOR ANY
! DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
! (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
! ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

cc for symmetric 	LCNTL(1)=.false. to save lower half matrix
cc for asymmetric 	LCNTL(1)=.false. to save full matrix
	program arrays
      implicit double precision (a-h,o-z)

      parameter(numnp=509,numel=2119,nen=4,ndm=3,ndf=1,
     $     neq=numnp*ndf,nquad=1,nsdm=4,nqdm=3,nst=nen*ndf,
     $     maxa=numel*100)

      dimension x(ndm,numnp),ix(nen,numel),id(ndf,numnp),
     $     u(ndf,numnp),f(ndf,numnp),a(maxa),ia(maxa),
     $     ja(maxa),AELT(numel*10),s(ndf*nen,ndf*nen,numel),
     $     IW(3*numnp),INFO(10)
      INTEGER ELTPTR(numel+1),ELTVAR(numel*nen),LIRN,LA,LP
      LOGICAL LCNTL(10)
	
        OPEN(10,FILE='coordinates.dat')
        OPEN(11,FILE='elements.dat')

	read(10,*)nn
	if (nn .ne. numnp) then
	   write(*,*)'incorrect number of nodes'
	   stop
	 end if

	do i=1,numnp
	   read(10,*)x(1,i),x(2,i),x(3,i)
	end do
	read(11,*)ne
	if (ne .ne. numel) then
	   write(*,*)'incorrect number of elements'
	   stop
	 end if
	DO I=1,numel
	   read(11,*)ix(1,I),ix(2,I),ix(3,I),ix(4,I)
	END DO

	close(10)
	close(11)

	LCNTL(1)=.false.
	LCNTL(2)=.true.
	LCNTL(3)=.false.
	LCNTL(4)=.false.
	k=0
	do i=1,numel
	   do j=1,nen
	      K=k+1
	      ELTVAR(k)=ix(j,i)
	   end do
	end do
	do i=1,numel
	   ELTPTR(i)=i*4-3
	end do
	k=1
	do n=1,numel*10
	   AELT(k)=0.5
	end do
	LIRN=100*numel
	LA=1
	LP=1

	call MC57AD(LCNTL,numnp,numel,ELTVAR,ELTPTR,AELT,N,
     $    LIRN,ja,ia,LA,A,IW,LP,INFO)

	write(*,*)info
        OPEN(21,FILE='ia.in')
        OPEN(22,FILE='ja.in')

	write(*,*)'size of A is',ia(numnp+1)-1
	do i=1,numnp+1
	   write(21,*)ia(i)
	end do
	do i=1,ia(numnp+1)-1
	   write(22,*)ja(i)
	end do
	close(21)
	close(22)


      stop
      end
C COPYRIGHT (c) 1999 Council for the Central Laboratory
*                    of the Research Councils
C Original date 30 November 1999
C +++ ++++++++++++++++++++++++++++++++++++++++++++++
C     Subroutine for assembling elements
C     Jennifer A. Scott
C +++ ++++++++++++++++++++++++++++++++++++++++++++++

C 12th July 2004 Version 1.0.0. Version numbering added.
C 25 August 2009 Version 1.0.1. Minor change to specification sheet.
C 20 April  2010 Version 1.1.0. Bug fixed so that IW(J) = J
C    on return if INFO(3) = N (this was not the case if
C    LCNTL(3) = .TRUE.)

      SUBROUTINE MC57AD(LCNTL,NMAX,NELT,ELTVAR,ELTPTR,AELT,N,LIRN,IRN,
     +                  IP,LA,A,IW,LP,INFO)

C Subroutine for assembling element matrices.
C The resulting matrix has a symmetric sparsity pattern but may have
C unsymmetric values. Optionally assembles only the sparsity pattern.
C
C Argument list (* indicates may be changed by rouitne)

C LCNTL : LOGICAL array of controlling parameters of length 10. The user
C         must set LCNTL(1) to LCNTL(4).
C LCNTL(1) : must be TRUE if symmetric and only lower
C            triangle wanted.
C LCNTL(2) : should be set to TRUE if user would like the entries
C            within each column to be ordered by rows.
C LCNTL(3) : should be set to TRUE if user wants null rows removed.
C LCNTL(4) : should be set to FALSE if only sparsity pattern wanted.
C         The other entries of LCNTL are not currently accessed.
C  NMAX:  INTEGER. On entry, must be at least as large as the
C         largest integer used to index
C         a variable.
C NELT :  INTEGER. Must be set to number of elements.
C ELTVAR :INTEGER array. On entry, leading
C         entries must hold element index lists. Unchanged on exit.
C ELTPTR :INTEGER array length NELT+1
C         On entry, must hold pointers for ELTVAR.
C AELT :  REAL (DP) array.  LCNTL(4) = TRUE, on entry,
C         element entries must be in AELT, with the entries for
C         each element stored column-by-column. If LCNTL(1) = TRUE ,
C         only entries in lower triangular part of element matrix
C         should be stored in packed form.
C         Not accessed if LCNTL(4) = FALSE.
C *N   :  INTEGER. Not set on entry.
C         On exit, holds the order of the assembled matrix
C LIRN  : INTEGER. length of array IRN. Note that, because of
C         duplicates, must be larger than number of entries in final
C         assembled matrix.
C *IRN :  INTEGER array length LIRN. Not set on entry. On exit
C         holds row indices of entries in assembled matrix, ordered
C         by columns. If SYM = .TRUE., only lower triangular part
C         is stored.
C *IP  :  INTEGER array length at least NMAX+1.
C         Not set on entry. On exit, first N+1 entries hold col.
C         pointers for assembled matrix.
C LA  :   INTEGER. Length of array A. LA must be set to 1 if
C         only pattern of A required.
C *A   :  REAL (DP) array length LA.
C         Not set on entry. If LCNTL(4) = TRUE,
C         on exit holds entries of assembled matrix (ordered by rows).
C         Not accessed if LA = 1.
C *IW   : INTEGER array length at least 2*N.
C         On exit, IW(J) is the original
C         variable index for row (or column) J in assembled matrix.
C LP    : INTEGER. Output stream for error messages.
C *INFO : INTEGER array length 10. INFO(1) used as error flag.
C         If INFO(1) = -3 on exit INFO(2) holds minimum value for LIRN
C         If INFO(1) = -4 on exit INFO(2) holds minimum value for LA
C         INFO(3) holds largest integer used to index a variable.
C         INFO(4) number of entries in assembled matrix.
C         INFO(5) holds first element with duplicated indices.
C         INFO(6) holds number of null rows/cols (removed if
C         LCNTL(3) = TRUE).
C         INFO(7) holds smallest integer used to index a variable.

C     .. Scalar Arguments ..
      INTEGER LA,LIRN,LP,N,NELT,NMAX
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(LA),AELT(*)
      INTEGER ELTPTR(NELT+1),ELTVAR(*),INFO(10),IP(NMAX+1),IRN(LIRN),
     +        IW(2*NMAX)
      LOGICAL LCNTL(10)
C     ..
C     .. Local Scalars ..
      INTEGER I,ICOL,IELT,IROW,IVAR,J,JMAX,JMIN,JSTOP,JSTRT,JVAR,K,L,
     +        NUMD,NUME,NVAR,NZ,NZJ
C     ..

C     .. External Subroutines ..
      EXTERNAL MC57BD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MAX,MIN
C     ..
      DO 10 I = 1,10
         INFO(I) = 0
   10 CONTINUE
C Check input data.
      IF (NMAX.LT.1 .OR. NELT.LT.1) THEN
         INFO(1) = -1
         IF (LP.GE.0) THEN
            WRITE (LP,FMT=9000) INFO(1)
            IF (NMAX.LT.1) WRITE (LP,FMT=9010) 'NMAX','NMAX',NMAX
            IF (NELT.LT.1) WRITE (LP,FMT=9010) 'NELT','NELT',NELT
         END IF
         GO TO 450
      END IF

C Initialise IP and IW
      DO 20 I = 1,NMAX
         IP(I) = 0
         IW(I) = 0
         IW(I+NMAX) = 0
   20 CONTINUE

C Count up number of entries in each column (duplicates will be
C removed later).
C Loop over the elements.
C IP(J) accumulates no. of entries in col. J (including duplicates).
      IF (LCNTL(1)) THEN
C Symmetric case. Only want lower triangular part.
C Note: we do not assume that the variables are ordered
C within the element variable lists
         JMAX = 1
         JMIN = NMAX
         DO 50 IELT = 1,NELT
            JSTRT = ELTPTR(IELT)
            NVAR = ELTPTR(IELT+1) - ELTPTR(IELT)
            JSTOP = JSTRT + NVAR - 1
            DO 40 J = JSTRT,JSTOP
               JVAR = ELTVAR(J)
C Store largest and smallest integers used to index a variable
               JMAX = MAX(JMAX,JVAR)
               JMIN = MIN(JMIN,JVAR)
C Check JVAR not out of range
               IF (JVAR.LT.1 .OR. JVAR.GT.NMAX) GO TO 40
C Check we don't have duplicates within an element
               IF (IW(JVAR).LT.IELT) THEN
                  IW(JVAR) = IELT
               ELSE
                  GO TO 80
               END IF
               DO 30 I = J,JSTOP
                  IVAR = ELTVAR(I)
                  IF (IVAR.GE.JVAR) THEN
C Entry will be in column JVAR
                     IP(JVAR) = IP(JVAR) + 1
                  ELSE
C Entry will be in column IVAR
                     IP(IVAR) = IP(IVAR) + 1
                  END IF
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
      ELSE
C Need upper and lower parts.
         JMAX = 1
         JMIN = NMAX
         DO 70 IELT = 1,NELT
            JSTRT = ELTPTR(IELT)
            NVAR = ELTPTR(IELT+1) - ELTPTR(IELT)
            JSTOP = JSTRT + NVAR - 1
            DO 60 J = JSTRT,JSTOP
               JVAR = ELTVAR(J)
C Store largest and smallest integers used to index a variable
               JMAX = MAX(JMAX,JVAR)
               JMIN = MIN(JMIN,JVAR)
C Check JVAR not out of range
               IF (JVAR.LT.1 .OR. JVAR.GT.NMAX) GO TO 60
C Check we don't have duplicates within an element
               IF (IW(JVAR).LT.IELT) THEN
                  IW(JVAR) = IELT
               ELSE
                  GO TO 80
               END IF
               IP(JVAR) = IP(JVAR) + NVAR
   60       CONTINUE
   70    CONTINUE
      END IF
      GO TO 90

   80 CONTINUE
      INFO(1) = -5
      INFO(5) = IELT
      IF (LP.GE.0) WRITE (LP,FMT=9000) INFO(1)
      IF (LP.GE.0) WRITE (LP,FMT=9030) JVAR,IELT
      GO TO 450

   90 CONTINUE
C Note that if IW(JVAR).ne.0 (1.le.jvar.le.jmax) then
C JVAR is used to index a variable (and if IW(JVAR)=0, then JVAR not
C used)
      INFO(3) = JMAX
      INFO(7) = JMIN
      IF (JMAX.GT.NMAX .OR. JMIN.LT.1) THEN
         INFO(1) = -2
         IF (LP.GE.0) THEN
            WRITE (LP,FMT=9000) INFO(1)
            IF (JMAX.GT.NMAX) WRITE (LP,FMT=9020) 'NMAX',JMAX
            IF (JMIN.LT.1) WRITE (LP,FMT=9040) JMIN
         END IF
         GO TO 450
      END IF

C Set IP(J) to first entry in column J+1
      IP(1) = IP(1) + 1
      DO 100 J = 2,JMAX
         IP(J) = IP(J) + IP(J-1)
  100 CONTINUE
      IP(JMAX+1) = IP(JMAX)
      INFO(2) = IP(JMAX) - 1

      IF (IP(JMAX)-1.GT.LIRN) THEN
         INFO(1) = -3
         IF (LP.GE.0) WRITE (LP,FMT=9000) INFO(1)
         IF (LP.GE.0) WRITE (LP,FMT=9020) 'LIRN',IP(JMAX) - 1
         GO TO 450
      END IF

      IF (LCNTL(4) .AND. IP(JMAX)-1.GT.LA) THEN
         INFO(1) = -4
         IF (LP.GE.0) WRITE (LP,FMT=9000) INFO(1)
         IF (LP.GE.0) WRITE (LP,FMT=9020) 'LA',IP(JMAX) - 1
         GO TO 450
      END IF

      IF (LCNTL(4)) GO TO 270

C Matrix pattern only.

C Fill-in the matrix entries.
C Loop over the elements.
      IF (LCNTL(1)) THEN
         DO 130 IELT = 1,NELT
            JSTRT = ELTPTR(IELT)
            NVAR = ELTPTR(IELT+1) - ELTPTR(IELT)
            JSTOP = JSTRT + NVAR - 1
            DO 120 J = JSTRT,JSTOP
               JVAR = ELTVAR(J)
               DO 110 I = J,JSTOP
                  IVAR = ELTVAR(I)
                  IF (IVAR.GE.JVAR) THEN
C Entry in lower triangular part
                     IP(JVAR) = IP(JVAR) - 1
                     K = IP(JVAR)
C Set row index
                     IRN(K) = IVAR
                  ELSE
C Entry in upper triangular part
                     IP(IVAR) = IP(IVAR) - 1
                     K = IP(IVAR)
C Set row index
                     IRN(K) = JVAR
                  END IF
  110          CONTINUE
  120       CONTINUE
  130    CONTINUE
      ELSE
         DO 160 IELT = 1,NELT
            JSTRT = ELTPTR(IELT)
            NVAR = ELTPTR(IELT+1) - ELTPTR(IELT)
            JSTOP = JSTRT + NVAR - 1
            DO 150 J = JSTRT,JSTOP
               JVAR = ELTVAR(J)
               DO 140 I = JSTRT,JSTOP
                  IVAR = ELTVAR(I)
                  IP(JVAR) = IP(JVAR) - 1
                  K = IP(JVAR)
C Set row index
                  IRN(K) = IVAR
  140          CONTINUE
  150       CONTINUE
  160    CONTINUE
      END IF

C IP(J) now points to the start of column J.
C We optionally remove null columns (these occur if largest integer
C used to index a variable is greater than order of assembled matrix).

C NUMD is number of null columns (= number of null rows).
      NUMD = 0
      IF (.NOT.LCNTL(3)) GO TO 230

      NZ = IP(JMAX+1) - 1
      DO 180 I = 1,NZ
         K = IRN(I)
C There is an entry in row K
         IW(NMAX+K) = 1
  180 CONTINUE

      DO 190 J = 1,JMAX
         IF (IP(J+1).EQ.IP(J)) THEN
C No entries in col J.
            NUMD = NUMD + 1
         ELSE
            IP(J-NUMD) = IP(J)
         END IF
  190 CONTINUE
      IP(JMAX-NUMD+1) = IP(JMAX+1)

      IF (NUMD.GT.0) THEN
C Null rows to remove by shifting all row indices
         NUME = 0
         DO 200 K = 1,JMAX
            IF (IW(NMAX+K).EQ.0) THEN
C No entries in row K
               NUME = NUME + 1
            ELSE
C Store number of null rows with index less than K
               IW(NMAX+K) = NUME
            END IF
  200    CONTINUE
C Loop over all row indices.
         DO 210 I = 1,NZ
            K = IRN(I)
C Shift row K by number of null rows with index less than K
            IRN(I) = K - IW(NMAX+K)
  210    CONTINUE

C Set IW so that IW(J) is the original
C variable index for row (or column) J in assembled matrix.
         DO 220 K = 1,JMAX
            IF (IW(K).NE.0) THEN
C Variable K used in element variable lists
               IF (IW(NMAX+K).EQ.0) THEN
                  IW(K) = K
               ELSE
                  J = K - IW(NMAX+K)
                  IW(J) = K
               END IF
            END IF
  220    CONTINUE
      END IF

  230 CONTINUE
      INFO(6) = NUMD

C Set N to order of system
      N = JMAX - NUMD
      IF (NUMD.EQ.0) THEN
         DO 170 J = 1,JMAX
            IW(J) = J
  170    CONTINUE
      END IF

C Remove duplicates
      DO 240 I = 1,N
         IW(N+I) = 0
  240 CONTINUE

C NZ is number of entries in assembled matrix (once duplicates removed)
      NZ = 0
      JSTRT = 1
C Loop over the columns
      DO 260 J = 1,N
         JSTOP = IP(J+1) - 1
         IP(J+1) = IP(J)
         DO 250 ICOL = JSTRT,JSTOP
            IROW = IRN(ICOL)
            IF (IW(N+IROW).LT.J) THEN
               NZ = NZ + 1
               IRN(NZ) = IROW
               IP(J+1) = IP(J+1) + 1
               IW(N+IROW) = J
            END IF
  250    CONTINUE
         JSTRT = JSTOP + 1
  260 CONTINUE

C We have now finished.
      INFO(4) = NZ

      GO TO 440

  270 CONTINUE

C Fill-in the matrix entries
C Loop over elements. At start of each loop, L
C points to first entry in AELT for element IELT.
      IF (LCNTL(1)) THEN
C Symmetric case. Only lower triangular part of each element
C is stored, in packed column form.
         L = 1
         DO 300 IELT = 1,NELT
            JSTRT = ELTPTR(IELT)
            NVAR = ELTPTR(IELT+1) - ELTPTR(IELT)
            JSTOP = JSTRT + NVAR - 1
            DO 290 J = JSTRT,JSTOP
               JVAR = ELTVAR(J)
C Set the row indices of the entries in col. JVAR.
C Remember we only have lower triangular part of element matrix
               DO 280 I = J,JSTOP
                  IVAR = ELTVAR(I)
                  IF (IVAR.GE.JVAR) THEN
C Belongs in column JVAR
                     IP(JVAR) = IP(JVAR) - 1
                     K = IP(JVAR)
                     IRN(K) = IVAR
                     A(K) = AELT(L)
                  ELSE
C Belongs in column IVAR
                     IP(IVAR) = IP(IVAR) - 1
                     K = IP(IVAR)
                     IRN(K) = JVAR
                     A(K) = AELT(L)
                  END IF
                  L = L + 1
  280          CONTINUE
  290       CONTINUE
  300    CONTINUE
      ELSE
         L = 1
         DO 330 IELT = 1,NELT
            JSTRT = ELTPTR(IELT)
            NVAR = ELTPTR(IELT+1) - ELTPTR(IELT)
            JSTOP = JSTRT + NVAR - 1
            DO 320 J = JSTRT,JSTOP
               JVAR = ELTVAR(J)
C Set the row indices of the entries in col. JVAR
               DO 310 I = JSTRT,JSTOP
                  IVAR = ELTVAR(I)
                  IP(JVAR) = IP(JVAR) - 1
                  K = IP(JVAR)
                  IRN(K) = IVAR
                  A(K) = AELT(L)
                  L = L + 1
  310          CONTINUE
  320       CONTINUE
  330    CONTINUE
      END IF

C IP(J) now points to the start of column J.
C Remove null columns if LCNTL(3) = TRUE.

      NUMD = 0
      IF (.NOT.LCNTL(3)) GO TO 400

      NZ = IP(JMAX+1) - 1
      DO 350 I = 1,NZ
         K = IRN(I)
         IW(NMAX+K) = 1
  350 CONTINUE

      DO 360 J = 1,JMAX
         IF (IP(J+1).EQ.IP(J)) THEN
            NUMD = NUMD + 1
         ELSE
            IP(J-NUMD) = IP(J)
         END IF
  360 CONTINUE
      IP(JMAX-NUMD+1) = IP(JMAX+1)

      IF (NUMD.GT.0) THEN
C Shift all indices
         NUME = 0
         DO 370 K = 1,JMAX
            IF (IW(NMAX+K).EQ.0) THEN
               NUME = NUME + 1
            ELSE
               IW(NMAX+K) = NUME
            END IF
  370    CONTINUE
         DO 380 I = 1,NZ
            K = IRN(I)
            IRN(I) = K - IW(NMAX+K)
  380    CONTINUE

C Set IW so that IW(J) is the original
C variable index for row (or column) J in assembled matrix.
         DO 390 K = 1,JMAX
            IF (IW(K).NE.0) THEN
C Variable K used in element variable lists
               IF (IW(NMAX+K).EQ.0) THEN
                  IW(K) = K
               ELSE
                  J = K - IW(NMAX+K)
                  IW(J) = K
               END IF
            END IF
  390    CONTINUE
      END IF

  400 CONTINUE
      INFO(6) = NUMD

C Set N to order of system
      N = JMAX - NUMD
      IF (NUMD.EQ.0) THEN
         DO 340 J = 1,JMAX
            IW(J) = J
  340    CONTINUE
      END IF

C Remove duplicates
      DO 410 I = 1,N
         IW(N+I) = 0
  410 CONTINUE

      NZ = 0
      JSTRT = 1
      NZJ = 0
C Loop over the columns
      DO 430 J = 1,N
         JSTOP = IP(J+1) - 1
         IP(J+1) = IP(J)
         DO 420 ICOL = JSTRT,JSTOP
            IROW = IRN(ICOL)
            IF (IW(N+IROW).LE.NZJ) THEN
               NZ = NZ + 1
               IRN(NZ) = IROW
               A(NZ) = A(ICOL)
               IP(J+1) = IP(J+1) + 1
               IW(N+IROW) = NZ
            ELSE
C We have a duplicate in col J
               A(IW(N+IROW)) = A(IW(N+IROW)) + A(ICOL)
            END IF
  420    CONTINUE
         JSTRT = JSTOP + 1
         NZJ = NZ
  430 CONTINUE
      INFO(4) = NZ

  440 CONTINUE
C Does the user want to order by rows?
      IF (LCNTL(2)) CALL MC57BD(N,LCNTL(4),LA,A,NZ,IRN,IP)

  450 CONTINUE
      RETURN
 9000 FORMAT (/' Error return . INFO(1) =  ',I2)
 9010 FORMAT (' Value of ',A,' out of range. ',A,' = ',I10)
 9020 FORMAT (' ',A,' must be at least ',I12)
 9030 FORMAT (' Duplicate occurrence of variable ',I10,' in element ',
     +       I6)
 9040 FORMAT (' Smallest integer used to index variable is ',I12)
      END

C*************************************************
      SUBROUTINE MC57BD(N,YESA,LA,A,NZ,IRN,IP)
C This is essentially MC20B/BD but allows for sparsity pattern only
C     . .
C     .. Scalar Arguments ..
      INTEGER LA,N,NZ
      LOGICAL YESA
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(LA)
      INTEGER IP(N),IRN(NZ)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ACE
      INTEGER ICE,IK,J,JJ,K,KDUMMY,KLO,KMAX,KOR
C     ..
      IF (.NOT.YESA) GO TO 60

      KMAX = NZ
      DO 50 JJ = 1,N
         J = N + 1 - JJ
         KLO = IP(J) + 1
         IF (KLO.GT.KMAX) GO TO 40
         KOR = KMAX
         DO 30 KDUMMY = KLO,KMAX
C Items KOR, KOR+1, .... ,KMAX are in order
            ACE = A(KOR-1)
            ICE = IRN(KOR-1)
            DO 10 K = KOR,KMAX
               IK = IRN(K)
               IF (ICE.LE.IK) GO TO 20
               IRN(K-1) = IK
               A(K-1) = A(K)
   10       CONTINUE
            K = KMAX + 1
   20       IRN(K-1) = ICE
            A(K-1) = ACE
            KOR = KOR - 1
   30    CONTINUE
C Next column
   40    KMAX = KLO - 2
   50 CONTINUE
      GO TO 120

   60 CONTINUE
C Sparsity pattern only
      KMAX = NZ
      DO 110 JJ = 1,N
         J = N + 1 - JJ
         KLO = IP(J) + 1
         IF (KLO.GT.KMAX) GO TO 100
         KOR = KMAX
         DO 90 KDUMMY = KLO,KMAX
C Items KOR, KOR+1, .... ,KMAX are in order
            ICE = IRN(KOR-1)
            DO 70 K = KOR,KMAX
               IK = IRN(K)
               IF (ICE.LE.IK) GO TO 80
               IRN(K-1) = IK
   70       CONTINUE
            K = KMAX + 1
   80       IRN(K-1) = ICE
            KOR = KOR - 1
   90    CONTINUE
C Next column
  100    KMAX = KLO - 2
  110 CONTINUE

  120 CONTINUE
      RETURN
      END
