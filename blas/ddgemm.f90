SUBROUTINE ddgemm (transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc&
&)
 use multifloats
 TYPE(FLOAT64X2) alpha, beta
 INTEGER k, lda, ldb, ldc, m, n
 CHARACTER transa, transb
 TYPE(FLOAT64X2) a(lda,*), b(ldb,*), c(ldc,*)
 LOGICAL lsame
 EXTERNAL :: lsame
 EXTERNAL :: xerbla
 ! INTRINSIC :: max
 TYPE(FLOAT64X2) temp
 INTEGER i, info, j, l, nrowa, nrowb
 LOGICAL nota, notb
 TYPE(FLOAT64X2) one, zero
 one=1.0d+0
 zero=0.0d+0
 nota = lsame(transa, "N")
 notb = lsame(transb, "N")
 IF (nota) THEN
  nrowa = m
 ELSE
  nrowa = k
 END IF
 IF (notb) THEN
  nrowb = k
 ELSE
  nrowb = n
 END IF
 info = 0
 IF ((.NOT.nota).AND.(.NOT.lsame(transa, "C")).AND.(.NOT.lsame(transa, "T"))) &
 &THEN
  info = 1
 ELSE IF ((.NOT.notb).AND.(.NOT.lsame(transb, "C")).AND.(.NOT.lsame(transb, "T&
 &"))) THEN
  info = 2
 ELSE IF (m<0) THEN
  info = 3
 ELSE IF (n<0) THEN
  info = 4
 ELSE IF (k<0) THEN
  info = 5
 ELSE IF (lda<max(1, nrowa)) THEN
  info = 8
 ELSE IF (ldb<max(1, nrowb)) THEN
  info = 10
 ELSE IF (ldc<max(1, m)) THEN
  info = 13
 END IF
 IF (info/=0) THEN
  CALL xerbla("DGEMM ", info)
  RETURN
 END IF
 IF ((m==0).OR.(n==0).OR.(((alpha==zero).OR.(k==0)).AND.(beta==one))) RETURN
 IF (alpha==zero) THEN
  IF (beta==zero) THEN
   DO 20 j=1,n
   DO 10 i=1,m
   c(i, j) = zero
   10 CONTINUE
   20 CONTINUE
  ELSE
   DO 40 j=1,n
   DO 30 i=1,m
   c(i, j) = beta*c(i, j)
   30 CONTINUE
   40 CONTINUE
  END IF
  RETURN
 END IF
 IF (notb) THEN
  IF (nota) THEN
   DO 90 j=1,n
   IF (beta==zero) THEN
    DO 50 i=1,m
    c(i, j) = zero
    50 CONTINUE
   ELSE IF (beta/=one) THEN
    DO 60 i=1,m
    c(i, j) = beta*c(i, j)
    60 CONTINUE
   END IF
   DO 80 l=1,k
   temp = alpha*b(l, j)
   DO 70 i=1,m
   c(i, j) = c(i, j)+temp*a(i, l)
   70 CONTINUE
   80 CONTINUE
   90 CONTINUE
  ELSE
   DO 120 j=1,n
   DO 110 i=1,m
   temp = zero
   DO 100 l=1,k
   temp = temp+a(l, i)*b(l, j)
   100 CONTINUE
   IF (beta==zero) THEN
    c(i, j) = alpha*temp
   ELSE
    c(i, j) = alpha*temp+beta*c(i, j)
   END IF
   110 CONTINUE
   120 CONTINUE
  END IF
 ELSE
  IF (nota) THEN
   DO 170 j=1,n
   IF (beta==zero) THEN
    DO 130 i=1,m
    c(i, j) = zero
    130 CONTINUE
   ELSE IF (beta/=one) THEN
    DO 140 i=1,m
    c(i, j) = beta*c(i, j)
    140 CONTINUE
   END IF
   DO 160 l=1,k
   temp = alpha*b(j, l)
   DO 150 i=1,m
   c(i, j) = c(i, j)+temp*a(i, l)
   150 CONTINUE
   160 CONTINUE
   170 CONTINUE
  ELSE
   DO 200 j=1,n
   DO 190 i=1,m
   temp = zero
   DO 180 l=1,k
   temp = temp+a(l, i)*b(j, l)
   180 CONTINUE
   IF (beta==zero) THEN
    c(i, j) = alpha*temp
   ELSE
    c(i, j) = alpha*temp+beta*c(i, j)
   END IF
   190 CONTINUE
   200 CONTINUE
  END IF
 END IF
 RETURN
END SUBROUTINE
