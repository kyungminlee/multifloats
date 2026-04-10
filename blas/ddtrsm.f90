SUBROUTINE ddtrsm (side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
 use multifloats
 TYPE(FLOAT64X2) alpha
 INTEGER lda, ldb, m, n
 CHARACTER diag, side, transa, uplo
 TYPE(FLOAT64X2) a(lda,*), b(ldb,*)
 LOGICAL lsame
 EXTERNAL :: lsame
 EXTERNAL :: xerbla
 ! INTRINSIC :: max
 TYPE(FLOAT64X2) temp
 INTEGER i, info, j, k, nrowa
 LOGICAL lside, nounit, upper
 TYPE(FLOAT64X2) one, zero
 one=1.0d+0
 zero=0.0d+0
 lside = lsame(side, "L")
 IF (lside) THEN
  nrowa = m
 ELSE
  nrowa = n
 END IF
 nounit = lsame(diag, "N")
 upper = lsame(uplo, "U")
 info = 0
 IF ((.NOT.lside).AND.(.NOT.lsame(side, "R"))) THEN
  info = 1
 ELSE IF ((.NOT.upper).AND.(.NOT.lsame(uplo, "L"))) THEN
  info = 2
 ELSE IF ((.NOT.lsame(transa, "N")).AND.(.NOT.lsame(transa, "T")).AND.(.NOT.ls&
 &ame(transa, "C"))) THEN
  info = 3
 ELSE IF ((.NOT.lsame(diag, "U")).AND.(.NOT.lsame(diag, "N"))) THEN
  info = 4
 ELSE IF (m<0) THEN
  info = 5
 ELSE IF (n<0) THEN
  info = 6
 ELSE IF (lda<max(1, nrowa)) THEN
  info = 9
 ELSE IF (ldb<max(1, m)) THEN
  info = 11
 END IF
 IF (info/=0) THEN
  CALL xerbla("DTRSM ", info)
  RETURN
 END IF
 IF (m==0.OR.n==0) RETURN
 IF (alpha==zero) THEN
  DO 20 j=1,n
  DO 10 i=1,m
  b(i, j) = zero
  10 CONTINUE
  20 CONTINUE
  RETURN
 END IF
 IF (lside) THEN
  IF (lsame(transa, "N")) THEN
   IF (upper) THEN
    DO 60 j=1,n
    IF (alpha/=one) THEN
     DO 30 i=1,m
     b(i, j) = alpha*b(i, j)
     30 CONTINUE
    END IF
    DO 50 k=m,1,-1
    IF (b(k, j)/=zero) THEN
     IF (nounit) b(k, j) = b(k, j)/a(k, k)
     DO 40 i=1,k-1
     b(i, j) = b(i, j)-b(k, j)*a(i, k)
     40 CONTINUE
    END IF
    50 CONTINUE
    60 CONTINUE
   ELSE
    DO 100 j=1,n
    IF (alpha/=one) THEN
     DO 70 i=1,m
     b(i, j) = alpha*b(i, j)
     70 CONTINUE
    END IF
    DO 90 k=1,m
    IF (b(k, j)/=zero) THEN
     IF (nounit) b(k, j) = b(k, j)/a(k, k)
     DO 80 i=k+1,m
     b(i, j) = b(i, j)-b(k, j)*a(i, k)
     80 CONTINUE
    END IF
    90 CONTINUE
    100 CONTINUE
   END IF
  ELSE
   IF (upper) THEN
    DO 130 j=1,n
    DO 120 i=1,m
    temp = alpha*b(i, j)
    DO 110 k=1,i-1
    temp = temp-a(k, i)*b(k, j)
    110 CONTINUE
    IF (nounit) temp = temp/a(i, i)
    b(i, j) = temp
    120 CONTINUE
    130 CONTINUE
   ELSE
    DO 160 j=1,n
    DO 150 i=m,1,-1
    temp = alpha*b(i, j)
    DO 140 k=i+1,m
    temp = temp-a(k, i)*b(k, j)
    140 CONTINUE
    IF (nounit) temp = temp/a(i, i)
    b(i, j) = temp
    150 CONTINUE
    160 CONTINUE
   END IF
  END IF
 ELSE
  IF (lsame(transa, "N")) THEN
   IF (upper) THEN
    DO 210 j=1,n
    IF (alpha/=one) THEN
     DO 170 i=1,m
     b(i, j) = alpha*b(i, j)
     170 CONTINUE
    END IF
    DO 190 k=1,j-1
    IF (a(k, j)/=zero) THEN
     DO 180 i=1,m
     b(i, j) = b(i, j)-a(k, j)*b(i, k)
     180 CONTINUE
    END IF
    190 CONTINUE
    IF (nounit) THEN
     temp = one/a(j, j)
     DO 200 i=1,m
     b(i, j) = temp*b(i, j)
     200 CONTINUE
    END IF
    210 CONTINUE
   ELSE
    DO 260 j=n,1,-1
    IF (alpha/=one) THEN
     DO 220 i=1,m
     b(i, j) = alpha*b(i, j)
     220 CONTINUE
    END IF
    DO 240 k=j+1,n
    IF (a(k, j)/=zero) THEN
     DO 230 i=1,m
     b(i, j) = b(i, j)-a(k, j)*b(i, k)
     230 CONTINUE
    END IF
    240 CONTINUE
    IF (nounit) THEN
     temp = one/a(j, j)
     DO 250 i=1,m
     b(i, j) = temp*b(i, j)
     250 CONTINUE
    END IF
    260 CONTINUE
   END IF
  ELSE
   IF (upper) THEN
    DO 310 k=n,1,-1
    IF (nounit) THEN
     temp = one/a(k, k)
     DO 270 i=1,m
     b(i, k) = temp*b(i, k)
     270 CONTINUE
    END IF
    DO 290 j=1,k-1
    IF (a(j, k)/=zero) THEN
     temp = a(j, k)
     DO 280 i=1,m
     b(i, j) = b(i, j)-temp*b(i, k)
     280 CONTINUE
    END IF
    290 CONTINUE
    IF (alpha/=one) THEN
     DO 300 i=1,m
     b(i, k) = alpha*b(i, k)
     300 CONTINUE
    END IF
    310 CONTINUE
   ELSE
    DO 360 k=1,n
    IF (nounit) THEN
     temp = one/a(k, k)
     DO 320 i=1,m
     b(i, k) = temp*b(i, k)
     320 CONTINUE
    END IF
    DO 340 j=k+1,n
    IF (a(j, k)/=zero) THEN
     temp = a(j, k)
     DO 330 i=1,m
     b(i, j) = b(i, j)-temp*b(i, k)
     330 CONTINUE
    END IF
    340 CONTINUE
    IF (alpha/=one) THEN
     DO 350 i=1,m
     b(i, k) = alpha*b(i, k)
     350 CONTINUE
    END IF
    360 CONTINUE
   END IF
  END IF
 END IF
 RETURN
END SUBROUTINE
