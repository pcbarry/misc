**************************************************************
***       File contains various interpolation codes        ***
**************************************************************
***      Code obtained from Alberto Accardi, Dec. 2008     ***
**************************************************************
*
* II  Polynomial function interpolation of given order
      subroutine pinterp(xa,ya,n,x,y,dy,order)
*     programmer: Alberto Accardi
*     date: 2/05/01
*
*  A. COMMENTARY
*
*     Performs an interpolation using a polynomial function
*     interpolation at a given order: given an x, it uses "order" points 
*     to its left and "order" to its right to perform the interpolation
*
*     xa(*) = (DP) array with tabulated abscissae (any dimension)
*     ya(*) = (DP) array with tabulated function  (any dimension)
*     n     = (I)  number of tabulated points  
*     x     = (DP) abscissa at which compute the interpolated function
*     y     = (DP) value of the function at x
*     dy    = (DP) estimated error (usually larger than real error)
*     order = (I)  order of the interpolation (see intro)  
*                  If order = 0 performs a linear interpolation
*                  between the nearest neighbours lattice point
*
*  B. DECLARATIONS
*
      implicit none

*    *** variables
      
      integer n, order

      double precision xa(*), ya(*), x, y, dy, tempx(n)
     :     , x1(2*order), y1(2*order), xmax, xmin, ymax, ymin  

      integer i, nlow, nmin

*
*  C. ACTION
*

      do i = 1, n
         tempx(i) = xa(i)
      end do
      call hunt(tempx,n,x,nlow)

      if (order.ge.1) then
         if (nlow.lt.order) then
            nmin = 0
         else if (nlow.le.n-order) then
            nmin = nlow-order
         else
            nmin = n-2*order
         end if
         do i = 1, 2*order
            x1(i) = xa(nmin+i) 
            y1(i) = ya(nmin+i) 
         end do
         call polintnum(x1,y1,2*order,x,y,dy)
      else
         ymax = ya(nlow+1)
         ymin = ya(nlow)
         xmax = xa(nlow+1)
         xmin = xa(nlow)
         y = ymin + (ymax-ymin)/(xmax-xmin) * (x-xmin)
      end if

      return
      end


************************************************************************
*
* III search in an ordered table 
      SUBROUTINE hunt(xx,n,x,jlo)
*     from "NUMERICAL RECIPES IN F77"
*
*  A. COMMENTARY
*
*     Given an array xx(1:n) and given a value x, returns a value j
*     suchthat x is between xx(j) and xx(j+1). xx(1:n) must be monotonic,
*     either decreasing or increasing. j=0 or j=n is returned to
*     indicate that x is out of range.
*
*  B. DECLARATIONS
*
      implicit none

*    *** variables
      
      INTEGER jlo,n

      double precision x,xx(n)
      INTEGER inc,jhi,jm
      LOGICAL ascnd

*
*  C. ACTION
*

      ascnd=xx(n).gt.xx(1)
      if(jlo.le.0.or.jlo.gt.n)then
        jlo=0
        jhi=n+1
        goto 3
      endif
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then
1       jhi=jlo+inc
        if(jhi.gt.n)then
          jhi=n+1
        else if(x.ge.xx(jhi).eqv.ascnd)then
          jlo=jhi
          inc=inc+inc
          goto 1
        endif
      else
        jhi=jlo
2       jlo=jhi-inc
        if(jlo.lt.1)then
          jlo=0
        else if(x.lt.xx(jlo).eqv.ascnd)then
          jhi=jlo 
          inc=inc+inc
          goto 2
        endif
      endif
3     if(jhi-jlo.eq.1)return
      jm=(jhi+jlo)/2
      if(x.gt.xx(jm).eqv.ascnd)then
        jlo=jm
      else
        jhi=jm
      endif
      goto 3
      END


************************************************************************
*
* IV  Polynomial interpolation and extrapolation
      SUBROUTINE polintnum(xa,ya,n,x,y,dy)
*     from "NUMERICAL RECIPES IN F77"
*
*  A. COMMENTARY
*
*     Given arrays xa and ya of length n, and given a value x, this
*     routine returns a value y and an error estimate dy. If P(x) is the
*     polynomial of degree N-1 such that P(xa_i) = ya_i, i=1,...,n
*     then the returned value y = P(x).
*
*  B. DECLARATIONS
*
      implicit none

*    *** variables
      
      INTEGER n,NMAX
      double precision dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      double precision den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)

*
*  C. ACTION
*

      if (n.gt.nmax) then
         print*, 'ERROR(polintnum): order larger than max', n,'>', nmax 
         stop
      end if
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.)pause 'failure in polintnum'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END

      
**************************************************************
****    (Code obtained from Peter Monaghan, Dec. 2008)     ***
****                (need to modify for use)               ***
**************************************************************
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C     Linear Interpolation Routine:
C
C       Calculates Y(X0) given array Y(X) assuming a linear
C       dependence:  Y = AX + B
C
C       This routine assumes the X values are arranged in
C       ascending order but does not assume a uniform X spacing
C       of the array.
C------------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION RINTERP(X,Y,X0,NBINS)
      IMPLICIT NONE
C
ccc      INCLUDE 'summary.cmn'
C
      DOUBLE PRECISION X(1),Y(1),X0,DET,A,B
      INTEGER NBINS,I1,I2,I,I_EXTRAP
C
      IF(X0 .LT. X(1) .OR. X0 .GT. X(NBINS))
CXXX #     TYPE *,' Extrapolating outside range:  X = ',X0
     #     I_EXTRAP = I_EXTRAP + 1
C
      IF(X0 .LE. X(1)) THEN
         I1 = 1
         I2 = 2
      ELSEIF(X0 .GE. X(NBINS-1)) THEN
         I1 = NBINS-1
         I2 = NBINS
      ELSE
         DO I=2,NBINS-1
            IF(X0 .LE. X(I)) THEN
               I1 = I-1
               I2 = I
               GOTO 1
            ENDIF
         ENDDO
      ENDIF
C
1     DET = X(I1)-X(I2)
      A = (Y(I1)-Y(I2))/DET
      B = (X(I1)*Y(I2)-X(I2)*Y(I1))/DET
C
      RINTERP = A*X0 + B
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C     Quadratic Interpolation Routine:
C
C       Calculates Y(X0) given array Y(X) assuming a quadratic
C       dependence:  Y = AX^2 + BX + C
C
C       This routine assumes the X values are arranged in
C       ascending order but does not assume a uniform X spacing
C       of the array.
C------------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION RINTERPQ(X,Y,X0,NBINS)
      IMPLICIT NONE
C
ccc      INCLUDE 'summary.cmn'
C
      DOUBLE PRECISION X(1),Y(1),X0,DET,A,B,C
      INTEGER NBINS,I1,I2,I3,I,I_EXTRAP
C
      IF(X0 .LT. X(1) .OR. X0 .GT. X(NBINS))
CXXX #     TYPE *,' Extrapolating outside range:  X = ',X0
     #     I_EXTRAP = I_EXTRAP + 1
C
      IF(X0 .LE. X(1)) THEN
         I1 = 1
         I2 = 2
         I3 = 3
      ELSEIF(X0 .GE. X(NBINS-1)) THEN
         I1 = NBINS-2
         I2 = NBINS-1
         I3 = NBINS
      ELSE
         DO I=2,NBINS-1
            IF(X0 .LE. X(I)) THEN
               I1 = I-1
               I2 = I
               I3 = I+1
               GOTO 1
            ENDIF
         ENDDO
      ENDIF
C
1     DET = (X(I2)-X(I3))*(X(I1)**2 - X(I1)*(X(I2)+X(I3)) + X(I2)*X(I3))
      A = ( Y(I1)*(X(I2)-X(I3)) - X(I1)*(Y(I2)-Y(I3)) + Y(I2)*X(I3)
     #        - X(I2)*Y(I3) )/DET
      B = ( X(I1)**2*(Y(I2)-Y(I3)) - Y(I1)*(X(I2)**2-X(I3)**2)
     #        + X(I2)**2*Y(I3) - X(I3)**2*Y(I2) )/DET
      C = ( X(I1)**2*(X(I2)*Y(I3)-X(I3)*Y(I2))
     #        - X(I1)*(X(I2)**2*Y(I3)-X(I3)**2*Y(I2))
     #        + Y(I1)*(X(I2)**2*X(I3)-X(I3)**2*X(I2)) )/DET
C
      RINTERPQ = A*X0**2 + B*X0 + C
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C     Exponential Interpolation Routine:
C
C       Calculates Y(X0) given array Y(X) assuming the exponential
C       form:  Y = EXP(AX^2 + BX + C)
C
C       This routine assumes the X values are arranged in
C       ascending order but does not assume a uniform X spacing
C       of the array.
C------------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION RINTERPEXP(X,Y,X0,NBINS)
      IMPLICIT NONE
C
ccc      INCLUDE 'summary.cmn'
C
      DOUBLE PRECISION X(1),Y(1),X0,DET,A,B,C,Y1,Y2,Y3
      INTEGER NBINS,I1,I2,I3,I,I_EXTRAP,I_INTERPEXP_NP
C
      IF(X0 .LT. X(1) .OR. X0 .GT. X(NBINS))
CXXX #     TYPE *,' Extrapolating outside range:  X = ',X0
     #     I_EXTRAP = I_EXTRAP + 1
C
      IF(X0 .LE. X(1)) THEN
         I1 = 1
         I2 = 2
         I3 = 3
      ELSEIF(X0 .GE. X(NBINS-1)) THEN
         I1 = NBINS-2
         I2 = NBINS-1
         I3 = NBINS
      ELSE
         DO I=2,NBINS-1
            IF(X0 .LE. X(I)) THEN
               I1 = I-1
               I2 = I
               I3 = I+1
               GOTO 1
            ENDIF
         ENDDO
      ENDIF
C
C ----------------------------------------------------------------------
C     If all three Y-values are > 0, perform quadratic interpolation
C     on their logarithms; otherwise return zero as the result.
C ----------------------------------------------------------------------
C
1     IF(Y(I1).GT.0.D0 .AND. Y(I2).GT.0.D0 .AND. Y(I3).GT.0.D0) THEN
         Y1 = LOG(Y(I1))
         Y2 = LOG(Y(I2))
         Y3 = LOG(Y(I3))
      ELSE
CXXX     TYPE *,' RINTERPEXP:  non-positive y-value; setting to zero '
         I_INTERPEXP_NP = I_INTERPEXP_NP + 1
         RINTERPEXP = 0.D0
         RETURN
      ENDIF
C
      DET = (X(I2)-X(I3))*(X(I1)**2 - X(I1)*(X(I2)+X(I3)) + X(I2)*X(I3))
      A = ( Y1*(X(I2)-X(I3)) - X(I1)*(Y2-Y3) + Y2*X(I3)
     #        - X(I2)*Y3 )/DET
      B = ( X(I1)**2*(Y2-Y3) - Y1*(X(I2)**2-X(I3)**2)
     #        + X(I2)**2*Y3 - X(I3)**2*Y2 )/DET
      C = ( X(I1)**2*(X(I2)*Y3-X(I3)*Y2)
     #        - X(I1)*(X(I2)**2*Y3-X(I3)**2*Y2)
     #        + Y1*(X(I2)**2*X(I3)-X(I3)**2*X(I2)) )/DET
C
      RINTERPEXP = EXP(A*X0**2 + B*X0 + C)
C
      RETURN
      END
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C     2-Dimensional Linear Interpolation Routine:
C
C       Calculates F(X0,Y0) given array F(X,Y) by two successive
C       interpolations, first in X and then in Y.
C
C       Assumes uniform spacing of array in X and Y.
C------------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION RINTERP2D(X,Y,F,X0,Y0,NX,NY,DELX,DELY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
ccc      INCLUDE 'summary.cmn'
C
      DIMENSION X(300),Y(600),F(300,600)
C
      I = DINT((X0-X(1))/DELX) + 1
      J = DINT((Y0-Y(1))/DELY) + 1
      IF((I+1.GT.NX).OR.(I.LT.1).OR.(J+1.GT.NY).OR.(J.LT.1))THEN
        I_EXTRAP = I_EXTRAP + 1
        RINTERP2D = 0.D0
        RETURN
      ENDIF
C
      RINTX1 = ((X0-X(I))/(X(I+1)-X(I)))*(F(I+1,J)-F(I,J))
     #                  + F(I,J)
      RINTX2 = ((X0-X(I))/(X(I+1)-X(I)))*(F(I+1,J+1)-F(I,J+1))
     #                  + F(I,J+1)
      RINTERP2D = ((Y0-Y(J))/(Y(J+1)-Y(J)))*(RINTX2-RINTX1) + RINTX1
C
      RETURN
      END
