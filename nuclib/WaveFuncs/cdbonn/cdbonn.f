C **********************************************************************
	SUBROUTINE CDBONN (q,u0,u2)
C
C  Deuteron wave function from CD-Bonn NN potential model.
C  q in 1/fm, u0,u2 in fm^3/2.
C
C  Normalization \int dq q^2 (u0^2+u2^2) = 1.
C
C  Sent by Charlotte Elster, April 8, 2009.
C **********************************************************************
        IMPLICIT NONE
        INTEGER id,nq
        PARAMETER (nq=95)
        REAL*8  q,u0,u2,
     &		qgrid(nq),uqgrid(nq),wqgrid(nq),weight(nq),dum
        LOGICAL init /.FALSE./
        REAL*8  pi,hc
	SAVE

Cf2py intent(in) q
Cf2py intent(out) u0
Cf2py intent(out) u2

        pi = 4*DATAN(1.D0)
        hc = 197.327D0		! MeV.fm conversion factor

        IF (init) GO TO 999     ! Data already read
C...Read data from file
	OPEN (10, FORM='FORMATTED',
     &		FILE='WaveFuncs/cdbonn/cdbn.qwave',
     &		STATUS='OLD')
        READ (10,100)
C...Momentum space [qgrid in MeV, uqgrid in MeV^-3/2]
        DO id=1,nq
          READ (10,*) qgrid(id), weight(id), uqgrid(id), wqgrid(id)
          qgrid(id) = qgrid(id) / hc            ! MeV => 1/fm
          uqgrid(id) = uqgrid(id) * hc**1.5D0   ! MeV^-3/2 => fm^3/2
          wqgrid(id) = wqgrid(id) * hc**1.5D0
        ENDDO
 100    FORMAT (2(1X/))
 101    FORMAT (3X,D13.6,3X,D20.13,3X,D20.13)
	PRINT *, '... CD-Bonn model read...'
	init = .TRUE.

C...Evaluate wave function
c 999	u0 = DQDVAL (q,nq,qgrid,uqgrid,.FALSE.)
c	u2 = DQDVAL (q,nq,qgrid,wqgrid,.FALSE.)
 999	CALL Pinterp (qgrid,uqgrid,nq,q,u0,dum,2)
	CALL Pinterp (qgrid,wqgrid,nq,q,u2,dum,2)

        RETURN
        END
