!************************************************************************

! File to read and re-write  He3 spectral function data

!************************************************************************
PROGRAM HE3SS
IMPLICIT NONE

INTEGER	nerg,nq,i,j

REAL*8 e(30),qmesh(30),qweigh(30)
REAL*8 r0t0(30,30),r1t0(30,30),r2t0(30,30)
REAL*8 r0t1(30,30),r1t1(30,30),r2t1(30,30)
real*8 factor,pi,hc
pi=2d0*dasin(1d0)
hc=197.327053D0
factor= hc
!factor=1d0


OPEN (UNIT=11,FILE='He3ss.dat',FORM='FORMATTED',STATUS='OLD')

  READ (11,*) nerg, nq
  READ (11,*) (e(i), i=1,nerg)
  READ (11,*) (qmesh(j), j=1,nq)
  READ (11,*) (qweigh(j),j=1,nq)
  READ (11,*) ((r0t0(i,j),r1t0(i,j),r2t0(i,j),j=1,nq),i=1,nerg)
  READ (11,*) ((r0t1(i,j),r1t1(i,j),r2t1(i,j),j=1,nq),i=1,nerg)

OPEN (UNIT=13,FILE='f0pd.dat',FORM='FORMATTED',STATUS='UNKNOWN')
OPEN (UNIT=14,FILE='f1pd.dat',FORM='FORMATTED',STATUS='UNKNOWN')
OPEN (UNIT=15,FILE='f2pd.dat',FORM='FORMATTED',STATUS='UNKNOWN')
DO j=1,nq
  WRITE (13,*) e(1),qmesh(j),r0t0(1,j)/2.D0 * factor
  WRITE (14,*) e(1),qmesh(j),r1t0(1,j)/2.D0 * factor
  WRITE (15,*) e(1),qmesh(j),r2t0(1,j)/2.D0 * factor
ENDDO

OPEN (UNIT=16,FILE='f0pc.dat',FORM='FORMATTED',STATUS='UNKNOWN')
OPEN (UNIT=17,FILE='f1pc.dat',FORM='FORMATTED',STATUS='UNKNOWN')
OPEN (UNIT=18,FILE='f2pc.dat',FORM='FORMATTED',STATUS='UNKNOWN')
DO i=2,nerg
   DO j=1,nq
      WRITE (16,*) e(i),qmesh(j),(r0t0(i,j)+r0t1(i,j))/2.D0 * factor
      WRITE (17,*) e(i),qmesh(j),(r1t0(i,j)+r1t1(i,j))/2.D0 * factor
      WRITE (18,*) e(i),qmesh(j),(r2t0(i,j)+r2t1(i,j))/2.D0 * factor
   ENDDO
ENDDO

OPEN (UNIT=19,FILE='f0nc.dat',FORM='FORMATTED',STATUS='UNKNOWN')
OPEN (UNIT=20,FILE='f1nc.dat',FORM='FORMATTED',STATUS='UNKNOWN')
OPEN (UNIT=21,FILE='f2nc.dat',FORM='FORMATTED',STATUS='UNKNOWN')
DO i=2,nerg
   DO j=1,nq
      WRITE (19,*) e(i),qmesh(j),2.D0*r0t1(i,j) * factor
      WRITE (20,*) e(i),qmesh(j),2.D0*r1t1(i,j) * factor
      WRITE (21,*) e(i),qmesh(j),2.D0*r2t1(i,j) * factor
   ENDDO
ENDDO

END PROGRAM
