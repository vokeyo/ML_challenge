      SUBROUTINE L2_GENERIC( PARAM_SET_NUMBER,HARDWARE,RESULT_FLAG,
     &  EXTRA_FLAG)
C----------------------------------------------------------------------
C-
C-   Purpose and Methods : Generic TOOL :
C-
C-   Inputs  : PARAM_SET_NUMBER : # of parameter set to use
C-             HARDWARE:          mask of set bits for LV1 trigger which started
C-                                  this filter.
C-   Outputs : RESULT_FLAG :      Flag set to TRUE when we want to pass tool
C-                                  under this PARAM_SET_NUMBER
C-             EXTRA_FLAG  :      Set to TRUE when we want to pass event and
C-                                  do no further filtering. (NOT IMPLEMENTED)
C-   Controls:
C-
C-   Created  23-Dec-1991   James T. Linnemann
C-   Updated   8-JAN-1992   James T. Linnemann  use ESUM bank
C-
C----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER PARAM_SET_NUMBER,HARDWARE
      LOGICAL EXTRA_FLAG,RESULT_FLAG
      INCLUDE 'D0$INC:ZEBCOM.INC'               ! zebra main store
      CHARACTER*80 MSG
      INTEGER NUM_FOUND
      REAL    ET_MIN
      REAL    CUT1,CUT2   ! generic cuts of type real
      INTEGER ICUT1,ICUT2 ! generic cuts of type integer
      LOGICAL LOG1,LOG2   ! generic parameters of type logical
      CHARACTER*64 CHAR1,CHAR2  ! generic parameters of type character
C
      INTEGER IP,NPARIN,IER
      LOGICAL EZERROR,OK
      INCLUDE 'D0$PARAMS:ESUM.PARAMS'
      INTEGER NGOT
      INTEGER NFOUND(ID_ALL:LAST_TYPE),ID,IFLAG
      REAL ET,ETA_PHYS,ETA_DET,THETA,PHI
      INTEGER I,J,N,NCHR
      REAL    E1,E2,PX1,PX2,PY1,PY2,PZ1,PZ2,MASS
      INTEGER NMAX
      PARAMETER(NMAX = 40)    !be generous
      INTEGER IORDER(NMAX),WORK(NMAX)
      REAL L2_VERT, THETA_FROM_ETA, ETA_FROM_THETA
      REAL RETA, RTHETA, CAL_TH
C----------------------------------------------------------------------
C: Statement functions
      THETA_FROM_ETA(RETA ) = 2*ATAN(EXP(-(RETA)))
      ETA_FROM_THETA(RTHETA) =  -ALOG(MAX(TAN((RTHETA)/2.),1.E-9))
C----------------------------------------------------------------------
      EXTRA_FLAG = .FALSE.
      RESULT_FLAG = .FALSE.
C
C...first, carefully retrieve cuts from RCP
      IP = PARAM_SET_NUMBER !I'm getting tired of typing this
      CALL EZPICK('L2_GENERIC') ! downloaded from configuration file (.FILTs)
      OK = .NOT.EZERROR(IER)
      IF (OK) THEN
C...is IP consistent with the number of sets which exist?
        IF (IER .EQ. 0) CALL EZGET('NUMBER_OF_SETS',NPARIN,IER)
        IF (IER.EQ.0) THEN
          IF ((IP.LE.0).OR.(IP.GT.NPARIN)) THEN
            WRITE(MSG,'(A,I5,A,I5)')' parameter set requested = ',
     &        PARAM_SET_NUMBER, ' but only had ',NPARIN
            CALL ERRMSG('L2_GENERIC','L2_GENERIC',MSG,'F')
            GO TO 999
          ENDIF
C...Get the parameters from the IPth set
          IF (IER.EQ.0) CALL EZGETA('NUM_FOUND',IP,IP,1,NUM_FOUND,IER)
          IF (IER.EQ.0) CALL EZGETA('ET_MIN',IP,IP,1,ET_MIN,IER)
          IF (IER.EQ.0) CALL EZGETA('ICUT1',IP,IP,1,ICUT1,IER)
          IF (IER.EQ.0) CALL EZGETA('ICUT2',IP,IP,1,ICUT2,IER)
          IF (IER.EQ.0) CALL EZGETA('CUT1',IP,IP,1,CUT1,IER)
          IF (IER.EQ.0) CALL EZGETA('CUT2',IP,IP,1,CUT2,IER)
          IF (IER.EQ.0) CALL EZGETS('CHAR1',IP,CHAR1,NCHR,IER)
          IF (IER.EQ.0) CALL EZGETS('CHAR2',IP,CHAR2,NCHR,IER)
          IF (IER.EQ.0) CALL EZGETA('LOG1',IP,IP,1,LOG1,IER)
          IF (IER.EQ.0) CALL EZGETA('LOG2',IP,IP,1,LOG2,IER)
        ENDIF
      ENDIF
      IF (.NOT.OK) THEN
        CALL ERRMSG('L2_GENERIC','L2_GENERIC','Couldn''t find bank','F')
      ELSEIF (IER .NE. 0) THEN      ! Error reading RCP
        CALL ERRMSG('L2_GENERIC','L2_GENERIC_PARAMETERS',
     &          'Couldn''t find parameter','F')
      ENDIF
C
C...now you can actually do some cutting
C
C...this is an example of cutting on electron pair mass
C...this example should be run after L2_EM with NUM_EM set to 2
C----------------------------------------------------------------------
      CALL GTESUM_COUNTS('FILT',NFOUND,IER)
      IF (IER.NE.0) THEN
        WRITE(MSG,'(A,I)')'IER From GTESUM_COUNTS = ',IER
        CALL ERRMSG('GTESUM_COUNTS','L2_GENERIC',MSG,'W')
      ENDIF
C...it's overkill, but I'll sort anyway; could do IORDER(k) -> k
      NGOT = NFOUND(ID_ELECTRON)
      CALL GTESUM_SORT('FILT',ID_ELECTRON,NMAX,IORDER,WORK,IER)
      IF (IER.NE.0) THEN
C...this should only happen if you didn't have enough space reserved (NMAX)
        WRITE(MSG,'(A,I)')'IER From GTESUM_SORT (FILT) = ',IER
        CALL ERRMSG('GTESUM_SORT','L2_GENERIC',MSG,'W')
      ENDIF
      N = 0
      DO I = 1,MIN(NGOT,NMAX)
        CALL GTESUM('FILT',ID_ELECTRON,IORDER(I),ET,ETA_PHYS,
     &    ETA_DET,PHI,IFLAG,IER)
        IF (IER.NE.0) THEN
C...this should only happen if you asked for more than was claimed
          WRITE(MSG,'(A,I)')'IER From GTESUM = ',IER
          CALL ERRMSG('GTESUM','L2_GENERIC',MSG,'W')
        ELSE
          THETA = THETA_FROM_ETA(ETA_PHYS)
          E1 = ET/SIN(THETA)
          PX1 = ET*COS(PHI)
          PY1 = ET*SIN(PHI)
          PZ1 = E1*COS(THETA)
          DO J = I+1,MIN(NGOT,NMAX)
            CALL GTESUM('FILT',ID_ELECTRON,IORDER(J),ET,ETA_PHYS,
     &        ETA_DET,PHI,IFLAG,IER)
            IF (IER.NE.0) THEN
C...this should only happen if you asked for more than was claimed
              WRITE(MSG,'(A,I)')'IER From GTESUM = ',IER
              CALL ERRMSG('GTESUM','L2_GENERIC',MSG,'W')
            ELSE
              THETA = THETA_FROM_ETA(ETA_PHYS)
              E2 = ET/SIN(THETA)
              PX2 = ET*COS(PHI)
              PY2 = ET*SIN(PHI)
              PZ2 = E2*COS(THETA)
              MASS = E1*E2 - PX1*PX2 - PY1*PY2 - PZ1*PZ2
              IF (MASS.GT.0) MASS = SQRT(2.*MASS)
              IF (MASS.GT.ET_MIN) THEN
                N = N + 1
                IF (N.GE.NUM_FOUND) GO TO 999
              ENDIF
            ENDIF
          ENDDO
        ENDIF
      ENDDO
  999 CONTINUE
      RESULT_FLAG = N.GE.NUM_FOUND
      IF (OK) CALL EZRSET
      RETURN
      END
