      subroutine bayes(noc, noe, nsteps, mode, er_mode, keep, ier)
      implicit none
      integer SIZE_C, SIZE_E
      PARAMETER (SIZE_C = 100)
      PARAMETER (SIZE_E = 100)
      integer noc, noe, nsteps, mode, er_mode, keep, ier
      
      include 'bayes_c.for'
*     Attention: ALL variables of the common are REAL !! 

***********************************************************************
*   B A Y E S : Unfolding program based on the Bayes Theorem
*               by Giulio D'Agostini, 19/12/93
*                                     28/12/93
*                                     12/01/94
*                                     08/02/94 ( including also 
*                                                comments from L. Feld )
*
*   Update :    21/7/94:   - faster calculation of Vc( , )
*                          - option OFFD to evaluate only diagonal 
*                            elements of Vc( , )
*
*                5/7/95    - Vc_u( , ) to calculate the covariance
*                            matrix taking into account also for the 
*                            overall absolute normalization due to the
*                            observed number of events (important for
*                            small numbers of observed events)
*                            ( the option is effective only if the 
*                              parameter IE of er_mode is 2 )
*
*               23/3/96    - Modification by Hartmut Rick (rick@dice2.desy.de)
*                            to speed up the calculation of the covariance 
*                            matrix of the results.
*               12/4/96    - Minor update of the Rick modification
*
*   Reference              NIM A362 (1995) 487.
***********************************************************************
*
*  I/O (*)
*   I   mode       : running mode: 
*                      1: - the smearing matrix is provided in terms of 
*                           probabilities pec(j,i), and they are assumed 
*                           without errors;  
*                      2: - the smearing matrix is provided giving the 
*                           number of MC events generated for each of the
*                           cause-cells nc_mc(i) and array of the number of 
*                           reconstructed events nec_mc(j,i). 
*                           pec(j,i) is calculated by the program.
*                           This mode allows to take into account the errors
*                           on the smearing matrix.  
*   I   er_mode    : options on the error calculation er_mode=OFFD*(10*IE+IM)
*                    ( Standard value should be "er_mode=22" ) 
*                    IE ( errors from the observed number of events ):
*                       0 : no error
*                       1 : Poisson approximation
*                       2 : multinomial distribution 
*                    IM ( errors from smearing matrix ):
*                       0 : no error ( even if "mode=2" )
*                       1 : ( only if "mode=2" ): Poisson approximation
*                       2 : ( only if "mode=2" ): multinomial distribution
*                  OFFD ( calculation of the offdiagonal elements of
*                         the covariance matrix ):
*                      +1 :  offdiagonal elements calculated;
*                      -1 :  only diagonal elements are evaluated.
*              Attention: Please be aware that the computation of the
*                         full covariance matrix can be extremely time
*                         consuming in case of a large number of cells.
*                         In this case it is convenient to request the
*                         calculation only at the very end.
*                         The parameter OFFD allows to calculate at least
*                         the variances in case the complete covariance
*                         matrix requires too much time. 
*                         (see also note added atthe end of the distribution
*                          file)
*   I   noc        : number of cells of the "true" physical distribution
*   I   noe        : number of cells of the measured distribution
*  I/O  nsteps     : I: max number of iterations :
*                       > 0  -> nsteps iterations are performed, UNLESS chi^2
*                               between two unfolded distributions becomes 
*                               smaller than "delta_chi2";
*                       = 0  -> nsteps = 9999 ( infinite ), UNLESS chi^2 
*                               between two unfolded distributions becomes 
*                               smaller than "delta_chi2";
*                       < 0  -> (-nsteps) iterations are forced independently
*                               of delta_chi2
*                  : O: actual number of steps performed
*   I   keep       : it should be 0 in the first call
*                             and 1 in the following ones ( it avoids to 
*                             perform again some inializations, but the 
*                             user should be carefull not to change the 
*                             common BAYES ! )  
*   O   ier        : returns 0 if OK
*   I   ne(j)      : number of events in the j-th cell of the measured values.
*   I   pec(j,i)   : P(E_j|C_i) : conditional probabilities of the
*                    j-th cell of the measured values to come from i-th 
*                    cell of the true values.
*                    ( only if mode = 1;
*                      if mode=2, they are calculated automatically ) 
*   I   nc_mc(i)   : number of MC events generated for each of the cause-cells.
*                    ( only if mode = 2)
*                    They DO NOT need to be genearated according to the 
*                    initial probabilities, and in fact it is sometime 
*                    preferable to generate them almost flat.
*   I   nec_mc(j,i): number of MC events generated in the i-th cause-cell
*                    and reconstructed in the j-th effect-cell. 
*                    ( only if mode = 2)
*   O   pce(i,j)   : P(C_i|E_j) : conditional probabilities of the
*                    i-th cell of the true values to come from j-th 
*                    cell of the measured values 
*  I/O  pc(i)      : I : initial probability of the i-th true value cell;
*                    O : final probability of the i-th true value cell;
*   O   nc(i)      : estimated number of events in the i-th cell of the 
*                    true values.
*   O   Vc0( , )   : contribution to the covariance matrix of nc(i) due
*                    to observed events and assuming the constraint 
*                    to the total number of observed events
*   O   Vc0_u( , ) : like Vc0( , ), but taking into account also the 
*                    overal normalization uncertainty due to the
*                    total number of observed events. 
*                    Notice: if IE=1 ( see er_mode ) the Vc0_u=Vc0 !!
*   O   Vc1( , )   : contribution to the covariance matrix of nc(i) due
*                    the smearing matrix 
*   O   Vc( , )    : covariance matrix of nc(i)  without overal 
*                    normalization uncertainty  ( = Vc0 + Vc1 )
*   O   Vc_u( , )  : like Vc( , ), but taking into account also the
*                    overal normalization uncertainty due to the
*                    total number of observed events. USE THIS 
*                    FOR THE GLOBAL UNCERTANTIES (expecially if the number
*                    of observed events is small)
*                    Notice: if IE=1 ( see er_mode ) the Vc_u=Vc !!
*   O   Vpc( , )   : covariance matrix on the final probabilities (it 
*                    may be usefull if one is interested to the shape
*                    of the distribution ( -> pc(i)  )
*   O   Vpc0( , )  : like Vc0, but referred to the final probabilities
*   O   Vpc1( , )  : like Vc1, but referred to the final probabilities
*                    ( notice: Vpc0_u( , ) and Vpc_u( , ) make no 
*                      sense, as the the total normalization is ininfluent )
* 
*                ---------------------------------------
*   (*) I: input; O: output
*       Notice: - the quantities of type "I" are not overwritten by
*                 the program;
*               - the quantities of kind "I/O" are overwritten. 
*              
*   ATTENTION:
*               1) the upper limits of the arrays are the parameters
*                       SIZE_C ( max # of causes ) 
*                       SIZE_E ( max # of effects ). 
*                  They must be set ALSO in the calling program;
*               2) The above description follows the convenction:
*                       - the index "i" runs from 1 to noc;
*                       - the index "j" runs from 1 to noe;
*               3) noc and noe may be different, but smaller respectivelly
*                  of SIZE_C e SIZE_E  
*
*
*   Please report problems and suggestions to the author:
*
*                 dagostini@vaxrom.roma1.infn.it
*                 dagostini@vxdesy.desy.de            
*   
***********************************************************************

      integer i,j,k,l,l0,u,ncycle,msg1,msg2,IE,IM,OFFD
      logical force
      real nc1(SIZE_C)
      real denom, chi2, delta_chi2, Nobs, Ntrue, ptot
      real invdenom
      real Cov_M_kilj
C
C-----  additional variables used for speedup of Vc1 calculation
      real nc_inv_mc(SIZE_C),npec_inv(SIZE_E,SIZE_C)
      real neff_inv,M_tmp(SIZE_E,SIZE_E)
            
c      delta_chi2 = 0.1
      delta_chi2 = noc/100. 
      msg1 = 0
      msg2 = 0

      if (keep .ne. 1) then 
        print *,' ### Bayes unfolding called; last update: 5/7/95 '
        print *,'     (no substantial differences wrt to 1st release)'
      endif
      force = .false.  
      if ( nsteps .eq. 0) then 
        nsteps = 9999
      else if( nsteps .lt. 0) then 
        nsteps = - nsteps
        force = .true. 
      endif
      
      ier = 0
      if (noc .lt. 1 .or. noc .gt. SIZE_C) ier = 1
      if (noe .lt. 1 .or. noe .gt. SIZE_E) ier = 1

      if (ier .ne. 0) return

      if (keep .ne. 1) then 
        if (mode .eq. 1) then 
          print *,' mode = 1: the program uses pec(j,i) '
          do i =1,noc
            eff(i) = 0.
            do j =1,noe
              eff(i) = eff(i) +  pec(j,i)
            enddo
            if (eff(i) .gt. 1.001) then
              print *,' eff(',i,') = ',eff(i),' > 1'
              ier = 7
              return
            endif
            if (eff(i) .eq. 0.) then
              print *,' BAYES: Warning: Cause nr. ',i,' has no Effect '
              print *,'    (it may be OK, but it could be a mistake !)'
            endif
          enddo  ! loop on noc          
        else if(mode .eq. 2) then
          print *,' mode = 2: the program uses MC events  '
          print *,'           ne_mc(i) and nec_mc(j,i)    '
          print *,'           and propagates their errors'
          print *,'           to the results'
       
          do i =1,noc
            eff(i) = 0.
            do j =1,noe
              if (nc_mc(i) .gt. 0) then 
                 pec(j,i) = nec_mc(j,i)/nc_mc(i)
                 if (i.eq.41) then
                    print *,'j=',j,', nec_mc(j,41)=',nec_mc(j,41),
     &                 ',nc_mc(i)=',nc_mc(i)
                 end if
              else
                 pec(j,i) = 0
              endif
              eff(i) = eff(i) +  pec(j,i)
            enddo
            if (eff(i) .gt. 1.000001) then 
              print *,' nec_mc(,i) and nc_mc(i) are inconsistent ( i = '
     &             ,i,' )'
              print *, 'eff(i)=',eff(i),'= nec_mc/nc_mc'
              ier = 5
              return
            endif 
            if (eff(i) .eq. 0.) then
              print *,' BAYES: Warning: Cause nr. ',i,' has no Effect '
              print *,'    (it may be OK, but it could be a mistake !)'
            endif
          enddo  ! loop over noc            
        else     ! ( end of mode=2 )
          print *,' unkown mode = ',mode
          ier = 6
          return
        endif

        if (ier .ne. 0) return

        Nobs =  0.
        do j = 1,noe
          Nobs = Nobs + ne(j)
        enddo
        eff_0 = 0.
        ptot = 0.
        do i = 1,noc
          ptot = ptot + pc(i)
          eff_0 = eff_0 + eff(i)*pc(i)
          nc(i) = pc(i)*Nobs
        enddo
        if (ptot .lt. 0.999 .or. ptot .gt. 1.001) then 
          print *,'   Unnormalized initial probabilities: ptot = ',ptot
          ier = 8
          return
        else 
          eff_0 = eff_0/ptot
          if (eff_0 .le. 0) then 
            print *,' initial efficiency = ', eff_0
            ier = 9
            return
          endif
          do i = 1,noc
            nc(i) = pc(i)*Nobs/eff_0
          enddo
        endif      

        if (Nobs .eq. 0) then 
          print *,' BAYES: # of observed events = 0 ! '
          ier = 4
        endif
      endif  ! keep = 1    
      
      ncycle = 0
    1 continue      
      ncycle = ncycle + 1

      do j=1,noe
        denom=0.
        do i=1,noc
          denom=denom+pec(j,i)*pc(i)
        enddo
        if (denom .eq. 0.) then
          invdenom = 0.
        else 
          invdenom = 1./denom
        endif
        do i=1,noc
          pce(i,j) = pec(j,i)*pc(i)*invdenom
        enddo
      enddo
      Ntrue = 0.
      do i=1,noc
        nc1(i) = 0.
        do j=1,noe
          if (eff(i) .ne. 0) then 
            M_unf(i,j) = pce(i,j)/eff(i)
          else 
            M_unf(i,j) = 0.
          endif
          nc1(i) = nc1(i) + ne(j)* M_unf(i,j)
        enddo
        Ntrue = Ntrue + nc1(i)
      enddo
***  compare the initial and the final distribution
      chi2 = 0.
      do i =1,noc
        if (nc(i)+nc1(i) .gt. 1) then
          chi2 = chi2+ (nc1(i)-nc(i))**2/(nc(i)+nc1(i)) 
        else 
          chi2 = chi2+ (nc1(i)-nc(i))**2 
        endif
      enddo
      do i = 1,noc
        nc(i) = nc1(i)
        pc(i) = nc1(i)/Ntrue
      enddo
      if ( ncycle .ge. nsteps ) goto 777
      if ( chi2 .lt. delta_chi2 .and. .not. force) goto 777
      goto 1

  777 nsteps = ncycle
 
      eff_true = Nobs/Ntrue    

      print *, ' ## Bayes unfold terminated ## '
      print *, '     Number of steps           = ',nsteps
      print *, '     Observed number of events = ',Nobs 
      print *, '     Unfolded number of events = ',Ntrue
      print *, '     Initial efficiency        = ',eff_0
      print *, '     final efficiency          = ',eff_true
       
c-------------  Error evaluation -----------
      OFFD = isign(1, er_mode)
      IE = iabs(er_mode)/10
      IM = mod(iabs(er_mode),10)
C  1. errors from observed number of events 
      do k = 1,noc
        do l = 1,noc
          Vc0(k,l) = 0.
        enddo
      enddo
      if(IE .ne. 0) then  
        do k = 1,noc
          l0 = 1
          if (OFFD .lt. 0) l0 = k          
          do l = l0,k
            do i =1,noe
              do j =1,noe
                if (i .eq. j) then 
                  if (IE .eq. 1) then 
                    Vc0(k,l) = Vc0(k,l) + M_unf(k,i)*M_unf(l,j)*ne(j)
                  else
                    Vc0(k,l) = Vc0(k,l) + M_unf(k,i)*M_unf(l,j)*
     &                                    ne(j) * (1. - ne(j)/Ntrue)  
                  endif
                else  ! (i .ne. i)
                  if (IE .eq. 2) then
                    Vc0(k,l) = Vc0(k,l) - M_unf(k,i)*M_unf(l,j)*
     &                                    ne(i) * ne(j)/Ntrue  
                  endif
                endif
              enddo
            enddo
            Vc0(l,k) = Vc0(k,l)
          enddo
        enddo
      endif

C  2. errors from smearing matrix ( only if mode = 2 )    
      do k = 1,noc
        do l = 1,noc
          Vc1(k,l) = 0.
        enddo
      enddo
      if (mode .eq. 2 .and. IM .ne. 0) then

C----- calculate first some temporary quantities to speed up
C      processing inside the main 4-fold nested do-loop
        do k = 1,noc
          if (nc_mc(k).ne.0) then
            nc_inv_mc(k)=1./nc_mc(k)
          else
            nc_inv_mc(k)=0.
          endif
          do i=1,noe
            if (pec(i,k).ne.0) then
              npec_inv(i,k)=nc_inv_mc(k)/pec(i,k)
            else
              npec_inv(i,k)=0.
            endif
          enddo
        enddo
        do i=1,noe
          M_tmp(i,i)=0
          do u=1,noc
            if (IM.eq.2) then
              M_tmp(i,i)=M_tmp(i,i)+M_unf(u,i)*M_unf(u,i)*eff(u)*eff(u)
     +                                *(npec_inv(i,u)-nc_inv_mc(u))
            else
C           (Poisson approximation)
              M_tmp(i,i)=M_tmp(i,i)+M_unf(u,i)*M_unf(u,i)*eff(u)*eff(u)
     +                                *npec_inv(i,u)
            endif
          enddo
          do j=i+1,noe
            M_tmp(i,j)=0
            if (IM.eq.2) then
              do u=1,noc
                M_tmp(i,j)=M_tmp(i,j)-M_unf(u,i)*M_unf(u,j)
     +                                *eff(u)*eff(u)*nc_inv_mc(u)
              enddo
            endif
            M_tmp(j,i)=M_tmp(i,j)
          enddo
        enddo
C-----
c
        do k = 1,noc
          if (eff(k).ne.0) then
            neff_inv=nc_inv_mc(k)/eff(k)
          else
            neff_inv=0.
          endif
          l0 = 1
          if (OFFD .lt. 0) l0 = k
          do l = l0,k
            do i = 1,noe
              do j = 1,noe
C----------   covariance matrix of unfolding matrix
                Cov_M_kilj = M_unf(l,i)*nc_inv_mc(l)
     +                     + M_unf(k,j)*nc_inv_mc(k)
     +                     + M_tmp(i,j)
                if (k.eq.l) then
                  Cov_M_kilj=Cov_M_kilj-neff_inv
                  if (i.eq.j) then
                    Cov_M_kilj=Cov_M_kilj+npec_inv(i,k)
                  endif
                endif
                if (i.eq.j) then
                  Cov_M_kilj=Cov_M_kilj-M_unf(l,i)*eff(l)*npec_inv(i,l)
     +                                 -M_unf(k,i)*eff(k)*npec_inv(i,k)
                endif
                Cov_M_kilj = Cov_M_kilj*M_unf(k,i)*M_unf(l,j)

C----------   covariance matrix of results 
                Vc1(k,l) = Vc1(k,l) + ne(i)*ne(j)*Cov_M_kilj 
              enddo ! j
            enddo ! i
            Vc1(l,k) = Vc1(k,l)
          enddo ! l
        enddo ! k
      endif

********  covariance matrices   

      do k = 1,noc
        do l = 1,noc
c adding the normalization uncertainty to the contribution from 
c the observed data 
          if (IE .eq. 2) then 
            Vc0_u(k,l) = Vc0(k,l) + nc(k)*nc(l)/Nobs         
          else 
            Vc0_u(k,l) = Vc0(k,l)
          endif 
c total COVARIANCE MATRIX (with constraint on the normalization) 
          Vc(k,l) = Vc0(k,l) + Vc1(k,l)
c total unconstrained COVARIANCE MATRIX 
          if (IE .eq. 2) then
            Vc_u(k,l) = Vc(k,l) + nc(k)*nc(l)/Nobs
          else 
            Vc_u(k,l) = Vc(k,l)
          endif 
c covariance matrices of final probabilities
          Vpc(k,l) = Vc(k,l)/Ntrue**2
          Vpc0(k,l) = Vc0(k,l)/Ntrue**2
          Vpc1(k,l) = Vc1(k,l)/Ntrue**2
        enddo  
      enddo

      print *,' ### Exit from Bayes unfolding '
      return
      end
