PROGRAM Two_Filter_Cross_Correlation

  IMPLICIT NONE

  ! Parameters in terms of decay rate gamma
  ! Atom decay rate
  REAL, PARAMETER :: gamma = 1.0
  ! Drive strength (Rabi Frequency)
  REAL, PARAMETER :: omega = 40.0
  ! Drive detuning from two-photon resonance, \omega_{d} - \omega_{gf}
  REAL, PARAMETER :: delta = 0.0
  ! Drive strength ratio of the two levels, |g> <-> |e> and |e> <-> |f>
  REAL, PARAMETER :: xi = 1.0
  ! Difference between two energy levels, \omega_{ef} - \omega_{ge}
  REAL, PARAMETER :: alpha = -120.0

  ! Eigenfrequencies for position of spectrum peaks for delta = 0
  ! eigen drive strength
  REAL, PARAMETER :: Omega_t = SQRT(((0.25*alpha)**2) + ((0.5*omega) &
                                     & ** 2) * (1 + (xi ** 2)))
  ! Positive eigenfrequency
  REAL, PARAMETER :: wp = -(0.25 * alpha) + Omega_t
  ! Negative eigenfrequency
  REAL, PARAMETER :: wm = -(0.25 * alpha) - Omega_t

  ! Filter parameter stuff
  ! Detuning of cavity "a" resonance frequency with drive frequency
  ! \Delta_{f} = \omega_{0} - \omega_{d}. \Delta_{f} = 0 is resonant with
  ! \omega_{gf} / 2 if \delta = 0.
  REAL, PARAMETER :: D_a = -wm
  ! Detuning of cavity "b" resonance frequency with drive frequency
  REAL, PARAMETER :: D_b = wm
  ! Cavity linewidth/transmission of cavity a
  REAL, PARAMETER :: kappa_a = 10.0
  ! Cavity linewidth/transmission of cavity b
  REAL, PARAMETER :: kappa_b = 10.0

  ! Quantum object stuff
  ! Hamiltonian from Wallraff et. al
  COMPLEX, DIMENSION(3,3) :: H
  ! Non-Hermitian Hamiltonian for continuous evolution
  ! Total Hamiltonian H_full + H_nH
  COMPLEX, DIMENSION(3,3) :: H_full
  ! Raising and Lowering Operators
  COMPLEX, DIMENSION(3,3) :: sigmap, sigmam
  ! Lowering operator for cascade system for cavity a
  ! -i * SQRT(\gamma\kappa_a) * \sqrt{na} * \Sigma_{-}
  COMPLEX, DIMENSION(3,3) :: cassig_a
  ! Lowering operator for cascade system for cavity b
  ! -i * SQRT(\gamma\kappa_b) * \sqrt{nb} * \Sigma_{-}
  COMPLEX, DIMENSION(3,3) :: cassig_b
  ! sigmap * sigmam
  COMPLEX, DIMENSION(3,3) :: sigmapm
  ! Hilbert Space - max number of photons in cavity 0 -> N
  INTEGER, PARAMETER :: N = 5
  ! Hilbert truncated dimension. 3 atomic states x (N+1)^2 Fock states for
  ! two cavities.
  INTEGER, PARAMETER :: d_Hilb = INT(((N + 1) ** 2) * 3)
  ! Dynamical state vector |\psi> = |Na> \otimes |Nb> \otimes |i> where |Na> is
  ! the Fock basis for cavity a,|Nb> is the Fock basis for cavity b, and
  ! i=g,e,f are the atomic states |g>, |e>, |f>.
  ! psi = |Na>|Nb>|i> = |0,0,g> + |0,0,e> + |0,0,f> + |0,1,g> + ...+|1,1,f> +...
  COMPLEX, DIMENSION(d_Hilb) :: psi
  ! A copy of the state vector for calculations
  COMPLEX, DIMENSION(d_Hilb) :: psi_clone
  ! A copy of the state vector for a cavity jump
  COMPLEX, DIMENSION(d_Hilb) :: psi_cav

  ! Time stuff
  ! Time step
  REAL, PARAMETER :: dt = 0.001
  ! Max time in units of \gamma \tau
  REAL, PARAMETER :: max_time = 100000
  ! Total steps
  INTEGER, PARAMETER :: steps = INT(max_time / dt, KIND=8)
  ! Max time for tau calculations
  REAL, PARAMETER :: tau = 10
  ! Max number of time steps for tau calculations
  INTEGER, PARAMETER :: tau_max = INT(tau / dt)


  ! Useful stuff
  ! Time counter
  INTEGER :: k
  ! Integer counter for atomic states m = {1,2,3} = {g,e,f}
  INTEGER :: m
  ! Integer counters for a photon number and b photon number
  INTEGER :: na, nb
  ! Other integer counters
  INTEGER :: j, l
  ! Integer place for photon number in psi vector (na + 1) * (nb + 1)
  INTEGER :: nplace, nplace_pm1
  ! Complex i = SQRT(-1)
  COMPLEX, PARAMETER :: i = CMPLX(0,1)
  ! Runge-Kutta 4th Order Vectors/Matrices
  COMPLEX, DIMENSION(d_Hilb) :: k1, k2, k3, k4
  ! 1 / 6 as a parameter to be calculated only once
  REAL, PARAMETER :: xis = 1.0 / 6.0
  ! Random number to be called fromm RANDOM_NUMBER()
  REAL :: rand
  ! Trace for normalisation
  REAL :: trace
  ! Square root of the trace
  REAL :: tracesq
  ! Temporal variables
  COMPLEX :: temp, temp1, temp2
  ! Temportal integer
  INTEGER :: tempint
  ! Jump counters
  INTEGER :: count_r_a, count_t_a, count_r_b, count_t_b

  ! Probability stuff
  ! Probability to jump for cavity a
  REAL :: prob_t_a, prob_r_a
  ! Probability to jump for cavity b
  REAL :: prob_t_b, prob_r_b
  ! Probability for an atomic decay gamma * dt * <\Sigma^{\dagger}\Sigma>
  REAL :: prob_atom_decay
  ! Total probability
  REAL :: prob_total
  ! Probability Array
  REAL, DIMENSION(4) :: P

  ! Data stuff
  ! State probabilities of the atom (|g>, |e>, |f>)
  REAL :: p_gg, p_ee, p_ff
  ! Photon number inside cavity a and b
  REAL :: photon_a, photon_b
  ! Steady state mean photon in cavity a and b
  REAL :: mean_photon_a, mean_photon_b
  ! Check to start writing correlation
  LOGICAL :: corr_write
  ! Correlation calculation
  REAL :: corr
  ! Integer for correlation array
  INTEGER :: tau_counter
  ! Recording counters
  INTEGER :: number_of_recordings
  ! Array to store correlation
  REAL, DIMENSION(0:tau_max) :: correlation
  ! Array to store cross correlation
  REAL, DIMENSION(0:tau_max) :: cross_correlation
  ! File name for saving correlation and cross correlation data
  CHARACTER(LEN=28) :: filename_correlation = "./data_files/correlation.txt"
  ! Filename for saving parameters and time
  CHARACTER(LEN=20) :: filename_time = "./data_files/tau.txt"

  ! Initialising matrices
  ! Raising and Lowering Operators
  sigmam = 0
  sigmam(1,2) = 1.0
  sigmam(2,3) = xi

  sigmap = 0
  sigmap(2,1) = sigmam(1,2)
  sigmap(3,2) = sigmam(2,3)

  ! S_+ x S_-
  sigmapm = 0
  sigmapm = MATMUL(sigmap, sigmam)

  ! Full Hamiltonian
  H = 0
  H(1,2) = 0.5 * omega
  H(2,1) = 0.5 * omega
  H(2,2) = -(0.5 * alpha + delta)
  H(2,3) = 0.5 * omega * xi
  H(3,2) = 0.5 * omega * xi
  H(3,3) = -2 * delta

  ! Total Hamiltonian for continuous evolution
  H = H - i * 0.5 * gamma * sigmapm

  ! Initialise data arrays
  correlation = 0
  cross_correlation = 0
  mean_photon_a = 0
  mean_photon_b = 0

  ! Initialise jump counters
  count_t_a = 0
  count_r_a = 0
  count_t_b = 0
  count_r_b = 0
  number_of_recordings = 0
  tau_counter = 0

  ! Initialise the seed for random number generation
  CALL init_random_seed

  ! Intialise Runge-Kutta vectors
  k1 = 0
  k2 = 0
  k3 = 0
  k4 = 0

  ! Atom is initial in the ground state |psi(0)> = |0,0,g> = (1,0,0,...,0)
  psi = 0.0
  psi(1) = 1.0
!##############################################################################!
!                     Section A: Start trajectory loop                         !
!##############################################################################!
  DO k=0,steps
    !##########################################################!
    !    Calculate mean photon number inside cavities          !
    !##########################################################!
    ! Mean photon number for cavity a
    photon_a = 0.0
    DO na=0,N
      DO nb=0,N
        ! photon number index
        nplace = (3 * (N + 1) * na) + (3 * nb)
        photon_a = photon_a + (na * (ABS(psi(nplace + 1)) ** 2))
        photon_a = photon_a + (na * (ABS(psi(nplace + 2)) ** 2))
        photon_a = photon_a + (na * (ABS(psi(nplace + 3)) ** 2))
      END DO
    END DO
    ! Steady state mean photon number in cavity a.
    ! Will be avaraged by number of time steps at the end of the loop.
    mean_photon_a = mean_photon_a + photon_a

    ! Mean photon number for cavity b
    photon_b = 0.0
    DO na=0,N
      DO nb=0,N
        ! photon number index
        nplace = (3 * (N + 1) * na) + (3 * nb)
        photon_b = photon_b + (nb * (ABS(psi(nplace + 1)) ** 2))
        photon_b = photon_b + (nb * (ABS(psi(nplace + 2)) ** 2))
        photon_b = photon_b + (nb * (ABS(psi(nplace + 3)) ** 2))
      END DO
    END DO
    ! Steady state mean photon number in cavity b.
    ! Will be avaraged by number of time steps at the end of the loop.
    mean_photon_b = mean_photon_b + photon_b
    !##########################################################!
    !    Calculate state probabilities |g>,|e>,|f>             !
    !##########################################################!
    p_gg = 0
    p_ee = 0
    p_ff = 0
    DO na=0,N
      DO nb=0,N
        ! photon number index
        nplace = (3 * (N + 1) * na) + (3 * nb)
        ! Probability to be in |g>
        p_gg = p_gg + ABS(psi(nplace + 1) ** 2)
        ! Probability to be in |e>
        p_ee = p_ee + ABS(psi(nplace + 2) ** 2)
        ! Probability to be in |f>
        p_ff = p_ff + ABS(psi(nplace + 3) ** 2)
      END DO
    END DO
!##############################################################################!
!                  Section B: Write correlation to array                       !
!##############################################################################!

    !###############################!
    !   Record correlation Values   !
    !###############################!

    ! As the transmission detector clicks, set that detection to tau=0 then
    ! record the probaility for a transmission jump to occur and write to the
    ! corr_array array for tau_max time steps.

    ! Set counter for correlation array. When a transmission jump occurs, JUMP
    ! is set to 1. We then set the corr_write value to 1 so we know to start
    ! recording values. If another transmission jump is detected while we are
    ! recording results it will have no effect on the current recording.
    ! IF (JUMP == 1) THEN
    !   corr_write = 1
    ! END IF

    ! If tau_counter increases beyond tau_max, we know we have recorded the
    ! correlation for as long as we planned to so we can reset the tau_counter
    ! and wait for another jump to occur.
    IF(tau_counter > tau_max) THEN
      corr_write = .FALSE.
      tau_counter = 0
    END IF

    ! If corr_write is set to one we know that a transmission jump has occured
    ! and we can start recording
    IF (corr_write .EQV. .TRUE.) THEN
      ! Update correlation of cavity a
      correlation(tau_counter) = correlation(tau_counter) + photon_a
      ! Update cross correlation of cavity b conditional on a detection from
      ! cavity a
      cross_correlation(tau_counter) = cross_correlation(tau_counter) + photon_b
      ! Now that a time step has been recorded, increase the tau_counter index
      ! by one time step.
      tau_counter = tau_counter + 1
    END IF
!##############################################################################!
!                        Calculate probabilities                               !
!##############################################################################!
    !#########################################!
    !      Atom decay jump probability        !
    !#########################################!
    ! \Sigma = |g><e| + \xi |e><f|
    ! Probability for an atom decay to occur is
    ! P = dt * gamma * \Sigma^{\dagger} + \Sigma
    !   = dt * gamma * (< |e><e| > + \xi^{2} < |f><f| >)
    prob_atom_decay = 0.0
    temp1 = 0.0
    temp2 = 0.0
    DO na=0,N
      DO nb=0,N
        ! photon number index
        nplace = (3 * (N + 1) * na) + (3 * nb)
        ! Probability to be in |e>
        temp1 = temp1 + ABS(psi(nplace + 2) ** 2)
        ! Probability to be in |f> \times \xi^{2}
        temp2 = temp2 + ((xi ** 2) * ABS(psi(nplace + 3) ** 2))
      END DO
    END DO
    prob_atom_decay = gamma * dt * REAL(temp1 + temp2)
    ! prob_atom_decay = gamma * dt * (p_ee + ((xi ** 2) * p_ff))

    !#########################################!
    !      Cavity a Probabilities             !
    !#########################################!
    ! Probability for a transmission jump through cavity a
    ! prob_t_a = kappa_a * dt * <\psi|a^{\dagger}a|\psi>
    prob_t_a = 0.0
    prob_t_a = kappa_a * dt * photon_a

    ! Calculate probability for an entangled jump to occur, < c^{\dagger}c >
    ! where c = \sqrt{\gamma}\Sigma + \sqrt{\kappa_a}a. The probability is
    ! then \gamma <\Sigma_{+}\Sigma_{-}> + \kappa_a <a^{\dagger}a> +
    ! \sqrt{\gamma\kappa_a}< \Sigma_{+}a > +
    ! \sqrt{\gamma\kappa_a}< \Sigma_{-}a^{\dagger} >. These last two terms
    ! are complex conjugates of eachother.

    ! First we add the probability of an atom decay and a transmission jump
    ! together
    prob_r_a = (0.5 * prob_atom_decay) + prob_t_a

    ! Now we calculate the imaginary parts
    ! dt * \sqrt{\gamma\kappa} <\psi| a \Sigma^{\dagger} | \psi >
    ! <\psi| a = (a^{\dagger}|\psi>)*
    !          = \sqrt{na} <na-1,nb,m|
    ! \Sigma^{\dagger} |\psi> = |na,nb,g> + \xi|na,nb,e>
    temp = 0.0
    DO na=1,N
      DO nb=0,N
        ! photon number index
        nplace = (3 * (N + 1) * na) + (3 * nb)
        ! photon number index for na - 1
        nplace_pm1 = (3 * (N + 1) * (na - 1)) + (3 * nb)
        ! <\Sigma^{\dagger} a>
        temp = temp + SQRT(1.0 * na) * &
             & CONJG(psi(nplace_pm1 + 2)) * psi(nplace + 1)
        ! <\Sigma a^{\dagger}>
        temp = temp + SQRT(1.0 * na) * xi * &
             & CONJG(psi(nplace_pm1 + 3)) * psi(nplace + 2)
      END DO
    END DO
    ! Update probability
    prob_r_a = prob_r_a + (dt * SQRT(gamma * kappa_a) * SQRT(2.0) * REAL(temp))
    !#########################################!
    !      Cavity b Probabilities             !
    !#########################################!
    ! Probability for a transmission jump through cavity b
    ! prob_t_b = kappa_b * dt * <\psi|b^{\dagger}b|\psi>
    prob_t_b = 0.0
    prob_t_b = kappa_b * dt * photon_b

    ! Same as before but for the b cavity operator
    ! First we add the probability of an atom decay and a transmission jump
    ! together
    prob_r_b = 0.0
    prob_r_b = (0.5 * prob_atom_decay) + prob_t_b

    ! Now we calculate the imaginary parts
    temp = 0.0
    DO na=0,N
      DO nb=1,N
        ! photon number index
        nplace = (3 * (N + 1) * na) + (3 * nb)
        ! photon number index for na - 1
        nplace_pm1 = (3 * (N + 1) * na) + (3 * (nb - 1))
        ! <\Sigma^{\dagger} a>
        temp = temp + SQRT(1.0 * nb) * &
             & CONJG(psi(nplace_pm1 + 2)) * psi(nplace + 1)
        ! <\Sigma a^{\dagger}>
        temp = temp + SQRT(1.0 * nb) * xi * &
             & CONJG(psi(nplace_pm1 + 3)) * psi(nplace + 2)
      END DO
    END DO
    ! Update probability
    prob_r_b = prob_r_b + (dt * SQRT(gamma * kappa_b) * SQRT(2.0) * REAL(temp))

    ! Add together for total probability
    prob_total = (prob_r_a + prob_t_a) + (prob_r_b + prob_t_b)
!##############################################################################!
!                      Section C: Evolve the system                            !
!   If rand > prob then evolve the system continuously with the Hamiltonian    !
!         using the Runge-Kutta 4th Order. Otherwise a jump occurs             !
!##############################################################################!
    ! Call random number
    CALL RANDOM_NUMBER(rand)

    ! If rand > prob then evolve the system continuously with the Hamiltonian
    IF (rand >= prob_total) THEN
      ! Evolve using Runge-Kutta 4th Order
      psi_clone = psi
      ! First loop (na) is for photon number in cavity a |0>, |1>, ... |N>
      ! Second loop (nb) is for photon number in cavity b |0>, |1>, ... |N>
      ! Third loop (m) is for atomic state 1=|g>, 2=|e>, 3=|f>

      !###################!
      !    Calculate k1   !
      !###################!
      DO na=0,N
        DO nb=0,N
          ! photon number index
          nplace = (3 * (N + 1) * na) + (3 * nb)
          ! Set the Hamiltonian to include cavity terms
          H_full = H
          H_full(1,1) = H_full(1,1) + (D_a - i * kappa_a) * na &
                                  & + (D_b - i * kappa_b) * nb
          H_full(2,2) = H_full(2,2) + (D_a - i * kappa_a) * na &
                                  & + (D_b - i * kappa_b) * nb
          H_full(3,3) = H_full(3,3) + (D_a - i * kappa_a) * na &
                                  & + (D_b - i * kappa_b) * nb
          ! Hamiltonian matrix for cascade system
          ! cavity a
          cassig_a = 0
          cassig_a = -i * SQRT(0.5 * gamma * kappa_a) * SQRT(1.0 * na) * sigmam
          ! cavity b
          cassig_b = 0
          cassig_b = -i * SQRT(0.5 * gamma * kappa_b) * SQRT(1.0 * nb) * sigmam
          ! Matrix multiplication
          DO m=1,3
            k1(nplace + m) = -i * dt * H_full(m,1) * psi_clone(nplace + 1) &
                            &-i * dt * H_full(m,2) * psi_clone(nplace + 2) &
                            &-i * dt * H_full(m,3) * psi_clone(nplace + 3)
            ! Only calculate the next part for photon number na > 0 as the
            ! a^{\dagger} operator couples the |N-1> state to the |N> state.
            temp1 = 0
            temp2 = 0
            IF (na /= 0) THEN
              ! cavity a
              ! photon number index for (na - 1)
              nplace_pm1 = (3 * (N + 1) * (na - 1)) + (3 * nb)
              temp1 = -i * dt * cassig_a(m,1) * psi_clone(nplace_pm1 + 1) &
                    & -i * dt * cassig_a(m,2) * psi_clone(nplace_pm1 + 2) &
                    & -i * dt * cassig_a(m,3) * psi_clone(nplace_pm1 + 3)
            END IF
            IF (nb /= 0) THEN
               ! cavity b
               ! photon number index for (nb - 1)
               nplace_pm1 = (3 * (N + 1) * na) + (3 * (nb - 1))
               temp2 = -i * dt * cassig_b(m,1) * psi_clone(nplace_pm1 + 1) &
                     & -i * dt * cassig_b(m,2) * psi_clone(nplace_pm1 + 2) &
                     & -i * dt * cassig_b(m,3) * psi_clone(nplace_pm1 + 3)
            END IF
            ! Update k vector
            k1(nplace + m) = k1(nplace + m) + (temp1 + temp2)
          END DO
        END DO
      END DO

      !###################!
      !    Calculate k2   !
      !###################!
      DO na=0,N
        DO nb=0,N
          ! photon number index
          nplace = (3 * (N + 1) * na) + (3 * nb)
          ! Set the Hamiltonian to include cavity terms
          H_full = H
          H_full(1,1) = H_full(1,1) + (D_a - i * kappa_a) * na &
                                  & + (D_b - i * kappa_b) * nb
          H_full(2,2) = H_full(2,2) + (D_a - i * kappa_a) * na &
                                  & + (D_b - i * kappa_b) * nb
          H_full(3,3) = H_full(3,3) + (D_a - i * kappa_a) * na &
                                  & + (D_b - i * kappa_b) * nb
          ! Hamiltonian matrix for cascade system
          ! cavity a
          cassig_a = 0
          cassig_a = -i * SQRT(0.5 * gamma * kappa_a) * SQRT(1.0 * na) * sigmam
          ! cavity b
          cassig_b = 0
          cassig_b = -i * SQRT(0.5 * gamma * kappa_b) * SQRT(1.0 * nb) * sigmam
          ! Matrix multiplication
          DO m=1,3
            k2(nplace + m) = -i * dt * H_full(m,1) * &
                           & (psi_clone(nplace + 1) + 0.5 * k1(nplace + 1))&
                           & -i * dt * H_full(m,2) * &
                           & (psi_clone(nplace + 2) + 0.5 * k1(nplace + 2))&
                           & -i * dt * H_full(m,3) * &
                           & (psi_clone(nplace + 3) + 0.5 * k1(nplace + 3))
           ! Only calculate the next part for photon number l > 0 as the
           ! a^{\dagger} operator couples the |N-1> state to the |N> state.
            temp1 = 0
            IF (na /= 0) THEN
              ! cavity a
              ! photon number index for (na - 1)
              nplace_pm1 = (3 * (N + 1) * (na - 1)) + (3 * nb)
              temp1 = -i * dt * cassig_a(m,1) * &
                    & (psi_clone(nplace_pm1 + 1) + 0.5 * k1(nplace_pm1 + 1))&
                    & -i * dt * cassig_a(m,2) * &
                    & (psi_clone(nplace_pm1 + 2) + 0.5 * k1(nplace_pm1 + 2))&
                    & -i * dt * cassig_a(m,3) * &
                    & (psi_clone(nplace_pm1 + 3) + 0.5 * k1(nplace_pm1 + 3))
             END IF
             temp2 = 0
             IF (nb /= 0) THEN
                ! cavity b
                ! photon number index for (nb - 1)
                nplace_pm1 = (3 * (N + 1) * na) + (3 * (nb - 1))
                temp2 = -i * dt * cassig_b(m,1) * &
                      & (psi_clone(nplace_pm1 + 1) + 0.5 * k1(nplace_pm1 + 1)) &
                      & -i * dt * cassig_b(m,2) * &
                      & (psi_clone(nplace_pm1 + 2) + 0.5 * k1(nplace_pm1 + 2)) &
                      & -i * dt * cassig_b(m,3) * &
                      & (psi_clone(nplace_pm1 + 3) + 0.5 * k1(nplace_pm1 + 3))
              END IF
            ! Update k vector
            k2(nplace + m) = k2(nplace + m) + (temp1 + temp2)
          END DO
        END DO
      END DO

      !###################!
      !    Calculate k3   !
      !###################!
      DO na=0,N
        DO nb=0,N
          ! photon number index
          nplace = (3 * (N + 1) * na) + (3 * nb)
          ! Set the Hamiltonian to include cavity terms
          H_full = H
          H_full(1,1) = H_full(1,1) + (D_a - i * kappa_a) * na &
                                  & + (D_b - i * kappa_b) * nb
          H_full(2,2) = H_full(2,2) + (D_a - i * kappa_a) * na &
                                  & + (D_b - i * kappa_b) * nb
          H_full(3,3) = H_full(3,3) + (D_a - i * kappa_a) * na &
                                  & + (D_b - i * kappa_b) * nb
          ! Hamiltonian matrix for cascade system
          ! cavity a
          cassig_a = 0
          cassig_a = -i * SQRT(0.5 * gamma * kappa_a) * SQRT(1.0 * na) * sigmam
          ! cavity b
          cassig_b = 0
          cassig_b = -i * SQRT(0.5 * gamma * kappa_b) * SQRT(1.0 * nb) * sigmam
          ! Matrix multiplication
          DO m=1,3
            k3(nplace + m) = -i * dt * H_full(m,1) * &
                           & (psi_clone(nplace + 1) + 0.5 * k2(nplace + 1))&
                           & -i * dt * H_full(m,2) * &
                           & (psi_clone(nplace + 2) + 0.5 * k2(nplace + 2))&
                           & -i * dt * H_full(m,3) * &
                           & (psi_clone(nplace + 3) + 0.5 * k2(nplace + 3))
           ! Only calculate the next part for photon number l > 0 as the
           ! a^{\dagger} operator couples the |N-1> state to the |N> state.
            temp1 = 0
            IF (na /= 0) THEN
              ! cavity a
              ! photon number index for (na - 1)
              nplace_pm1 = (3 * (N + 1) * (na - 1)) + (3 * nb)
              temp1 = -i * dt * cassig_a(m,1) * &
                    & (psi_clone(nplace_pm1 + 1) + 0.5 * k2(nplace_pm1 + 1))&
                    & -i * dt * cassig_a(m,2) * &
                    & (psi_clone(nplace_pm1 + 2) + 0.5 * k2(nplace_pm1 + 2))&
                    & -i * dt * cassig_a(m,3) * &
                    & (psi_clone(nplace_pm1 + 3) + 0.5 * k2(nplace_pm1 + 3))
            END IF
            temp2 = 0
            IF (nb /= 0) THEN
               ! cavity b
               ! photon number index for (nb - 1)
               nplace_pm1 = (3 * (N + 1) * na) + (3 * (nb - 1))
               temp2 = -i * dt * cassig_b(m,1) * &
                     & (psi_clone(nplace_pm1 + 1) + 0.5 * k2(nplace_pm1 + 1)) &
                     & -i * dt * cassig_b(m,2) * &
                     & (psi_clone(nplace_pm1 + 2) + 0.5 * k2(nplace_pm1 + 2)) &
                     & -i * dt * cassig_b(m,3) * &
                     & (psi_clone(nplace_pm1 + 3) + 0.5 * k2(nplace_pm1 + 3))
             END IF
            ! Update k vector
            k3(nplace + m) = k3(nplace + m) + (temp1 + temp2)
          END DO
        END DO
      END DO

      !###################!
      !    Calculate k4   !
      !###################!
      DO na=0,N
        DO nb=0,N
          ! photon number index
          nplace = (3 * (N + 1) * na) + (3 * nb)
          ! Set the Hamiltonian to include cavity terms
          H_full = H
          H_full(1,1) = H_full(1,1) + (D_a - i * kappa_a) * na &
                                  & + (D_b - i * kappa_b) * nb
          H_full(2,2) = H_full(2,2) + (D_a - i * kappa_a) * na &
                                  & + (D_b - i * kappa_b) * nb
          H_full(3,3) = H_full(3,3) + (D_a - i * kappa_a) * na &
                                  & + (D_b - i * kappa_b) * nb
          ! Hamiltonian matrix for cascade system
          ! cavity a
          cassig_a = 0
          cassig_a = -i * SQRT(0.5 * gamma * kappa_a) * SQRT(1.0 * na) * sigmam
          ! cavity b
          cassig_b = 0
          cassig_b = -i * SQRT(0.5 * gamma * kappa_b) * SQRT(1.0 * nb) * sigmam
          ! Matrix multiplication
          DO m=1,3
            k4(nplace + m) = -i * dt * H_full(m,1) * &
                           & (psi_clone(nplace + 1) + k3(nplace + 1))&
                           & -i * dt * H_full(m,2) * &
                           & (psi_clone(nplace + 2) + k3(nplace + 2))&
                           & -i * dt * H_full(m,3) * &
                           & (psi_clone(nplace + 3) + k3(nplace + 3))
           ! Only calculate the next part for photon number l > 0 as the
           ! a^{\dagger} operator couples the |N-1> state to the |N> state.
            temp1 = 0
            IF (na /= 0) THEN
              ! cavity a
              ! photon number index for (na - 1)
              nplace_pm1 = (3 * (N + 1) * (na - 1)) + (3 * nb)
              temp1 = -i * dt * cassig_a(m,1) * &
                    & (psi_clone(nplace_pm1 + 1) + k3(nplace_pm1 + 1))&
                    & -i * dt * cassig_a(m,2) * &
                    & (psi_clone(nplace_pm1 + 2) + k3(nplace_pm1 + 2))&
                    & -i * dt * cassig_a(m,3) * &
                    & (psi_clone(nplace_pm1 + 3) + k3(nplace_pm1 + 3))
            END IF
            temp2 = 0
            IF (nb /= 0) THEN
               ! cavity b
               ! photon number index for (nb - 1)
               nplace_pm1 = (3 * (N + 1) * na) + (3 * (nb - 1))
               temp2 = -i * dt * cassig_b(m,1) * &
                     & (psi_clone(nplace_pm1 + 1) + k3(nplace_pm1 + 1)) &
                     & -i * dt * cassig_b(m,2) * &
                     & (psi_clone(nplace_pm1 + 2) + k3(nplace_pm1 + 2)) &
                     & -i * dt * cassig_b(m,3) * &
                     & (psi_clone(nplace_pm1 + 3) + k3(nplace_pm1 + 3))
             END IF
            ! Update k vector
            k4(nplace + m) = k4(nplace + m) + (temp1 + temp2)
          END DO
        END DO
      END DO

      psi_clone = psi_clone + xis * (k1 + 2.0*(k2 + k3) + k4)

      ! Normalise the state
      trace = 0
      DO l=1,d_Hilb
        trace = trace + (ABS(psi_clone(l))) ** 2
      END DO
      tracesq = SQRT(trace)
      psi = psi_clone / tracesq

    ELSE IF (rand <= prob_total) THEN
!##############################################################################!
!                      Section D: A Jump Occurs!!!                             !
! First we sort all the probabilities into an array normalised by the total    !
! probabilty. A new random number is called and compared to the values inside  !
! the probability array. Depending on where it falls in the "steps" will decide!
! which jump occurs.                                                           !
!##############################################################################!
      CALL RANDOM_NUMBER(rand)

      ! Sort jump probabilities into a step array
      P = 0.0
      P(1) = prob_t_a
      P(2) = P(1) + prob_r_a
      P(3) = P(2) + prob_t_b
      P(4) = P(3) + prob_r_b

      ! Normalise by total probability
      P = P / prob_total

      ! Compare probabilities with a new random number
      IF (rand <= P(1)) THEN
        !############################################!
        !     A cavity emission from cavity a        !
        !############################################!
        count_t_a = count_t_a + 1
        ! Set logical for correlation.
        IF (tau_counter == 0) THEN
          corr_write = .TRUE.
          number_of_recordings = number_of_recordings + 1
        END IF
        ! Calculate the state for a cavity decay from a
        psi_cav = 0.0
        DO na=0,N-1
          DO nb=0,N
            ! photon number index
            nplace = (3 * (N + 1) * na) + (3 * nb)
            ! photon number index for (na + 1)
            nplace_pm1 = (3 * (N + 1) * (na + 1)) + (3 * nb)
            psi_cav(nplace + 1) = SQRT(kappa_a) * SQRT(1.0 * (na + 1)) * &
                                & psi(nplace_pm1 + 1)
            psi_cav(nplace + 2) = SQRT(kappa_a) * SQRT(1.0 * (na + 1)) * &
                                & psi(nplace_pm1 + 2)
            psi_cav(nplace + 3) = SQRT(kappa_a) * SQRT(1.0 * (na + 1)) * &
                                & psi(nplace_pm1 + 3)
          END DO
        END DO
        psi_clone = psi_cav

        ! Normalise the state
        trace = 0
        DO l=1,d_Hilb
          trace = trace + (ABS(psi_clone(l))) ** 2
        END DO
        tracesq = SQRT(trace)
        psi = psi_clone / tracesq

      ELSE IF (rand > P(1) .AND. rand <= P(2)) THEN
        !############################################!
        !     A reflection jump from cavity a        !
        !############################################!
        count_r_a = count_r_a + 1

        ! Calculate the state for a cavity jump from a
        psi_cav = 0.0
        DO na=0,N-1
          DO nb=0,N
            ! photon number index
            nplace = (3 * (N + 1) * na) + (3 * nb)
            ! photon number index for (na + 1)
            nplace_pm1 = (3 * (N + 1) * (na + 1)) + (3 * nb)
            psi_cav(nplace + 1) = SQRT(kappa_a) * SQRT(1.0 * (na + 1)) * &
                                & psi(nplace_pm1 + 1)
            psi_cav(nplace + 2) = SQRT(kappa_a) * SQRT(1.0 * (na + 1)) * &
                                & psi(nplace_pm1 + 2)
            psi_cav(nplace + 3) = SQRT(kappa_a) * SQRT(1.0 * (na + 1)) * &
                                & psi(nplace_pm1 + 3)
          END DO
        END DO

        ! Calculate the state after an atom decay
        psi_clone = 0.0
        DO na=0,N
          DO nb=0,N
            ! photon number index
            nplace = (3 * (N + 1) * na) + (3 * nb)
            ! |e> -> |g>
            psi_clone(nplace + 1) = SQRT(0.5 * gamma) * psi(nplace + 2)
            ! |f> -> \xi |e>
            psi_clone(nplace + 2) = SQRT(0.5 * gamma) * xi * psi(nplace + 3)
          END DO
        END DO

        ! Add this state with a cavity jump state for a superposition
        psi_clone = psi_clone + psi_cav

        ! Normalise the state
        trace = 0
        DO l=1,d_Hilb
          trace = trace + (ABS(psi_clone(l))) ** 2
        END DO
        tracesq = SQRT(trace)
        psi = psi_clone / tracesq

      ELSE IF (rand > P(2) .AND. rand <= P(3)) THEN
        !############################################!
        !     A cavity emission from cavity b        !
        !############################################!
        count_t_b = count_t_b + 1

        ! Calculate the state for a cavity decay from b
        psi_cav = 0.0
        DO na=0,N
          DO nb=0,N-1
            ! photon number index
            nplace = (3 * (N + 1) * na) + (3 * nb)
            ! photon number index for (nb + 1)
            nplace_pm1 = (3 * (N + 1) * na) + (3 * (nb + 1))
            psi_cav(nplace + 1) = SQRT(kappa_b) * SQRT(1.0 * (nb + 1)) * &
                                & psi(nplace_pm1 + 1)
            psi_cav(nplace + 2) = SQRT(kappa_b) * SQRT(1.0 * (nb + 1)) * &
                                & psi(nplace_pm1 + 2)
            psi_cav(nplace + 3) = SQRT(kappa_b) * SQRT(1.0 * (nb + 1)) * &
                                & psi(nplace_pm1 + 3)
          END DO
        END DO
        psi_clone = psi_cav

        ! Normalise the state
        trace = 0
        DO l=1,d_Hilb
          trace = trace + (ABS(psi_clone(l))) ** 2
        END DO
        tracesq = SQRT(trace)
        psi = psi_clone / tracesq
      ELSE IF (rand > P(3) .AND. rand <= P(4)) THEN
        !############################################!
        !     A reflection jump from cavity b        !
        !############################################!
        count_r_b = count_r_b + 1

        ! Calculate the state for a cavity jump from a
        psi_cav = 0.0
        DO na=0,N
          DO nb=0,N-1
            ! photon number index
            nplace = (3 * (N + 1) * na) + (3 * nb)
            ! photon number index for (nb + 1)
            nplace_pm1 = (3 * (N + 1) * na) + (3 * (nb + 1))
            psi_cav(nplace + 1) = SQRT(kappa_b) * SQRT(1.0 * (nb + 1)) * &
                                & psi(nplace_pm1 + 1)
            psi_cav(nplace + 2) = SQRT(kappa_b) * SQRT(1.0 * (nb + 1)) * &
                                & psi(nplace_pm1 + 2)
            psi_cav(nplace + 3) = SQRT(kappa_b) * SQRT(1.0 * (nb + 1)) * &
                                & psi(nplace_pm1 + 3)
          END DO
        END DO

        ! Calculate the state after an atom decay
        psi_clone = 0.0
        DO na=0,N
          DO nb=0,N
            ! photon number index
            nplace = (3 * (N + 1) * na) + (3 * nb)
            ! |e> -> |g>
            psi_clone(nplace + 1) = SQRT(0.5 * gamma) * psi(nplace + 2)
            ! |f> -> \xi |e>
            psi_clone(nplace + 2) = SQRT(0.5 * gamma) * xi * psi(nplace + 3)
          END DO
        END DO

        ! Add this state with a cavity jump state for a superposition
        psi_clone = psi_clone + psi_cav

        ! Normalise the state
        trace = 0
        DO l=1,d_Hilb
          trace = trace + (ABS(psi_clone(l))) ** 2
        END DO
        tracesq = SQRT(trace)
        psi = psi_clone / tracesq

      ! Close IF statement for deciding on which jump
      END IF
    ! Close IF statement for Evolve/Jump decision
    END IF

    ! Print completion percentage of program
    IF (MOD(1.0 * k, 1E3) == 0.0) THEN
      PRINT*, 100.0 * (1.0 * k) / (1.0 * steps), "% complete"
    END IF

  ! Close time loop
  END DO
!##############################################################################!
!                     Section E: Write results to file                         !
!##############################################################################!
  ! Average steady state mean photon number by number of time steps
  mean_photon_a = mean_photon_a / steps
  mean_photon_b = mean_photon_b / steps
  ! Average corr_array by number of samples taken
  correlation = correlation / (1.0 * number_of_recordings)
  cross_correlation = cross_correlation / (1.0 * number_of_recordings)
  ! Normalise by steady state mean photon number in cavity
  correlation = correlation / mean_photon_a
  cross_correlation = cross_correlation / mean_photon_b

  ! Open file for paramters, time and mean photon number to be written to
  OPEN(UNIT=1, file=filename_time, STATUS='replace', ACTION='write')
  ! Write the values of the parameters to the first 9 lines. The time values
  ! will start on the 10th line
  WRITE(1,*) omega
  WRITE(1,*) delta
  WRITE(1,*) xi
  WRITE(1,*) alpha
  WRITE(1,*) D_a
  WRITE(1,*) kappa_a
  WRITE(1,*) D_b
  WRITE(1,*) kappa_b
  WRITE(1,*) ' '

  ! Open files for state probabilities to be written to
  OPEN(UNIT=2, file=filename_correlation, STATUS='replace', ACTION='write')

  DO k=0,tau_max
    WRITE(1,*) k * dt
    WRITE(2,*) correlation(k), cross_correlation(k)
  END DO

  ! Print number of jumps that occured

  PRINT*, 100.0 * (count_t_a + count_r_a) / steps, "% of steps were jumps"

  PRINT*, count_t_a, "transmission jumps from a"
  PRINT*, count_r_a, "reflection jumps from a"

  PRINT*, count_t_b, "transmission jumps from b"
  PRINT*, count_r_b, "reflection jumps from b"

END PROGRAM Two_Filter_Cross_Correlation

! Subroutine for creating a random seed based on the computer clock.
SUBROUTINE init_random_seed()
  IMPLICIT NONE
  INTEGER :: l, n, clock
  INTEGER, DIMENSION(:), ALLOCATABLE :: seed
  CALL RANDOM_SEED(size = n)
  ALLOCATE(seed(n))
  CALL SYSTEM_CLOCK(COUNT=clock)
  seed = clock + 37 * (/ (l - 1, l = 1, n) /)
  CALL RANDOM_SEED(PUT = seed)
  DEALLOCATE(seed)
END SUBROUTINE init_random_seed
