!==============================================================================
! MOLECULAR DYNAMICS SIMULATION PROGRAM
! This program simulates the dynamics of a system of atoms using the Verlet 
! algorithm and Lennard-Jones potential for the ACT subject of the TCCM
! 
! This program has been done by: 
! - Francisco Javier Patiño López
! - Alejandro Barrios
! - Maeva Rochefort
!
! Key features:
! - Reads atomic coordinates from an XYZ file
! - Computes interatomic distances
! - Calculates Lennard-Jones potential energy
! - Implements velocity Verlet algorithm for time integration
! - Outputs trajectory and energies to a file
!==============================================================================

program mdhw
  implicit none
  ! Variable declarations for the program 
  integer :: Natoms          ! total number of atoms in the system
  integer :: i_stat          ! status for memory allocation
  integer :: io_stat         ! status for input/output operations
  integer :: i, j, k         ! loop variables
  ! Atomic properties
  character(len=2), allocatable :: sat(:)
  double precision, allocatable :: coord(:,:)
  double precision, allocatable :: mass(:)
  double precision, allocatable :: distance(:,:)
  double precision, allocatable :: velocity(:,:)
  double precision, allocatable :: acceleration(:,:)
  ! Other variables
  character(len=256) :: input_file
  double precision :: sigma, eps, lj_potential, timestep, kin_energy, tot_energy
  ! Lennard Jones parameters are initialized
  eps = 0.0661
  sigma = 0.3345
  
  ! Asks the user which is the name of the input file
  write(*,*) "Write the name of your input file"
  read(*,*) input_file
  
  ! Read the number of atoms and are printed
  Natoms = read_Natoms(input_file)
  print *, "The number of atoms is equal to:", Natoms

  ! ----------------------------------------------------------
  !                 ALLOCATION PART
  !  Manage all the allocations of the arrays used during the
  !  simulation
  ! ----------------------------------------------------------

  allocate(coord(Natoms, 3), stat=i_stat)
  if (i_stat /= 0) then
        print *, "Memory allocation failed on the coordinates!"
        stop
  end if

  allocate(mass(Natoms), stat=i_stat)
  if(i_stat/=0) then
          print *, "Memory allocation failed on the mass!"
          stop
  end if

  allocate(sat(Natoms), stat=i_stat)  
  if(i_stat/=0) then
          print *, "Memory allocation failed on the atomic symbol!"
          stop
  end if
  
  allocate(distance(Natoms, Natoms), stat=i_stat)
  if (i_stat/=0) then 
          print *, "Memory allocation failed on the distance!"
          stop 
  end if 

  allocate(velocity(Natoms,3), stat=i_stat)
  if (i_stat/=0) then 
          print*, "Memory allocation failed on the velocity!"
          stop 
  end if 

  allocate(acceleration(Natoms, 3), stat=i_stat) 
  if (i_stat/=0) then 
          print *, "Memory allocation failed on the acceleration!"
          stop 
  end if 
  ! --------------------------------------------------------------
  !            CALLING ALL THE SUBROUTINES
  ! --------------------------------------------------------------

  ! read the structure of the molecule
  call read_molecule(input_file, Natoms, coord, mass, sat)
  do i=1, Natoms
      print *, sat(i), coord(i,:)
  end do
  
  ! Compute initial distances and accelerations
  call compute_distances(Natoms, coord, distance)
  call compute_acc(Natoms, coord, mass, distance, acceleration)

  ! Run the MD simulation with the verlet algorithm
  call verlet_algorithm(coord, velocity, timestep, acceleration)
  
  ! --------------------------------------------------------------
  !            DEALLOCATION OF THE ARRAYS
  ! --------------------------------------------------------------
  deallocate(mass)
  deallocate(coord)
  deallocate(sat)
  deallocate(distance)
  deallocate(velocity)
  deallocate(acceleration)

contains
  ! --------------------------------------------------------------
  !            DESCRIBING THE ATOMS
  ! --------------------------------------------------------------
  

  ! --------------------------------------------------------------
  !       THE FUNCTION BELOW READS THE NUMBER OF ATOMS:
  ! Reads the number of the atoms from the first line of the 
  ! .xyz file
  ! --------------------------------------------------------------

  integer function read_Natoms(input_file) result(Natoms)
    implicit none
    integer :: io_stat
    character(*), intent(in) :: input_file

    open(10, file=input_file, status='old', action='read', iostat=io_stat)
    if (io_stat /= 0) then
      print *, "Error opening file:", input_file
      stop
    else
      print *, "File opened successfully: ", input_file
    end if

    read(10, *, iostat=io_stat) Natoms
    if (io_stat /= 0) then
        print *, "Error on reading the number of atoms in the file", input_file
        print *, "Please verify that the format of the file is the correct one"
    else
        print *, "The format of the file is correct"
    end if

    close(10)
  end function read_Natoms
  ! --------------------------------------------------------------
  !       THE FUNCTION BELOW OPENS AND READ THE .XYZ FILE:
  ! Read atomic coordinates and masses from the .xyz file and 
  ! assign the atomic symbol to the mass
  ! --------------------------------------------------------------

  subroutine read_molecule(input_file, Natoms, coord, mass, sat)
          implicit none
          character(len=*), intent(in) :: input_file
          character(len=2), intent(out) :: sat(Natoms)
          integer :: io_stat, i
          integer, intent(in) :: Natoms
          double precision, intent(out) :: coord(Natoms, 3)
          double precision, intent(out) :: mass(Natoms)

          open(10, file=input_file, status='old', action='read', iostat=io_stat)
          if (io_stat /= 0) then
              print *, "Error on opening the input"
              stop
          end if

          read(10,*)
          if (io_stat /= 0) then
              print *, "Error on skipping the first line of the input"
              stop
          end if 
          
          do i=1,Natoms
               read(10,*) coord(i,:), mass(i)
               if (io_stat /= 0) then
                   print *, "Error on reading the atom", i, "in the file:", input_file
                   stop
               end if
               if (abs(mass(i) - 39.948d0) < 1.0d-6) then 
                       sat(i) = 'Ar'
               else 
                       print *, "Element not recognized for atom", i, "with mass", mass(i)
               end if              
          end do
          close(10)
  end subroutine read_molecule

  ! --------------------------------------------------------------
  !  THE SUBROUTINE BELOW COMPUTE THE DISTANCES BETWEEN THE ATOMS:
  ! Calculates the distances between the atoms and stores the 
  ! results on a distance matrix
  ! --------------------------------------------------------------

  subroutine compute_distances(Natoms, coord, distance)
          implicit none 
          integer :: io_stat, xdist, i, j
          integer, intent(in) :: Natoms
          double precision, intent(in) :: coord(Natoms, 3)
          double precision, intent(out) :: distance(Natoms, Natoms)

          do i = 1, Natoms
                do j = 1, Natoms
                        if (i == j) then 
                                distance(i,j) = 0
                        else 
                                distance(i,j) = sqrt((coord(i, 1) - coord(j, 1))**2 + &
                                    (coord(i, 2) - coord(j, 2))**2 + &
                                    (coord(i, 3) - coord(j, 3))**2) 
                        end if 
                 end do 
          end do 
  end subroutine compute_distances 

  ! --------------------------------------------------------------
  !    THE SUBROUTINE BELOW COMPUTES THE LENNARD-JONES POTENTIAL
  ! Computes the Lennard-Jones potential energy: 
  ! V(r) = 4*eps*((sigma/r)^12 - (sigma/r)^6) 
  ! --------------------------------------------------------------

  double precision function V(eps, sigma, Natoms, distance)
          implicit none
          double precision, intent(in) :: eps
          double precision, intent(in) :: sigma
          integer, intent(in) :: Natoms
          double precision, intent(inout) :: distance(Natoms, Natoms)
          integer :: i, j
          
          V = 0.0
          do i = 1, Natoms
                do j = i + 1, Natoms
                        if (i==j) then
                                CYCLE
                        else
                                V = V + 4 * eps * ((sigma / distance(i,j))**12 - (sigma / distance(i,j))**6) 
                        end if
                end do
          end do 
  end function V

  ! --------------------------------------------------------------
  !    THE SUBROUTINE BELOW COMPUTES THE KINETIC ENERGY:
  ! Computes the kinetic energy: 
  ! T = Σ(1/2 * m * v^2)
  ! --------------------------------------------------------------

  double precision function T(Natoms, velocity, mass)
         implicit none
         integer, intent(in) :: Natoms
         double precision, intent(out) :: velocity(Natoms, 3)
         double precision, intent(in) :: mass(Natoms)
         integer :: i
         
         T = 0.0  ! Initialize T
         do i = 1, Natoms
                T = T + 0.5*mass(i)*(velocity(i,1)**2 + velocity(i,2)**2 + velocity(i,3)**2) 
         end do       
         
  end function T

  ! --------------------------------------------------------------
  !    THE SUBROUTINE BELOW COMPUTES THE TOTAL ENERGY:
  ! Computes the total energy adding the kinetic and the potential: 
  ! E = T + V
  ! --------------------------------------------------------------

  double precision function E(T, V)
         implicit none
         double precision, intent(in) :: T
         double precision, intent(in) :: V
         E = T + V
  end function E

  ! --------------------------------------------------------------
  !    THE SUBROUTINE BELOW COMPUTES THE ACCELERATION:
  !  Computes the acceleration from forces derived from LJ potential: 
  ! F = -∇V, a = F/m
  ! --------------------------------------------------------------

  subroutine compute_acc(Natoms, coord, mass, distance, acceleration) 
         implicit none
         integer, intent(in) :: Natoms
         double precision, intent(in) :: coord(Natoms, 3)
         double precision, intent(in) :: mass(Natoms)
         double precision, intent(in) :: distance(Natoms, Natoms)
         double precision, intent(out) :: acceleration(Natoms, 3) 
         integer :: i, j, k
         double precision :: eps, sigma
         
         eps = 0.0661
         sigma = 0.3345


         do i=1, Natoms
                 acceleration(i, 1) = 0.0
                 acceleration(i, 2) = 0.0
                 acceleration(i, 3) = 0.0
                do j=1, Natoms
                        do k = 1, 3
                        if (i==j) then 
                              CYCLE  
                        else
                               acceleration(i, k) = acceleration(i,k) - 24 * eps / mass(i) * &
                     ((sigma / distance(i, j))**6 - 2 * (sigma / distance(i, j))**12) * &
                     (coord(i, k) - coord(j, k)) / distance(i, j)
                        end if
                        end do 
                end do 
          end do 
  end subroutine compute_acc

  ! --------------------------------------------------------------
  !  IN THE SUBROUTINE BELOW IS IMPLEMENTED THE VERLET ALGORITHM
  !  AND THE OUTPUT FILE IS CREATED, FORMATTED AND WRITTEN:
  ! Computes the Verlet Algorithm that consists on the following
  ! steps: 
  ! 1. Update positions: r(t+dt) = r(t) + v(t)*dt + 0.5*a(t)*dt^2
  ! 2. Update velocities (first part): v(t+dt/2) = v(t) + 0.5*a(t)*dt
  ! 3. Compute new accelerations from updated positions
  ! 4. Update velocities (second part): v(t+dt) = v(t+dt/2) + 0.5*a(t+dt)*dt
  !
  ! Also manage the trajectory output.
  ! --------------------------------------------------------------

  subroutine verlet_algorithm(coord, velocity, timestep, acceleration)
          implicit none
          double precision, intent(out) :: coord(Natoms, 3)
          double precision, intent(out) :: acceleration(Natoms, 3)
          double precision, intent(out) :: velocity(Natoms, 3)
          double precision :: timestep
          integer :: steps, iteration, X, i, j
          logical :: file_exists
          double precision :: lj_potential, kin_energy, tot_energy
          
          timestep = 0.2 
          steps = 1000
          velocity(:,:) = 0.0

          ! Check if trajectory file exists and delete it at the start
          inquire(file='trajectory.dat', exist=file_exists)
          if (file_exists) then
              open(10, file='trajectory.dat', status='old')
              close(10, status='delete')
          end if
          
          write(*,*) "Each how many iterations you want to write info:"
          read(*,*) X
          
          do iteration = 1, steps 
                
                ! Update positions
                do i=1, Natoms
                      do j=1, 3
                              coord(i,j) = coord(i,j) + velocity(i,j) * timestep + &
                                         0.5 * acceleration(i,j) * (timestep**2)
                      end do 
                end do
                
                ! First part of velocity verlet
                do i = 1, Natoms
                      do j = 1, 3
                              velocity(i, j) = velocity(i, j) + 0.5 * acceleration(i, j) * timestep
                      end do
                end do

                ! Update distances and accelerations
                call compute_distances(Natoms, coord, distance)
                call compute_acc(Natoms, coord, mass, distance, acceleration)

                ! Second part of velocity verlet
                do i = 1, Natoms
                      do j = 1, 3
                              velocity(i, j) = velocity(i, j) + 0.5 * acceleration(i, j) * timestep
                      end do
                end do
           
                lj_potential = V(eps, sigma, Natoms, distance)
                kin_energy = T(Natoms, velocity, mass)
                tot_energy = E(kin_energy, lj_potential)

                if(mod(iteration, X)==0) then
                        ! Open file in append mode, creating it if it doesn't exist
                        open(10, file='trajectory.xyz', status='unknown', position='append', action='write')
                        
                        write(10,*) Natoms
                        write(10,*) "STEP:", iteration, "KIN:", kin_energy, &
                                   "POT:", lj_potential, "TOT:", tot_energy
                        do i = 1, Natoms
                              write(10,'(A2,1X,3F12.6)') sat(i), coord(i,1), coord(i,2), coord(i,3)
                        end do 
                        
                        close(10)
                end if
          end do
          print *, "PROGRAM HAS FINISHED SUCCESFULLY THE RESULT ARE STORAGED ON TRAJECTORY.XYZ FILE"
   end subroutine verlet_algorithm
          
end program mdhw
