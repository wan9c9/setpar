       include "tpar_mod.f90" 
	program main_program
        use tpar_mod
        implicit none
        
	call test_asycov_tpar
	!call test_asycov
        !call simulation_threshold
        !call simulation
        !call rd_tpar
        end program main_program

