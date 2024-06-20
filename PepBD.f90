Program PepBD

use datatypes
use sys_vars
use utilities
use input
use pdbfile
use database
use energy_calculation
use advanced_function
use optimization_techniques
use omp_lib

implicit none
integer                                         :: step, i, sub_circle
real                                            :: ran2, rmsd=0.0
double precision                                :: cur_energies(6)
integer                                         :: SA_cycle = 1, net_step = 1
character*100                                   :: pdb_name
integer, allocatable                            :: cg_flags(:)
type(databackup), allocatable                   :: backbone_backup(:) 
type(groupdetails), allocatable                 :: group_pep(:), group_rec(:), original_pep(:)
!type(groupdetails), allocatable                 :: SA_best_pep(:) !stores best peptide found during a simulated annealing interval, which is passed along to next round
type(energyparameters), allocatable             :: pep_para(:), rec_para(:) !, SA_best_para(:)
     
!!!! print cool ASCII art (doesn't print correctly at the moment, so it is commented out) !!!!
!print ("(A)"),"   &
!&.______    _______ .______   .______    _______&
!&|   _  \  |   ____||   _  \  |   _  \  |       \&
!&|  |_)  | |  |__   |  |_)  | |  |_)  | |  .--.  |&
!&|   ___/  |   __|  |   ___/  |   _  <  |  |  |  |&
!&|  |      |  |____ |  |      |  |_)  | |  '--'  |&
!&| _|      |_______|| _|      |______/  |_______/"
!
!print ("(A)"), "  &
!&  ______     ___      .______        ______    __          __  ___     __    __       ___       __       __          __          ___      .______&
!& /      |   /   \     |   _  \      /  __  \  |  |        |  |/  /    |  |  |  |     /   \     |  |     |  |        |  |        /   \     |   _  \&
!&|  ,----'  /  ^  \    |  |_)  |    |  |  |  | |  |        |  '  /     |  |__|  |    /  ^  \    |  |     |  |        |  |       /  ^  \    |  |_)  |&
!&|  |      /  /_\  \   |      /     |  |  |  | |  |        |    <      |   __   |   /  /_\  \   |  |     |  |        |  |      /  /_\  \   |   _  <& 
!&|  `----./  _____  \  |  |\  \----.|  `--'  | |  `----.   |  .  \     |  |  |  |  /  _____  \  |  `----.|  `----.   |  `----./  _____  \  |  |_)  |&
!& \______/__/     \__\ | _| `._____| \______/  |_______|   |__|\__\    |__|  |__| /__/     \__\ |_______||_______|   |_______/__/     \__\ |______/" 
!
!print ("(A)"), "  &
!&.__   __.   ______         _______.___________.    ___   .___________. _______&
!&|  \ |  |  /      |       /       |           |   /   \  |           ||   ____|&
!&|   \|  | |  ,----'      |   (----`---|  |----`  /  ^  \ `---|  |----`|  |__&
!&|  . `  | |  |            \   \       |  |      /  /_\  \    |  |     |   __|&
!&|  |\   | |  `----.   .----)   |      |  |     /  _____  \   |  |     |  |____ &
!&|__| \__|  \______|   |_______/       |__|    /__/     \__\  |__|     |_______|"              
                                                                                                                                        
call inputfile !read in the input file
call calc_system_size !determine size of system

!seed random number generator
if (ran_seed(1).eq.-1) then
    call RANDOM_SEED()
else
    call RANDOM_SEED(put=[ran_seed,ran_seed,ran_seed,ran_seed,ran_seed,ran_seed,ran_seed,ran_seed])
endif

!#### BIG BLOCK OF MEMORY ALLOCATIONS #####

!peptide parameters
allocate(backbone_backup(pep_res))
allocate(group_pep(pep_res))
allocate(Tgroup_pep(pep_res))
allocate(pep_para(pep_res))
!allocate(SA_best_pep(pep_res))
!allocate(SA_best_para(pep_res))
allocate(residuesite(pep_res))
allocate(res_starts(pep_res))
allocate(res_ends(pep_res))
allocate(sc_starts(pep_res))
allocate(sc_ends(pep_res))
allocate(original_pep(pep_res))
allocate(atom(atom_per_pep_res, pep_res))
allocate(numex_pep(atom_per_pep_res, pep_res))
allocate(numex4_pep(atom_per_pep_res, pep_res))
allocate(inb_pep(2,20,atom_per_pep_res, pep_res))
allocate(inb4_pep(2,60,atom_per_pep_res, pep_res))
allocate(atomid_pep(atom_per_pep_res, pep_res))
allocate(mdcrd_pep(3,atom_per_pep_res, pep_res))
allocate(charge_pep(atom_per_pep_res, pep_res)) 
allocate(epsion_pep(atom_per_pep_res, pep_res))
allocate(r_pep(atom_per_pep_res, pep_res))
allocate(rborn_pep(atom_per_pep_res, pep_res))
allocate(fs_pep(atom_per_pep_res, pep_res))
allocate(dielecons_pep(atom_per_pep_res, pep_res))
allocate(phis(pep_res))
allocate(psis(pep_res))
allocate(energy_per_res(pep_res, 4))
allocate(res_coords_flags(pep_res))
allocate(res_params_flags(pep_res))
allocate(aa_groups(40,pep_res))
allocate(cg_flags(pep_res))
allocate(atoms_per_res(pep_res))

!Receptor parameters
allocate(group_rec(rec_res))
allocate(rec_para(rec_res))
allocate(mdcrd_rec(3, atom_num_rec))
allocate(charge_rec(atom_num_rec)) 
allocate(epsion_rec(atom_num_rec))
allocate(r_rec(atom_num_rec))
allocate(rborn_rec(atom_num_rec))
allocate(fs_rec(atom_num_rec))
allocate(dielecons_rec(atom_num_rec))

!GB Solvation parameters
if (solvation_flag.eq.1) then
    allocate(alpha_pep_sys(atom_per_pep_res, pep_res))
    allocate(alpha_pep_iso(atom_per_pep_res, pep_res))
    allocate(integrals_rec_iso(atom_num_rec))
    allocate(alpha_rec_sys(atom_num_rec))
    allocate(alpha_rec_iso(atom_num_rec))
    allocate(integrals_comp_rec_from_pep(atom_num_rec, pep_res)) 
    allocate(integrals_comp_pep_from_rec(atom_per_pep_res, pep_res))
    allocate(GB_area_calc_flags(pep_res)) 
endif

!pairwise distances
allocate(rij_matrix_pep(atom_per_pep_res, atom_per_pep_res, pep_res, pep_res))
allocate(rij_matrix_sys(atom_num_rec, atom_per_pep_res, pep_res))
allocate(rij_matrix_rec(atom_num_rec, atom_num_rec))

!General simulation parameters
allocate(SEQ_PerRes_attempts(pep_res))
allocate(SEQ_PerRes_success(pep_res))
allocate(lbres(atom_num))
if (NP_flag .eq. 1) allocate(LCPO_indices(atom_num))
SEQ_PerRes_attempts = 0; SEQ_PerRes_success = 0

!############################################ 

call prep_output !prepare output files
call load_system(group_pep, group_rec, original_pep, backbone_backup) !read pdb file for system
call PH_checking(group_pep, group_rec) !check amino acid names in pdb file match with pH of design
call energy_parameter(group_pep, pep_para, group_rec, rec_para) !load the force field parameters for the system
!if (NP_flag .eq. 1) call load_LCPO_params() !load parameters for calculating solvent accessible surface area (SASA)
call rotamerlib !load rotamers for all amino acids
call ramachandranmap !load ramachandran map

!store receptor coordinates and parameters in array format
call load_coords(group_rec, 2)
call load_params(group_rec, rec_para, 2)
atom_links_flag = 1

!calculate GB and NP energy of receptor (note that NP energy is currently not calculated)
if (solvation_flag.eq.1) then
    call calc_rec_iso_GB_energy
    GB_area_calc_flags = 1
endif

cg_flags = 0 !turn off coarse-graining of sidechains for all design (not implemented into PepBD yet) 
best_energies = 10000

!if not restarting design, then modify peptide sequence so it follows the criteria specified by the user
if(restart_flag.eq.0) then
    if (check_init_structure.eq.1) then !verify there are no steric overlaps in the starting structure
        call check_backbone_overlap(group_pep, pep_para, i)
        if (i.eq.0) then
            write(5,"(A)") "Starting structure has steric overlap with receptor!"
            write(5,"(A)") "Please fix the input structure and restart design."
            write(5,"(A)") "Terminating program."
            stop
        endif
    endif
    call fit_to_hydration_range(group_pep, pep_para, backbone_backup, sub_circle) !adjust hydration properties of peptide to match input file specifications
    call replace_extra_tryptophans(group_pep, pep_para, backbone_backup) !replace tryptophans with valine if above the tryptophan limit
    call replace_excluded_AAs(group_pep, pep_para, backbone_backup) !replace all amino acids excluded in design
    write(pdb_name, "(I0,A,I0,A,I0)") 0, "_", 0, "_", 0
    call writepdb_sys(group_pep, group_rec, pdb_name)
endif

!if doing simulated annealing, then initialize best group and parameters with the starting peptide 
!if (SA_flag.eq.1) then
!    SA_best_pep = group_pep
!    SA_best_para = pep_para
!endif

!*********************************
! THE INITIAL RANDOMIZATION STAGE
!*********************************

if (initial_randomize.eq.1.and.restart_flag.eq.0) then
    !do a short design before randomizing peptide sequence to get baseline for a "good" peptide binder

    call bindingenergy_full(group_pep, pep_para, cur_energies)
    best_energies=cur_energies
    open(5, file="energyprofile.txt", access="append") !creates output file for energy calculations
        write(5,6) 0, 0, cur_energies(1), cur_energies(2), cur_energies(3), cur_energies(4), cur_energies(5), cur_energies(6)
    close(5)

    !If we are using simulated annealing for design, then use the final temperature for this design stage
    if (SA_flag.eq.1) then
        ekt_seq = SA_T_end
        ekt_backbone = SA_T_end
    endif

    step=1
    do while (step.le.steps_before_randomize)
        call sequence_optimization(group_pep,               group_rec, pep_para, backbone_backup, step, 0, cur_energies)

        !write results to output files
        open(5, file="energyprofile.txt", access="append")
        write(5,6) 0, step, cur_energies(1), cur_energies(2), cur_energies(3), cur_energies(4), cur_energies(5), cur_energies(6)
        close(5)
        step = step+1
    enddo

    write(pdb_name, "(I0,A,I0,A,I0)") 0, "_", step+1, "_", 0
    call writepdb_sys(group_pep, group_rec, pdb_name)

    !store coordinates of backbone 
    open(5, file="backup4backbone.txt", status="replace")
        do i=1, pep_res
        write(5,704) backbone_backup(i)%coo(1), backbone_backup(i)%coo(2), backbone_backup(i)%coo(3)
        enddo
        write(5,"(i1)") flag4conformer
        write(5,"(f6.3)") best_energies(1)
        write(5,"(a6,i4)") "cycle=", SA_cycle
        write(5,"(a12, i4)") "final step=", step
    close(5)
endif

!*********************************
!    THE MAIN DESIGN PROCESS
!*********************************

!if we doing simulated annealing, then specify starting annealing temperature
if (SA_flag.eq.1) then
    ekt_seq = SA_T_start
    ekt_backbone = SA_T_start
endif

!if peptide designe being started from intermediate point, then reduce the temperature based on how far annealing progressed in previous design
if (restart_flag.eq.1) then
    ekt_seq = ekt_seq * (SA_fraction)**(SA_start_cycle-1)
    ekt_backbone = ekt_seq
endif

!The main design loop. Very pretty.
do SA_cycle=SA_start_cycle, SA_num_cycles
    !first randomize the sequence to exit out of local minimum
    if (SA_flag.eq.0.or.SA_cycle.eq.1) then !if doing simulated annealing, only want to randomize at the very beginning of design
        
        call sequence_randomization(group_pep, pep_para, backbone_backup)   

        !initialize SA storage for the current peptide
        !SA_best_pep=group_pep
        !SA_best_para=pep_para
        !SA_best_hyd = cur_hyd
        !SA_above_hyd = above_hyd
        !SA_below_hyd = below_hyd
        !SA_trp_count = trp_count
    endif

    atom_links_flag = 1
    call bindingenergy_full(group_pep, pep_para, cur_energies) !get binding energy of initial complex

    !write data of peptide at start of this cycle to output file
    open(5, file="energyprofile.txt", access="append") 
        write(5,6) SA_cycle, 0, cur_energies(1), cur_energies(2), cur_energies(3), cur_energies(4), cur_energies(5), cur_energies(6)
    close(5)

    !write pdbfile for peptide before starting design at current system temperature
    write(pdb_name, "(I0,A,I0,A,I0)") SA_cycle, "_", 0, "_", 0
    call writepdb_sys(group_pep, group_rec, pdb_name)
    
    !determine number of steps, depending on whether using simulated annealing or constant temperature design
    if (SA_flag.eq.0) then
        step = nstep_start
    else
        step = 1
    endif

    !design peptide for N steps at current system temperature!
    do while (step.le.steps_per_cycle)
        if(step.le.steps_before_conf_change) then ! for first few steps, don't do conformation changes since bound state is not stable (conformation changes not helpful)
            call sequence_optimization(group_pep, group_rec, pep_para, backbone_backup, step, SA_cycle, cur_energies)
        else
            !choose one of the move types to make randomly 
            call RANDOM_NUMBER(ran2)
            if(ran2.lt.seq_switch) then !change peptide sequence
                call sequence_optimization(group_pep,               group_rec, pep_para, backbone_backup, step, SA_cycle, cur_energies)
            !elseif (ran2.le.backbone_switch) then !change peptide backbone
            else
                call backbone_optimization(group_pep, original_pep, group_rec, pep_para, backbone_backup, step, SA_cycle, cur_energies)
            !else !rotate or translate peptide
            !    call rigid_body_optimization(group_pep, group_rec, pep_para, SA_best_pep,             backbone_backup, step, 0, binding_energy_old, vdw_old, ele_old, sgb_old, SASA_old, pep_stability_old)
            endif
        endif

        !write energy of peptide to output file
        open(5, file="energyprofile.txt", access="append")
            write(5,6) SA_cycle, step, cur_energies(1), cur_energies(2), cur_energies(3), cur_energies(4), cur_energies(5), cur_energies(6)
        close(5)

        !write RMSD of peptide to output file
        open(5, file="rmsd.txt", access="append")
            call calc_rmsd(group_pep, original_pep, rmsd)
            rmsd_flag = 0
            write(5,703) SA_cycle, step, rmsd
        close(5)

        !periodically write the system into a pdb file
        if(mod(step,PDB_gen_freq).eq.0) then
            write(pdb_name, "(I0,A,I0,A,I0)") SA_cycle, "_", step+1, "_", 0
            call writepdb_sys(group_pep, group_rec, pdb_name)
            open(5, file="backup4backbone.txt", status="replace")
                do i=1, pep_res
                    write(5,704) backbone_backup(i)%coo(1), backbone_backup(i)%coo(2), backbone_backup(i)%coo(3)
                enddo
                write(5,"(i1)") flag4conformer
                write(5,"(f6.3)") best_energies(1)
                write(5,"(a6,i4)") "cycle=", SA_cycle
                write(5,"(a12, i4)") "final step=", step
            close(5)
        endif
        net_step = net_step + 1
        step = step + 1
    enddo

    !do simulated annealing bookkeeping
    if (SA_flag.eq.1) then
        ekt_seq = ekt_seq * SA_fraction
        ekt_backbone = ekt_backbone * SA_fraction
        ekt_rigid = ekt_rigid * SA_fraction
        open(5, file="output.txt", access="append")
            write(5,"(A,F5.3)") "Starting next cycle, reference kBT value is now ", ekt_seq
        close(5)

        !reload the best peptide found in previous cycle
        !group_pep = SA_best_pep
        !pep_para = SA_best_para 
        !cur_hyd = SA_best_hyd 
        !above_hyd = SA_above_hyd
        !below_hyd = SA_below_hyd
        !trp_count = SA_trp_count
        !res_coords_flags = 1
        !res_params_flags = 1
        !atom_links_flag = 1
        !GB_area_calc_flags = 1
    endif
enddo

!write out probability data
open(102, file='output.txt', access='append')
write(102,*); write(102,"(A)") "********************"
write(102, "(A)") "Probalities of Accepting Peptide Changes"
if (SEQ_SWAP_attempts.eq.0) SEQ_SWAP_attempts = 1
write(102,705) "Swap Mutations", real(SEQ_SWAP_success)*1./real(SEQ_Swap_attempts)
if (SEQ_POINT_attempts.eq.0) SEQ_POINT_attempts = 1
write(102,705) "Point Mutations", real(SEQ_POINT_success)*1./real(SEQ_POINT_attempts)
if (CONF_N2C_attempts.eq.0) CONF_N2C_attempts = 1
write(102,705) "N and C termini CONROT moves", real(CONF_N2C_success)*1./real(CONF_N2C_attempts)
if (CONF_N_attempts.eq.0) CONF_N_attempts = 1
write(102,705) "N terminus CONROT moves", real(CONF_N_success)*1./real(CONF_N_attempts)
if (CONF_C_attempts.eq.0) CONF_C_attempts = 1
write(102,705) "C terminus CONROT moves", real(CONF_C_success)*1./real(CONF_C_attempts)
if (CONF_Internal_attempts.eq.0) CONF_Internal_attempts = 1
write(102,705) "Internal CONROT moves", real(CONF_Internal_success)*1./real(CONF_Internal_attempts)
write(102, *); write(102,"(A)") "Per Residue probabilities of successful sequence changes"
do i=1, pep_res
    if (SEQ_PerRes_attempts(i).eq.0) SEQ_PerRes_attempts(i) = 1
    write(102,706) "Residue", i, real(SEQ_PerRes_success(i))*1./real(SEQ_PerRes_attempts(i))
enddo
close(102)

!deallocate memory
deallocate(backbone_backup)
deallocate(group_pep)
deallocate(Tgroup_pep)
deallocate(pep_para)
!deallocate(SA_best_pep)
!deallocate(SA_best_para)
deallocate(residuesite)
deallocate(res_starts)
deallocate(res_ends)
deallocate(sc_starts)
deallocate(sc_ends)
deallocate(original_pep)
deallocate(atom)
deallocate(numex_pep)
deallocate(numex4_pep)
deallocate(inb_pep)
deallocate(inb4_pep)
deallocate(atomid_pep)
deallocate(mdcrd_pep)
deallocate(charge_pep) 
deallocate(epsion_pep)
deallocate(r_pep)
deallocate(rborn_pep)
deallocate(fs_pep)
deallocate(dielecons_pep)
deallocate(phis)
deallocate(psis)
deallocate(energy_per_res)
deallocate(res_coords_flags)
deallocate(res_params_flags)
deallocate(aa_groups) 
deallocate(cg_flags)

!Receptor parameters
deallocate(group_rec)
deallocate(rec_para)
deallocate(mdcrd_rec)
deallocate(charge_rec) 
deallocate(epsion_rec)
deallocate(r_rec)
deallocate(rborn_rec)
deallocate(fs_rec)
deallocate(dielecons_rec)

!pairwise distances
deallocate(rij_matrix_pep)
deallocate(rij_matrix_sys)
deallocate(rij_matrix_rec)

!General simulation parameters
deallocate(SEQ_PerRes_attempts)
deallocate(SEQ_PerRes_success)
if (NP_flag .eq. 1) deallocate(LCPO_indices)
deallocate(lbres)

!GB Solvation parameters
if (solvation_flag.eq.1) then
    deallocate(alpha_pep_sys)
    deallocate(alpha_pep_iso)
    deallocate(alpha_rec_sys)
    deallocate(alpha_rec_iso)
    deallocate(integrals_rec_iso)
    deallocate(integrals_comp_rec_from_pep)
    deallocate(integrals_comp_pep_from_rec)
    deallocate(GB_area_calc_flags)
endif

!string formats
6   format(i6, i6, 7f12.3)
703 format(i3, i5, f7.3)
704 format(3f10.3)
705 format(a30, f8.3)
706 format(a8, i3, f8.3)

end program PepBD