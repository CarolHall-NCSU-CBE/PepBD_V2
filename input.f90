module input

use datatypes
use sys_vars

contains

function file_exists(filename) result(res)
implicit none
character(len=*),intent(in) :: filename
logical                     :: res

! Check if the file exists
inquire( file=trim(filename), exist=res )
end function

subroutine inputfile
implicit none
integer               :: i=0, status, str_status, check_hyd, l_index, r_index, len
character*40          :: name, val
character*80          :: line
logical               :: max_hyd_found

nstep_start = -1
nstep_terminal = -1
pep_res =-1
receptor_name='temp'
ph_value = -1
rmsd_max = -1
rmsd_min = -1
min_hyd(1:6) = (/-1, -1, -1, -1, -1, -1/) !order is glycine, hydrophobic, positive, hydrophilic, other, and negative
max_hyd(1:6) = (/-1, -1, -1, -1, -1, -1/) !order is glycine, hydrophobic, positive, hydrophilic, other, and negative
nstep_start = -1
nstep_terminal = -1
ekt_seq = -1
ekt_backbone = -1

!read input file
open(10, file="input.txt", status="old")
do while (.true.)
    read(10, '(A)', iostat=status) line
    if(status.ne.0) exit !end of file has been reached

    !split entry into the two parts (I do this instead of a list read because LIB_PATH variable may have a "/" character, which is not read properly in a list read)
    len = LEN_TRIM(line) !get length of input, without trailing spaces
    l_index = INDEX(line(1:len), ' ') !get location if first space
    if (l_index.eq.0) then !the values in the line not separated by spaces, so try tabs
        l_index = INDEX(line(1:len), achar(9)) !achar(9) is apparently the ASCII table entry for the tab character
        r_index = INDEX(line(1:len), achar(9), .true.) 
    else
        r_index = INDEX(line(1:len), ' ', .true.) !get location of last space, not considering trailing spaces
    endif

    !make sure the line was properly split
    if ((l_index.eq.0 .or. r_index.eq.0) .and. len.gt.0) then
        write(*,*) 'Unable to parse data in input file!'
        write(*,*) 'Make sure each line has exactly two values separated by either a space or tab'
        write(*,*) 'Terminating program'
    endif

    name = line(1:l_index-1) !get first entry in line
    val = line(r_index+1:len) !get second entry in line

    !BASIC SYSTEM PARAMETERS
    if(name=="SYSTEM_PDB") then
        system_pdb = val   
    elseif(name=="RESTARTPDB") then
        restart_pdb = val
    elseif(name=='LIB_PATH') then
        lib_path = val
    elseif(name=='PDB_PATH') then
        pdb_path = val
    elseif(name=="RESTART_FLAG") then
        read(val,*,iostat=str_status) restart_flag
    elseif (name.eq."RANDOM_SEED") then
        read(val,*,iostat=str_status) ran_seed(1)
    elseif(name=="PEP_RES") then
        read(val,*,iostat=str_status)  pep_res
    elseif(name=="RECEPTOR_NAME") then
        read(val,*,iostat=str_status) receptor_name
    elseif(name=="PH") then
        read(val,*,iostat=str_status) ph_value
    !elseif(name=="BACKBONE_CHANGE_PROB") then
    !    read(val,*,iostat=str_status) backbone_prob

    !basic design parameters
    elseif(name=="SEQ_CHANGE_PROB") then
        read(val,*,iostat=str_status) seq_prob
    elseif(name=="RES_MUTATE_PROB") then
        read(val,*,iostat=str_status) res_mutate_prob
    elseif(name=="MAX_RMSD") then
        read(val,*,iostat=str_status) rmsd_max
    !elseif(name=="MIN_RMSD") then !this variable not used
    !    read(val,*,iostat=str_status) rmsd_min
    elseif(name=="PDB_GEN_FREQ") then
        read(val,*,iostat=str_status) PDB_gen_freq
    elseif(name=="PEP_WEIGHT_FACTOR") then
        read(val,*,iostat=str_status) weighting_factor
    elseif(name=="CONROT_ONE_TERM") then
        read(val,*,iostat=str_status) CONROT_one_term
    elseif(name=="CONROT_TWO_TERM") then
        read(val,*,iostat=str_status) CONROT_two_term
    elseif(name=="INITIAL_RANDOMIZE") then
        read(val,*,iostat=str_status) initial_randomize 
    elseif(name=="STEPS_BEFORE_CONF_CHANGE") then
        read(val,*,iostat=str_status) steps_before_conf_change 
    elseif(name=="END_POINT_FLAG") then
        read(val,*,iostat=str_status) end_point_flag
    !elseif(name=="VDW_CUTOFF") then
    !    read(val,*,iostat=str_status) vdw_cutoff
    !elseif(name=="MIN_RIJ2_RATIO") then
    !    read(val,*,iostat=str_status) min_rij2_ratio
    !elseif (name.eq."OVERLAP_DIST") then
    !    read(val,*,iostat=str_status) overlap_dist
    elseif (name.eq."OVERLAP_TYPE") then
        read(val,*,iostat=str_status) overlap_type
    elseif(name=="DEBUG_FLAG") then
        read(val,*,iostat=str_status) debug_flag
    !elseif(name=="E_SIDECHAIN_OPT_CUT") then
    !    read(val,*,iostat=str_status) e_sidechain_opt_cut
    elseif(name=="SOLVATION_FLAG") then
        read(val,*,iostat=str_status) solvation_flag
    elseif(name=="NEIGHBOR_FLAG") then
        read(val,*,iostat=str_status) neighbor_flag
    !elseif(name=="NEIGHBOR_CUTOFF") then
    !    read(val,*,iostat=str_status) neighbor_cutoff
    elseif(name=="CHECK_INIT_STRUCT") then
        read(val,*,iostat=str_status) check_init_structure

    !excluding certain amino acids from design
    elseif(name=="EXCLUDED_AA") then
        num_ex_AAs = num_ex_AAs + 1
        if (num_ex_AAs.gt.20) then
            !open(20, file="error.txt", access="append")
            write(*,"(A)") "Error in input file!"
            write(*,"(A)") "At most, only 20 amino acids can be excluded in design, but more than 20 excluded AAs are in the input file!"
            write(*,"(A)") "Please fix input file and run again"
            !close(20)
            stop
        endif
        excluded_AAs(num_ex_AAs) = trim(val)

    !parameters for non-polar solvation energy
    elseif(name=="NP_FLAG") then
        read(val,*,iostat=str_status) NP_flag
    elseif(name=="PROBE_RADIUS") then
        read(val,*,iostat=str_status) probe_r
    elseif(name=="SURFACE_TENSION_SOLV") then
        read(val,*,iostat=str_status) surftens_solv
    elseif(name=="SURFACE_TENSION_REC") then
        read(val,*,iostat=str_status) surftens_rec
    
    !parameters for either simulated annealing or metropolis monte carlo sampling
    elseif(name=="SA_FLAG") then
        read(val,*,iostat=str_status) SA_flag

    !simulated annealing parameters
    elseif(name=="SA_INTERVAL") then
        read(val,*,iostat=str_status) SA_interval    
    elseif(name=="SA_FRACTION") then
        read(val,*,iostat=str_status) SA_fraction
    elseif(name=="SA_T_START") then
        read(val,*,iostat=str_status) SA_T_start
    elseif(name=="SA_T_END") then
        read(val,*,iostat=str_status) SA_T_end
    elseif(name=="SA_START_CYCLE") then
        read(val,*,iostat=str_status) SA_start_cycle
    
    !basic monte carlo parameters
    elseif(name=="STEP_START") then
        read(val,*,iostat=str_status) nstep_start
    elseif(name=="STEP_END") then
        read(val,*,iostat=str_status) nstep_terminal
    !elseif(name=="STEPS_PER_CYCLE") then
    !    read(val,*,iostat=str_status) steps_per_cycle
    elseif(name=="KT_SEQUENCE") then
        read(val,*,iostat=str_status) ekt_seq
    elseif(name=="KT_BACKBONE") then
        read(val,*,iostat=str_status) ekt_backbone
    !elseif(name=="KT_RIGID_RATIO") then
    !    read(val,*,iostat=str_status) ekt_rigid_ratio

    !hydration property parameters
    elseif(name=="NUM_GLY_MIN") then
        read(val,*,iostat=str_status) min_hyd(1)
    elseif(name=="NUM_PHOBIC_MIN") then
        read(val,*,iostat=str_status) min_hyd(2)
    elseif(name=="NUM_POS_MIN") then
        read(val,*,iostat=str_status) min_hyd(3)
    elseif(name=="NUM_POLAR_MIN") then
        read(val,*,iostat=str_status) min_hyd(4)
    elseif(name=="NUM_OTHER_MIN") then
        read(val,*,iostat=str_status) min_hyd(5)
    elseif(name=="NUM_NEG_MIN") then
        read(val,*,iostat=str_status) min_hyd(6)
    elseif(name=="NUM_GLY_MAX") then
        read(val,*,iostat=str_status) max_hyd(1)
    elseif(name=="NUM_PHOBIC_MAX") then
        read(val,*,iostat=str_status) max_hyd(2)
    elseif(name=="NUM_POS_MAX") then
        read(val,*,iostat=str_status) max_hyd(3)
    elseif(name=="NUM_POLAR_MAX") then
        read(val,*,iostat=str_status) max_hyd(4)
    elseif(name=="NUM_OTHER_MAX") then
        read(val,*,iostat=str_status) max_hyd(5)
    elseif(name=="NUM_NEG_MAX") then
        read(val,*,iostat=str_status) max_hyd(6)
    elseif(name=="TRP_LIMIT") then
        read(val,*,iostat=str_status) trp_limit

    !structural parameters
    elseif(name=="NUM_HELIX") then
        read(val,*,iostat=str_status) helix_num
        if(helix_num.gt.maximum_helix_num) then
            !open(20, file="error.txt", access="append")
                write(*,"(A)") "Error in input file!"
                write(*,"(A,I0,A,I0)") "helix_num=", helix_num, "is more than the maximum helix number", maximum_helix_num
                write(*,"(A)") "Please adjust the maximum_helix_num in the code!"
            !close(20)
            stop
        endif
    elseif(name=="HELIX_START") then
        read(val,*,iostat=str_status) helix_start(i)
    elseif(name=="HELIX_END") then
        read(val,*,iostat=str_status) helix_end(i)
        i = i + 1
  
    !stuff nno longer used or used in other programs associated with PepBD
    
    !elseif(name=="PEP_PDB") then
    !    pep_pdb = val
    !elseif(name=="REC_PDB") then
    !    rec_pdb = val    
    !elseif(name=="NUM_RESTARTS") then
    !    read(val,*,iostat=str_status) num_restarts
    !elseif(name=="RANDOMIZE_LENGTH") then
    !    read(val,*,iostat=str_status) randomize_length   
    !elseif(name=="STEPS_BEFORE_RANDOMIZE") then
    !    read(val,*,iostat=str_status) steps_before_randomize    
    ! elseif(name=="SMART_RESTART_FLAG") then
    !    read(val,*,iostat=str_status) smart_restart_flag
    !elseif(name=="SMART_RESTART_LEN") then
    !    read(val,*,iostat=str_status) smart_restart_len    
    !elseif(name=="SMART_RESTART_OFFSET") then
    !    read(val,*,iostat=str_status) smart_restart_offset
    !elseif(name=="SMART_RESTART_WAIT_LEN") then
    !    read(val,*,iostat=str_status) smart_restart_wait_len
    !elseif(name=="SMART_RESTART_SCORE_METHOD") then
    !    read(val,*,iostat=str_status) smart_restart_score_method
    !elseif(name=="DECOMP_MUTATE_FLAG") then
    !    read(val,*,iostat=str_status) decomp_mutate_flag
    !elseif(name=="DECOMP_SCORE_TYPE") then
    !    read(val,*,iostat=str_status) decomp_score_type
    !elseif(name=="ROTATION_SWITCH") then
    !    read(val,*,iostat=str_status) rotation_switch
    !elseif(name=="TRANSLATE_Z_FLAG") then
    !    read(val,*,iostat=str_status) translate_z_flag
    !elseif(name=="MIN_ZDIST_FROM_SURFACE") then
    !    read(val,*,iostat=str_status) min_zdist_from_surface
    !elseif(name=="MAX_ZDIST_FROM_SURFACE") then
    !    read(val,*,iostat=str_status) max_zdist_from_surface
    !elseif(name=="TABU_FLAG") then
    !    read(val,*,iostat=str_status) TABU_flag
    !elseif(name=="TABU_BAN_LENGTH") then
    !    read(val,*,iostat=str_status) TABU_ban_length
    !elseif(name=="MOVE_PEP_OVER_SURFACE") then
    !    read(val,*,iostat=str_status) move_pep_over_surface_flag
    !elseif(name=="PEP_OVER_SURFACE_Z") then
    !    read(val,*,iostat=str_status) target_pep_distance

    !amino acids to exclude from design

    !elseif(name=="INSERT_PEP_RIGOR") then
    !    read(val,*,iostat=str_status) insert_pep_rigor
    !elseif(name=="NUM_BACKBONES") then
    !    read(val,*,iostat=str_status) num_backbones

    !backbone sampling inputs (not used by main PepBD program)
    !elseif (name.eq."X_BUFFER") then
    !    read(val,*,iostat=str_status) x_buffer
    !else if (name.eq."Y_BUFFER") then
    !    read(val,*,iostat=str_status) y_buffer
    !elseif (name.eq."ITERATIONS") then
    !elseif (name.eq."Z_MIN") then
    !    read(val,*,iostat=str_status) z_min 
    !elseif (name.eq."Z_MAX") then
    !    read(val,*,iostat=str_status) z_max
    !elseif(name=="EVAL_TWOBODY") then
    !    read(val,*,iostat=str_status) eval_twobody_flag
    !elseif(name=="EVAL_ONEBODY") then
    !    read(val,*,iostat=str_status) eval_onebody_flag

    else
        !open(20, file="error.txt", access="append")
        write(*,"(A)") "Error in input file!"
        write(*,"(A,A,A)") "Unrecognized parameter ", name, " in input file, terminating program."
        !close(20)
        stop
    endif
enddo

!!!!!!check the input values don't have problems !!!!!!!!!!!
!open(20, file="error.txt", access="append")

!check start/stop positions were given for each helix
if (i.ne.helix_num) then
    write(*,"(A)") "Error in input file!"
    write(*,"(A)") "Incorrect number of helix start/stop positions given, terminating program."
    write(*,"(A,I0,A,I0,A)") "i is ", i, " while there are ", helix_num, "helices"
    !close(20)
    stop
endif

!check the input pdbs exist
!
!INQUIRE(file=rec_pdb, exist=rec_found)
!INQUIRE(file=system_pdb, exist=system_found)

!if (system_found.eqv..FALSE. .and. (pep_found.eqv..FALSE. .or. rec_found.eqv..FALSE.)) then
INQUIRE(FILE = system_pdb, exist=system_found)
if (system_found) then
    write(*,"(A)") "Found pdb file for the complex!"
else
    write(*,"(A)") "Error in input file!"
    write(*,"(A)") "Could not find system pdb file - please check file path."
    write(*,"(A)") "Terminating program."
    !close(20)
    stop
endif

!if (system_found.eqv..TRUE. .and. pep_found.eqv..TRUE. .or. rec_found.eqv..TRUE.) then
!    write(*,"(A)") "Error in input file!"
!    write(*,"(A)") "PDB files found for both complex and isolated peptide/receptor"
!    write(*,"(A)") "Please only provide one option and restart design"
!    write(*,"(A)") "Terminating program."
!    !close(20)
!    stop
!endif

!check the receptor type specified
if (receptor_name == 'temp') then
    write(*,"(A)") "Error in input file!"
    write(*,"(A)") "Receptor type not specified in input file (parameter: RECEPTOR_NAME). Please add to input file and restart"
    write(*,"(A)") "Terminating program."
    !close(20)
    stop
endif

!check number of peptide residues is specified (only needed if complex pdb provided)
if (pep_res.eq.-1) then
    write(*,"(A)") "Error in input file!"
    write(*,"(A)") "Please specify number of peptide residues in the input file and restart."
    write(*,"(A)") "Terminating program."
    !close(20)
    stop
endif

!check the pH value
if (pH_value.eq.-1) then
    write(*,"(A)") "Possible issue in input file!"
    write(*,"(A)") "pH of system not specified, setting pH of system to 7."
    !close(20)
    pH_value=7
endif

!check the rmsd_max
if (rmsd_max.eq.-1) then
    write(*,"(A)") "Possible issue in input file!"
    write(*,"(A)") "Max RMSD not specified, will set max RMSD to 1000"
    !close(20)
    rmsd_max=1000
endif

!check the hydration property specifications
max_hyd_found = .TRUE.
do i=1, 6
    if (min_hyd(i).eq.-1 .or. max_hyd(i).eq.-1) then
        write(*,"(A)") "Possible issue in input file!"
        write(*,"(A)") "Minimum/maximum hydration constraints not set for all amino acid types"
        write(*,"(A)") "All values not provided will be set the max limits (0 or minima, peptide length for maxima)"
        !close(20)
        max_hyd_found= .FALSE.
        exit
    endif
enddo

if (max_hyd_found.eqv..TRUE.) then
    check_hyd = max_hyd(1) + max_hyd(2) + max_hyd(3) + max_hyd(4) + max_hyd(5) + max_hyd(6) 
    if (check_hyd < pep_res) then
        write(*,"(A)") "Error in input file!"
        write(*,"(A,I0,A,I0,A)") "The sum of the maximum hydration properties, ", check_hyd, "is less than the total number of peptide amino acids, ", pep_res, ", so design is not possible!"
        write(*,"(A)") "Fix the maximum hydration values so it at least equals the number of amino acids in the peptide"
        !close(20)
        stop
    endif
endif 

!check the restart parameters make sense
!if (num_restarts.eq.-1) then 
!    if (nstep_terminal.eq.-1.and.SA_flag.eq.0) then
!        write(*,"(A)") "Error in input file!"
!        write(*,"(A)") "Please either provide restart parameters (NUM_RESTARTS and STEPS_PER_CYCLE)."
!        write(*,"(A)") "or length of single run (STEP_END). Ending program"
!        !close(20)
!        stop
!    else
!        num_restarts=1 
!        steps_per_cycle = nstep_terminal - nstep_start
!    endif
!endif

!if (num_restarts.gt.1.and.steps_per_cycle.eq.-1) then
!    write(*,"(A)") "Error in input file!"
!    write(*,"(A)") "Trying to use restart method, but number of steps before restarting not specified!"
!    write(*,"(A)") "Add STEPS_PER_CYCLE parameter to input.txt file"
!    !close(20)
!    stop
!endif

!check the move probabilities are usable
if (CONROT_two_term.ge.CONROT_one_term) then
    write(*,"(A)") "Error in input file!"
    write(*,"(A)") "CONROT_two_term probability is larger than CONROT_one_term! Please fix input file so CONROT_one_term is larger"
    !close(20)
    stop
endif

!check mutation probabilities make sense
if(seq_prob.lt.0.0 .or. seq_prob.gt.1) then
    write(*,"(A)") "Error in input file!"
    write(*,"(A)") "SEQ_PROB must be between 0 and 1 (inclusive). Please fix the input file."
    stop
else
    seq_switch = seq_prob
    backbone_switch=1. 
endif
!if (backbone_prob + seq_prob > 1.00001) then
!    write(*,"(A)") "Error in input file!"
!    write(*,"(A)") "Sum of BACKBONE_PROB and SEQ_PROB cannot be greater than 1! Please fix the input file"
!    !close(20)
!    stop
!elseif(backbone_prob.lt.0.0.or.seq_prob.lt.0.0) then
!    write(*,"(A)") "Error in input file!"
!    write(*,"(A)") "Neither BACKBONE_PROB nor SEQ_PROB can be less than 0! Please fix the input file"
!    !close(20)
!    stop
!else
!    seq_switch = 1.0-seq_prob
!    backbone_switch = seq_switch - backbone_prob
!endif

!check the smart restart inputs
!if (smart_restart_flag.eq.1.and.smart_restart_len.gt.steps_per_cycle) then
!    write(*,"(A)") "Error in input file!"
!    write(*,"(A)") "value of SMART_RESTART_LEN cannot be greater than STEPS_PER_CYCLE! Please update input file and rerun"
!    !close(20)
!    stop
!endif
!
!if (smart_restart_score_method.lt.1.or.smart_restart_score_method.gt.3) then
!    write(*,"(A)") "Error in input file!"
!    write(*,"(A)") "value of SMART_RESTART_SCORE_METHOD not usable! This value can only be 1, 2, or 3. Please update input file and rerun"
!    !close(20)
!    stop
!endif

!check the simulated annealing inputs
if (SA_flag.eq.1) then
    if (SA_T_start.le.SA_T_end) then  
        write(*,"(A)") "Error in input file!"      
        write(*,"(A)") "Starting temperature for simulated annealing is smaller than end temperature! Please update input file and rerun"
        !close(20)
        stop
    endif
    if (smart_restart_flag.eq.1) then
        write(*,"(A)") "Error in input file!"
        write(*,"(A)") "Currently do not allow both simulated annealing and random restart."
        write(*,"(A)") "Please update input file to remove one of these options and rerun"
        !close(20)
        stop
    end if

    if (SA_interval.le.10) then
        write(*,"(A)") "Error in input file!"
        write(*,"(A)") "Simulated annealing interval length should be greater than 10 to get reasonable results."
        write(*,"(A)") "Please update input file and rerun"
        !close(20)
        stop
    endif

    if (SA_fraction.ge.1.or.SA_fraction.lt.0.5) then
        write(*,"(A)") "Error in input file!"
        write(*,"(A)") "Simulated annealing fraction should be between 0.5 and 1 to get reasonable results"
        write(*,"(A)") "Please update input file and rerun"
        !close(20)
        stop
    endif

    !calculate number of intervals needed for annealing
    SA_num_cycles = int((log(SA_T_end) - log(SA_T_start))/log(SA_fraction)) + 1
    steps_per_cycle = SA_interval
else
    !check the parameters for regular metropolis Monte Carlo
    if (nstep_start.eq.-1) then
        write(*,"(A)") "STEP_START not specified in input file, setting to 1"
        nstep_start = 1
    endif
    if (nstep_terminal.eq.-1) then
        write(*,"(A)") "STEP_STOP not specified in input file, setting to 1000"
        nstep_terminal = 1000
    endif
    if (ekt_seq.eq.-1) then
        write(*,"(A)") "EKT_SEQ not specified in input file, setting to 1"
        ekt_seq = 1
    endif
    if (ekt_backbone.eq.-1) then
        write(*,"(A)") "EKT_BACKBONE not specified in input file, setting to 1"
        ekt_backbone = 1
    endif

    !to work with current PepBD.f90 file, set SA values based on input 
    SA_start_cycle = 1
    SA_num_cycles = 1
    steps_per_cycle = nstep_terminal - nstep_start + 1
endif

!other miscellaneous checks
!if (translate_z_flag.ne.0.and.translate_z_flag.ne.1) then
!    write(*,"(A)") "Error in input file!"
!    write(*,"(A)") "TRANSLATE_Z_FLAG can only be 0 or 1. Please update input file and rerun."
!    !close(20)
!    stop
!endif

!if (min_zdist_from_surface.lt.0) then
!    write(*,"(A)") "Error in input file!"
!    write(*,"(A)") "MIN_ZDIST_FROM_SURFACE should not be less than 0, as this causes steric overlaps and results in many translations being rejected."
!    write(*,"(A)") "Please update input file to increase value above 0 and rerun."
!    !close(20)
!    stop
!endif

if (ph_value.lt.6.0.or.ph_value.gt.8.3) then
    write(*,"(A)") "Error in input file!"
    write(*,"(A)") "Input pH is outside the range of 6.0 - 8.3, and the program may not work in this region."
    write(*,"(A)") "This will be fixed in the future. Terminating program."
    stop
endif

!if (decomp_mutate_flag.eq.1.and.(decomp_score_type.ne."REG".and.decomp_score_type.ne."INV")) then
!    write(*,"(A)") "Error in input file!"
!    write(*,"(A)") "When using per-residue energy composition for mutations, only valid options for DECOMP_SCORE_TYPE are INV or REG."
!    write(*,"(A)") "Input value of ", decomp_score_type, " not allowed!"
!    write(*,"(A)") "Please update input file and rerun."
!    !close(20)
!    stop
!endif

!if (TABU_flag.eq.1.and.decomp_mutate_flag.eq.1) then
!    write(*,"(A)") "Error in input file!"
!    write(*,"(A)") "Currently cannot use both TABU search and select residues for mutation based on energy decomposition, as residues will not be selected properly"
!    write(*,"(A)") "Please update input file to use only one of these options and rerun."
!    !close(20)
!    stop
!endif

if(end_point_flag.ne.0 .and. end_point_flag.ne.1) then
    write(*,"(A)") "Error in input file!"
    write(*,"(A,I1)") "end_point_flag needs to be either 0 or 1, but input value is ", end_point_flag
    write(*,"(A)") "Please update input file and rerun."
    !close(20)
    stop
endif

close(20)

!neighbor cutoff terms
!neighbor_cutoff2 = neighbor_cutoff * neighbor_cutoff

!setup reference temperatures
ekt_backbone = ekt_seq * ekt_backbone_ratio
ekt_rigid = ekt_seq * ekt_rigid_ratio

!get amino acid names based on pH of design and excluded amino acid list
call get_AA_list()

!setup the TABU lists
!if (TABU_flag.eq.1) then
!    call initialize_TABU(1)
!endif

return
end subroutine inputfile

!subroutine initialize_TABU(flag)
!implicit none
!integer    :: i,flag
!
!if (flag.eq.1) then !only need to allocate memory the first time we initialize TABU list
!    allocate(TABU_ban_lengths(pep_res))
!    allocate(TABU_allowed_residues(pep_res))
!endif
!
!TABU_ban_lengths = 0
!do i=1,pep_res
!    TABU_allowed_residues(i) = i
!enddo
!
!TABU_allowed_length = pep_res
!end subroutine initialize_TABU

subroutine fix_hyd_props
implicit none
integer    :: i

do i=1, 6
    if (min_hyd(i).eq.-1) min_hyd(i) = 0
    if (max_hyd(i).eq.-1) max_hyd(i) = pep_res
enddo

end subroutine fix_hyd_props

subroutine get_AA_list()
implicit none
integer        :: i, j, k, start_limit
character*3    :: swap

!order of amino acid types: glycine, hydrophobic, positive, hydrophilic, other, and negative
!create initial array of available amino acids, if all can be used during design
AAs_available(1,1) = 'GLY'
AAs_available(1:7,2) = (/'ILE', 'LEU', 'MET', 'PHE', 'TRP', 'TYR', 'VAL'/)
AAs_available(1:2,3) = (/'ARG', 'LYS'/)
AAs_available(1:5,4) = (/'SER', 'THR', 'ASN', 'GLN', 'HIE'/)
AAs_available(1:3,5) = (/'ALA', 'CYS', 'PRO'/)
AAs_available(1:2,6) = (/'ASP', 'GLU'/)
num_AA_types = (/1, 7, 2, 5, 3, 2/)
num_AAs_available(1:6) = (/1, 7, 2, 5, 3, 2/)

!now go through the exclusion list and remove them from the available amino acids
do i=1,num_ex_AAs
    do j=1, 6
        start_limit = num_AAs_available(j)
        do k=1, start_limit
            if (trim(excluded_AAs(i)).eq.trim(AAs_available(k,j))) then
                swap = AAs_available(k,j)
                AAs_available(k,j) = AAs_available(num_AAs_available(j), j)
                AAs_available(num_AAs_available(j), j) = swap
                num_AAs_available(j) = num_AAs_available(j) - 1
                exit
            endif
        enddo
    enddo
enddo

!move tryptophan to end of hydrophobic amino acid list, to make easier to maintain limits on tryptophan
do i=1, num_AAs_available(2)
    if (AAs_available(i,2).eq."TRP") then
        AAs_available(i,2) = AAs_available(num_AAs_available(2),2)
        AAs_available(num_AAs_available(2),2) = "TRP"
        exit
    endif
enddo

!Lastly, update names of residues based on the pH of the system
!!NOTE: This will be done in the future; not hard to do, but not a concern at the moment

end subroutine get_AA_list

subroutine prep_output
implicit none
integer     :: stat, i

!Delete files from previous runs, if they exist and are not restarting design
if (restart_flag.eq.0) then
    open(5, file="energydetails.txt", status="old", iostat=stat)
    if (stat == 0) close(5, status="delete")

    open(5, file="minimum_energy.txt",status="old", iostat=stat)
    if (stat == 0) close(5, status="delete")

    open(5, file="randomnumber.txt", status="old", iostat=stat)
    if (stat == 0) close(5, status="delete")

    open(5, file="rmsd.txt", status="old", iostat=stat)
    if (stat == 0) close(5, status="delete")

        open(5, file="error.txt", status="old", iostat=stat)
    if (stat == 0) close(5, status="delete")

    open(5, file="output.txt", status="old", iostat=stat)
    if (stat == 0) close(5, status="delete")

    open(5, file="energydecomp.txt", status="old", iostat=stat)
    if (stat == 0) close(5, status="delete")

    !Add header info to a few output files
    open(5, file="energyprofile.txt", status="replace")
            write(5,"(A)") "Cycle      Step    E_Net      E_VDW        E_ELE       E_SGB       E_SASA      E_Pep"
    close(5)

    !open(5, file="energydecomp.txt", status="replace")
    !    write(5,"(A)") "Resname   VDW    ELE    EGB    NET"
    !close(5)
endif

!Write parameters used for peptide design
open(5, file='output.txt', access='append')

write(5,"(A)") "######## PEPBD DESIGN PARAMETERS ########"
write(5,"(A,I0)") 'gnum: ', gnum
write(5,"(A,I0)") 'atom_num_pep: ', atom_num_pep
write(5,"(A,I0)") 'atom_num_rec: ', atom_num_rec
write(5,"(A,I0)") 'nstep_start: ', nstep_start
write(5,"(A,I0)") 'nstep_terminal: ', nstep_terminal
write(5,"(A,I0)") 'pep_res: ', pep_res
write(5,"(A,I0)") 'fragmentnum: ', fragmentnum
write(5,"(A,I0)") 'helix_num: ', helix_num
do i=1,helix_num
    write(5,"(A,I0,A,I0)") 'helix_start: ', helix_start(i), 'helix_end: ', helix_end(i)
enddo
write(5,"(A,I0)") 'steps_before_conf_change: ', steps_before_conf_change
write(5,*) 'ekt_seq: ', ekt_seq
write(5,*) 'ekt_backbone_ratio: ', ekt_backbone_ratio 
write(5,*) 'ekt_backbone: ', ekt_backbone
write(5,*) 'ekt_rigid_ratio: ', ekt_rigid_ratio
write(5,*) 'ekt_rigid: ', ekt_rigid
write(5,*) 'ekt_scmf: ', ekt_scmf
write(5,*) 'ph_value: ', ph_value 
write(5,*) 'lampda: ', lampda
write(5,*) 'surftens_rec: ', surftens_rec
write(5,*) 'surftens_solv: ', surftens_solv
write(5,*) 'SASA_offset: ', SASA_offset
write(5,*) 'weighting_factor: ', weighting_factor
write(5,*) 'dihedral_weighting_factor: ', dihedral_weighting_factor
write(5,*) 'vdw14_coeff: ', vdw14_coeff
write(5,*) 'ele14_coeff: ', ele14_coeff
write(5,"(A,I0)") 'PDB_gen_freq: ', PDB_gen_freq
write(5,"(A,I0)") 'rmsd_flag: ', rmsd_flag
write(5,*) 'delta_phi_min: ', delta_phi_min 
write(5,*) 'delta_phi_max: ', delta_phi_max
write(5,*) 'delta_phi_step: ', delta_phi_step
write(5,*) 'translate_max: ', translate_max
write(5,*) 'translate_min: ', translate_min
write(5,*) 'rotate_max: ', rotate_max
write(5,*) 'max_dev_from_center: ', max_dev_from_center
!write(5,*) 'num_top_layer:', num_top_layer
!write(5,*) 'translate_z_flag:', translate_z_flag
write(5,*) 'min_zdist_from_surface: ', min_zdist_from_surface
write(5,*) 'max_zdist_from_surface: ', max_zdist_from_surface
write(5,*) 'backbone_prob: ', backbone_prob
write(5,*) 'seq_prob: ', seq_prob
write(5,*) 'rigid_body_prob: ', 1 - backbone_prob - seq_prob
write(5,*) 'res_mutate_prob: ', res_mutate_prob
write(5,*) 'CONROT_two_term: ', CONROT_two_term
write(5,*) 'CONROT_one_term: ', CONROT_one_term
write(5,*) 'rotation_switch: ', rotation_switch
write(5,"(A,I0)") 'trp_limit: ', trp_limit
write(5,"(A,I0)") 'glycine_min: ', min_hyd(1)
write(5,"(A,I0)") 'hydrophobic_min: ', min_hyd(2)
write(5,"(A,I0)") 'positive_min: ', min_hyd(3)
write(5,"(A,I0)") 'hydrophilic_min: ', min_hyd(4)
write(5,"(A,I0)") 'other_min: ', min_hyd(5)
write(5,"(A,I0)") 'negative_min: ', min_hyd(6)
write(5,"(A,I0)") 'glycine_max: ', max_hyd(1)
write(5,"(A,I0)") 'hydrophobic_max: ', max_hyd(2)
write(5,"(A,I0)") 'positive_max: ', max_hyd(3)
write(5,"(A,I0)") 'hydrophilic_max: ', max_hyd(4)
write(5,"(A,I0)") 'other_max: ', max_hyd(5)
write(5,"(A,I0)") 'negative_max: ', max_hyd(6)
write(5,"(A,I0)") 'negative_max: ', max_hyd(6)
!write(5,*) 'num_restart:', num_restarts
!if (num_restarts.gt.1) then
!    write(5,*) 'smart_restart_flag:', smart_restart_flag
!    if (smart_restart_flag.eq.1) then
!        write(5,*) 'smart_restart_len:', smart_restart_len
!        write(5,*) 'smart_restart_offset:', smart_restart_offset
!        write(5,*) 'smart_restart_wait_len:', smart_restart_wait_len
!        write(5,*) 'smart_restart_score_method:', smart_restart_score_method
!    endif
!    write(5,*) 'steps_per_cycle:', steps_per_cycle
!    write(5,*) 'randomize_length:', randomize_length
!endif
write(5,"(A,I0)") 'SA_flag:', SA_flag
if (SA_flag.eq.1) then
    write(5,"(A,I0)") 'SA_interval: ', SA_interval
    write(5,*) 'SA_fraction: ', SA_fraction
    write(5,*) 'SA_T_start: ', SA_T_start
    write(5,*) 'SA_T_end: ', SA_T_end
    write(5,"(A,I0)") 'SA_num_cycles: ', SA_num_cycles
endif
!write(5,"(A,I0)") 'decomp_mutate_flag:', decomp_mutate_flag
!if (decomp_mutate_flag.eq.1) then
!    write(5,*) 'decomp_score_type:', decomp_score_type
!endif
write(5,*) 'probe_r:', probe_r
write(5,"(A,I0)") 'max_neighbors: ', max_neighbors
write(5,"(A,I0)") 'neighbor_flag: ', neighbor_flag
!if (neighbor_flag.eq.1) then 
!    write(5,*) 'neighbor_cutoff:', neighbor_cutoff
!    write(5,*) 'neighbor_length:', neighbor_length
!endif
!write(5,"(A,I0)") 'TABU_flag :', TABU_flag 
!if (TABU_flag.eq.1) then
!    write(5,*) 'TABU_ban_length:', TABU_ban_length
!endif
write(5,"(A,I0)") 'restart_flag: ', restart_flag
!write(5,"(A,I0)") 'move_pep_over_surface_flag:', move_pep_over_surface_flag
!if (move_pep_over_surface_flag.eq.1) then
!    write(5,*) 'target_pep_distance:', target_pep_distance
!endif
write(5,*) 'receptor_name: ', receptor_name
write(5,*) 'system_pdb: ', system_pdb
!if (system_found.eqv..TRUE.) then 
    !write(5,*) 'system_pdb:', system_pdb
!else
!    write(5,*) 'pep_pdb:', pep_pdb
!    write(5,*) 'rec_pdb:', rec_pdb
!endif
if (restart_flag.eq.1) write(5,*) 'restart_pdb: ', restart_pdb
write(5,"(A,I0)") 'end_point_flag: ', end_point_flag
write(5,"(A,I0)") 'num_ex_AAs: ', num_ex_AAs
do i=1, num_ex_AAs
    write(5,*) 'Excluded AA: ', excluded_AAs(i)
enddo
write(5,*) 'vdw_cutoff: ', vdw_cutoff
write(5,*) 'min_rij2_ratio: ', min_rij2_ratio
write(5,"(A,I0)") 'check_init_structure: ', check_init_structure
write(5,*)
write(5,*) "######## STARTING PEPTIDE DESIGN ##########"
write(5,*)

end subroutine prep_output
  
end module input