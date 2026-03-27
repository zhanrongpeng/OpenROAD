source /root/OpenROAD/test/helpers.tcl
read_liberty Nangate45/Nangate45_typ.lib
read_lef Nangate45/Nangate45.lef
read_def s35932.def
read_sdc s35932.sdc
source Nangate45/Nangate45.rc
set_wire_rc -layer metal3
estimate_parasitics -placement 

report_worst_slack

report_design_area

set tiehi "LOGIC1_X1/Z"
set tielo "LOGIC0_X1/Z"

ord::set_thread_count "3"
restructure -liberty_file Nangate45/Nangate45_typ.lib -target area -abc_logfile results/abc_rcon.log  -tielo_port $tielo -tiehi_port $tiehi -work_dir ./results 

report_worst_slack
report_design_area
