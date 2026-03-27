source helpers.tcl
read_liberty Nangate45/Nangate45_typ.lib
read_lef Nangate45/Nangate45.lef
read_def gcd_placed.def
read_sdc gcd.sdc

set db [ord::get_db]
set block [[$db getChip] getBlock]
set inst [lindex [$block getInsts] 0]
puts "AFTER read_def - Instance: [$inst getName]"
foreach iterm [$inst getITerms] {
    puts "  - [$iterm getName] connected=[$iterm isConnected] net=[$iterm getNet]"
}
exit
