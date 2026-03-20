read_lib -v Nangate45/Nangate45_typ.lib
read_blif -n ./results/gcd_crit_path.blif
strash
read_coords ./results/gcd_coords.txt
propagate_coords
write_blif ./results/gcd_test_out.blif
quit
