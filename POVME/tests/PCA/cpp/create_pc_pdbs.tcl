#!/usr/bin/tclsh

mol new common_points.xyz

for {set i 0} {$i<100} {incr i} {
	set filename "principal_component"
	append filename $i ".vect"
	set fp [open $filename r]
	set file_data [read $fp]
	close $fp
	set sel [atomselect 0 all]
	$sel set beta $file_data
	set output_filename "principal_component"
	append output_filename $i ".pdb"
	$sel writepdb $output_filename
}

exit
