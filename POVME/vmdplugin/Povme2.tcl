#!/usr/bin/tcl

#Name:
# POVME GUI
#Synopsis:
# A GUI for the pocket volume program POVME2
#Version 
# 1.0
#

###################################
#package require Tcl 8.2
#package require struct::matrix 2.0.1

package provide povme2 1.0
package require Tk

set PI 3.14159265

variable w
namespace eval ::povme2:: {
  variable w
  variable mainwin
  variable sep "/"
  variable help_url "http://nbcr.ucsd.edu/data/sw/hosted/POVME/"
  
  variable load_into_vmd 1 ;# whether the volume map is loaded into VMD at the end of the POVME run
  variable plot_volumes 1 ;# whether the volume results are plotted after running a trajectory on VMD
  variable arglist {}
  array set descarray {} ;# a description of all arguments
  array set min_array {}
  array set max_array {}
  
  variable workdir "./" ;# the directory where we execute from 
  variable tempdir $workdir
  
  variable pdb_filename ""
  variable temp_pdb "${tempdir}povme2_tmp.pdb"
  
  variable inclusion_list {}
  variable exclusion_list {}
  
  # povme2 options & temporary variables
  # IMPORTANT NOTE ABOUT VARIABLES BELOW:
  #  Any variable that ends with _temp, _format, or _caption are variables that are changed by the user
  #  _temp stores the value that the user has inputted, where if they click "OK", then it saves it into
  # the variable itself, but if they click "cancel", that the _temp value is discarded. The _format 
  # variable can take values "file", "checkbox", "list", "slider", etc. The _format variable tells the 
  # TK interpreter what type of widget to place there. Finally, the _caption variable is merely what
  # is displayed in a label above the input.
  #  A blank or nonexistent value for _format defaults to a text entry. 
  
  # Files menu properties
  #variable workdir "./" ;# declaration above
  variable workdir_temp "$workdir"
  variable workdir_format "file"
  variable workdir_caption "Working Directory"
  #variable tempdir $workdir ;# declaration above
  variable tempdir_temp $workdir ;# this is not in the files menu but has to actually be here
  
  variable povme2_directory_format "file"
  variable povme2_directory_caption "Temporary Directory"
  
  variable povme2_directory "./POVME2.py" ;# depending on whether the user has 
  variable povme2_directory_temp ""
  variable povme2_directory_format "file"
  variable povme2_directory_caption "POVME2 location"
  
  variable python_executable "python"
  variable python_executable_temp ""
  #variable python_executable_format "file"
  variable python_executable_caption "Python Executable"
  # Point menu Properties
  variable grid_point_spacing "1.0"
  variable grid_point_spacing_temp ""
  variable grid_point_spacing_caption "Grid Point Spacing:"
  
  variable distance_cutoff "1.09"
  variable distance_cutoff_temp ""
  variable distance_cutoff_caption "Distance Cutoff:"
  
  variable make_point_field_pdb "0"
  variable make_point_field_pdb_temp ""
  variable make_point_field_pdb_format "checkbox"
  variable make_point_field_pdb_caption "Make Point-Field PDB (No Volume Calc):"
  # Contiguous Point/Convex-Hull Exclusion menu properties
  variable use_contiguous_points "1"
  variable use_contiguous_points_temp ""
  variable use_contiguous_points_format "checkbox"
  variable use_contiguous_points_caption "Use Contiguous Points:"
  
  variable exclude_points_outside_convex_hull "1"
  variable exclude_points_outside_convex_hull_temp ""
  variable exclude_points_outside_convex_hull_format "checkbox"
  variable exclude_points_outside_convex_hull_caption "Exclude Points Outside Convex Hull:"
  #  need to be able to select contiguous points location
  variable contiguous_point_criteria "3"
  variable contiguous_point_criteria_temp ""
  variable contiguous_point_criteria_caption "Contiguous Point Criteria:"
  # Output menu properties
  variable output_filename_prefix "./"
  variable output_filename_prefix_temp ""
  variable output_filename_prefix_format "file"
  variable output_filename_prefix_caption "Output Filename Prefix:"
  
  variable compress_output "0"
  variable compress_output_temp ""
  variable compress_output_format "checkbox"
  variable compress_output_caption "Compress Output:"
  
  variable separate_volume_pdbs "0"
  variable separate_volume_pdbs_temp ""
  variable separate_volume_pdbs_format "checkbox"
  variable separate_volume_pdbs_caption "Separate Volume PDBs:"
  
  variable volume_trajectory "0"
  variable volume_trajectory_temp ""
  variable volume_trajectory_format "checkbox"
  variable volume_trajectory_caption "Volume Trajectory:"
  
  variable equal_num_points_per_frame "0"
  variable equal_num_points_per_frame_temp ""
  variable equal_num_points_per_frame_format "checkbox"
  variable equal_num_points_per_frame_caption "Equal # of Points per Frame:"
  
  variable tabbed_volume_file "1"
  variable tabbed_volume_file_temp ""
  variable tabbed_volume_file_temp_format "checkbox"
  variable tabbed_volume_file_temp_caption "Tabbed Volume File"
  
  variable volumetric_density_file "1"
  variable volumetric_density_file_temp ""
  variable volumetric_density_file_format "checkbox"
  variable volumetric_density_file_caption "Volumetric Density File"
  
  variable colored_density_files "1"
  variable colored_density_files_temp ""
  variable colored_density_files_format "checkbox"
  variable colored_density_files_caption "Colored Density Files"
  
  variable num_processors "1"
  variable num_processors_temp ""
  variable num_processors_caption "Number of Processors:"
  
  variable disk_instead_of_memory "0"
  variable disk_instead_of_memory_temp ""
  variable disk_instead_of_memory_format "checkbox"
  variable disk_instead_of_memory_caption "Disk Instead of Memory:"
  
  variable python_executable "python"
  variable python_executable_temp ""
  variable python_executable_caption "Python Executable:"
  
  # NOTE about arglist variables: Each of these variables represents a category that appears in the drop down menus
  # above the main window. Notice that the each of the user input options that we declared in the above variables
  # appears in one of these lists below. If you want to add a new user input to one of these items, you need to 
  # declare the variable above, as well as its associated _temp, _format, and _caption variables. Then you need
  # to add it to one of the lists below. Finally, you will need to make sure that the variable is written to the
  # make_input_file function below.
  variable arglist_files "povme2_directory python_executable tempdir workdir"
  variable arglist_point_properties "grid_point_spacing distance_cutoff make_point_field_pdb"
  variable arglist_contiguous_points "use_contiguous_points exclude_points_outside_convex_hull contiguous_point_criteria"
  variable arglist_output "output_filename_prefix compress_output separate_volume_pdbs volume_trajectory equal_num_points_per_frame tabbed_volume_file volumetric_density_file colored_density_files num_processors disk_instead_of_memory"
  variable arglist_other ""
  
  # User inputs on the main window; involving the selection of molecules and atomselections
  variable moltxt "(none)"
  variable which_mol -1
  variable candidate_mol -1
  variable selection "all"
  
  # Variables concerning the specification of inclusion and exclusion shapes. Only spheres so far: no prisms yet
  variable mousetrack_mode 0
  variable prev_label_list ""
  
  variable shape_type "sphere"
  variable shape_x "0.0"
  variable shape_y "0.0"
  variable shape_z "0.0"
  variable shape_r "10.0"
  variable shape_width "10.0"
  variable shape_height "10.0"
  variable shape_depth "10.0"
  variable last_action "" ;# inclusion, exclusion, contig
  variable last_shape "" ;# the last type of shape we just added
  variable last_x ""; variable last_y ""; variable last_z ""; variable last_r ""
  variable last_width ""; variable last_height ""; variable last_depth ""
  variable inclusion_list ""
  variable exclusion_list ""
  variable contig_list ""
  variable shape_mode ""
  variable highlight_num ""
  variable highlight_list ""
  variable graphics_mol ""
  variable graphics_resolution "20"
  variable graphics_highlight_color "yellow"
  variable graphics_inclusion_color "blue"
  #variable graphics_inclusion_material "Transparent"
  variable graphics_exclusion_color "red"
  #variable graphics_inclusion_material "Transparent"
  variable graphics_contig_color "ochre"
  variable graphics_contig_shape_size "5.0"
  variable graphics_material "Transparent"
  variable graphics_index_list ""
  variable cylinder_resolution 20
  
  variable big_incr "10.0"
  variable little_incr "1.0"
  variable python_arguments "-W ignore" ;# arguments to run python with. '-W ignore' suppresses warnings.
  
  variable suppress_GLSL_message 0
  variable automatically_set_mouse_mode 1
  
  variable plot_program "multiplot"
  variable plot_data ""
  variable plot_title "No data"
}
set ::povme2::auto_path $auto_path

# maybe also search the $PATH and $PYTHONPATH variables
proc ::povme2::UpdateMolecule {args} {
	# args is meaningless, but necessary
	set mainwin $::povme2::mainwin
	variable w
    	global vmd_molecule

    	# Update the molecule browser
    	set mollist [molinfo list]
    	$w.mollist.m configure -state disabled
        $w.mollist.m.menu delete 0 end
        set ::povme2::moltxt "(none)"
        #puts "Mainwin: $mainwin.mollist.m"
    	#if { [llength $mollist] > 0 } {
       # 	#$w.foot configure -state normal
       # 	$mainwin.mollist.m configure -state normal 
       # 	#puts $mollist
       # 	foreach id $mollist {
       #     		$mainwin.mollist.m insert end "$id - [molinfo $id get name]"            		
       # 	}
#	}
    if { [llength $mollist] > 0 } {
        #$w.foot configure -state normal
        $w.mollist.m configure -state normal 
        #puts $mollist
        foreach id $mollist {
            if { [molinfo $id get name] == "{POVME2 Graphics}" } { continue }
            $w.mollist.m.menu add radiobutton -value $id \
                -command {  global vmd_molecule ; variable w
                            # CHECK whether something exists in the shape list and warn the user
                            if { ($::povme2::which_mol != $::povme2::candidate_mol) && (([llength $::povme2::inclusion_list] != 0) || ([llength $::povme2::exclusion_list] != 0)) } {
                              set answer [tk_messageBox -title "Changing molecule" -message "Are you sure you want to change your molecule? All the shapes will be lost." -type yesno -icon warning] ;# warn the user
                              #puts "answer: $answer"
                              if { $answer == no} {  ;# then dont actually switch the molecule
                                #puts "clicked No"
                                #$w.mollist.m.menu activate $::povme2::which_mol ;# reset the listbox to select the right one
                                set ::povme2::candidate_mol $::povme2::which_mol
                                return 
                              } else {
                                #puts "clicked Yes"
                                set ::povme2::exclusion_list "" ;# delete shape list etc.
                                set ::povme2::inclusion_list ""
                                ::povme2::Update_shape_list inclusion
                                ::povme2::Update_shape_list exclusion
                                ::povme2::delete_highlight_graphics $::povme2::which_mol
                                ::povme2::draw_shape_graphics $::povme2::which_mol
                              }
                            } 
                            set ::povme2::which_mol $::povme2::candidate_mol ;# now actually set the selected mol to what they want
                            if {[info exists vmd_molecule($::povme2::which_mol)]} {
                              set ::povme2::moltxt "${::povme2::which_mol}: [molinfo $::povme2::which_mol get name]"
                            } else {
                              set ::povme2::moltxt "(none)" ; set ::povme2::which_mol -1
                            } 
                          } \
                -label "$id [molinfo $id get name]" \
                -variable ::povme2::candidate_mol
            if {$id == $::povme2::which_mol} {
              puts "id: $id"
              puts "which_mol: $::povme2::which_mol"
              puts "exists: [info exists vmd_molecule($::povme2::which_mol)]"
                if {[info exists vmd_molecule($::povme2::which_mol)]} {
                    set ::povme2::moltxt "$::povme2::which_mol:[molinfo $::povme2::which_mol get name]"  
                } else {
                    set ::povme2::moltxt "(none)"
                    set ::povme2::which_mol -1
                }
            }
        }
    }
}

proc draw_all_selections { w } {
  set ::povme2::highlight_list ""
  foreach shape_index [$w curselection] {
    set shape_string [$w get $shape_index]
    #puts "shape_string: $shape_string"
    lappend ::povme2::highlight_list $shape_string
  }
  ::povme2::draw_shape_graphics $::povme2::which_mol ;#$::povme2::graphics_mol
}

proc ::povme2::povme2_mainwin {} {
  variable w
  if { [winfo exists .povme2] } {
    #destroy $::povme2::w
    wm deiconify $w
    #::povme2::UpdateMolecule $w
    return
  }
  set w [toplevel ".povme2"]
  wm title $w "POVME2"
  wm resizable $w 0 0
  set ::povme2::mainwin $w
	
  wm protocol $w WM_DELETE_WINDOW {
    grab release $::povme2::mainwin
    after idle destroy $::povme2::mainwin
    set ::povme2::exclusion_list ""
    set ::povme2::inclusion_list ""
    ::povme2::delete_highlight_graphics $::povme2::which_mol
    ::povme2::draw_shape_graphics $::povme2::which_mol ;# delete graphics stuff
  }
	##
	## make menu bar
	##
  set ::povme2::w $w
  frame $w.menubar -relief raised -bd 2 ;# frame for menubar
  pack $w.menubar -padx 1 -fill x
	
  menubutton $w.menubar.help -text Help -underline 0 -menu $w.menubar.help.menu
  menubutton $w.menubar.settings -text Settings -underline 0 -menu $w.menubar.settings.menu
  menubutton $w.menubar.data -text Data -underline 0 -menu $w.menubar.data.menu
  menu $w.menubar.settings.menu -tearoff no
  $w.menubar.settings.menu add command -label "Files..." -command {::povme2::povme2_settings files}
  $w.menubar.settings.menu add command -label "Point Properties..." -command {::povme2::povme2_settings point_properties}
  $w.menubar.settings.menu add command -label "Contiguous Points/Convex-Hull Exclusion..." -command {::povme2::povme2_settings contiguous_points}
  #$w.menubar.settings.menu add command -label "Path Search..." -command {::povme2::povme2_settings pathsearching} ;# not currently needed
  $w.menubar.settings.menu add command -label "Output..." -command {::povme2::povme2_settings output}
  #$w.menubar.settings.menu add command -label "Graphics..." -command {::povme2::povme2_settings graphics}
  #$w.menubar.settings.menu add command -label "Advanced..." -command {::povme2::povme2_settings advanced}
  #puts "arglist_other: $::povme2::arglist_other"
  if {$::povme2::arglist_other != ""} { ;# if there are any subsequent settings, put them here
    $w.menubar.settings.menu add command -label "Other..." -command {::povme2::povme2_settings other}
  }
  
  menu $w.menubar.data.menu -tearoff no
  $w.menubar.data.menu add command -label "Load Volume Data" -command {
    set which_tk_getopen [tk_getOpenFile]
    if {$which_tk_getopen != ""} {
      set ::povme2::plot_data [::povme2::read_volume_data $which_tk_getopen]
      set ::povme2::plot_title $which_tk_getopen
      if { $::povme2::plot_data != "" } {
        ::povme2::do_plot [lindex $::povme2::plot_data 0] [lindex $::povme2::plot_data 1] $::povme2::plot_title
      } else {
        tk_messageBox -type ok -title "Data not loaded" -icon error -message "There was a problem with the volume data file loaded. Data could not be displayed."
      }
    }
    
  }
  $w.menubar.data.menu add command -label "Plot Volume Data" -command {
    if { $::povme2::plot_data != "" } {
      ::povme2::do_plot [lindex $::povme2::plot_data 0] [lindex $::povme2::plot_data 1] $::povme2::plot_title
    } else {
      tk_messageBox -type ok -title "No data loaded" -icon error -message "Please load POVME2 volume data or run a new POVME2 calculation."
    }
  }
  
  menu $w.menubar.help.menu -tearoff no
  $w.menubar.help.menu add command -label "About.." -command {tk_messageBox -type ok -title "About POVME2" -message "VMD GUI plugin for running pocket volume calculations using POVME2.\n\nDeveloped in the:\n\tAmaro Laboratory\n\tUniversity of California, San Diego\n\thttp://amarolab.ucsd.edu\n\nDevelopers:\n\tJacob Durrant\n\tLane Votapka\n\tJeff Wagner\n\nPlease see README in the installation directory for additional information" -icon info}
  $w.menubar.help.menu add command -label "Help..." -command "vmd_open_url $::povme2::help_url"
  # XXX - set menubutton width to avoid truncation in OS X
  $w.menubar.help config -width 5
  $w.menubar.settings config -width 5
  $w.menubar.data config -width 5
	
  pack $w.menubar.help -side right
  pack $w.menubar.settings -side left
  pack $w.menubar.data -side left
  # Mol List
  frame $w.mollist
  pack [label $w.mollist.lbl -text "Select Molecule:"] -side top -anchor w
  #listbox $w.mollist.m -relief raised -selectmode single ;# molecule listbox
  #$w.mollist.m
  #scrollbar $w.mollist.scrollbar -command [list $w.mollist.m yview] ;# molecule listbox scrollbar
  #$w.mollist.m configure -yscrollcommand [list $w.mollist.scrollbar set]
  #pack $w.mollist.m -side left -fill x -expand 1  ;# pack the listbox
  #pack $w.mollist.scrollbar -side left -fill y  ;# pack the scrollbar
  #
  menubutton $w.mollist.m -relief raised -bd 2 -direction flush \
  	-textvariable ::povme2::moltxt \
	-menu $w.mollist.m.menu
  menu $w.mollist.m.menu -tearoff no
  pack $w.mollist.m -side left -fill x -expand yes
  pack $w.mollist -side top -fill x -padx 10 -pady 0
  
  frame $w.sel
  pack [label $w.sel.lbl -text "Selection:"] -side top
  pack [entry $w.sel.selection  -width 40 -textvariable ::povme2::selection] -side top
  pack $w.sel -side top -anchor n
  
  frame $w.shapes
  # INCLUSION SHAPES
  frame $w.shapes.inclusion -relief groove -bd 2
  pack [label $w.shapes.inclusion.lbl -text "Inclusion Shapes:"] -side top -anchor w
  listbox $w.shapes.inclusion.m -relief raised -selectmode multiple -height 5 -width 40 ;# molecule listbox
  bind $w.shapes.inclusion.m <<ListboxSelect>> {draw_all_selections %W}
  
  scrollbar $w.shapes.inclusion.scrollbar -command [list $w.shapes.inclusion.m yview] ;# molecule listbox scrollbar
  $w.shapes.inclusion.m configure -yscrollcommand [list $w.shapes.inclusion.scrollbar set]
  pack $w.shapes.inclusion.m -side left -fill x -expand 1  ;# pack the listbox
  pack $w.shapes.inclusion.scrollbar -side left -fill y  ;# pack the scrollbar
  
  pack $w.shapes.inclusion -side top -anchor n -padx 10 -pady 5
  frame $w.shapes.inclusionbuttons
  button $w.shapes.inclusionbuttons.add -text "Add new shape..." -width 12 -command ::povme2::add_new_shape
  #button $w.shapes.inclusionbuttons.edit -text "Edit shape..." -width 10 -command ::povme2::edit_shape
  button $w.shapes.inclusionbuttons.delete -text "Delete shape" -width 10 -command ::povme2::delete_shape
  pack $w.shapes.inclusionbuttons.add -side left
  #pack $w.shapes.inclusionbuttons.edit -side left ;# Not working properly
  pack $w.shapes.inclusionbuttons.delete -side left
  pack $w.shapes.inclusionbuttons -side top -anchor n -pady 5
  # EXCLUSION SHAPES
  frame $w.shapes.exclusion -relief groove -bd 2
  pack [label $w.shapes.exclusion.lbl -text "Exclusion Shapes:"] -side top -anchor w
  listbox $w.shapes.exclusion.m -relief raised -selectmode multiple -height 5 -width 40 ;# molecule listbox
  bind $w.shapes.exclusion.m <<ListboxSelect>> {draw_all_selections %W}
  scrollbar $w.shapes.exclusion.scrollbar -command [list $w.shapes.exclusion.m yview] ;# molecule listbox scrollbar
  $w.shapes.exclusion.m configure -yscrollcommand [list $w.shapes.exclusion.scrollbar set]
  pack $w.shapes.exclusion.m -side left -fill x -expand 1  ;# pack the listbox
  pack $w.shapes.exclusion.scrollbar -side left -fill y  ;# pack the scrollbar
  
  pack $w.shapes.exclusion -side top -anchor n -padx 10 -pady 5
  frame $w.shapes.exclusionbuttons
  button $w.shapes.exclusionbuttons.add -text "Add new shape..." -width 12 -command {::povme2::add_new_shape exclusion}
  #button $w.shapes.exclusionbuttons.edit -text "Edit shape..." -width 10 -command {::povme2::edit_shape exclusion}
  button $w.shapes.exclusionbuttons.delete -text "Delete shape" -width 10 -command {::povme2::delete_shape exclusion}
  pack $w.shapes.exclusionbuttons.add -side left
  #pack $w.shapes.exclusionbuttons.edit -side left ;# not working properly
  pack $w.shapes.exclusionbuttons.delete -side left
  pack $w.shapes.exclusionbuttons -side top -anchor n
  # CONTIGUOUS SHAPES
  frame $w.shapes.contig -relief groove -bd 2
  pack [label $w.shapes.contig.lbl -text "Contiguous Pocket Seeds:"] -side top -anchor w
  listbox $w.shapes.contig.m -relief raised -selectmode multiple -height 5 -width 40 ;# molecule listbox
  bind $w.shapes.contig.m <<ListboxSelect>> {draw_all_selections %W}
  
  scrollbar $w.shapes.contig.scrollbar -command [list $w.shapes.contig.m yview] ;# molecule listbox scrollbar
  $w.shapes.contig.m configure -yscrollcommand [list $w.shapes.contig.scrollbar set]
  pack $w.shapes.contig.m -side left -fill x -expand 1  ;# pack the listbox
  pack $w.shapes.contig.scrollbar -side left -fill y  ;# pack the scrollbar
  
  pack $w.shapes.contig -side top -anchor n -padx 10 -pady 5
  frame $w.shapes.contigbuttons
  button $w.shapes.contigbuttons.add -text "Add new seed..." -width 12 -command {::povme2::add_new_shape contig}
  
  button $w.shapes.contigbuttons.delete -text "Delete seed" -width 10 -command {::povme2::delete_shape contig}
  pack $w.shapes.contigbuttons.add -side left
  pack $w.shapes.contigbuttons.delete -side left
  pack $w.shapes.contigbuttons -side top -anchor n -pady 5
  
  pack $w.shapes -side top -anchor n -padx 10 -pady 10
  
  frame $w.options
  pack [checkbutton $w.options.load_into_vmd -text "Load volumes into VMD when finished" -variable ::povme2::load_into_vmd ] -side top
  pack [checkbutton $w.options.plot_volumes -text "Plot volume results after run" -variable ::povme2::plot_volumes ] -side top
  #pack [label $w.src_sink.lbl0 -text "Source selection:"] -side top
  #pack [entry $w.src_sink.src_selection -textvariable ::povme2::src_selection] -side top
  #pack [label $w.src_sink.lbl1 -text "Sink selection:"] -side top
  #pack [entry $w.src_sink.sink_selection -textvariable ::povme2::sink_selection] -side top
  pack $w.options -side top -anchor n -padx 10 -pady 10
  # num_paths
  #frame $w.num_paths -relief groove -bd 2
  #pack [label $w.num_paths.lbl1 -text "Desired Number of Paths:"] -side top
  #pack [entry $w.num_paths.entry -textvariable ::povme2::desired_number_of_paths] -side top
  #pack $w.num_paths -side top -anchor n -padx 10 -pady 10
  # buttons
  frame $w.buttons
  button $w.buttons.run -text "Run POVME2" -width 10 -command ::povme2::run_povme2
  pack $w.buttons.run -side left
  pack $w.buttons -side top -anchor n -padx 10 -pady 10
  ::povme2::UpdateMolecule $w
  global vmd_molecule
  trace variable vmd_molecule w ::povme2::UpdateMolecule;
}

proc ::povme2::Update_shape_list {{mode "inclusion"}} {
  # args is meaningless, but necessary
  #set mainwin $::povme2::mainwin
  #variable w
  #variable moltxt
  #variable molid
  set w ".povme2"
    	
  $w.shapes.$mode.m delete 0 end
  $w.shapes.$mode.m configure -state disabled
  set ourname ::povme2::${mode}_list
  eval "set ourlist \"\$$ourname\""
  puts "ourlist: $ourlist"
  if { "[llength \$$ourlist]" > 0 } {
    $w.shapes.$mode.m configure -state normal 
    foreach id $ourlist {
      $w.shapes.$mode.m insert end "$id"     }
  }
}

proc ::povme2::draw_highlight_graphics { mol { shape "sphere" } { highlight_index "" } } {
  # draw the shape that is highlighted
  graphics $mol color $::povme2::graphics_highlight_color
  graphics $mol material $::povme2::graphics_material
  if { $shape == "sphere" } {
    if { $highlight_index != "" } {
      if { [llength $highlight_index] == 1 } {
        graphics $mol replace $highlight_index
      } else { ;# this its greater
        foreach graphic_id $highlight_index {
          graphics $mol delete $graphic_id
        }
      }
    }
    set ::povme2::highlight_num [ graphics $mol sphere "$::povme2::shape_x $::povme2::shape_y $::povme2::shape_z" radius $::povme2::shape_r resolution $::povme2::graphics_resolution ]
    
  } elseif { $shape == "prism" } {
    # PROBLEM: how to replace the graphics here... allowing highlight_index to be a list of graphics to replace
    
    set ::povme2::highlight_num [ ::povme2::draw_prism $mol $::povme2::shape_x $::povme2::shape_y $::povme2::shape_z $::povme2::shape_width $::povme2::shape_height $::povme2::shape_depth $highlight_index]
  } elseif { $shape == "cylinder" } {
    
    set ::povme2::highlight_num [ ::povme2::draw_cylinder $mol $::povme2::shape_x $::povme2::shape_y $::povme2::shape_z $::povme2::shape_width $::povme2::shape_height $::povme2::shape_depth $::povme2::shape_r $highlight_index]
  }
  
}

proc ::povme2::delete_highlight_graphics { mol } {
  if { [llength $::povme2::highlight_num] > 1 } {
    foreach index $::povme2::highlight_num {
      graphics $mol delete $index
    }
  } elseif { $::povme2::highlight_num > -1 } {
    graphics $mol delete $::povme2::highlight_num
  }
  set ::povme2::highlight_list ""
  
}

proc ::povme2::draw_shape_graphics { mol } {
  # delete all old shapes
  foreach shape $::povme2::graphics_index_list {
    foreach index $shape {
      graphics $mol delete $index
    }
  }
  set ::povme2::graphics_index_list "" ;# reset the graphics index list to be empty
  
  # inclusion shapes
  graphics $mol material $::povme2::graphics_material
  foreach mode {inclusion exclusion contig} { ;# maybe also contig
    set ourcolorname ::povme2::graphics_${mode}_color
    eval "set ourcolor \"\$$ourcolorname\""
    set ourlistname ::povme2::${mode}_list
    eval "set ourlist \"\$$ourlistname\""
    
    foreach entry $ourlist {
      if {[lsearch $::povme2::highlight_list $entry] == -1 } {
        graphics $mol color "$ourcolor"
      } else {
        graphics $mol color $::povme2::graphics_highlight_color
      }
      #parse the shape info
      set shape [lindex $entry 0]
      set x [lindex $entry 2]; set y [lindex $entry 3]; set z [lindex $entry 4]; 
      if { $shape == "sphere" } {
        set radius [lindex $entry 6]; 
        set graphics_index [graphics $mol sphere "$x $y $z" radius $radius resolution $::povme2::graphics_resolution]
      } elseif { $shape == "prism" } {
        set width [lindex $entry 6]; set height [lindex $entry 7]; set depth [lindex $entry 8]
        set graphics_index [::povme2::draw_prism $mol $x $y $z $width $height $depth]
      }  elseif { $shape == "cylinder" } {
        set width [lindex $entry 6]; set height [lindex $entry 7]; set depth [lindex $entry 8]; set radius [lindex $entry 10]
        set graphics_index [::povme2::draw_cylinder $mol $x $y $z $width $height $depth $radius]
      }
      lappend ::povme2::graphics_index_list $graphics_index
    }
  }    
  # exclusion shapes
  
}

proc ::povme2::enable_entries { our_entry args } {
  $our_entry.muchless configure -state normal
  $our_entry.less configure -state normal
  $our_entry.textbox configure -state normal
  $our_entry.more configure -state normal
  $our_entry.muchmore configure -state normal
}

proc ::povme2::disable_entries { our_entry args } {
  $our_entry.muchless configure -state disabled
  $our_entry.less configure -state disabled
  $our_entry.textbox configure -state disabled
  $our_entry.more configure -state disabled
  $our_entry.muchmore configure -state disabled
}

proc ::povme2::povme2_addshape {color {mode inclusion} {shape_str ""} } {
  variable w
  variable addshape_win
  set ::povme2::shape_type "sphere"
  set win_name $w.addshape
  if { [winfo exists $win_name] } {
    wm deiconify $addshape_win
    return
  }
  set addshape_win [toplevel "$w.addshape"]
  wm transient $addshape_win $w
  wm protocol $addshape_win WM_DELETE_WINDOW {
    grab release $::povme2::addshape_win
    after idle destroy $::povme2::addshape_win
    trace remove variable ::vmd_pick_atom_silent write ::povme2::track_mouse
    trace remove variable ::vmd_pick_atom write ::povme2::mouse_click
  }
  wm title $addshape_win "Add shape for $mode"
  wm resizable $addshape_win no no
  
  # let user know that they should activate GLSL
  if { !($::povme2::suppress_GLSL_message) && ([display get rendermode] != "GLSL") } { ;# then display a message box asking the user to switch to GLSL
    set answer [tk_messageBox -message "It is recommended that you enable GLSL display render mode to visualize POVME shapes. Do this now?" -type yesno -icon question]
    if {$answer == "yes"} {display rendermode "GLSL"}
  }
  if { $::povme2::automatically_set_mouse_mode } {
    mouse mode labelatom
  }
  
  # Add description telling option for using mouse if requisite variables exist in vmd
  
  #option for Sphere/Prism
  if {$::povme2::which_mol == -1} {
    #set ::povme2::shape_x "0.0"
    #set ::povme2::shape_y "0.0"
    #set ::povme2::shape_z "0.0"
    #set ::povme2::shape_r "10.0"
    grab release $::povme2::addshape_win
    after idle destroy $::povme2::addshape_win 
    tk_messageBox -type ok -title "No molecule selected" -icon error -message "Please select a molecule before adding shapes."
    
    return
  } elseif { $shape_str == "" } {
    set coords [measure center [atomselect $::povme2::which_mol "all"]]
      set min [lindex [measure minmax [atomselect $::povme2::which_mol "all"]] 0 ]
    if { $mode == "contig" } {
      set size $::povme2::graphics_contig_shape_size
      if { ($::povme2::last_x != "") && ($::povme2::last_y != "") && ($::povme2::last_z != "")} {
        set coords "$::povme2::last_x $::povme2::last_y $::povme2::last_z"
      }
    } else {
      set size [vecsub $min $coords]
    }
    set ::povme2::shape_x [format {%0.3f} [lindex $coords 0]] ;# find the center of the entire molecule
    set ::povme2::shape_y [format {%0.3f} [lindex $coords 1]]
    set ::povme2::shape_z [format {%0.3f} [lindex $coords 2]]
    #set ::povme2::shape_r [format {%0.3f} [veclength $size]]
    set ::povme2::shape_r [format {%0.3f} 10 ]
    set ::povme2::shape_width [format {%0.3f} [veclength $size]]
    set ::povme2::shape_height [format {%0.3f} [veclength $size]]
    set ::povme2::shape_depth [format {%0.3f} [veclength $size]]
  } else {
    set ::povme2::shape_x [lindex $shape_str 2]
    set ::povme2::shape_y [lindex $shape_str 3]
    set ::povme2::shape_z [lindex $shape_str 4]
    if {[llength $shape_str] == 7} { ;# then its a sphere
      set ::povme2::shape_r [lindex $shape_str 6]
    } elseif { [lindex $shape_str 0] == "prism" } { ;# then it must be a prism
      set ::povme2::shape_width [lindex $shape_str 6]
      set ::povme2::shape_height [lindex $shape_str 7]
      set ::povme2::shape_depth [lindex $shape_str 8]
    } elseif { [lindex $shape_str 0] == "cylinder" } { ;# then it must be a prism
      set ::povme2::shape_width [lindex $shape_str 6]
      set ::povme2::shape_height [lindex $shape_str 7]
      set ::povme2::shape_depth [lindex $shape_str 8]
      set ::povme2::shape_r [lindex $shape_str 10]
    }
  }
  #set highlight_index [::povme2::draw_highlight_graphics $::povme2::graphics_mol "sphere" ""]
  set highlight_index [::povme2::draw_highlight_graphics $::povme2::which_mol "sphere" ""]
  #::povme2::update_highlight_graphics
  
  set ::povme2::shape_mode $mode
  pack [label $addshape_win.lblshape -text "Shape:"] -side top -anchor w
  pack [radiobutton $addshape_win.option_sphere -text "Sphere" -value "sphere" -variable ::povme2::shape_type -command "::povme2::delete_highlight_graphics $::povme2::which_mol; ::povme2::draw_highlight_graphics $::povme2::which_mol sphere $::povme2::highlight_num; ::povme2::enable_entries $addshape_win.r; ::povme2::disable_entries $addshape_win.width; ::povme2::disable_entries $addshape_win.height; ::povme2::disable_entries $addshape_win.depth"] -side top -anchor w ;
  pack [radiobutton $addshape_win.option_prism -text "Prism" -value "prism" -variable ::povme2::shape_type -command " ::povme2::delete_highlight_graphics $::povme2::which_mol; ::povme2::draw_highlight_graphics $::povme2::which_mol prism $::povme2::highlight_num; ::povme2::disable_entries $addshape_win.r; ::povme2::enable_entries $addshape_win.width; ::povme2::enable_entries $addshape_win.height; ::povme2::enable_entries $addshape_win.depth"] -side top -anchor w
  pack [radiobutton $addshape_win.option_cylinder -text "Cylinder" -value "cylinder" -variable ::povme2::shape_type -command " ::povme2::delete_highlight_graphics $::povme2::which_mol; ::povme2::draw_highlight_graphics $::povme2::which_mol cylinder $::povme2::highlight_num; ::povme2::enable_entries $addshape_win.r; ::povme2::enable_entries $addshape_win.width; ::povme2::enable_entries $addshape_win.height; ::povme2::enable_entries $addshape_win.depth"] -side top -anchor w
  
  pack [label $addshape_win.lblx -text "Center X-Coordinate:"] -side top -anchor w
  frame $addshape_win.x
  pack [button $addshape_win.x.muchless -text "<<" -command {set ::povme2::shape_x [format {%0.3f} [expr "$::povme2::shape_x - $::povme2::big_incr"]]; ::povme2::draw_highlight_graphics $::povme2::which_mol $::povme2::shape_type $::povme2::highlight_num } ] -side left -anchor w
  pack [button $addshape_win.x.less -text "<" -command {set ::povme2::shape_x [format {%0.3f} [expr "$::povme2::shape_x - $::povme2::little_incr"]]; ::povme2::draw_highlight_graphics $::povme2::which_mol $::povme2::shape_type $::povme2::highlight_num } ] -side left -anchor w
  pack [entry $addshape_win.x.textbox -width 10 -textvariable ::povme2::shape_x] -side left -anchor w
  pack [button $addshape_win.x.more -text ">" -command {set ::povme2::shape_x [format {%0.3f} [expr "$::povme2::shape_x + $::povme2::little_incr"]]; ::povme2::draw_highlight_graphics $::povme2::which_mol $::povme2::shape_type $::povme2::highlight_num } ] -side left -anchor w
  pack [button $addshape_win.x.muchmore -text ">>" -command {set ::povme2::shape_x [format {%0.3f} [expr "$::povme2::shape_x + $::povme2::big_incr"]]; ::povme2::draw_highlight_graphics $::povme2::which_mol $::povme2::shape_type $::povme2::highlight_num } ] -side left -anchor w
  pack $addshape_win.x -side top -anchor n
  
  pack [label $addshape_win.lbly -text "Center Y-Coordinate:"] -side top -anchor w
  frame $addshape_win.y
  pack [button $addshape_win.y.muchless -text "<<" -command {set ::povme2::shape_y [format {%0.3f} [expr "$::povme2::shape_y - $::povme2::big_incr"]]; ::povme2::draw_highlight_graphics $::povme2::which_mol $::povme2::shape_type $::povme2::highlight_num } ] -side left -anchor w
  pack [button $addshape_win.y.less -text "<" -command {set ::povme2::shape_y [format {%0.3f} [expr "$::povme2::shape_y - $::povme2::little_incr"]]; ::povme2::draw_highlight_graphics $::povme2::which_mol $::povme2::shape_type $::povme2::highlight_num } ] -side left -anchor w
  pack [entry $addshape_win.y.textbox -width 10 -textvariable ::povme2::shape_y] -side left -anchor w
  pack [button $addshape_win.y.more -text ">" -command {set ::povme2::shape_y [format {%0.3f} [expr "$::povme2::shape_y + $::povme2::little_incr"]]; ::povme2::draw_highlight_graphics $::povme2::which_mol $::povme2::shape_type $::povme2::highlight_num } ] -side left -anchor w
  pack [button $addshape_win.y.muchmore -text ">>" -command {set ::povme2::shape_y [format {%0.3f} [expr "$::povme2::shape_y + $::povme2::big_incr"]]; ::povme2::draw_highlight_graphics $::povme2::which_mol $::povme2::shape_type $::povme2::highlight_num } ] -side left -anchor w
  pack $addshape_win.y -side top -anchor n
  
  pack [label $addshape_win.lblz -text "Center Z-Coordinate:"] -side top -anchor w
  frame $addshape_win.z
  pack [button $addshape_win.z.muchless -text "<<" -command {set ::povme2::shape_z [format {%0.3f} [expr "$::povme2::shape_z - $::povme2::big_incr"]]; ::povme2::draw_highlight_graphics $::povme2::which_mol $::povme2::shape_type $::povme2::highlight_num } ] -side left -anchor w
  pack [button $addshape_win.z.less -text "<" -command {set ::povme2::shape_z [format {%0.3f} [expr "$::povme2::shape_z - $::povme2::little_incr"]]; ::povme2::draw_highlight_graphics $::povme2::which_mol $::povme2::shape_type $::povme2::highlight_num } ] -side left -anchor w
  pack [entry $addshape_win.z.textbox -width 10 -textvariable ::povme2::shape_z] -side left -anchor w
  pack [button $addshape_win.z.more -text ">" -command {set ::povme2::shape_z [format {%0.3f} [expr "$::povme2::shape_z + $::povme2::little_incr"]]; ::povme2::draw_highlight_graphics $::povme2::which_mol $::povme2::shape_type $::povme2::highlight_num } ] -side left -anchor w
  pack [button $addshape_win.z.muchmore -text ">>" -command {set ::povme2::shape_z [format {%0.3f} [expr "$::povme2::shape_z + $::povme2::big_incr"]]; ::povme2::draw_highlight_graphics $::povme2::which_mol $::povme2::shape_type $::povme2::highlight_num } ] -side left -anchor w
  pack $addshape_win.z -side top -anchor n
  
  # add frame for prism/cylinder options
   frame $addshape_win.width
  pack [label $addshape_win.width.lblr -text "Axis X-Coordinate:"] -side top -anchor w
  pack [button $addshape_win.width.muchless -state disabled -text "<<" -command { if {$::povme2::shape_width >= $::povme2::big_incr} {set ::povme2::shape_width [format {%0.3f} [expr "$::povme2::shape_width - $::povme2::big_incr"]]; ::povme2::draw_highlight_graphics $::povme2::which_mol $::povme2::shape_type $::povme2::highlight_num }} ] -side left -anchor w
  pack [button $addshape_win.width.less -state disabled -text "<" -command { if {$::povme2::shape_width >= $::povme2::little_incr} {set ::povme2::shape_width [format {%0.3f} [expr "$::povme2::shape_width - $::povme2::little_incr"]]; ::povme2::draw_highlight_graphics $::povme2::which_mol $::povme2::shape_type $::povme2::highlight_num }} ] -side left -anchor w
  pack [entry $addshape_win.width.textbox -state disabled -width 10 -textvariable ::povme2::shape_width] -side left -anchor w
  pack [button $addshape_win.width.more -state disabled -text ">" -command {set ::povme2::shape_width [format {%0.3f} [expr "$::povme2::shape_width + $::povme2::little_incr"]]; ::povme2::draw_highlight_graphics $::povme2::which_mol $::povme2::shape_type $::povme2::highlight_num } ] -side left -anchor w
  pack [button $addshape_win.width.muchmore -state disabled -text ">>" -command {set ::povme2::shape_width [format {%0.3f} [expr "$::povme2::shape_width + $::povme2::big_incr"]]; ::povme2::draw_highlight_graphics $::povme2::which_mol $::povme2::shape_type $::povme2::highlight_num } ] -side left -anchor w
  pack $addshape_win.width -side top -anchor n
  
   frame $addshape_win.height
  pack [label $addshape_win.height.lblr -text "Axis Y-Coordinate:"] -side top -anchor w
  pack [button $addshape_win.height.muchless -state disabled -text "<<" -command { if {$::povme2::shape_height >= $::povme2::big_incr} {set ::povme2::shape_height [format {%0.3f} [expr "$::povme2::shape_height - $::povme2::big_incr"]]; ::povme2::draw_highlight_graphics $::povme2::which_mol $::povme2::shape_type $::povme2::highlight_num }} ] -side left -anchor w
  pack [button $addshape_win.height.less -state disabled -text "<" -command { if {$::povme2::shape_height >= $::povme2::little_incr} {set ::povme2::shape_height [format {%0.3f} [expr "$::povme2::shape_height - $::povme2::little_incr"]]; ::povme2::draw_highlight_graphics $::povme2::which_mol $::povme2::shape_type $::povme2::highlight_num }} ] -side left -anchor w
  pack [entry $addshape_win.height.textbox -state disabled -width 10 -textvariable ::povme2::shape_height] -side left -anchor w
  pack [button $addshape_win.height.more -state disabled -text ">" -command {set ::povme2::shape_height [format {%0.3f} [expr "$::povme2::shape_height + $::povme2::little_incr"]]; ::povme2::draw_highlight_graphics $::povme2::which_mol $::povme2::shape_type $::povme2::highlight_num } ] -side left -anchor w
  pack [button $addshape_win.height.muchmore -state disabled -text ">>" -command {set ::povme2::shape_height [format {%0.3f} [expr "$::povme2::shape_height + $::povme2::big_incr"]]; ::povme2::draw_highlight_graphics $::povme2::which_mol $::povme2::shape_type $::povme2::highlight_num } ] -side left -anchor w
  pack $addshape_win.height -side top -anchor n
  
   frame $addshape_win.depth
  pack [label $addshape_win.depth.lblr -text "Axis Z-Coordinate:"] -side top -anchor w
  pack [button $addshape_win.depth.muchless -state disabled -text "<<" -command { if {$::povme2::shape_depth >= $::povme2::big_incr} {set ::povme2::shape_depth [format {%0.3f} [expr "$::povme2::shape_depth - $::povme2::big_incr"]]; ::povme2::draw_highlight_graphics $::povme2::which_mol $::povme2::shape_type $::povme2::highlight_num }} ] -side left -anchor w
  pack [button $addshape_win.depth.less -state disabled -text "<" -command { if {$::povme2::shape_depth >= $::povme2::little_incr} {set ::povme2::shape_depth [format {%0.3f} [expr "$::povme2::shape_depth - $::povme2::little_incr"]]; ::povme2::draw_highlight_graphics $::povme2::which_mol $::povme2::shape_type $::povme2::highlight_num }} ] -side left -anchor w
  pack [entry $addshape_win.depth.textbox -state disabled -width 10 -textvariable ::povme2::shape_depth] -side left -anchor w
  pack [button $addshape_win.depth.more -state disabled -text ">" -command {set ::povme2::shape_depth [format {%0.3f} [expr "$::povme2::shape_depth + $::povme2::little_incr"]]; ::povme2::draw_highlight_graphics $::povme2::which_mol $::povme2::shape_type $::povme2::highlight_num } ] -side left -anchor w
  pack [button $addshape_win.depth.muchmore -state disabled -text ">>" -command {set ::povme2::shape_depth [format {%0.3f} [expr "$::povme2::shape_depth + $::povme2::big_incr"]]; ::povme2::draw_highlight_graphics $::povme2::which_mol $::povme2::shape_type $::povme2::highlight_num } ] -side left -anchor w
  pack $addshape_win.depth -side top -anchor n
  
  # add frame here for spheres; that hides when prisms are selected
  frame $addshape_win.r
  pack [label $addshape_win.r.lblr -text "Radius:"] -side top -anchor w
  pack [button $addshape_win.r.muchless -text "<<" -command { if {$::povme2::shape_r >= $::povme2::big_incr} {set ::povme2::shape_r [format {%0.3f} [expr "$::povme2::shape_r - $::povme2::big_incr"]]; ::povme2::draw_highlight_graphics $::povme2::which_mol $::povme2::shape_type $::povme2::highlight_num }} ] -side left -anchor w
  pack [button $addshape_win.r.less -text "<" -command { if {$::povme2::shape_r >= $::povme2::little_incr} {set ::povme2::shape_r [format {%0.3f} [expr "$::povme2::shape_r - $::povme2::little_incr"]]; ::povme2::draw_highlight_graphics $::povme2::which_mol $::povme2::shape_type $::povme2::highlight_num }} ] -side left -anchor w
  pack [entry $addshape_win.r.textbox -width 10 -textvariable ::povme2::shape_r] -side left -anchor w
  pack [button $addshape_win.r.more -text ">" -command {set ::povme2::shape_r [format {%0.3f} [expr "$::povme2::shape_r + $::povme2::little_incr"]]; ::povme2::draw_highlight_graphics $::povme2::which_mol $::povme2::shape_type $::povme2::highlight_num } ] -side left -anchor w
  pack [button $addshape_win.r.muchmore -text ">>" -command {set ::povme2::shape_r [format {%0.3f} [expr "$::povme2::shape_r + $::povme2::big_incr"]]; ::povme2::draw_highlight_graphics $::povme2::which_mol $::povme2::shape_type $::povme2::highlight_num } ] -side left -anchor w
  pack $addshape_win.r -side top -anchor n
  
  
  
  
  frame $addshape_win.okaycancel
  button $addshape_win.okaycancel.okay -text OK -width 6 \
	  -command {
	  if {$::povme2::shape_type == "sphere"} {
	    lappend ::povme2::${::povme2::shape_mode}_list "$::povme2::shape_type coords: [format {%0.3f} $::povme2::shape_x] [format {%0.3f} $::povme2::shape_y] [format {%0.3f} $::povme2::shape_z] radius: [format {%0.3f} ${::povme2::shape_r}]"
	  } elseif {$::povme2::shape_type == "prism"} {
	    lappend ::povme2::${::povme2::shape_mode}_list "$::povme2::shape_type coords: [format {%0.3f} $::povme2::shape_x] [format {%0.3f} $::povme2::shape_y] [format {%0.3f} $::povme2::shape_z] corner: [format {%0.3f} ${::povme2::shape_width}] [format {%0.3f} ${::povme2::shape_height}] [format {%0.3f} ${::povme2::shape_depth}]"
	  } elseif {$::povme2::shape_type == "cylinder"} {
	    lappend ::povme2::${::povme2::shape_mode}_list "$::povme2::shape_type coords: [format {%0.3f} $::povme2::shape_x] [format {%0.3f} $::povme2::shape_y] [format {%0.3f} $::povme2::shape_z] corner: [format {%0.3f} ${::povme2::shape_width}] [format {%0.3f} ${::povme2::shape_height}] [format {%0.3f} ${::povme2::shape_depth}] radius: [format {%0.3f} ${::povme2::shape_r}]"
	  }
	  ::povme2::Update_shape_list ${::povme2::shape_mode}
	  if { $::povme2::shape_mode == "inclusion" } {
	    set ::povme2::last_action ${::povme2::shape_mode}
	    set ::povme2::last_shape ${::povme2::shape_type}
	    set ::povme2::last_x [format {%0.3f} $::povme2::shape_x]; set ::povme2::last_y [format {%0.3f} $::povme2::shape_y]; set ::povme2::last_z [format {%0.3f} $::povme2::shape_z] 
	    set ::povme2::last_width [format {%0.3f} ${::povme2::shape_width}]; set ::povme2::last_height [format {%0.3f} ${::povme2::shape_height}]; set ::povme2::last_depth [format {%0.3f} ${::povme2::shape_depth}]; set ::povme2::last_r [format {%0.3f} ${::povme2::shape_r}]
	  }
	  grab release $::povme2::addshape_win
	  after idle destroy $::povme2::addshape_win 
	  #::povme2::delete_highlight_graphics $::povme2::graphics_mol 
	  ::povme2::delete_highlight_graphics $::povme2::which_mol
	  #::povme2::draw_shape_graphics $::povme2::graphics_mol 
	  ::povme2::draw_shape_graphics $::povme2::which_mol 
	  trace remove variable ::vmd_pick_atom_silent write ::povme2::track_mouse
          trace remove variable ::vmd_pick_atom write ::povme2::mouse_click
          if { $::povme2::automatically_set_mouse_mode } {
            mouse mode rotate
          }
	}
  button $addshape_win.okaycancel.cancel -text Cancel -width 6 \
	  -command {
	  	grab release $::povme2::addshape_win
	  	after idle destroy $::povme2::addshape_win
	  	#::povme2::delete_highlight_graphics $::povme2::graphics_mol
	  	::povme2::delete_highlight_graphics $::povme2::which_mol
	  	#::povme2::draw_shape_graphics $::povme2::graphics_mol 
	        trace remove variable ::vmd_pick_atom_silent write ::povme2::track_mouse
                trace remove variable ::vmd_pick_atom write ::povme2::mouse_click
                if { $::povme2::automatically_set_mouse_mode } {
                  mouse mode rotate
                }
	}
  pack $addshape_win.okaycancel.okay $addshape_win.okaycancel.cancel -side left
  pack $addshape_win.okaycancel -side top -anchor n
  #return "sphere"
  
  # track the mouse
  set ::povme2::prev_label_list []
  mouse callback on ;# make it so that we can track which atom the mouse is near
  trace variable ::vmd_pick_atom_silent w ::povme2::track_mouse
  trace variable ::vmd_pick_atom w ::povme2::mouse_click
}

proc ::povme2::track_mouse { varname args } {
  set mol $::vmd_pick_mol_silent
  set atom $::vmd_pick_atom_silent
  #puts "POVME detected the mouse moving: mol: $mol, atom: $atom"
  set center [measure center [atomselect $mol "index $atom"]]
  if {$::povme2::mousetrack_mode == 0} { ;# shift is not held down
    set ::povme2::shape_x [format {%0.3f} [lindex $center 0]]
    set ::povme2::shape_y [format {%0.3f} [lindex $center 1]]
    set ::povme2::shape_z [format {%0.3f} [lindex $center 2]]
    ::povme2::draw_highlight_graphics $::povme2::which_mol $::povme2::shape_type $::povme2::highlight_num
  } elseif {$::povme2::mousetrack_mode == 2} { ;# shift is held down
    set diffvec [vecsub $center "$::povme2::shape_x $::povme2::shape_y $::povme2::shape_z"]
    set ::povme2::shape_r [ format {%0.3f} [veclength $diffvec]]
    set ::povme2::shape_width [ format {%0.3f} [lindex $diffvec 0]]
    set ::povme2::shape_height [ format {%0.3f} [lindex $diffvec 1]]
    set ::povme2::shape_depth [ format {%0.3f} [lindex $diffvec 2]]
    ::povme2::draw_highlight_graphics $::povme2::which_mol $::povme2::shape_type $::povme2::highlight_num
  }
}

proc ::povme2::mouse_click { varname args } {
  set mol $::vmd_pick_mol
  set atom $::vmd_pick_atom
  set shift $::vmd_pick_shift_state
  
  
  set lastindex [expr "[llength [label list Atoms]]-1"]
  label delete Atoms $lastindex
  #puts "POVME detected the mouse clicking: mol: $mol, atom: $atom, shift: $shift, mousemode: $::povme2::mousetrack_mode"
  if {$::povme2::mousetrack_mode == 0} {
    if {$shift == 0} { ;# then the shift is not pressed
      set ::povme2::mousetrack_mode 1
    } else {
      set ::povme2::mousetrack_mode 2
    }
  } elseif {$::povme2::mousetrack_mode == 1} {
    if {$shift == 0} { ;# then the shift is not pressed
      set ::povme2::mousetrack_mode 0
    } else {
      set ::povme2::mousetrack_mode 2
    }
  } elseif {$::povme2::mousetrack_mode == 2} {
    set ::povme2::mousetrack_mode 1
  }
  
}

proc ::povme2::add_new_shape {{mode "inclusion"} {shape_str ""} } {
  set color "blue"
  puts "mode: $mode"
  set new_shape [::povme2::povme2_addshape $color $mode $shape_str]
  
  # add shape to $w.shapes.inclusion.m
}

proc ::povme2::edit_shape {{mode "inclusion"}} {
  # retrieve the selected shape's information
  if { [llength $::povme2::highlight_list] == 1 } {
    set shape_str [lindex $::povme2::highlight_list 0]
    ::povme2::delete_shape $mode
    ::povme2::delete_highlight_graphics $::povme2::which_mol ;#$::povme2::graphics_mol
    ::povme2::add_new_shape $mode $shape_str
  } elseif { [llength $::povme2::highlight_list] == 0 } {
    # error dialog
    tk_messageBox -icon error -title "No shape selected" -message "An $mode shape must be selected in order to edit"
  } else {
    # error dialog
    tk_messageBox -icon error -title "Too many shapes" -message "Only one shape may be selected when editing"
  }
}

proc ::povme2::delete_shape {{mode "inclusion"}} {
  #if {$mode == "inclusion"} {
  #::povme2::delete_highlight_graphics $::povme2::graphics_mol
  ::povme2::delete_highlight_graphics $::povme2::which_mol
  set ourname ::povme2::mainwin.shapes.$mode.m
  eval "set ourname2 \"\$$ourname\""
    set selected_indeces [$ourname2 curselection]
    foreach index [lreverse $selected_indeces] { ;# we have to reverse it, so we're not deleting members of a list before we can reference them
      set ourname ::povme2::${mode}_list
      eval "set ourlist \"\$$ourname\""
      set ::povme2::${mode}_list [lreplace $ourlist $index $index]
      ::povme2::Update_shape_list $mode
    }
  #}
  #::povme2::draw_shape_graphics $::povme2::graphics_mol
  ::povme2::draw_shape_graphics $::povme2::which_mol
}


proc ::povme2::run_povme2 {} {
  #set selected_mol_index [$::povme2::mainwin.mollist.m curselection]
  if {$::povme2::which_mol == -1} {
    tk_messageBox -type ok -title "Error" -icon error -message "Error: Please select a molecule."; return
  }
  if {[file exists $::povme2::povme2_directory] == 0} {	;# povme2 is not found
	#give an error message
	tk_messageBox -type ok -title "POVME2 not found" -icon error -message "Error: The POVME2 program location does not exist at the specified location: $::povme2::povme2_directory. Go to Settings > Files to change location"
	return
  }
  if {[file readable $::povme2::povme2_directory] == 0} {	;# povme2 is not readable
	#give an error message
	tk_messageBox -type ok -title "POVME2 not readable" -icon error -message "Error: The specified POVME2 program does not have read permissions."
	return
  }
  set selected_mol $::povme2::which_mol ;#[lindex [molinfo list] $::povme2::which_mol]
  
  
  puts "Now writing POVME2 input file..."
  

  set dialog [toplevel .wait]
  pack [frame $dialog.frame -borderwidth 2 -relief raised] -fill both
  pack [message $dialog.frame.msg -text "Please wait\nPOVME2 is running..." -width 100] -padx 100 -pady 100
  update
  grab $dialog
  
  set input_filename [::povme2::make_input_file "${::povme2::workdir}${::povme2::sep}povme_input.ini"] ;# PROBLEM HERE
  
  #animate write pdb $::povme2::temp_pdb sel $pdb_struct $selected_mol  ;# this writes the pdb trajectory explicitly
  
  #set pdb_trajectory_filename $::povme2::temp_pdb ;#[molinfo $selected_mol get filename]
  set starttime [clock seconds] ;# get the time of the start of the job
  puts "starttime: $starttime" ;# TEMP
  set timename [clock format [clock seconds] -format "POVME2_output__%h_%d_%Y__%H_%M_%p"]
  #if {$::povme2::output_directory == ""} {
  #  set outdir $timename
  #} else {
  #  set outdir [file join $::povme2::output_directory $timename]
  #}  
  set curdir [pwd]
  cd $::povme2::workdir
  set command "$::povme2::python_executable $::povme2::python_arguments $::povme2::povme2_directory $input_filename"
  puts "running command: $command"
  update
  
  set error_code [catch {exec {*}$command} results options]
  #set results "TESTING"; set error "0"
  destroy $dialog
  puts "POVME2 std. output: $results"
  if {$error_code} {tk_messageBox -type ok -icon warning -title "Failure" -message "Alert: POVME2 failed to run properly. See VMD standard output for more detailed information. Error code: $error_code"}
  #set error [catch {exec rm $::povme2::temp_pdb}] ;# delete the temporary trajectory
  if {$::povme2::load_into_vmd == 1} {
    #cd $outdir ;# cd to the directory where everything is
  #  source "visualize.tcl"
  set ::povme2::iso_value 0.050000
  # sets iso value for volumetric and color density map visualizations

    if { $::povme2::volumetric_density_file == 1 } { ;# load volumetric density
      set vol_dx_filename "${::povme2::output_filename_prefix}volumetric_density.dx"
      if {([file exists $vol_dx_filename] == 1) } { ;# make sure it was modified since running povme2
        if { ([file mtime $vol_dx_filename] >= $starttime)} {
          mol representation Isosurface $::povme2::iso_value 0 0 1 1 1
          mol new $vol_dx_filename type "dx"
        } else {
          puts "Alert: POVME2 failed to load volumetric density file. The file was timestamped later than expected; indicating that POVME2 may have failed with errors: $vol_dx_filename."
        }
      } else {
        puts "Alert: POVME2 failed to load volumetric density file. The file was not located where expected: $vol_dx_filename."
        #tk_messageBox -type ok -icon warning -title "Failure to load volumetric density" -message "Alert: POVME2 failed to load volumetric density file. The file was not located where expected."
      }
    }
    if { $::povme2::colored_density_files == 1 } { ;# load color density files
      set ::povme2::feature_dict [dict create hbondAcceptor 1 hbondDonor 0 aromatic 3 hydrophobic 13 hydrophilic 7 adjacency 8 occupancy 8] 
      set ::povme2::file_append "_volumetric_density.dx"
      set ::povme2::dict_keys [dict keys $::povme2::feature_dict]
      foreach term $::povme2::dict_keys {
        set col_dx_filename "${::povme2::output_filename_prefix}$term$::povme2::file_append"

        if {[file exists $col_dx_filename] == 1 } {
          set color_num [dict get $::povme2::feature_dict $term]
          mol color ColorID $color_num
          mol representation Isosurface $::povme2::iso_value 0 0 1 1 1 #set default isoval, wireframe and no box
          mol new $col_dx_filename type "dx"
          if {$term in [list adjacency occupancy]} {
            set index [molinfo top get id]
            mol off $index
          }
        } else {
          puts "Alert: POVME2 failed to load $term colored density files. The files were not located where expected."
          #tk_messageBox -type ok -icon warning -title "Failure to load volumetric density" -message "Alert: POVME2 failed to load volumetric density file. The file was not located where expected."
        }
      }
    }
    
    if { $::povme2::volume_trajectory == 1 } { ;# load pocket structure file
      set vol_traj_filename "${output_filename_prefix}volume_trajectory.pdb"
      if {([file exists $vol_traj_filename] == 1) } { ;# see if the file exists and make sure it was modified since running povme2
        if { ([file mtime $vol_traj_filename] >= $starttime) } {
          mol new $vol_traj_filename type "PDB"
        } else {
          puts "Alert: POVME2 failed to load volume trajectory file. The file was timestamped later than expected; indicating that POVME2 may have failed with errors: $vol_traj_filename."
        }
      } else {
        puts "Alert: POVME2 failed to load volume trajectory file. The file was not located where expected: $vol_traj_filename"
        #tk_messageBox -type ok -icon warning -title "Failure to load volumetric density" -message "Alert: POVME2 failed to load volumetric density file. The file was not located where expected."
      }
    }
  } ;# load the povme2 tcl file for visualization
  if { $::povme2::plot_volumes == 1} {
    set vol_tabbed_filename "${::povme2::output_filename_prefix}volumes.tabbed.txt"
    if {([file exists $vol_dx_filename] == 1)  } { ;# see if the file exists 
      if {([file mtime $vol_dx_filename] >= $starttime) } { ;# and make sure it was modified since running povme2
        set vol_data [::povme2::read_volume_data $vol_tabbed_filename]
        #if { [llength [lindex $vol_data 0]] > 1 } { ;# if there's more than one datapoint
          set ::povme2::plot_data $vol_data
          set ::povme2::plot_title $vol_tabbed_filename ;# save this info for later in case the user wants to load it
          ::povme2::do_plot [lindex $vol_data 0] [lindex $vol_data 1] $vol_tabbed_filename
        #}
      } else {
        puts "Alert: POVME2 failed to load tabbed volume file. The file was timestamped earlier than expected; indicating that POVME2 may have failed with errors. Unable to plot volume data. Location: $vol_tabbed_filename"
      }
    } else {
      puts "Alert: POVME2 failed to load tabbed volume file. The file was not located where expected. Unable to plot volume data. Location: $vol_tabbed_filename"
    }
  }
  cd $curdir
  
  ::povme2::UpdateMolecule
  return
}

proc ::povme2::read_volume_data {volume_filename} {
  set vol_data ""
  set fp [open $volume_filename r] ;# open the file
  set file_data [read $fp]
  set lines [split $file_data "\n"]
  set counter 1
  set indeces ""
  set volumes ""
  foreach line $lines { ;# process each of the lines
    set index [lindex $line 0]
    set volume [lindex $line 1]
    if {$index != ""} {lappend indeces $index}
    if {$volume != ""} {lappend volumes $volume}
    incr counter
  }
  lappend vol_data $indeces
  lappend vol_data $volumes
  return $vol_data
}

proc ::povme2::parse_helpfile {{help_command "-help"}} {
  set full_help_command "python $::povme2::povme2_directory $help_command"
  set helpstring [exec {*}$full_help_command]
  set helplist [split $helpstring "\n"]
  set ::povme2::desclist {}
  set append_desc False ;# whether we are appending the description contents to a particular argument
  set ::povme2::arglist {}
  set descvar ""
  set arg ""
  foreach line $helplist {
    set firstchar [lindex [split $line ""] 0] ;# get the first character in a line
    if {($firstchar == " ") || ($firstchar == "\t")} { ;# then the first character is a whitespace
      if {$append_desc == True} { ;# then we are appending this to the description list
        set descvar "${descvar}[string range $line 3 end]"
      }
      continue
    } elseif {([scan $firstchar %c] >= 97) && ([scan $firstchar %c] <= 122)} { ;# if the line begins with a lowercase letter
      if {$arg != ""} {set ::povme2::descarray($arg) $descvar}
      set arg [lindex [split $line ":"] 0]
      lappend ::povme2::arglist $arg
      #puts "descvar: $descvar, arg: $arg"
      set append_desc True
      set descvar [lindex [split $line ":"] 1]
    } else {
      set append_desc False
    }
  }
  #puts $::povme2::arglist
  #puts $::povme2::desclist
  foreach arg $::povme2::arglist {
    if {!([info exists "::povme2::$arg"])} { ;# if the variable doesn't already exist, create it
      set ::povme2::$arg ""
      lappend ::povme2::arglist_other $arg
    } 
  }
  
}

proc ::povme2::find_py {} {
  # attempts to find the povme2.py python script
  if {[info exists ::povme2::povme2_directory] && [file exists $::povme2::povme2_directory]} {return True} ;# search the defined location
  foreach path $::povme2::auto_path {
    if {[file exists [file join $path POVME POVME2.py]]} { ;# then we've found it
      set ::povme2::povme2_directory [file join $path POVME POVME2.py]
      return True
    } elseif {[file exists [file join $path POVME2.py]]} { ;# then we've found it
      set ::povme2::povme2_directory [file join $path POVME2.py]
      return True
    }
  }
  set ::povme2::povme2_directory "./POVME2.py" ;# then we couldn't find it
  return False
  
}

proc ::povme2::default_args {} { ;# in case there is some sort of problem with parsing the helpfile, this will assign the default arguments
  set ::povme2::arglist {}
  #lappend ::povme2::arglist output_directory
  lappend ::povme2::arglist povme2_directory
  lappend ::povme2::arglist grid_point_spacing
  lappend ::povme2::arglist distance_cutoff
  lappend ::povme2::arglist make_point_field_pdb
  lappend ::povme2::arglist use_contiguous_points
  lappend ::povme2::arglist exclude_points_outside_convex_hull
  lappend ::povme2::arglist contiguous_point_criteria
  lappend ::povme2::arglist output_filename_prefix
  lappend ::povme2::arglist compress_output
  lappend ::povme2::arglist separate_volume_pdbs
  lappend ::povme2::arglist volume_trajectory
  lappend ::povme2::arglist equal_num_points_per_frame
  lappend ::povme2::arglist tabbed_volume_file
  lappend ::povme2::arglist volumetric_density_file
  lappend ::povme2::argList colored_density_files
  lappend ::povme2::arglist num_processors
  lappend ::povme2::arglist disk_instead_of_memory
  lappend ::povme2::arglist python_executable

}

proc ::povme2::num_to_bool {num} { ;# returns lowercase boolean based on number
  if {$num == "0"} { return "false" 
  } elseif {$num == "1"} { return "true"
  } else { error "cannot convert value $num into a boolean" }
}

proc ::povme2::make_input_file {{filename "povme_input.ini"}} {
  #set olddir [pwd] ;# find where we are
  #cd $::povme2::workdir ;# go to the working directory
  set infile [open $filename w] ;# open a new file for writing
  
  set pdbfile [molinfo $::povme2::which_mol get filename]
  
  set pdb_selection [atomselect $::povme2::which_mol $::povme2::selection] 
  
  if {([file exists $pdbfile]) && ([molinfo $::povme2::which_mol get filetype] == "pdb") && ([$pdb_selection num] == [molinfo $::povme2::which_mol get numatoms])} { ;# then the file exists in PDB format and the selection has the same number of atoms
    puts $infile "PDBFileName\t$pdbfile"
  } else { ;# we will need to write the file itself in PDB format
    puts "Now writing pdb trajectory..."
    animate write pdb $::povme2::temp_pdb sel [atomselect $::povme2::which_mol $::povme2::selection] ;# this writes the pdb trajectory explicitly
    set pdb_trajectory_filename $::povme2::temp_pdb ;#[molinfo $selected_mol get filename]
    puts $infile "PDBFileName\t$::povme2::temp_pdb"
  }
  puts $infile "GridSpacing\t$::povme2::grid_point_spacing"
  foreach shape_str $::povme2::inclusion_list {
    set shape [lindex $shape_str 0]
    set shape_x [lindex $shape_str 2]
    set shape_y [lindex $shape_str 3]
    set shape_z [lindex $shape_str 4]
    if { $shape == "sphere" } {
      set shape_r [lindex $shape_str 6]
      puts $infile "InclusionSphere\t$shape_x $shape_y $shape_z $shape_r"
    } elseif { $shape == "prism" } {
      set shape_width [lindex $shape_str 6]
      set shape_height [lindex $shape_str 7]
      set shape_depth [lindex $shape_str 8]
      if {$shape_width < 0} {set shape_width [expr "-1.0 * $shape_width"]} ;# can't have negative widths, etc.
      if {$shape_height < 0} {set shape_height [expr "-1.0 * $shape_height"]} ;# can't have negative widths, etc.
      if {$shape_depth < 0} {set shape_depth [expr "-1.0 * $shape_depth"]} ;# can't have negative widths, etc.
      #puts $infile "PointsInclusionBox\t[expr "$shape_x - $shape_width"] [expr "$shape_y - $shape_height"] [expr "$shape_z - $shape_depth"] [expr 2.0 * $shape_width] [expr 2.0 * $shape_height] [expr 2.0 * $shape_depth]"
      puts $infile "InclusionBox\t$shape_x $shape_y $shape_z [expr 2.0 * $shape_width] [expr 2.0 * $shape_height] [expr 2.0 * $shape_depth]"
    } elseif { $shape == "cylinder" } {
      set shape_width [lindex $shape_str 6]
      set shape_height [lindex $shape_str 7]
      set shape_depth [lindex $shape_str 8]
      set shape_r [lindex $shape_str 10]
      set point1x [expr "$shape_x - $shape_width"]; set point1y [expr "$shape_y - $shape_height"]; set point1z [expr "$shape_z - $shape_depth"]
      set point2x [expr "$shape_x + $shape_width"]; set point2y [expr "$shape_y + $shape_height"]; set point2z [expr "$shape_z + $shape_depth"]
      puts $infile "InclusionCylinder\t$point1x $point1y $point1z $point2x $point2y $point2z $shape_r"
    }
  }
  foreach shape_str $::povme2::exclusion_list {
    set shape [lindex $shape_str 0]
    set shape_x [lindex $shape_str 2]
    set shape_y [lindex $shape_str 3]
    set shape_z [lindex $shape_str 4]
    if { $shape == "sphere" } {
      set shape_r [lindex $shape_str 6]
      puts $infile "ExclusionSphere\t$shape_x $shape_y $shape_z $shape_r"
    } elseif { $shape == "prism" } {
      set shape_width [lindex $shape_str 6]
      set shape_height [lindex $shape_str 7]
      set shape_depth [lindex $shape_str 8]
      if {$shape_width < 0} {set shape_width [expr "-1.0 * $shape_width"]} ;# can't have negative widths, etc.
      if {$shape_height < 0} {set shape_height [expr "-1.0 * $shape_height"]} ;# can't have negative widths, etc.
      if {$shape_depth < 0} {set shape_depth [expr "-1.0 * $shape_depth"]} ;# can't have negative widths, etc.
      #puts $infile "PointsExclusionBox\t[expr "$shape_x - $shape_width"] [expr "$shape_y - $shape_height"] [expr "$shape_z - $shape_depth"] [expr 2.0 * $shape_width] [expr 2.0 * $shape_height] [expr 2.0 * $shape_depth]"
      puts $infile "ExclusionBox\t$shape_x $shape_y $shape_z [expr 2.0 * $shape_width] [expr 2.0 * $shape_height] [expr 2.0 * $shape_depth]"
    } elseif { $shape == "cylinder" } {
      set shape_width [lindex $shape_str 6]
      set shape_height [lindex $shape_str 7]
      set shape_depth [lindex $shape_str 8]
      set shape_r [lindex $shape_str 10]
      set point1x [expr "$shape_x - $shape_width"]; set point1y [expr "$shape_y - $shape_height"]; set point1z [expr "$shape_z - $shape_depth"]
      set point2x [expr "$shape_x + $shape_width"]; set point2y [expr "$shape_y + $shape_height"]; set point2z [expr "$shape_z + $shape_depth"]
      puts $infile "ExclusionCylinder\t$point1x $point1y $point1z $point2x $point2y $point2z $shape_r"
    }
  }
  foreach shape_str $::povme2::contig_list {
    set shape [lindex $shape_str 0]
    set shape_x [lindex $shape_str 2]
    set shape_y [lindex $shape_str 3]
    set shape_z [lindex $shape_str 4]
    if { $shape == "sphere" } {
      set shape_r [lindex $shape_str 6]
      puts $infile "SeedSphere\t$shape_x $shape_y $shape_z $shape_r"
    } elseif { $shape == "prism" } {
      set shape_width [lindex $shape_str 6]
      set shape_height [lindex $shape_str 7]
      set shape_depth [lindex $shape_str 8]
      if {$shape_width < 0} {set shape_width [expr "-1.0 * $shape_width"]} ;# can't have negative widths, etc.
      if {$shape_height < 0} {set shape_height [expr "-1.0 * $shape_height"]} ;# can't have negative widths, etc.
      if {$shape_depth < 0} {set shape_depth [expr "-1.0 * $shape_depth"]} ;# can't have negative widths, etc.
      #puts $infile "ContinuousPocketSeedBox\t[expr "$shape_x - $shape_width"] [expr "$shape_y - $shape_height"] [expr "$shape_z - $shape_depth"] [expr 2.0 * $shape_width] [expr 2.0 * $shape_height] [expr 2.0 * $shape_depth]"
      puts $infile "SeedBox\t$shape_x $shape_y $shape_z [expr 2.0 * $shape_width] [expr 2.0 * $shape_height] [expr 2.0 * $shape_depth]"
    } elseif { $shape == "cylinder" } {
      set shape_width [lindex $shape_str 6]
      set shape_height [lindex $shape_str 7]
      set shape_depth [lindex $shape_str 8]
      set shape_r [lindex $shape_str 10]
      set point1x [expr "$shape_x - $shape_width"]; set point1y [expr "$shape_y - $shape_height"]; set point1z [expr "$shape_z - $shape_depth"]
      set point2x [expr "$shape_x + $shape_width"]; set point2y [expr "$shape_y + $shape_height"]; set point2z [expr "$shape_z + $shape_depth"]
      puts $infile "SeedCylinder\t$point1x $point1y $point1z $point2x $point2y $point2z $shape_r"
    }
  }
  puts $infile "DistanceCutoff\t$::povme2::distance_cutoff"
  #puts $infile "SavePoints\t[::povme2::num_to_bool $::povme2::make_point_field_pdb]"
  puts $infile "ConvexHullExclusion\t[::povme2::num_to_bool $::povme2::exclude_points_outside_convex_hull]"
  if {$::povme2::use_contiguous_points == "True"} { ;# handle contiguous points
    # seed sphere
    # seed box
    # points criteria
  }
  puts $infile "OutputFilenamePrefix\t$::povme2::output_filename_prefix"
  puts $infile "CompressOutput\t[::povme2::num_to_bool $::povme2::compress_output]"
  #puts $infile "SaveIndividualPocketVolumes\t[::povme2::num_to_bool $::povme2::separate_volume_pdbs]"
  #puts $infile "SavePocketVolumesTrajectory\t[::povme2::num_to_bool $::povme2::volume_trajectory]"
  #puts $infile "OutputEqualNumPointsPerFrame\t[::povme2::num_to_bool $::povme2::equal_num_points_per_frame]"
  #puts $infile "SaveTabbedVolumeFile\t[::povme2::num_to_bool $::povme2::tabbed_volume_file]"
  #puts $infile "SaveVolumetricDensityDX\t[::povme2::num_to_bool $::povme2::volumetric_density_file]"
  #puts $infile "SaveColoredMap\t[::povme2::num_to_bool $::povme2::colored_density_files]"
  puts $infile "NumProcessors\t$::povme2::num_processors"
  #puts $infile "UseDiskNotMemory\t[::povme2::num_to_bool $::povme2::disk_instead_of_memory]"
  close $infile  
  #cd $olddir
  return $filename
}

proc ::povme2::povme2_volumes_tabbed { volumes_data } { ;# The window that displays the volumes after the POVME run is complete
  # $volumes_data needs to be a Nx2 list of lists
  variable w
  variable volumes_win
  set win_name $w.volumes
  if { [winfo exists $win_name] } {
    wm deiconify $volumes_win
    return
  }
  set volumes_win [toplevel "$w.volumes"]
  wm transient $volumes_win $w
  wm protocol $volumes_win WM_DELETE_WINDOW { ;# make sure that if a user closes a window, that it goes away the way we want it to
    grab release $::povme2::volumes_win
    after idle destroy $::povme2::volumes_win
  }
  wm title $volumes_win "POVME2 Volumes" ;# title of the window
  wm resizable $volumes_win yes yes ;# let the users resize the height and width of the window
  # now put all the widgets in place
  frame $volumes_win.header
  frame $volumes_win.data
  frame $volumes_win.buttons
  # pack the labels in the header
  grid [label $volumes_win.header.frame_number -text "Frame Number"] -column 0 -row 0
  grid [label $volumes_win.header.volume -text "Volume"] -column 1 -row 0 
  #listbox $w.mollist.m -relief raised -selectmode single ;# molecule listbox
  #$w.mollist.m
  #scrollbar $w.mollist.scrollbar -command [list $w.mollist.m yview] ;# molecule listbox scrollbar
  #$w.mollist.m configure -yscrollcommand [list $w.mollist.scrollbar set]
  #pack $w.mollist.m -side left -fill x -expand 1  ;# pack the listbox
  #pack $w.mollist.scrollbar -side left -fill y  ;# pack the scrollbar
  #
  # display the data
  frame $volumes_win.data.table
  set N [llength $volumes_data]
  for {set row 0} {$row < $N} {incr row} {
    for {set col 0} {$col < 2} {incr col} {
      set datum [lindex [lindex $volumes_data $row] $col]
      set label_name "row${row}_col${col}" ;# assign a name for this particular entry
      label $volumes_win.data.table.$label_name -text $datum -relief ridge
      grid $volumes_win.data.table.$label_name -column $col -row [expr "$row + 1"] ;# place them all into a regular grid
      grid configure $volumes_win.data.table.$label_name
    }
  }
  pack $volumes_win.data.table -side left -anchor w
  scrollbar $volumes_win.data.scrollbar -command [list $w.shapes.inclusion.m yview] ;# molecule listbox scrollbar
  $w.shapes.inclusion.m configure -yscrollcommand [list $w.shapes.inclusion.scrollbar set]
  pack $volumes_win.data.scrollbar -side left -anchor w
  # pack the buttons, and all the frames
  button $volumes_win.buttons.close_it -text Close -width 6 \
    -command {
    grab release $::povme2::volumes_win
    after idle destroy $::povme2::volumes_win
  }

  pack $volumes_win.buttons.close_it -side top -anchor n ;# $settings_win.okaycancel.default
  pack $volumes_win.header -side top -anchor n ;
  pack $volumes_win.data -side top -anchor n ;
  pack $volumes_win.buttons -side top -anchor n ;
  return
}

proc ::povme2::do_plot {x y {title "Volume plot" }} {
  if {$::povme2::plot_program == "multiplot"} {
    if [catch {package require multiplot} msg] {
      showMessage "Plotting in Multiplot not available: package multiplot not installed!\nDo you have the latest VMD version?"
      return
    }
    
    
      #set title $title
      set xlab "Frame"
    
    set ylab "Volume (A^3)"
    set plothandle [multiplot -title $title -xlabel $xlab -ylabel $ylab -nostats]
    
    set k 0
    
    set color black ;#[index2rgb $coln]
    set iname "Volume" ;# "[molinfo $i get name] ($i)"
      
    if {[llength $y] == 1} {
      $plothandle add $x $y -marker circle -radius 4 -nolines -fillcolor $color -linecolor $color -nostats -legend $iname
    } else {
      $plothandle add $x $y -marker point -radius 2 -fillcolor $color -linecolor $color -nostats -legend $iname
    }
    
    $plothandle replot
  }
}

proc ::povme2::povme2_settings {{mode point_properties}} { ;# mode can be graphics, advanced
	if {$mode == "files"} { ;# figure out which class of settings we are modifying. This determines the window's title and list of options
	  set arglist $::povme2::arglist_files
	  set window_title "POVME2 Files"
	} elseif {$mode == "point_properties"} {
	  set arglist $::povme2::arglist_point_properties
	  set window_title "POVME2 Point Properties"
	} elseif { $mode == "contiguous_points" } {
	  set arglist $::povme2::arglist_contiguous_points
	  set window_title "POVME2 Contiguous Points/Convex-Hull Exclusion"
	} elseif { $mode == "output" } {
	  set arglist $::povme2::arglist_output
	  set window_title "POVME2 Output"
	
	} elseif { $mode == "other"} {
	  set arglist $::povme2::arglist_other
	  set window_title "POVME2 Other Settings"
	}
	# specify window information, like for every window
	variable w
	variable settings_win
	set win_name $w.settings_$mode
	if { [winfo exists $win_name] } {
		wm deiconify $settings_win
		return
	}
	set settings_win [toplevel "$w.settings_$mode"]
	wm transient $settings_win $w
	wm protocol $settings_win WM_DELETE_WINDOW {
		grab release $::povme2::settings_win
		after idle destroy $::povme2::settings_win
	}
	wm title $settings_win "$window_title"
	wm resizable $settings_win no no
	
	# temp variables
	#set ::povme2::workdir_temp $::povme2::workdir
	#set ::povme2::povme2dir_temp $::povme2::povme2dir
	
	# menu widgets
	
	
	foreach arg $arglist {
          eval "set ::povme2::${arg}_temp \$::povme2::$arg" ;# set all temporary argument variables
        }
	#puts "shortest_path_opacity: $::povme2::shortest_path_opacity_temp"
	set len_arglist [llength $arglist] ;# get the length of the argument list
	set half_len_arglist [expr "([llength $arglist] / 2)"] ;# get half the length of the argument list for the two columns in the settings window
	
	frame $settings_win.leftcol ;# first do the left column
	frame $settings_win.rightcol
	
	
	#frame $settings_win.leftcol.temp_dir
	#pack [label $settings_win.leftcol.temp_dir.caption -text "Working Directory"] -side top -anchor w ; balloon $settings_win.leftcol.temp_dir.caption "description"
	#pack [entry $settings_win.leftcol.temp_dir.textbox -width 20 -textvariable ::povme2::workdir_temp] -side top -anchor w
	#pack $settings_win.leftcol.temp_dir -side top -anchor w -pady 10
#	
	#frame $settings_win.leftcol.povme2_dir
	#pack [label $settings_win.leftcol.povme2_dir.caption -text "povme2 location"] -side top -anchor w
	#pack [entry $settings_win.leftcol.povme2_dir.textbox -width 20 -textvariable ::povme2::povme2dir_temp] -side top -anchor w
	#pack $settings_win.leftcol.povme2_dir -side top -anchor w -pady 10
	
	#puts "len_arglist: $len_arglist"
	
	for {set i 0} {$i < $len_arglist} {incr i} {
	  set arg [lindex $arglist $i]
	  if {([info exists "::povme2::${arg}_format"])} {
	    eval "set arg_format_list \$::povme2::${arg}_format"
	  } else {
	    set arg_format_list "normal"
	  }
	  if {([info exists "::povme2::${arg}_caption"])} {
	    eval "set arg_caption \$::povme2::${arg}_caption"
	  } else {
	    set arg_caption "${arg}"
	  }
	  set arg_format [lindex $arg_format_list 0] ;# get the first element of the format list, because the rest of the stuff are parameters
	  if {[llength $arg_format_list] > 2} { ;# then we have more parameters for this argument
	    set arg_min [lindex $arg_format_list 1] ;# the min value this argument can take
	    set ::povme2::min_array($arg) $arg_min ;# set the array that holds minimum allowed values
	    set arg_max [lindex $arg_format_list 2] ;# the max value this argument can take
	    set ::povme2::max_array($arg) $arg_max ;# set the array that holds maximum allowed values
	  } elseif {[llength $arg_format_list] == 2} { ;# then its a list
	    set arg_list_options [lindex $arg_format_list 1]
	  }
	  if {$arg == "pdb_trajectory_filename" } { ;# or many other things...
	    continue ;# skip it, we already have it covered in the main menu
	  }
	  #if {($i < $half_len_arglist || $len_arglist < 6) } { ;# then we are on the left side
	    set column $settings_win.leftcol
	  #} 
	  if {$i >= $half_len_arglist && $len_arglist > 6} { ;# then we are on the right side
	    if {$arg_format != "little_slider_middle"} {
	      if {$i >= $half_len_arglist} {set column $settings_win.rightcol} ;# all this is very messy, but I'm trying to keep the little_sliders (which require the same frame) together on the same column
	    } elseif {$arg_format != "little_slider_end"} {
	      if {$i > [expr "$half_len_arglist + 1"]} {set column $settings_win.rightcol}
	    } else {
	      set column $settings_win.rightcol ;# if I'm not dealing with little_sliders, then no need to worry
	    }
	  }
	  set our_widget "$column.$arg"
	  #puts "now creating widget: $our_widget"
	  if {$arg_format == "little_slider_begin"}  {
	    set last_str_index [expr "[string length $arg] - 3"] ;# get the index of the third to last character in the string for the title
	    set widget_name [string range $arg 0 $last_str_index] ;# need a generic name for three subsequent widgets
	    set our_widget $column.$widget_name 
	    frame $our_widget ;# need to create the frame for the title and entries
	    pack [label $our_widget.caption -text "$widget_name (RGB)"] -side top -anchor w ; balloon $our_widget.caption ::povme2::descarray($arg) ;# title
	    frame $our_widget.subframe ;# a subframe for the three entries
	    #pack [entry $our_widget.subframe.textbox_begin -width 6 -textvariable ::povme2::${arg}_temp] -side left -anchor w; balloon $our_widget.subframe.textbox_begin ::povme2::descarray($arg)
	    #::povme2::make_scale "$our_widget.subframe.slider_begin" $arg $arg_max $arg_min
	    eval "scale $our_widget.subframe.slider_begin -to $arg_max -from $arg_min -orient horizontal -digits 3 -length 50 -sliderlength 20 -showvalue true -resolution 0.01 -command {set ::povme2::${arg}_temp }"
	    eval "$our_widget.subframe.slider_begin set \$::povme2::${arg}_temp"
	    pack $our_widget.subframe.slider_begin -side left -anchor w ; balloon $our_widget.subframe.slider_begin ::povme2::descarray($arg)
	  } elseif {$arg_format == "little_slider_middle"} { ;# for the middle field in a set of three small entries
	    set last_str_index [expr "[string length $arg] - 3"] ;# get the index of the third to last character in the string for the title
	    set widget_name [string range $arg 0 $last_str_index]
	    set our_widget $column.$widget_name
	    #pack [entry $our_widget.subframe.textbox_middle -width 6 -textvariable ::povme2::${arg}_temp] -side left -anchor w; balloon $our_widget.subframe.textbox_middle ::povme2::descarray($arg)
	    eval "scale $our_widget.subframe.slider_middle -to $arg_max -from $arg_min -orient horizontal -digits 3 -length 50 -sliderlength 20 -showvalue true -resolution 0.01 -command {set ::povme2::${arg}_temp }"
	    eval "$our_widget.subframe.slider_middle set \$::povme2::${arg}_temp"
	    pack $our_widget.subframe.slider_middle -side left -anchor w ; balloon $our_widget.subframe.slider_middle ::povme2::descarray($arg)
	  } elseif {$arg_format == "little_slider_end"} { ;# for the last field in a set of three small entries
	    set last_str_index [expr "[string length $arg] - 3"] ;# get the index of the third to last character in the string for the title
	    set widget_name [string range $arg 0 $last_str_index]
	    set our_widget $column.$widget_name
	    #pack [entry $our_widget.subframe.textbox_end -width 6 -textvariable ::povme2::${arg}_temp] -side left -anchor w; balloon $our_widget.subframe.textbox_end ::povme2::descarray($arg)
	    eval "scale $our_widget.subframe.slider_end -to $arg_max -from $arg_min -orient horizontal -digits 3 -length 50 -sliderlength 20 -showvalue true -resolution 0.01 -command {set ::povme2::${arg}_temp }"
	    eval "$our_widget.subframe.slider_end set \$::povme2::${arg}_temp"
	    pack $our_widget.subframe.slider_end -side left -anchor w ; balloon $our_widget.subframe.slider_end ::povme2::descarray($arg)
	    pack $our_widget.subframe -side top -anchor w -fill x -expand 1 ;# close the entry subframe
	    pack $our_widget -side top -anchor w -pady 10 -padx 10 ;# close the entire frame
	  } elseif {$arg_format == "slider"} { ;# a slider widget
	    frame $our_widget ;# balloon $our_widget "Test Help"
	    pack [label $our_widget.caption -text "$arg_caption"] -side top -anchor w ; balloon $our_widget.caption ::povme2::descarray($arg)
	    frame $our_widget.cutoffslider
	    #pack [entry $our_widget.cutoffslider.textbox -width 6 -textvariable ::povme2::${arg}_temp] -side left -fill x
	    
	    if {${arg} == "node_sphere_opacity"} { ;# this means that this argument is node_sphere opacity, and that a structure is loaded
	      eval "scale $our_widget.cutoffslider.slider -to $arg_max -from $arg_min -orient horizontal -digits 3 -length 150 -showvalue true -resolution 0.01 -command {if {[lsearch [material list] node_spheres] != -1} {material change opacity node_spheres \$::povme2::${arg}_temp}; set ::povme2::${arg}_temp }"
	    } elseif {${arg} == "shortest_path_opacity" || $arg == "longest_path_opacity"} {
	      ::povme2::make_scale $our_widget.cutoffslider.slider $arg $arg_max $arg_min
	    
	    } else { ;# otherwise, just treat it normally
	      eval "scale $our_widget.cutoffslider.slider -to $arg_max -from $arg_min -orient horizontal -digits 3 -length 150 -showvalue true -resolution 0.01 -command {set ::povme2::${arg}_temp }"
	    }
	    eval "$our_widget.cutoffslider.slider set \$::povme2::${arg}_temp"
	    pack $our_widget.cutoffslider.slider -side left -fill x -expand 1
	    pack $our_widget.cutoffslider -side left -fill x -expand 1
	    pack $our_widget -side top -anchor w -pady 10 -padx 10
	  } elseif {$arg_format == "file" || $arg_format == "directory"} {
	    if {$arg_format == "file"} { set which_tk_getopen "tk_getOpenFile" 
	    } elseif {$arg_format == "directory"} { set which_tk_getopen "tk_chooseDirectory" } ;# two different things the user might choose: files or directories
	    frame $our_widget -relief groove -bd 2  ;# balloon $our_widget "Test Help"
	    pack [label $our_widget.caption -text "$arg_caption"] -side top -anchor w ; balloon $our_widget.caption ::povme2::descarray($arg)
	    eval "pack \[entry $our_widget.textbox -width 25 -textvariable ::povme2::${arg}_temp -validate all -validatecommand {::povme2::isFile %P $our_widget.textbox}\] -side top -anchor w; balloon $our_widget.textbox ::povme2::descarray($arg)" ;# NOTE: will need to make an array in place of desclist
	    eval "pack \[button $our_widget.browse -width 10 -text \"Browse...\" -command {set ::povme2::${arg}_temp \[$which_tk_getopen\]}\] -side top -anchor w; balloon $our_widget.browse ::povme2::descarray($arg)" ;# NOTE: will need to make an array in place of desclist
	    pack $our_widget -side top -anchor w -pady 10 -padx 10 -fill x
	  } elseif {$arg_format == "list"} {
	    frame $our_widget -relief groove -bd 2 ;# balloon $our_widget "Test Help"
	    pack [label $our_widget.caption -text "$arg_caption"] -side top -anchor w ; balloon $our_widget.caption ::povme2::descarray($arg)
	    #pack [entry $our_widget.textbox -width 20 -textvariable ::povme2::${arg}_temp] -side top -anchor w; balloon $our_widget.textbox ::povme2::descarray($arg);# NOTE: will need to make an array in place of desclist
	    foreach option $arg_list_options {
	      pack [radiobutton $our_widget.option_$option -text "$option" -value $option -variable ::povme2::${arg}_temp] -side top -anchor w ;
	    }
	    pack $our_widget -side top -anchor w -pady 10 -padx 10 -fill x
	  } elseif {$arg_format == "checkbox"} {
	    frame $our_widget -relief groove -bd 2 ;# balloon $our_widget "Test Help"
	    pack [label $our_widget.caption -text "$arg_caption"] -side top -anchor w ; balloon $our_widget.caption ::povme2::descarray($arg)
	    #pack [entry $our_widget.textbox -width 20 -textvariable ::povme2::${arg}_temp] -side top -anchor w; balloon $our_widget.textbox ::povme2::descarray($arg);# NOTE: will need to make an array in place of desclist
	    pack [checkbutton $our_widget.checkbox -text "" -variable ::povme2::${arg}_temp] -side top -anchor w ;
	    pack $our_widget -side top -anchor w -pady 10 -padx 10 -fill x
	  } else {
	    frame $our_widget ;# balloon $our_widget "Test Help"
	    pack [label $our_widget.caption -text "$arg_caption"] -side top -anchor w ; balloon $our_widget.caption ::povme2::descarray($arg)
	    pack [entry $our_widget.textbox -width 20 -textvariable ::povme2::${arg}_temp] -side top -anchor w; balloon $our_widget.textbox ::povme2::descarray($arg);# NOTE: will need to make an array in place of desclist
	    pack $our_widget -side top -anchor w -pady 10 -padx 10
	  }
	}
	
	if { $mode == "graphics" } { ;# then place the contrast bar
	  set our_widget "$column.path_contrast"
	  set arg_max [lindex $::povme2::path_contrast_format 2]
	  set arg_min [lindex $::povme2::path_contrast_format 1]
	  frame $our_widget ;# balloon $our_widget "Test Help"
	  pack [label $our_widget.caption -text "path contrast"] -side top -anchor w ; balloon $our_widget.caption "Adjusts the contrast between paths of varying lengths"
	  frame $our_widget.cutoffslider
	  eval "scale $our_widget.cutoffslider.slider -to $arg_max -from $arg_min -orient horizontal -digits 3 -length 150 -showvalue true -resolution 0.1 -command {::povme2::adjust_path_contrast }"
	  # balloon?
	  $our_widget.cutoffslider.slider set $::povme2::path_contrast
	  pack $our_widget.cutoffslider.slider -side left -fill x -expand 1
	  pack $our_widget.cutoffslider -side left -fill x -expand 1
	  pack $our_widget -side top -anchor w -pady 10 -padx 10
	}

	grid $settings_win.leftcol -column 0 -row 0
	grid $settings_win.rightcol -column 1 -row 0

#frame $settings_win.leftcol.checkboxes
	#pack [checkbutton $settings_win.leftcol.checkboxes.sovereignty -text "No intrastructural clustering" -variable ::povme2::sovereignty_temp -onvalue 1 -offvalue 0] -side top -anchor w
	#pack [checkbutton $settings_win.leftcol.checkboxes.hydrophobic -text "Hydrophobic calculations" -variable ::povme2::find_hydrophobic_temp -onvalue 1 -offvalue 0] -side top -anchor w
	
	
	frame $settings_win.okaycancel
	# need to use an eval so that the $arglist doesn't have to be global
	eval "button \$settings_win.okaycancel.okay -text OK -width 6 \
	  -command {;\
	  	foreach arg [list $arglist] {\
	  	  eval \"set value \\\$::povme2::\${arg}_temp\";\
	  	  if {[info exists ::povme2::min_array(\$arg)] && [info exists ::povme2::max_array(\$arg)]} {\
	  	    if {\$value < \$::povme2::min_array(\$arg)} { \
	  	      tk_messageBox -type ok -title \"Error\" -icon error -message \"Error: the value you entered for \$arg is too low. It must be greater than \$::povme2::min_array(\$arg).\";\
	  	      return;\
	  	    } ;\
	  	    if {\$value > \$::povme2::max_array(\$arg)} {\
	  	      tk_messageBox -type ok -title \"Error\" -icon error -message \"Error: the value you entered for \$arg is too high. It must be less than \$::povme2::max_array(\$arg).\";\
	  	      return;\
	  	    };\
	  	  };\
     		  set ::povme2::\$arg \$value;\
     		} ;\
	  	grab release $::povme2::settings_win;\
	  	after idle destroy $::povme2::settings_win;\
	}"
	button $settings_win.okaycancel.cancel -text Cancel -width 6 \
	  -command {
	  	grab release $::povme2::settings_win
	  	after idle destroy $::povme2::settings_win
	}
	#button $settings_win.okaycancel.default -text Default -width 6 \
	#  -command {
	#  	set ::povme2::workdir_temp $::povme2::workdir
	#  	#...
	#	
	#}
	pack $settings_win.okaycancel.okay $settings_win.okaycancel.cancel -side left ;# $settings_win.okaycancel.default
	grid $settings_win.okaycancel -column 0 -row 1
	
}

proc ::povme2::draw_cylinder { mol centerx centery centerz width height depth radius {replace_list ""}} {
  if {[llength $replace_list] >= 3} { ;# then it must be a sphere that just needs to be deleted
    foreach graphic_id $replace_list {
      graphics $mol delete $graphic_id
    }
    set replace_list "" ;# then just remove all the graphic ids
  } elseif { [llength $replace_list] == 1 } { ;# this is making the assumption that the list will have exactly 12 graphic ids
    graphics $mol replace [lindex $replace_list 0]
  }
  set point1x [expr "$centerx - $width"]
  set point1y [expr "$centery - $height"]
  set point1z [expr "$centerz - $depth"]
  set point2x [expr "$centerx + $width"]
  set point2y [expr "$centery + $height"]
  set point2z [expr "$centerz + $depth"]
  set graphics_id [graphics $mol cylinder "$point1x $point1y $point1z" "$point2x $point2y $point2z" radius $radius resolution $::povme2::cylinder_resolution filled yes]
  return $graphics_id
}

proc ::povme2::draw_prism { mol centerx centery centerz width height depth {replace_list ""}} {
  # Draws a box (or prism) in the molecule "mol" that consists of 12 triangle graphics
  set graphics_id_list {} ;# keep track of all graphics indeces
  #set width [expr "$a - $centerx"]
  #set height [expr "$b - $centery"] ;# dimensions of the box
  #set depth [expr "$c - $centerz"]
  set center "$centerx $centery $centerz"
  set widths "$width $height $depth"
  set faces {}
  set triangles {}
  set counter 0
  if {[llength $replace_list] >= 3} { ;# then it must be a sphere that just needs to be deleted
    foreach graphic_id $replace_list {
      graphics $mol delete $graphic_id
    }
    set replace_list "" ;# then just remove all the graphic ids
  }
  foreach index {0 1 2} { ;# for each of the x y z - facing surfaces
    foreach sign {-1.0 1.0} { ;# each needs a parallel face, making 6 total
      set face $center
      lset face $index [expr "[lindex $face $index] + $sign * [lindex $widths $index]"]
      lappend faces $face
      set index1 [expr "($index + 1) % 3"] ;# the index to assign the y coordinate
      set index2 [expr "($index + 2) % 3"] ;# the index to assign the z coordinate
      
      set side1 [lindex $widths $index1]
      set side2 [lindex $widths $index2]
      
      set triangle1 "{$face} {$face} {$face}"
      lset triangle1 0 $index1 [expr "[lindex $face $index1] + [lindex $widths $index1]"]
      lset triangle1 1 $index1 [expr "[lindex $face $index1] - [lindex $widths $index1]"]
      lset triangle1 2 $index1 [expr "[lindex $face $index1] - [lindex $widths $index1]"]
      
      lset triangle1 0 $index2 [expr "[lindex $face $index2] + [lindex $widths $index2]"]
      lset triangle1 1 $index2 [expr "[lindex $face $index2] + [lindex $widths $index2]"]
      lset triangle1 2 $index2 [expr "[lindex $face $index2] - [lindex $widths $index2]"]
      
      set triangle2 $triangle1 ;# but each quadrilateral is made of 2 triangles
      lset triangle2 1 $index1 [expr "[lindex $face $index1] + [lindex $widths $index1]"]
      lset triangle2 1 $index2 [expr "[lindex $face $index2] - [lindex $widths $index2]"]
      
      
      # 1st triangle
      if { [llength $replace_list] == 12 } { ;# this is making the assumption that the list will have exactly 12 graphic ids
          graphics $mol replace [lindex $replace_list $counter]
      }
      lappend graphics_id_list [graphics $mol triangle [lindex $triangle1 0] [lindex $triangle1 1] [lindex $triangle1 2]] ;# draw the triangle
      incr counter
      # 2nd triangle
      if { [llength $replace_list] == 12 } { ;# this is making the assumption that the list will have exactly 12 graphic ids
          graphics $mol replace [lindex $replace_list $counter]
      }
      lappend graphics_id_list [graphics $mol triangle [lindex $triangle2 0] [lindex $triangle2 1] [lindex $triangle2 2]] ;# draw the triangle
      incr counter
    }
  }
  
  return $graphics_id_list ;# return the list of graphics indeces
}

proc ::povme2::make_scale {name arg arg_max arg_min} {
  eval "scale $name -to $arg_max -from $arg_min -orient horizontal -digits 3 -length 150 -sliderlength 40 -showvalue true -resolution 0.01 -command {::povme2::adjust_wispmaterial ::povme2::shortest_path_opacity_temp ::povme2::longest_path_opacity_temp ::povme2::shortest_path_r ::povme2::longest_path_r ::povme2::num_paths; set ::povme2::${arg}_temp }"
  #eval "$name set \$::povme2::${arg}_temp"
  
}


proc ::povme2::isFile {f w} {
    if { [file exists $f] && ([file isfile $f] || [file isdirectory $f]) } {
        $w configure -fg black
    } else {
        $w configure -fg red
    }
    return 1;
}

proc ::povme2::povme2_init {} {
  ;# the correct OS directory separator must be defined
  switch [vmdinfo arch] { 
    WIN64 -
    WIN32 {
      set ::povme2::sep "\\"
      set ::povme2:tempdir ""
    } default {
      set ::povme2::sep "/"
    }
  }
  #set ::povme2::povme2_directory "[file dirname $argv0]${sep}povme2.py"
  #if { [::povme2::find_py] == True } { ::povme2::parse_helpfile } ;# search for the python script
  if { $::povme2::graphics_mol == "" } {
    #set ::povme2::graphics_mol [mol load graphics "povme2_graphics"] ;# creating the graphics mol for povme
    #mol rename $::povme2::graphics_mol "POVME2 Graphics" ;# rename the graphics molecule
    #graphics $::povme2::graphics_mol material $::povme2::graphics_material
  }
  ::povme2::povme2_mainwin ;# call the main window
}

### Balloon help ###
proc balloon {w help} {
    if {![info exists $help]} {return}
    eval "bind $w <Any-Enter> \"after 1000 [list balloon:show %W [list \$$help]]\""
    bind $w <Any-Leave> "destroy %W.balloon"
}
proc balloon:show {w arg} {
    if {[eval winfo containing  [winfo pointerxy .]]!=$w} {return}
    set top $w.balloon
    catch {destroy $top}
    toplevel $top -bd 1 -bg black
    wm overrideredirect $top 1
    if {[string equal [tk windowingsystem] aqua]}  {
        ::tk::unsupported::MacWindowStyle style $top help none
    }   
    pack [message $top.txt -aspect 200 -bg lightyellow \
            -font fixed -text $arg]
    set wmx [winfo rootx $w]
    set wmy [expr [winfo rooty $w]+[winfo height $w]]
    wm geometry $top \
      [winfo reqwidth $top.txt]x[winfo reqheight $top.txt]+$wmx+$wmy
    raise $top
}
 

#this gets called by VMD the first time the menu is opened
proc povme2_tk_cb {} {
  ::povme2::povme2_init	;# start the tool
  return $::povme2::w
}
