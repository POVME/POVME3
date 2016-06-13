package require Tk

### POVME dialogue window ###
#set dialog [toplevel .wait]
#pack [message $dialog.frame.msg -text "Please wait\nPOVME2 is running..." -width 100] -padx 100 -pady 100
#pack [frame $dialog.frame -borderwidth 2 -relief raised] -fill both
#pack configure progress {W}

proc progress {W args} {
 
  array set map [list \
    -bd -borderwidth \
    -bg -background \
  ]
 
  array set arg [list \
    -activebackground blue \
    -borderwidth 1 \
    -from 0 \
    -to 100 \
    -orient horizontal \
    -sliderrelief flat \
    -sliderlength 0 \
    -troughcolor #AAAAAA \
    -showvalue 0 \
    -label 0 \
    -state active \
  ]
 
  foreach {option value} $args {
    if { [info exists map($option)] } { set option $map($option) }
    set arg($option) $value
  }

  eval [linsert [array get arg] 0 scale $W]

  bind $W <Enter>  {break}
  bind $W <Leave>  {break}
  bind $W <Motion> {break}
  bind $W <1>      {break}
  bind $W <ButtonRelease-1> {break}

  bind $W <Configure> [list [namespace current]::progress:redraw %W]
 
  return $W
 
}
 
proc progress:redraw {W} {
  set value [$W cget -label]
  set bd    [$W cget -bd]
  set ht    [$W cget -highlightthickness]
  set from  [$W cget -from]
  set to    [$W cget -to]
  set w [winfo width $W]
  set tw [expr {$w - (4 * $bd) - (2 * $ht)}]
  set range [expr {$to - $from}]
  set pc [expr {($value - $from) * 1.0 / $range}]
  set sl [expr {round($pc * $tw)}]
  $W configure -sliderlength $sl
  return
}

proc progress:set {W value} {
  $W configure -label $value
  progress:redraw $W
  update
  return
}

proc go {W value} {
  progress:set $W $value
  set outputcheck_running 1
    puts "Starting go"

  while {$outputcheck_running == 1} {
    ### change to location of output.txt ###
    set fptr [open "../output.txt" r]
    set contents [read -nonewline $fptr]
    close $fptr
    set splitCont [split $contents "\n"]
  
    set regexMatch 0 
    set lineSearch 0
    set totalnum 250
    #set totalnum {molinfo top get numframes}

    while {$regexMatch == 0} {
      set checkLine [lindex $splitCont end-$lineSearch] 
      set matchedProcessingLine [regexp {Further processing frame ([0-9]*)} $checkLine -> processingFrame] 

      if {$matchedProcessingLine == 1} {
        # Calculating scale value out of 100 (0-40)
        set value [expr double($processingFrame - 1)*50/($totalnum)]
        puts $value
        #after 50 [list go $W $value]
        set regexMatch 1 
      }

      set matchedAnalyzingLine [regexp {Frame ([0-9]*)} $checkLine -> analyzingFrame]
      if {$matchedAnalyzingLine == 1} {
        # Calculating scale value out of 100 (50-90)
        set value [expr double($analyzingFrame - 1)*50/($totalnum) + 50]
        puts $value
        #after 50 [list go $W $value]
        set regexMatch 1
      }

      set matchedEndLine "Execution time"
      if {[string match $matchedEndLine* $checkLine]} {
        set value 100
        puts $value
        #after 50 [list go $W $value]
        set regexMatch 1
        set outputcheck_running 0 
      }
      set lineSearch [expr {$lineSearch + 1}]
    }
    progress:set $W $value
    #update
    #set outputcheck_running 0
    after 1000
  }
  return
}

  if { [info exists argv0] && [string equal [info script] $argv0] } {
    progress .sc

    button .go -text go -default active \
      -command [list [namespace current]::go .sc 0]

    pack .sc -side top    -expand 1 -padx 40 -fill both
    pack .go -side bottom -expand 0 -padx 40 -fill none -anchor se
  }
