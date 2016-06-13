#!/usr/bin/python

from Tkinter import *
#------------------------------------------------------------------------------#
#                                                                              #
#                                  POVME2_GUI                                  #
#                                                                              #
#------------------------------------------------------------------------------#
class POVME2_GUI(Frame):
    def __init__(self,Master=None,**kw):
        #
        #Your code here
        #

        apply(Frame.__init__,(self,Master),kw)
        self.__Frame2 = Frame(self)
        self.__Frame2.pack(side='top')
        self.__Label1 = Label(self.__Frame2,background='#000000'
            ,foreground='White',text='Point Field',width=620)
        self.__Label1.pack(side='top')
        self.__Frame36 = Frame(self)
        self.__Frame36.pack(pady=1,side='top')
        self.__Frame1 = Frame(self,width=900)
        self.__Frame1.pack(side='top')
        self.__Frame4 = Frame(self,width=900)
        self.__Frame4.pack(side='top')
        self.__Frame3 = Frame(self)
        self.__Frame3.pack(side='top')
        self.__PointRegionsList = Listbox(self.__Frame3,height=5,width=40)
        self.__PointRegionsList.pack(side='top')
        self.__DeletePointRegion = Button(self.__Frame3,text='Delete Item')
        self.__DeletePointRegion.pack(side='bottom')
        self.__DeletePointRegion.bind('<ButtonRelease-1>' \
            ,self.__on_DeletePointRegion_ButRel_1)
        self.__Frame18 = Frame(self)
        self.__Frame18.pack(pady=1,side='top')
        self.__Label10 = Label(self.__Frame18,background='#000000'
            ,foreground='White',text='Point Properties',width=620)
        self.__Label10.pack(side='top')
        self.__Frame19 = Frame(self)
        self.__Frame19.pack(side='top')
        self.__Frame22 = Frame(self)
        self.__Frame22.pack(side='top')
        self.__Label13 = Label(self.__Frame22,background='#000000'
            ,foreground='White',text='PDB Filename',width=620)
        self.__Label13.pack(expand='yes',fill='both',pady=1,side='top')
        self.__Frame17 = Frame(self)
        self.__Frame17.pack(side='top')
        self.pdb_filename = StringVar()
        self.__PDBFilename = Label(self.__Frame17,text='Select PDB File...'
            ,textvariable=self.pdb_filename)
        self.__PDBFilename.pack(side='left')
        self.__FindPDBFilename = Button(self.__Frame17,text='Find PDB File')
        self.__FindPDBFilename.pack(padx=5,side='left')
        self.__FindPDBFilename.bind('<ButtonRelease-1>' \
            ,self.__on_FindPDBFilename_ButRel_1)
        self.__Frame23 = Frame(self)
        self.__Frame23.pack(side='top')
        self.__Label34 = Label(self.__Frame23,background='#000000'
            ,foreground='White',text='Contiguous Points/Convex-Hull Exclusion'
            ,width=720)
        self.__Label34.pack(expand='yes',fill='both',pady=1,side='top')
        self.__Frame6 = Frame(self)
        self.__Frame6.pack(side='top')
        self.__Frame24 = Frame(self)
        self.__Frame24.pack(side='top')
        self.__Frame33 = Frame(self,height=5)
        self.__Frame33.pack(side='top')
        self.__Frame41 = Frame(self)
        self.__Frame41.pack(side='top')
        self.__ContiguousRegionsList = Listbox(self.__Frame41,height=5,width=40)
        self.__ContiguousRegionsList.pack(side='top')
        self.__DeleteContiguousRegion = Button(self.__Frame41,state='disabled'
            ,text='Delete Item')
        self.__DeleteContiguousRegion.pack(side='bottom')
        self.__DeleteContiguousRegion.bind('<ButtonRelease-1>' \
            ,self.__on_DeleteContiguousRegion_ButRel_1)
        self.__Frame37 = Frame(self)
        self.__Frame37.pack(side='top')
        self.contig_points_criteria = StringVar()
        self.__ContigPointsCriteria = Entry(self.__Frame37,state='disabled'
            ,textvariable=self.contig_points_criteria,width=5)
        self.__ContigPointsCriteria.pack(side='right')
        self.__ContigLabel11 = Label(self.__Frame37,state='disabled'
            ,text='Contiguous Point Criteria:')
        self.__ContigLabel11.pack(side='right')
        self.__Frame28 = Frame(self)
        self.__Frame28.pack(pady=1,side='top')
        self.__Label25 = Label(self.__Frame28,background='#000000'
            ,foreground='White',text='Output',width=900)
        self.__Label25.pack(side='top')
        self.__Frame12 = Frame(self)
        self.__Frame12.pack(side='top')
        self.__Frame29 = Frame(self)
        self.__Frame29.pack(side='top')
        self.__Frame45 = Frame(self)
        self.__Frame45.pack(side='top')
        self.__Frame30 = Frame(self)
        self.__Frame30.pack(side='top')
        self.__Frame13 = Frame(self)
        self.__Frame13.pack(side='top')
        self.__Frame11 = Frame(self.__Frame1,height=30,width=410)
        self.__Frame11.pack(expand='yes',fill='both',side='left')
        self.__Label32 = Label(self.__Frame11)
        self.__Label32.pack(padx=40,side='left')
        self.__Label2 = Label(self.__Frame11,text='Center X/Y/Z:')
        self.__Label2.pack(side='left')
        self.center_x_1 = StringVar()
        self.__CenterX1 = Entry(self.__Frame11,textvariable=self.center_x_1
            ,width=5)
        self.__CenterX1.pack(anchor='w',side='left')
        self.__Label3 = Label(self.__Frame11,text='/')
        self.__Label3.pack(side='left')
        self.center_y_1 = StringVar()
        self.__CenterY1 = Entry(self.__Frame11,textvariable=self.center_y_1
            ,width=5)
        self.__CenterY1.pack(side='left')
        self.__Label4 = Label(self.__Frame11,text='/')
        self.__Label4.pack(side='left')
        self.center_z_1 = StringVar()
        self.__CenterZ1 = Entry(self.__Frame11,textvariable=self.center_z_1
            ,width=5)
        self.__CenterZ1.pack(side='left')
        self.__Frame5 = Frame(self.__Frame1,width=110)
        self.__Frame5.pack(expand='yes',fill='both',side='left')
        self.__AddSphere1 = Button(self.__Frame5,text='Add Sphere')
        self.__AddSphere1.pack(anchor='e',side='top')
        self.__AddSphere1.bind('<ButtonRelease-1>' \
            ,self.__on_AddSphere1_ButRel_1)
        self.__Frame9 = Frame(self.__Frame1,width=140)
        self.__Frame9.pack(expand='yes',fill='both',side='left')
        self.__ExcludeSphere1 = Button(self.__Frame9,text='Exclude Sphere')
        self.__ExcludeSphere1.pack(anchor='w',side='top')
        self.__ExcludeSphere1.bind('<ButtonRelease-1>' \
            ,self.__on_ExcludeSphere1_ButRel_1)
        self.__Frame7 = Frame(self.__Frame4,height=30,width=100)
        self.__Frame7.pack(expand='yes',fill='both',side='left')
        self.__Label5 = Label(self.__Frame7,text='Radius:')
        self.__Label5.pack(side='left')
        self.radius_1 = StringVar()
        self.__Radius1 = Entry(self.__Frame7,textvariable=self.radius_1,width=5)
        self.__Radius1.pack(side='left')
        self.__Frame14 = Frame(self.__Frame4,width=30)
        self.__Frame14.pack(expand='yes',fill='both',side='left')
        self.__Label6 = Label(self.__Frame14,text='OR')
        self.__Label6.pack(expand='yes',fill='both',side='top')
        self.__Frame16 = Frame(self.__Frame4,width=280)
        self.__Frame16.pack(expand='yes',fill='both',side='left')
        self.__Label31 = Label(self.__Frame16)
        self.__Label31.pack(padx=2,side='left')
        self.__Label9 = Label(self.__Frame16,text='Length X/Y/Z')
        self.__Label9.pack(side='left')
        self.length_x_1 = StringVar()
        self.__LengthX1 = Entry(self.__Frame16,textvariable=self.length_x_1
            ,width=5)
        self.__LengthX1.pack(side='left')
        self.__Label7 = Label(self.__Frame16,text='/')
        self.__Label7.pack(side='left')
        self.length_y_1 = StringVar()
        self.__LengthY1 = Entry(self.__Frame16,textvariable=self.length_y_1
            ,width=5)
        self.__LengthY1.pack(side='left')
        self.__Label8 = Label(self.__Frame16,text='/')
        self.__Label8.pack(side='left')
        self.length_z_1 = StringVar()
        self.__LengthZ1 = Entry(self.__Frame16,textvariable=self.length_z_1
            ,width=5)
        self.__LengthZ1.pack(side='left')
        self.__Frame10 = Frame(self.__Frame4,width=110)
        self.__Frame10.pack(expand='yes',fill='both',side='left')
        self.__AddPrism1 = Button(self.__Frame10,text='Add Prism')
        self.__AddPrism1.pack(anchor='e',side='top')
        self.__AddPrism1.bind('<ButtonRelease-1>',self.__on_AddPrism1_ButRel_1)
        self.__Frame8 = Frame(self.__Frame4,width=140)
        self.__Frame8.pack(expand='yes',fill='both',side='left')
        self.__ExcludePrism1 = Button(self.__Frame8,text='Exclude Prism')
        self.__ExcludePrism1.pack(anchor='w',side='top')
        self.__ExcludePrism1.bind('<ButtonRelease-1>' \
            ,self.__on_ExcludePrism1_ButRel_1)
        self.__Frame20 = Frame(self.__Frame19)
        self.__Frame20.pack(side='left')
        self.__Label11 = Label(self.__Frame20,text='Grid Point Spacing:')
        self.__Label11.pack(side='left')
        self.point_spacing = StringVar()
        self.__PointSpacing = Entry(self.__Frame20
            ,textvariable=self.point_spacing,width=5)
        self.__PointSpacing.pack(side='left')
        self.__Frame21 = Frame(self.__Frame19)
        self.__Frame21.pack(side='left')
        self.__Label33 = Label(self.__Frame21)
        self.__Label33.pack(padx=10,side='left')
        self.__Label12 = Label(self.__Frame21,text='Distance Cutoff:')
        self.__Label12.pack(side='left')
        self.padding = StringVar()
        self.__Padding = Entry(self.__Frame21,textvariable=self.padding,width=5)
        self.__Padding.pack(side='left')
        self.__Frame54 = Frame(self.__Frame19)
        self.__Frame54.pack(side='left')
        self.__Label38 = Label(self.__Frame54)
        self.__Label38.pack(padx=10,side='left')
        self.__Label39 = Label(self.__Frame54
            ,text='Make Point-Field PDB (No Volume Calc)')
        self.__Label39.pack(side='left')
        self.point_field_only = StringVar()
        self.__Checkbutton1 = Checkbutton(self.__Frame54
            ,variable=self.point_field_only)
        self.__Checkbutton1.pack(side='left')
        self.__Frame43 = Frame(self.__Frame6)
        self.__Frame43.pack(side='left')
        self.__Label16 = Label(self.__Frame43)
        self.__Label16.pack(padx=20,side='right')
        self.use_contig_points = StringVar()
        self.__CheckUseContigPoints = Checkbutton(self.__Frame43
            ,variable=self.use_contig_points, command=self.__on_CheckUseContigPoints_ButRel_1)
        self.__CheckUseContigPoints.pack(side='right')
        
        # Note, when the below was used, the checkbutton didn't function properly
        # on Mac OSX, so I used the command=callback technique above.
        #self.__CheckUseContigPoints.bind('<ButtonRelease-1>' \
        #    ,self.__on_CheckUseContigPoints_ButRel_1)
        
        self.__Label15 = Label(self.__Frame43,text='Use Contiguous Points')
        self.__Label15.pack(side='right')
        self.__Frame42 = Frame(self.__Frame6)
        self.__Frame42.pack(side='left')
        self.__Label14 = Label(self.__Frame42
            ,text='Exclude Points Outside Convex Hull')
        self.__Label14.pack(side='left')
        self.convex_hull = StringVar()
        self.__CheckConvexHull = Checkbutton(self.__Frame42
            ,variable=self.convex_hull)
        self.__CheckConvexHull.pack(side='left')
        self.__Frame25 = Frame(self.__Frame24)
        self.__Frame25.pack(side='left')
        self.__Frame26 = Frame(self.__Frame24,width=110)
        self.__Frame26.pack(expand='yes',fill='both',side='left')
        self.__AddSphere2 = Button(self.__Frame26,state='disabled'
            ,text='Add Sphere')
        self.__AddSphere2.pack(anchor='e',side='top')
        self.__AddSphere2.bind('<ButtonRelease-1>' \
            ,self.__on_AddSphere2_ButRel_1)
        self.__Frame27 = Frame(self.__Frame24,width=110)
        self.__Frame27.pack(expand='yes',fill='both',side='left')
        self.__AddPrism2 = Button(self.__Frame27,state='disabled'
            ,text='Add Prism')
        self.__AddPrism2.pack(anchor='w',side='top')
        self.__AddPrism2.bind('<ButtonRelease-1>',self.__on_AddPrism2_ButRel_1)
        self.__Frame48 = Frame(self.__Frame12)
        self.__Frame48.pack(side='left')
        self.points_files_prefix = StringVar()
        self.__PointsFilesPrefix = Entry(self.__Frame48
            ,textvariable=self.points_files_prefix,width=10)
        self.__PointsFilesPrefix.pack(side='right')
        self.__Label26 = Label(self.__Frame48,text='Output Filename Prefix:')
        self.__Label26.pack(side='right')
        self.__Frame47 = Frame(self.__Frame12)
        self.__Frame47.pack(side='left')
        self.__Label37 = Label(self.__Frame47,state='disabled')
        self.__Label37.pack(padx=10,side='left')
        self.__Label27 = Label(self.__Frame47,text='Compress Output')
        self.__Label27.pack(side='left')
        self.compress_points_files = StringVar()
        self.__CompressPointsFiles = Checkbutton(self.__Frame47
            ,variable=self.compress_points_files)
        self.__CompressPointsFiles.pack(side='left')
        self.__Frame49 = Frame(self.__Frame29)
        self.__Frame49.pack(pady=1,side='left')
        self.save_individual_pocket_volumes = StringVar()
        self.__CheckSaveIndividualPocketVolumes = Checkbutton(self.__Frame49
            ,variable=self.save_individual_pocket_volumes, command=self.__on_CheckSaveIndividualPocketVolumes_ButRel_1)
        self.__CheckSaveIndividualPocketVolumes.pack(side='right')
        
        # Note, when the below was used, the checkbutton didn't function properly
        # on Mac OSX, so I used the command=callback technique above.
        #self.__CheckSaveIndividualPocketVolumes.bind('<ButtonRelease-1>' \
        #    ,self.__on_CheckSaveIndividualPocketVolumes_ButRel_1)
        
        self.__Label19 = Label(self.__Frame49,text='Separate Volume PDBs')
        self.__Label19.pack(side='right')
        self.__Frame46 = Frame(self.__Frame29)
        self.__Frame46.pack(pady=1,side='left')
        self.__Label20 = Label(self.__Frame46)
        self.__Label20.pack(padx=15,side='left')
        self.__Label21 = Label(self.__Frame46,text='Volume Trajectory')
        self.__Label21.pack(side='left')
        self.save_pocket_volumes_trajectory = StringVar()
        self.__CheckSavePocketVolumesTrajectory = Checkbutton(self.__Frame46
            ,variable=self.save_pocket_volumes_trajectory, command=self.__on_CheckSavePocketVolumesTrajectory_ButRel_1)
        self.__CheckSavePocketVolumesTrajectory.pack(side='left')
        
        # Note, when the below was used, the checkbutton didn't function properly
        # on Mac OSX, so I used the command=callback technique above.
        #self.__CheckSavePocketVolumesTrajectory.bind('<ButtonRelease-1>' \
        #    ,self.__on_CheckSavePocketVolumesTrajectory_ButRel_1)
        
        self.__Frame50 = Frame(self.__Frame29)
        self.__Frame50.pack(pady=1,side='left')
        self.__Label22 = Label(self.__Frame50)
        self.__Label22.pack(padx=15,side='left')
        self.__Label23 = Label(self.__Frame50,state='disabled'
            ,text='Equal # of Points per Frame')
        self.__Label23.pack(side='left')
        self.output_equal_num_points_per_frame = StringVar()
        self.__CheckOutputEqualNumPointsPerFrame = Checkbutton(self.__Frame50
            ,state='disabled',variable=self.output_equal_num_points_per_frame)
        self.__CheckOutputEqualNumPointsPerFrame.pack(side='left')
        self.__Frame52 = Frame(self.__Frame45)
        self.__Frame52.pack(pady=1,side='left')
        self.save_tabbed_volume_file = StringVar()
        self.__CheckSaveTabbedVolumeFile = Checkbutton(self.__Frame52
            ,variable=self.save_tabbed_volume_file)
        self.__CheckSaveTabbedVolumeFile.pack(side='right')
        self.__Label24 = Label(self.__Frame52,text='Tabbed Volume File')
        self.__Label24.pack(side='right')
        self.__Frame51 = Frame(self.__Frame45)
        self.__Frame51.pack(pady=1,side='left')
        self.__Label28 = Label(self.__Frame51)
        self.__Label28.pack(padx=15,side='left')
        self.__Label30 = Label(self.__Frame51,text='Volumetric Density File')
        self.__Label30.pack(side='left')
        self.save_volumetric_density_map = StringVar()
        self.__CheckSaveVolumetricDensityMap = Checkbutton(self.__Frame51
            ,variable=self.save_volumetric_density_map)
        self.__CheckSaveVolumetricDensityMap.pack(side='left')
        self.__Frame44 = Frame(self.__Frame30)
        self.__Frame44.pack(pady=1,side='left')
        self.__Label18 = Label(self.__Frame44)
        self.__Label18.pack(padx=15,side='right')
        self.num_processors = StringVar()
        self.__Entry19 = Entry(self.__Frame44,textvariable=self.num_processors
            ,width=5)
        self.__Entry19.pack(side='right')
        self.__Label29 = Label(self.__Frame44,text='Num Processors:')
        self.__Label29.pack(side='right')
        self.__Frame53 = Frame(self.__Frame30)
        self.__Frame53.pack(side='left')
        self.__Label17 = Label(self.__Frame53,text='Disk Instead of Memory')
        self.__Label17.pack(side='left')
        self.use_disk = StringVar()
        self.__CheckbuttonUseDisk = Checkbutton(self.__Frame53
            ,variable=self.use_disk)
        self.__CheckbuttonUseDisk.pack(side='left')
        self.__Frame35 = Frame(self.__Frame30)
        self.__Frame35.pack(pady=1,side='left')
        self.__Label35 = Label(self.__Frame35)
        self.__Label35.pack(padx=15,side='left')
        self.__Label36 = Label(self.__Frame35,text='Python Executable:')
        self.__Label36.pack(side='left')
        self.python_exec = StringVar()
        self.__Entry1 = Entry(self.__Frame35,textvariable=self.python_exec
            ,width=10)
        self.__Entry1.pack(side='left')
        self.__Frame15 = Frame(self.__Frame13)
        self.__Frame15.pack(pady=1,side='left')
        self.__RunProgram = Button(self.__Frame15,text='Run Program!')
        self.__RunProgram.pack(side='top')
        self.__RunProgram.bind('<ButtonRelease-1>' \
            ,self.__on_RunProgram_ButRel_1)
        self.__Frame31 = Frame(self.__Frame13)
        self.__Frame31.pack(pady=1,side='left')
        self.__ButtonExit = Button(self.__Frame31,text='Exit')
        self.__ButtonExit.pack(side='top')
        self.__ButtonExit.bind('<ButtonRelease-1>' \
            ,self.__on_ButtonExit_ButRel_1)
        self.__Frame55 = Frame(self.__Frame13)
        self.__Frame55.pack(side='left')
        self.__ButtonHelp = Button(self.__Frame55,text='Help!')
        self.__ButtonHelp.pack(side='left')
        self.__ButtonHelp.bind('<ButtonRelease-1>' \
            ,self.__on_ButtonHelp_ButRel_1)
        self.__Frame34 = Frame(self.__Frame25,height=30,width=380)
        self.__Frame34.pack(side='top')
        self.__ContigLabel1 = Label(self.__Frame34,state='disabled')
        self.__ContigLabel1.pack(padx=40,side='left')
        self.__ContigLabel2 = Label(self.__Frame34,state='disabled'
            ,text='Center X/Y/Z:')
        self.__ContigLabel2.pack(side='left')
        self.center_x_2 = StringVar()
        self.__CenterX2 = Entry(self.__Frame34,state='disabled'
            ,textvariable=self.center_x_2,width=5)
        self.__CenterX2.pack(side='left')
        self.__ContigLabel3 = Label(self.__Frame34,state='disabled',text='/')
        self.__ContigLabel3.pack(side='left')
        self.center_y_2 = StringVar()
        self.__CenterY2 = Entry(self.__Frame34,state='disabled'
            ,textvariable=self.center_y_2,width=5)
        self.__CenterY2.pack(side='left')
        self.__ContigLabel4 = Label(self.__Frame34,state='disabled',text='/')
        self.__ContigLabel4.pack(side='left')
        self.center_z_2 = StringVar()
        self.__CenterZ2 = Entry(self.__Frame34,state='disabled'
            ,textvariable=self.center_z_2,width=5)
        self.__CenterZ2.pack(side='left')
        self.__Frame32 = Frame(self.__Frame25)
        self.__Frame32.pack(side='top')
        self.__Frame38 = Frame(self.__Frame32,height=30,width=100)
        self.__Frame38.pack(side='left')
        self.__ContigLabel5 = Label(self.__Frame38,state='disabled'
            ,text='Radius:')
        self.__ContigLabel5.pack(side='left')
        self.radius_2 = StringVar()
        self.__Radius2 = Entry(self.__Frame38,state='disabled'
            ,textvariable=self.radius_2,width=5)
        self.__Radius2.pack(side='left')
        self.__Frame40 = Frame(self.__Frame32,width=30)
        self.__Frame40.pack(side='left')
        self.__ContigLabel6 = Label(self.__Frame40,state='disabled',text='OR')
        self.__ContigLabel6.pack(side='top')
        self.__Frame39 = Frame(self.__Frame32,width=250)
        self.__Frame39.pack(side='left')
        self.__ContigLabel7 = Label(self.__Frame39,state='disabled')
        self.__ContigLabel7.pack(padx=2,side='left')
        self.__ContigLabel8 = Label(self.__Frame39,state='disabled'
            ,text='Length X/Y/Z:')
        self.__ContigLabel8.pack(side='left')
        self.length_x_2 = StringVar()
        self.__LengthX2 = Entry(self.__Frame39,state='disabled'
            ,textvariable=self.length_x_2,width=5)
        self.__LengthX2.pack(side='left')
        self.__ContigLabel9 = Label(self.__Frame39,state='disabled',text='/')
        self.__ContigLabel9.pack(side='left')
        self.length_y_2 = StringVar()
        self.__LengthY2 = Entry(self.__Frame39,state='disabled'
            ,textvariable=self.length_y_2,width=5)
        self.__LengthY2.pack(side='left')
        self.__ContigLabel10 = Label(self.__Frame39,state='disabled',text='/')
        self.__ContigLabel10.pack(side='left')
        self.length_z_2 = StringVar()
        self.__LengthZ2 = Entry(self.__Frame39,state='disabled'
            ,textvariable=self.length_z_2,width=5)
        self.__LengthZ2.pack(side='left')
        #
        #Your code here
        #
        
        self.__Frame11.pack_propagate(0)
        self.__Frame5.pack_propagate(0)
        self.__Frame9.pack_propagate(0)

        self.__Frame7.pack_propagate(0)
        self.__Frame14.pack_propagate(0)
        
        self.__Frame8.pack_propagate(0)
        self.__Frame10.pack_propagate(0)
        self.__Frame16.pack_propagate(0)
        self.__Frame26.pack_propagate(0)
        self.__Frame27.pack_propagate(0)
        self.__Frame34.pack_propagate(0)
        self.__Frame33.pack_propagate(0)
        
        self.pdb_filename.set("Select PDB File...")
        self.num_processors.set("1")
        self.point_spacing.set("1.0")
        self.padding.set("1.09")
        self.use_contig_points.set("0")
        self.contig_points_criteria.set("3")
        self.convex_hull.set("0")
        self.compress_points_files.set("0")
        self.use_disk.set("0")
        self.save_individual_pocket_volumes.set("0")
        self.save_pocket_volumes_trajectory.set("0")
        self.output_equal_num_points_per_frame.set("0")
        self.save_tabbed_volume_file.set("0")
        self.save_volumetric_density_map.set("0")
        self.points_files_prefix.set("." + os.sep)
        self.python_exec.set("python")
        self.point_field_only.set("0")

    #
    #Start of event handler methods
    #

    def __on_AddPrism1_ButRel_1(self,Event=None):
        x = self.center_x_1.get()
        y = self.center_y_1.get()
        z = self.center_z_1.get()
        lengthx = self.length_x_1.get()
        lengthy = self.length_y_1.get()
        lengthz = self.length_z_1.get()
        if not self.is_float(x) == True: 
            tkMessageBox.showerror("Error!", "The CenterX value must be a float.")
            return
        if not self.is_float(y) == True:
            tkMessageBox.showerror("Error!", "The CenterY value must be a float.")
            return
        if not self.is_float(z) == True:
            tkMessageBox.showerror("Error!", "The CenterZ value must be a float.")
            return        
        if not self.is_float(lengthx) == True:
            tkMessageBox.showerror("Error!", "The LengthX value must be a float.")
            return
        if not self.is_float(lengthy) == True:
            tkMessageBox.showerror("Error!", "The LengthY value must be a float.")
            return
        if not self.is_float(lengthz) == True:
            tkMessageBox.showerror("Error!", "The LengthZ value must be a float.")
            return
        self.__PointRegionsList.insert(END,"PointsInclusionBox " + x + " " + y + " " + z + " " + lengthx + " " + lengthy + " " + lengthz)

        self.center_x_1.set("")
        self.center_y_1.set("")
        self.center_z_1.set("")
        self.radius_1.set("")
        self.length_x_1.set("")
        self.length_y_1.set("")
        self.length_z_1.set("")

    def __on_AddPrism2_ButRel_1(self,Event=None):
        x = self.center_x_2.get()
        y = self.center_y_2.get()
        z = self.center_z_2.get()
        lengthx = self.length_x_2.get()
        lengthy = self.length_y_2.get()
        lengthz = self.length_z_2.get()
        if not self.is_float(x) == True: 
            tkMessageBox.showerror("Error!", "The CenterX value must be a float.")
            return
        if not self.is_float(y) == True:
            tkMessageBox.showerror("Error!", "The CenterY value must be a float.")
            return
        if not self.is_float(z) == True:
            tkMessageBox.showerror("Error!", "The CenterZ value must be a float.")
            return        
        if not self.is_float(lengthx) == True:
            tkMessageBox.showerror("Error!", "The LengthX value must be a float.")
            return
        if not self.is_float(lengthy) == True:
            tkMessageBox.showerror("Error!", "The LengthY value must be a float.")
            return
        if not self.is_float(lengthz) == True:
            tkMessageBox.showerror("Error!", "The LengthZ value must be a float.")
            return
        self.__ContiguousRegionsList.insert(END,"ContiguousPocketSeedBox " + x + " " + y + " " + z + " " + lengthx + " " + lengthy + " " + lengthz)
        self.center_x_2.set("")
        self.center_y_2.set("")
        self.center_z_2.set("")
        self.radius_2.set("")
        self.length_x_2.set("")
        self.length_y_2.set("")
        self.length_z_2.set("")

    def __on_AddSphere1_ButRel_1(self,Event=None):
        x = self.center_x_1.get()
        y = self.center_y_1.get()
        z = self.center_z_1.get()
        radius = self.radius_1.get()
        if not self.is_float(x) == True: 
            tkMessageBox.showerror("Error!", "The CenterX value must be a float.")
            return
        if not self.is_float(y) == True:
            tkMessageBox.showerror("Error!", "The CenterY value must be a float.")
            return
        if not self.is_float(z) == True:
            tkMessageBox.showerror("Error!", "The CenterZ value must be a float.")
            return        
        if not self.is_float(radius) == True:
            tkMessageBox.showerror("Error!", "The radius value must be a float.")
            return
        self.__PointRegionsList.insert(END,"PointsInclusionSphere " + x + " " + y + " " + z + " " + radius)
        self.center_x_1.set("")
        self.center_y_1.set("")
        self.center_z_1.set("")
        self.radius_1.set("")
        self.length_x_1.set("")
        self.length_y_1.set("")
        self.length_z_1.set("")   

    def __on_AddSphere2_ButRel_1(self,Event=None):
        x = self.center_x_2.get()
        y = self.center_y_2.get()
        z = self.center_z_2.get()
        radius = self.radius_2.get()
        if not self.is_float(x) == True: 
            tkMessageBox.showerror("Error!", "The CenterX value must be a float.")
            return
        if not self.is_float(y) == True:
            tkMessageBox.showerror("Error!", "The CenterY value must be a float.")
            return
        if not self.is_float(z) == True:
            tkMessageBox.showerror("Error!", "The CenterZ value must be a float.")
            return        
        if not self.is_float(radius) == True:
            tkMessageBox.showerror("Error!", "The radius value must be a float.")
            return
        self.__ContiguousRegionsList.insert(END,"ContiguousPocketSeedSphere " + x + " " + y + " " + z + " " + radius)
        self.center_x_2.set("")
        self.center_y_2.set("")
        self.center_z_2.set("")
        self.radius_2.set("")
        self.length_x_2.set("")
        self.length_y_2.set("")
        self.length_z_2.set("")

    def __on_ButtonExit_ButRel_1(self,Event=None):
        sys.exit(0)

    def __on_ButtonHelp_ButRel_1(self,Event=None):
        tkMessageBox.showerror("Get Help Online", "A helpful video can be found at: http://MY.VIDEO.HERE.COM")      

    def __on_CheckSaveIndividualPocketVolumes_ButRel_1(self,Event=None):
        self.set_output_equal_num_points_per_frame()

    def __on_CheckSavePocketVolumesTrajectory_ButRel_1(self,Event=None):
        self.set_output_equal_num_points_per_frame()

    def set_output_equal_num_points_per_frame(self):
       if self.save_individual_pocket_volumes.get() == "1" or self.save_pocket_volumes_trajectory.get() == "1": todo = 'normal'
       else: todo = 'disabled'
       
       self.__Label23['state'] = todo
       self.__CheckOutputEqualNumPointsPerFrame['state'] = todo        

    def __on_CheckUseContigPoints_ButRel_1(self,Event=None):
       try: 
            if self.use_contig_points.get() == "1": todo = 'normal'
            else: todo = 'disabled'
            
            self.__ContigLabel1['state'] = todo
            self.__ContigLabel2['state'] = todo
            self.__ContigLabel3['state'] = todo
            self.__ContigLabel4['state'] = todo
            self.__ContigLabel5['state'] = todo
            self.__ContigLabel6['state'] = todo
            self.__ContigLabel7['state'] = todo
            self.__ContigLabel8['state'] = todo
            self.__ContigLabel9['state'] = todo
            self.__ContigLabel10['state'] = todo
            self.__ContigLabel11['state'] = todo
            self.__CenterX2['state'] = todo
            self.__CenterY2['state'] = todo
            self.__CenterZ2['state'] = todo
            self.__Radius2['state'] = todo
            self.__LengthX2['state'] = todo
            self.__LengthY2['state'] = todo
            self.__LengthZ2['state'] = todo
            self.__AddSphere2['state'] = todo
            self.__AddPrism2['state'] = todo
            self.__DeleteContiguousRegion['state'] = todo
            self.__ContigPointsCriteria['state'] = todo
            self.__ContiguousRegionsList['state'] = todo
       except: pass

    def __on_DeleteContiguousRegion_ButRel_1(self,Event=None):
        try:
            self.__ContiguousRegionsList.delete(self.__ContiguousRegionsList.curselection())
        except:
            tkMessageBox.showerror("Error!", "Please select the item to delete from the list above.")    

    def __on_DeletePointRegion_ButRel_1(self,Event=None):
        try:
            self.__PointRegionsList.delete(self.__PointRegionsList.curselection())
        except:
             tkMessageBox.showerror("Error!", "Please select the item to delete from the list above.")   

    def __on_ExcludePrism1_ButRel_1(self,Event=None):
        x = self.center_x_1.get()
        y = self.center_y_1.get()
        z = self.center_z_1.get()
        lengthx = self.length_x_1.get()
        lengthy = self.length_y_1.get()
        lengthz = self.length_z_1.get()
        if not self.is_float(x) == True: 
            tkMessageBox.showerror("Error!", "The CenterX value must be a float.")
            return
        if not self.is_float(y) == True:
            tkMessageBox.showerror("Error!", "The CenterY value must be a float.")
            return
        if not self.is_float(z) == True:
            tkMessageBox.showerror("Error!", "The CenterZ value must be a float.")
            return        
        if not self.is_float(lengthx) == True:
            tkMessageBox.showerror("Error!", "The LengthX value must be a float.")
            return
        if not self.is_float(lengthy) == True:
            tkMessageBox.showerror("Error!", "The LengthY value must be a float.")
            return
        if not self.is_float(lengthz) == True:
            tkMessageBox.showerror("Error!", "The LengthZ value must be a float.")
            return
        self.__PointRegionsList.insert(END,"PointsExclusionBox " + x + " " + y + " " + z + " " + lengthx + " " + lengthy + " " + lengthz)

        self.center_x_1.set("")
        self.center_y_1.set("")
        self.center_z_1.set("")
        self.radius_1.set("")
        self.length_x_1.set("")
        self.length_y_1.set("")
        self.length_z_1.set("")   

    def __on_ExcludeSphere1_ButRel_1(self,Event=None):
        x = self.center_x_1.get()
        y = self.center_y_1.get()
        z = self.center_z_1.get()
        radius = self.radius_1.get()
        if not self.is_float(x) == True: 
            tkMessageBox.showerror("Error!", "The CenterX value must be a float.")
            return
        if not self.is_float(y) == True:
            tkMessageBox.showerror("Error!", "The CenterY value must be a float.")
            return
        if not self.is_float(z) == True:
            tkMessageBox.showerror("Error!", "The CenterZ value must be a float.")
            return        
        if not self.is_float(radius) == True:
            tkMessageBox.showerror("Error!", "The radius value must be a float.")
            return
        self.__PointRegionsList.insert(END,"PointsExclusionSphere " + x + " " + y + " " + z + " " + radius)
        self.center_x_1.set("")
        self.center_y_1.set("")
        self.center_z_1.set("")
        self.radius_1.set("")
        self.length_x_1.set("")
        self.length_y_1.set("")
        self.length_z_1.set("")     

    def __on_FindPDBFilename_ButRel_1(self,Event=None):
         afile = tkFileDialog.askopenfilename(parent=self.master,title='Choose a PDB File')
         if afile != "" and not afile == () and os.path.exists(afile):
             self.pdb_filename_full = afile
             self.pdb_filename.set(os.path.basename(afile))

    def __on_RunProgram_ButRel_1(self,Event=None):        
        # first, make sure all required paramters have been specified
        if not self.is_float(self.point_spacing.get()) == True: 
            tkMessageBox.showerror("Error!", "The Point Spacing value must be a number.")
            return
        if not self.is_float(self.padding.get()) == True: 
            tkMessageBox.showerror("Error!", "The Padding value must be a number.")
            return
        
        if self.point_field_only.get() == "0": # so the user hasn't asked for just a point-field PDB file.
            try:
                if not os.path.exists(self.pdb_filename_full):
                    tkMessageBox.showerror("Error!", "The specified PDB file does not exist.")
                    return
            except:
                    tkMessageBox.showerror("Error!", "The specified PDB file does not exist.")
                    return
                
            if self.use_contig_points.get() == "1" and not self.is_float(self.contig_points_criteria.get()):
                tkMessageBox.showerror("Error!", "If you choose to use contiguous points, the Contiguous Points Criteria must be a number.")
                return

            if not self.is_int(self.num_processors.get()):
                tkMessageBox.showerror("Error!", "The Number of Processors must be an number.")
                return

        # make any directories if required
        if os.sep in self.points_files_prefix.get():
            output_dirname = os.path.dirname(self.points_files_prefix.get())
            # if os.path.exists(output_dirname): shutil.rmtree(output_dirname) # So delete the directory if it already exists.
            try: os.mkdir(output_dirname)
            except: pass

        # create the input file
        f = open(self.points_files_prefix.get() + "POVME.in", 'w')
        
        f.write("GridSpacing " + self.point_spacing.get() + "\n")
        for item in list(self.__PointRegionsList.get(0, END)): f.write(item + "\n")
        f.write("DistanceCutoff " + self.padding.get() + "\n")
        f.write("OutputFilenamePrefix " + self.points_files_prefix.get() + "\n")

        if self.compress_points_files.get() == "1": f.write("CompressOutput true\n")           
        else: f.write("CompressOutput false\n")

        if self.point_field_only.get() == "0": # so the user hasn't asked for just a point-field PDB file.

            f.write("PDBFileName " + self.pdb_filename_full + "\n")

            if self.convex_hull.get() == "1": f.write("ConvexHullExclusion true\n")
            else: f.write("ConvexHullExclusion false\n")

            if self.use_contig_points.get() == "1":
                for item in list(self.__ContiguousRegionsList.get(0, END)): f.write(item + "\n")
                f.write("ContiguousPointsCriteria " + self.contig_points_criteria.get() + "\n")

            f.write("NumProcessors " + self.num_processors.get() + "\n")
        
            if self.use_disk.get() == "1": f.write("UseDiskNotMemory true\n")
            else: f.write("UseDiskNotMemory false\n")

            if self.save_individual_pocket_volumes.get() == "1": f.write("SaveIndividualPocketVolumes true\n")
            else: f.write("SaveIndividualPocketVolumes false\n")

            if self.save_pocket_volumes_trajectory.get() == "1": f.write("SavePocketVolumesTrajectory true\n")
            else: f.write("SavePocketVolumesTrajectory false\n")

            if self.output_equal_num_points_per_frame.get() == "1": f.write("OutputEqualNumPointsPerFrame true\n")
            else: f.write("OutputEqualNumPointsPerFrame false\n")

            if self.save_tabbed_volume_file.get() == "1": f.write("SaveTabbedVolumeFile true\n")
            else: f.write("SaveTabbedVolumeFile false\n")

            if self.save_volumetric_density_map.get() == "1": f.write("SaveVolumetricDensityMap true\n")
            else: f.write("SaveVolumetricDensityMap false\n")

        f.close()
        
        msg = 'The POVME calculation will now run. For your records, the output will be saved to "' + self.points_files_prefix.get() + '*".'        
        tkMessageBox.showinfo("Message", msg)
        
        if os.path.exists("POVME2.py"): # so you're already in the same directory as POVME2.py
            torun = self.python_exec.get() + " POVME2.py " + os.path.abspath(self.points_files_prefix.get())
        else:
            torun = self.python_exec.get() + " " + os.path.dirname(sys.argv[0]) + os.sep + "POVME2.py " + os.path.abspath(self.points_files_prefix.get())

        if self.points_files_prefix.get()[-1] == os.sep: torun = torun + os.sep
        torun = torun + "POVME.in" # > " + self.points_files_prefix.get() + "POVME.out" ### No need to direct output to a file. POVME does this automatically.
        os.system(torun)
        
        tkMessageBox.showinfo("Message", "The POVME calculation is done!")
        
    #
    #Start of non-Rapyd user code
    #

    def is_float(self, str):
        try:
            float(str)
            return True
        except ValueError:
            return False
            
    def is_int(self, str):
        try:
            int(str)
            return True
        except ValueError:
            return False

import tkMessageBox
import tkFileDialog
import os
#import shutil
import sys

if __name__ == '__main__':

    Root = Tk()
    App = POVME2_GUI(Root)
    App.pack(expand='yes',fill='both')

    Root.geometry('680x725+10+10')
    Root.resizable(0,0)
    Root.title('POVME 2.0 GUI')
    Root.mainloop()
