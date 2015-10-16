#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
Scripting to handle DNA quantification, normalisation and nextera processes for team214's single cell genomics pipeline.

-------------------------------------------------------------------------------
Team 214 Single Cell Genomics Core Facility
Author : Andrew Sparkes
Date   : Oct 2015
-------------------------------------------------------------------------------
Versions:
Initial version 								     Andrew Sparkes    Oct 2015

===============================================================================
Modes

-------------------------------------------------------------------------------
quantsetup - Quantification set up process
-------------------------------------------------------------------------------
GUI-based mode to handle the creation of the rundef and echo csv files required 
for the quantification process.

Inputs:

1. LIMS or manually created file holding details about the selected grouping of
   plates to be quantified (1->n plates), containing a grouping or job id, and 
   the following for each plate: barcode, sample type, concentration thresholds,
   and layout (location of Samples and relevant control wells).

2. Standards file that holds the layout of the relevant standards plate (which 
   file depends on the source plate sample type e.g. cDNA, gDNA etc.)

3. Template RunDef file into which to insert run-specific parameters, some of
   which depend upon the numbers of plates.

Outputs:

1. Quantification RunDef file for the Labcyte Tempo software.

2. Echo transfer csv files specific to source DNA plates and standards plate.

3. Experiment run log to detail what was done.

4. Input LIMS file moved to experiment run directory to detail what was 
   requested.

-------------------------------------------------------------------------------
quantcalc - Quantification calculations
-------------------------------------------------------------------------------
GUI-based mode to handle the calculation of quantification results from the BMG
files.
[Auto-timeout feature so run automatically continues?]

Inputs:

Outputs:

===============================================================================





Standards setup file
--------------------
Local standards layout file containing the following details:
*  Version number of this file
*  
*  DNA Standards ladder locations and concentrations by sample type (N.B. may be one set for cDNA, another for gDNA 
   because may expect different sample concentration ranges for these types).
   i.e. lists of locations and concentrations:
   gDNA
   		A23, 500 ng/ul
   		B23, 300 ng/ul
   		...
   		H23, 0 ng/ul
   cDNA
   		A23, 200 ng/ul
   		B23, 150 ng/ul
   		...
   		H23, 0 ng/ul
*  Starting locations where DNA plate pools should be put (max 16 sources?)
*  Number of pools to make e.g. 3






'''

import argparse # to parse command line arguments
import sys # for sys.exit
import os # for file directory selection
import json # for reading json files
import inspect # for getting method name

from configobj 		import ConfigObj 		# for reading config files
from Tkinter        import * 				# for the GUI interfaces
import ttk 									# for the GUI interface widgets
from tkMessageBox 	import askyesno 		# for pop-up message boxes
from tkFileDialog 	import askopenfilename 	# for file selection


from pprint import pprint # for pretty printing e.g. lists and dictionaries

# -----------------------------------------------------------------------------
# Variables
# -----------------------------------------------------------------------------
config_filename 			= 'config/access_system.cfg' # configuration filename
valid_modes     			= ['quantsetup'] # valid program modes
args     					= {} # stores parsed command line arguments
settings 					= {} # stores parsed configuration settings
valid_quant_standards  		= ['SS2'] # list of the valid standards types
quant_standards				= {} # stores parsed quantification standards

# GUI 
bg_colour 					= 'white'
fg_colour 					= 'black'
entry_field_colour  		= 'yellow'
title_colour  				= 'blue'
fg_ok_colour 				= '#32701e' # green
fg_error_colour 			= 'red'

font_arial_normal 			= 'Arial 16'
font_courier_normal			= 'Courier 12'
font_arial_large_bold 		= 'Arial 24 bold'

# -----------------------------------------------------------------------------
# Initialisation Methods
# -----------------------------------------------------------------------------
def parse_command_line_arguments():
	'''Parse the command line arguments.

	Parses in the various command line arguments.
	Most critical is the mode argument that tells us which functionality to use for this call to the program.
	'''
	parser = argparse.ArgumentParser()
	parser.add_argument("--debug", 			help="Debug mode for development, default is off.", 	action="store_true")
	parser.add_argument("-v", "--verbose", 	help="increase log output verbosity", 					action="store_true")
	parser.add_argument("-m", "--mode", 	help="Mode to run script in i.e. quantsetup or quant.", type=str)
	global args
	args = parser.parse_args()
	
	if(args.debug == True):
		print("DEBUG mode is active")		
		if(args.verbose == True):
			print("Verbose mode is active")

	# check if mode is valid
	if args.mode in valid_modes:
		if(args.debug == True):
			print("Mode is valid: " + args.mode)
	else:
		sys.exit("Access System Script: ERROR: Chosen mode is not recognised <" + str(args.mode)+ ">, cannot continue")

	return

def parse_access_system_config_file():
	'''Parse the settings from the access system configuration file.

	This configuration file contains directories and filepaths for the various input and outputs.
	The configuration file is presumed to be called access_system.cfg and to be found in the same directory as the script.

	Contains the following details:
	*  Tempo inbox directory location
	*  ECHO ECP files directory (or full filepaths for each ECP??)
	*  Experiment runs root directory (create run-specific directory in here once LIMS file validated and user Ok 
	   to hold the ECHO csv files, rundef file copy and log file)
	*  LIMS file network directory (copy from here to temp directory for validation, and if accepted paste into 
	   the experiment runs directory???)
	*  Transfer volumes and other fixed parameters
	*  Flag for debug mode
	*  Standards layout file filepath
	*  Maximum number of source plates allowed (check in validation of LIMS file received)
	'''

	if(args.debug == True):
		print("In %s" % inspect.currentframe().f_code.co_name)

	# Read the configuration file
	config = ConfigObj(config_filename)

	# Verify that the options imported from the file match with an expected list maintained here
	options_list = {
		'Version':{'version_number':'flt'},
		'Common':{'tempo_inbox_dir':'str', 
				'ecp_armadillo_filepath':'str',
				'ecp_greiner_black_filepath':'str',
				'ecp_third_type_filepath':'str',
				'expt_run_root_dir':'str',
				'gui_width':'int',
				'gui_height':'int',
				'gui_x_posn':'int',
				'gui_y_posn':'int'},
		'Quantification':{'quant_max_source_plates_allowed':'int', 
				'quant_lims_file_network_dir':'str', 
				'quant_standards_dir':'str',
				'quant_filename_standards_ss2':'str',
				'quant_cdna_min_conc':'flt',
				'quant_cdna_max_conc':'flt',
				'quant_cdna_min_percent_wells_ok':'flt',
				'quant_gdna_min_conc':'flt',
				'quant_gdna_max_conc':'flt',
				'quant_gdna_min_percent_wells_ok':'flt'
		}
	}

	global settings

	# copy the options into the global settings list so that they are available throughout the program
	for opt_group, opt_dict in options_list.iteritems():
		settings[opt_group] = {}
		for opt, opt_type in opt_dict.iteritems():
			if(config.has_key(opt_group)):
				if(config[opt_group].has_key(opt)):
					if(opt_type == 'str'):
						settings[opt_group][opt] = str(config[opt_group][opt])
					elif(opt_type == 'int'):
						settings[opt_group][opt] = int(config[opt_group][opt])
					elif(opt_type == 'flt'):
						settings[opt_group][opt] = float(config[opt_group][opt])
					elif(opt_type == 'bool'):
						settings[opt_group][opt] = bool(config[opt_group][opt])
						print(".")
					else:
						sys.exit("Access System Script: ERROR: Configuration file field type not understood for option group <" + str(opt_group)+ "> and option <" + str(opt) + ">, cannot continue")
				else:
					sys.exit("Access System Script: ERROR: Configuration file format option group <" + str(opt_group)+ "> is missing option <" + str(opt) + ">, cannot continue")
			else:
				sys.exit("Access System Script: ERROR: Configuration file format missing option group <" + str(opt_group)+ ">, cannot continue")

	# print out settings if in debug mode
	if(args.debug == True):
		print('-' * 80)
		print("Settings from config file are as follows:")
		print('-' * 80)
		pprint(settings)
		print('-' * 80)

	return


def parse_quantification_standards_config_file(stnd_type):
	'''Parse the settings from the selected standards configuration file.

	There will be a version-controlled standards file for each type of standards plate.
	These standards files are located in a directory specified in the main config file,
	and have filenames specified in the main config file.
	'''

	if(args.debug == True):
		print("In %s" % inspect.currentframe().f_code.co_name)
		print("Input standards type: %s" % stnd_type)

	# Read the relevant configuration file
	standards_config_filename = ""
	if(stnd_type == 'SS2'):
		standards_config_filename = settings.get('Quantification').get('quant_filename_standards_ss2')

	if(args.debug == True):
		print("Standards config filename : %s" % standards_config_filename)	

	standards_config_filepath = os.path.join(settings.get('Quantification').get('quant_standards_dir'), standards_config_filename)

	if(args.debug == True):
		print("Standards config filepath : %s" % standards_config_filepath)

	config = ConfigObj(standards_config_filepath)

	# Verify that the options imported from the file match with an expected list maintained here
	options_list = {
		'Version':{'version_number':'flt'},
		'Information':{'kit_name':'str',
				'manufacturer':'str',
				'order_number':'str',
				'volume_of_stds_to_transfer_nl':'flt'},
		'Ladder':{'number_of_ladder_wells':'int',
				'sheared_size_of_dna_kb':'flt'},
		'Pools':{'number_of_pools':'int',
				'max_num_sources':'int'}
	}

	global quant_standards
	quant_standards[stnd_type] = {}

	# copy the standards into the global standards list so that they are available throughout the program
	for opt_group, opt_dict in options_list.iteritems():
		quant_standards[stnd_type][opt_group] = {}
		for opt, opt_type in opt_dict.iteritems():
			if(config.has_key(opt_group)):
				if(config[opt_group].has_key(opt)):
					if(opt_type == 'str'):
						quant_standards[stnd_type][opt_group][opt] = str(config[opt_group][opt])
					elif(opt_type == 'int'):
						quant_standards[stnd_type][opt_group][opt] = int(config[opt_group][opt])
					elif(opt_type == 'flt'):
						quant_standards[stnd_type][opt_group][opt] = float(config[opt_group][opt])
					elif(opt_type == 'bool'):
						quant_standards[stnd_type][opt_group][opt] = bool(config[opt_group][opt])
						print(".")
					else:
						sys.exit("Access System Script: ERROR: Standards file field type not understood for option group <" + str(opt_group)+ "> and option <" + str(opt) + ">, cannot continue")
				else:
					sys.exit("Access System Script: ERROR: Standards file format option group <" + str(opt_group)+ "> is missing option <" + str(opt) + ">, cannot continue")
			else:
				sys.exit("Access System Script: ERROR: Standards file format missing option group <" + str(opt_group)+ ">, cannot continue")

	# retrieve sub-section information for the ladder
	lad_idx 			= 1
	num_ladder_wells 	= int(quant_standards[stnd_type]['Ladder']['number_of_ladder_wells'])

	quant_standards[stnd_type]['Ladder']['wells'] = {}

	while lad_idx <= num_ladder_wells:
		s_lad_idx 		= str(lad_idx)
		quant_standards[stnd_type]['Ladder']['wells'][s_lad_idx] 						= {}
		quant_standards[stnd_type]['Ladder']['wells'][s_lad_idx]['well_posn'] 			= str(config['Ladder'][s_lad_idx]['well_posn'])
		quant_standards[stnd_type]['Ladder']['wells'][s_lad_idx]['concentration_ng_ul']	= float(config['Ladder'][s_lad_idx]['concentration_ng_ul'])
		quant_standards[stnd_type]['Ladder']['wells'][s_lad_idx]['vol_to_dispense_nl'] 	= float(config['Ladder'][s_lad_idx]['vol_to_dispense_nl'])
		lad_idx 		+= 1

	# retrieve sub-section information for the pools
	src_idx 			= 1
	max_num_sources 	= int(quant_standards[stnd_type]['Pools']['max_num_sources'])

	quant_standards[stnd_type]['Pools']['sources'] = {}

	while src_idx <= max_num_sources:
		s_src_idx 		= str(src_idx)
		print("index : %s" % s_src_idx)
		pprint(config['Pools'][s_src_idx])
		quant_standards[stnd_type]['Pools']['sources'][s_src_idx] 						= {}
		quant_standards[stnd_type]['Pools']['sources'][s_src_idx]['initial_pool_posn'] 	= str(config['Pools'][s_src_idx]['initial_pool_posn'])
		src_idx 		+= 1

	# print out quant_standards if in debug mode
	if(args.debug == True):
		print('-' * 80)
		print("Standards are as follows:")
		print('-' * 80)
		pprint(quant_standards)
		print('-' * 80)

	return

def read_library_prep_config_file():
	'''Parse the settings from the library preparation configuration file.


	'''

	if(args.debug == True):
		print("In %s" % inspect.currentframe().f_code.co_name)

	return

# -----------------------------------------------------------------------------
# Processing Methods
# -----------------------------------------------------------------------------
def process_quantsetup():
	'''Entry method for processing the quantsetup mode'''

	if(args.debug == True):
		print("In %s" % inspect.currentframe().f_code.co_name)

	# build the user interface to ask the user to the select files required 
	display_gui_quantsetup()

	return

def process_quantcalc():
	'''Entry method for processing the quantcalc mode'''

	if(args.debug == True):
		print("In %s" % inspect.currentframe().f_code.co_name)

	return

# -----------------------------------------------------------------------------
# GUI Screen Classes
# -----------------------------------------------------------------------------
class QuantSetupGUI:
	'''GUI class for quantification setup'''

	# -----------------------------------------------------------------------------
	# Variables
	# -----------------------------------------------------------------------------
	selected_lims_file 		= "" # holds selected lims file filepath
	data_from_lims_file 	= {} # holds selected file data once read
	data_summary 			= {} # holds summarised data

	def __init__(self, master):
		'''Initialise a new frame for the GUI'''

		if(args.debug == True):
			print("In QuantSetupGUI.%s" % inspect.currentframe().f_code.co_name)

		# setup any common styles and themes for the screens
		setup_styles_and_themes()

		# containing frame for the contents
		frame 					= Frame(master, bg = "", colormap = "new")
		frame.pack(fill = BOTH, expand = 1) # expand the frame to fill the root window

		# give columns equal weighting
		col_num = 0
		while (col_num < 4):
			frame.grid_columnconfigure(col_num, weight=1)
			col_num += 1

		# define the GUI widget elements by row
		# ---------------------------------------------------------------------
		# Row 0 - Team name and Sanger logo
		# ---------------------------------------------------------------------
		frame.grid_rowconfigure(0, weight=1)
		lbl_scgcf_params 		= {'widg_text':"Single Cell Genomics Core Facility", 
									'widg_fg':'#18631e',
									'widg_font':'Verdana 12 italic',
									'grid_row':0,
									'grid_col':0,
									'grid_cols':2,
									'grid_sticky':NW}
		self.lbl_scgcf 			= create_widget_label(frame, lbl_scgcf_params)

		self.logo_image			= PhotoImage(file='/Users/as28/Documents/PYTHON/team214/images/sanger_logo.ppm')
		lbl_logo_params			= {'widg_image':self.logo_image,
									'widg_anchor':CENTER,
									'grid_row':0,
									'grid_col':2,
									'grid_cols':2,
									'grid_sticky':NE}
		self.lbl_logo 			= create_widget_label(frame, lbl_logo_params)

		# ---------------------------------------------------------------------
		# Row 1 - Title label
		# ---------------------------------------------------------------------
		frame.grid_rowconfigure(1, weight=1)
		lbl_title_params 		= {'widg_text':"Quantification Setup", 
									'widg_fg':title_colour,
									'widg_font':font_arial_large_bold,
									'widg_justify':CENTER,
									'widg_anchor':CENTER,
									'grid_row':1,
									'grid_col':0,
									'grid_cols':4,
									'grid_sticky':W+E}
		self.lbl_title 			= create_widget_label(frame, lbl_title_params)

		# ---------------------------------------------------------------------
		# Row 2 - Instructions for user text area
		# ---------------------------------------------------------------------
		frame.grid_rowconfigure(2, weight=1)
		self.instr_text  		= "Click 'Select File' to choose the Quantification plate grouping file downloaded " \
		"from the LIMS (or the equivalent manually created file). Check the summary of the information from the file, "\
		"confirm how many black plates are loaded into stack 4 and then press 'Create Access Files' to create the RunDef "\
		"and ECHO csv files."
		txt_instr_params 		= {'widg_text':self.instr_text, 
									'widg_height':3,
									'widg_state':DISABLED,
									'grid_row':2,
									'grid_col':0,
									'grid_cols':4,
									'grid_sticky':NW+NE}
		self.txt_instr 			= create_widget_text(frame, txt_instr_params)

		# ---------------------------------------------------------------------
		# Row 3 - LIMS filepath and button
		# ---------------------------------------------------------------------
		frame.grid_rowconfigure(3, weight=1)
		lbl_lims_fp_params 		= {'widg_bg':'yellow',
									'widg_font':font_courier_normal,
									'widg_width':90,
									'grid_row':3,
									'grid_col':0,
									'grid_cols':3,
									'grid_sticky':N+W+E+S,
									'grid_has_border':True}
		self.lbl_lims_fp 		= create_widget_label(frame, lbl_lims_fp_params)

		btn_sel_lims_file_params = {'widg_text':"Select LIMS File",
									'widg_width':16,
									'widg_command':self.lims_file_open_callback,
									'grid_row':3,
									'grid_col':3,
									'grid_cols':1,
									'grid_sticky':E}
		self.btn_sel_lims_file 	= create_widget_button(frame, btn_sel_lims_file_params)

		# ---------------------------------------------------------------------
		# Row 4 - Summary of plates text area
		# ---------------------------------------------------------------------
		frame.grid_rowconfigure(4, weight=1)
		txt_summary_params		= {'widg_height':22,
									'grid_row':4,
									'grid_col':0,
									'grid_cols':4,
									'grid_sticky':N+W+E+S,
									'grid_has_border':True}
		self.txt_summary 		= create_widget_text(frame, txt_summary_params)

		# ---------------------------------------------------------------------
		# Row 5 - Black plates required labels
		# ---------------------------------------------------------------------
		#TODO:
		# Validate num blk plates on pressing create files button.

		frame.grid_rowconfigure(5, weight=1)
		lbl_blk_plts_req_params = {'widg_text':"Number of black plates needed (incl. for standards) : ", 
									'grid_row':5,
									'grid_col':2,
									'grid_cols':1,
									'grid_sticky':NW+NE}
		self.lbl_blk_plts_req 	= create_widget_label(frame, lbl_blk_plts_req_params)

		self.var_num_blk_plts_reqd 	= IntVar(frame)
		self.var_num_blk_plts_reqd.set(0)
		lbl_num_blk_plts_reqd_params = {'widg_txt_var':self.var_num_blk_plts_reqd,
									'grid_row':5,
									'grid_col':3,
									'grid_cols':1,
									'grid_sticky':NW+NE}
		self.lbl_num_blk_plts_reqd 	= create_widget_label(frame, lbl_num_blk_plts_reqd_params)

		# ---------------------------------------------------------------------
		# Row 6 - Actual black plates in system label and combobox
		# ---------------------------------------------------------------------
		frame.grid_rowconfigure(6, weight=1)
		lbl_num_blk_plts_deck_params = {'widg_text':"Actual number of black plates loaded in stack four : ", 
									'grid_row':6,
									'grid_col':2,
									'grid_cols':1,
									'grid_sticky':NW+NE}
		self.lbl_num_blk_plts_deck 	= create_widget_label(frame, lbl_num_blk_plts_deck_params)


		self.lst_stk_plate_count	= list(range(0,41)) # create list of numbers
		self.num_blk_plates_deck 	= IntVar(frame) # integer variable to hold selected value
		self.num_blk_plates_deck.set(0) # initial value set at zero to make user check stack

		cmbbx_num_blk_plts_params 	= {'widg_values':self.lst_stk_plate_count,
									'widg_state':'readonly',
									'widg_width':3,
									'widg_txt_var':self.num_blk_plates_deck,
									'grid_row':6,
									'grid_col':3,
									'grid_cols':1,
									'grid_sticky':NW+NE}

		self.cmbbx_num_blk_plts = create_widget_combobox(frame, cmbbx_num_blk_plts_params)

		# ---------------------------------------------------------------------
		# Row 7 - Message panel
		# ---------------------------------------------------------------------
		frame.grid_rowconfigure(7, weight=1)
		txt_msg_panel_params		= {'widg_height':3,
									'grid_row':7,
									'grid_col':0,
									'grid_cols':4,
									'grid_sticky':N+W+E+S,
									'grid_has_border':True}
		self.txt_msg_panel 		= create_widget_text(frame, txt_msg_panel_params)

		# ---------------------------------------------------------------------
		# Row 8 - Quit and Create Access Files buttons
		# ---------------------------------------------------------------------
		frame.grid_rowconfigure(8, weight=1)
		btn_quit_params 		= {'widg_text':"Quit",
									'widg_width':12,
									'widg_command':self.quit_button_callback,
									'grid_row':8,
									'grid_col':0,
									'grid_cols':1,
									'grid_sticky':SW}
		self.btn_quit 			= create_widget_button(frame, btn_quit_params)

		btn_create_files_params = {'widg_text':"Create Access Files",
									'widg_width':16,
									'widg_command':self.create_access_files_button_callback,
									'widg_state':DISABLED,
									'grid_row':8,
									'grid_col':3,
									'grid_cols':1,
									'grid_sticky':SE}
		self.btn_create_files 	= create_widget_button(frame, btn_create_files_params)


	def lims_file_open_callback(self):
		'''Triggered when user has pressed the Select LIMS file button'''

		if(args.debug == True):
			print("In QuantSetupGUI.%s" % inspect.currentframe().f_code.co_name)

		# clear error msg line in screen, summary text and disable create files button
		self.clear_screen()

		# normpath fixes path for different OSs
		limsdir = os.path.normpath(settings.get('Quantification').get('quant_lims_file_network_dir'))

		if(args.debug == True):
			print("LIMS file directory: %s" % limsdir)

		# open file dialog and save selection N.B. if path not recognised it opens anyway
		self.selected_lims_file = askopenfilename(title  = "Select LIMS plate grouping file", initialdir = limsdir)

		# update the lable to display the filepath
		self.lbl_lims_fp.configure(text = self.selected_lims_file)

		if((len(self.selected_lims_file) > 0) and (os.path.isfile(self.selected_lims_file))):
			if(args.debug == True):
				print("File chosen: %s" % self.selected_lims_file)

			# read file into memory
			self.read_json_file()

			# validate data and create summary data
			if(self.validate_data_from_lims_file() == True):

				# at this point we need to know which standards we are using e.g. SS2 and load the relevant file
				parse_quantification_standards_config_file(self.data_summary["standards_type"])

				# TODO: need validation check to make sure standards layout can handle number of plates in this run

				# display a summary of the data in the file (one line per plate) and ask user to confirm generation of rundef and csv files
				self.display_summary_of_plates()

				# TODO: we should also be able to calculate and summarise amount of reagents reqd
				# maybe have links to display what the standards plate looks like, incl where pools will go

				# Enable the Create Files Button
				self.btn_create_files.config(state = NORMAL)
		else:
			if(args.debug == True):
				print("No file was selected!")

			# red error msg line in screen
			self.display_message(True, "No file was selected, please try again.")		

		return

	def create_access_files_button_callback(self):
		'''Triggered when user has pressed the Check File button'''

		if(args.debug == True):
			print("In QuantSetupGUI.%s" % inspect.currentframe().f_code.co_name)

		# Clear the message panel
		self.txt_msg_panel.delete('1.0', END)
		
		# validate that there are enough black plates in the stack
		num_blk_deck = self.num_blk_plates_deck.get()
		num_blk_reqd = self.var_num_blk_plts_reqd.get()

		if(args.debug == True):
			print("Num plates required = %s" % str(num_blk_reqd))
			print("Num plates in stack = %s" % str(num_blk_deck))

		if(num_blk_reqd == 0):
			# error should be > 0
			self.display_message(True, "ERROR: Number of black plates required is zero, should be 1 or more.")
			return
		elif(num_blk_deck == 0):
			# error user should set to num in stack
			self.display_message(True, "ERROR: Please load sufficient black plates (%s or more required) into stack four "\
				"and use the dropdown entry field to set how many black plates there are now in the stack." % str(num_blk_reqd))
			return
		elif(num_blk_deck < num_blk_reqd):
			# error user needs to load at least reqd plates to stack
			self.display_message(True, "ERROR: Insufficient black plates (%s or more required), please add more and use the "\
				"dropdown entry field to set how many black plates there are now in the stack." % str(num_blk_reqd))
			return
		if askyesno('Verify', 'This will generate the Access System RunDef and Echo files. Are you sure?'):
			self.display_message(False, "Starting file creation process, please wait...")


		# TODO:
		# create expt dir using RUN_ID
		# generate rundef and Echo csvs in expt dir
		# copy rundef to INBOX dir
		# move LIMS file to expt dir

		return

	def display_message(self, is_error, message):
		'''Display message on the screen'''

		if(args.debug == True):
			print("In QuantSetupGUI.%s" % inspect.currentframe().f_code.co_name)
			if(is_error):
				print("Input is Error message")
			else:
				print("Input is Standard message")
			print("Input message: %s" % message)			

		# delete current message content
		self.txt_msg_panel.delete('1.0', END)

		# set up tags
		self.txt_msg_panel.tag_configure('msg_error', 	 font=font_arial_normal, foreground=fg_error_colour)
		self.txt_msg_panel.tag_configure('msg_standard', font=font_arial_normal, foreground=fg_ok_colour)

		# display message
		if(is_error):
			print("Error message : %s" % message)
			self.txt_msg_panel.insert('1.0', message, ('msg_error'))
		else:
			print("Standard message : %s" % message)
			self.txt_msg_panel.insert('1.0', message, ('msg_standard'))

		return

	def quit_button_callback(self):
		'''Triggered when user has pressed the Quit button'''

		if(args.debug == True):
			print("In QuantSetupGUI.%s" % inspect.currentframe().f_code.co_name)

		if askyesno('Verify', 'Are you sure you want to Quit?'):
			sys.exit()

		return

	def read_json_file(self):
		'''Reads the data from a json file'''

		if(args.debug == True):
			print("In QuantSetupGUI.%s" % inspect.currentframe().f_code.co_name)

		# read the file N.B. the with clause closes the file as soon as we are out of the with context.
		# don't print or do anything in the with open clause as that keeps the file open longer.
		with open(self.selected_lims_file) as data_file:    
			self.data_from_lims_file = json.load(data_file)
		
		if(args.debug == True):
			print("LIMS file data:")
			print('-' * 80)
			pprint(self.data_from_lims_file)
			print('-' * 80)

		return

	def validate_data_from_lims_file(self):
		'''Validates the data extracted from the LIMS json file

		This file is created either by the LIMS or manually. It holds the details about the grouping of plates
		about to be processed on the Labcyte Access system. This information is used to generate the Access RunDef
		and Echo transfer csv files required for the quantification process.

		Method returns a boolean indicating validation success or failure.
		'''

		if(args.debug == True):
			print("In QuantSetupGUI.%s" % inspect.currentframe().f_code.co_name)

		# validation check: there should be certain fields in the json file
		expected_fields = ['RUN_ID', 'PLATES']
		for expected_field in expected_fields:
			if self.data_from_lims_file.has_key(expected_field):
				if(args.debug == True):
					print("Field located : %s" % expected_field)
			else:
				self.display_message(True, "ERROR: Key field <%s> missing from this file. Cannot continue." % expected_field)
				return False

		self.data_summary["run_id"] 	= self.data_from_lims_file["RUN_ID"]
		self.data_summary["num_plates"] = len(self.data_from_lims_file["PLATES"])
		self.data_summary["plates"] 	= {}

		if(args.debug == True):
			print("Run ID = %s" % self.data_summary["run_id"])
			print("Num plates found = %s" % self.data_summary["num_plates"])

		# validation check: there should be at least one plate
		if(self.data_summary["num_plates"] == 0):
			self.display_message(True, "ERROR: Unable to identify any plates in this file. Cannot continue.")
			return False

		# set black plates required (n + 1 for standards intermediate plate)
		self.num_black_plates_reqd = self.data_summary["num_plates"] + 1

		# cycle through the plates in the file and extract key information
		plate_index 					= 1
		num_standards_types 			= 0
		standards_types                 = {}

		while plate_index <= self.data_summary["num_plates"]:
			s_plt_idx 					= str(plate_index)
			
			# store plate summary details
			self.data_summary["plates"][s_plt_idx] 							= {}
			self.data_summary["plates"][s_plt_idx]["barcode"] 				= self.data_from_lims_file["PLATES"][s_plt_idx]["BARCODE"]
			self.data_summary["plates"][s_plt_idx]["library_prep_params"] 	= self.data_from_lims_file["PLATES"][s_plt_idx]["LIBRARY_PREP_PARAMS"]
			self.data_summary["plates"][s_plt_idx]["standards_params"] 		= self.data_from_lims_file["PLATES"][s_plt_idx]["STANDARDS_PARAMS"]

			if self.data_from_lims_file["PLATES"][s_plt_idx]["STANDARDS_PARAMS"] not in valid_quant_standards:
				self.display_message(True, "ERROR: plate with barcode <%s> has a standards type of <%s> which is not currently supported. Cannot continue." % (self.data_from_lims_file["PLATES"][s_plt_idx]["BARCODE"], self.data_from_lims_file["PLATES"][s_plt_idx]["STANDARDS_PARAMS"]))
				return False

			# count number of types of standards plates, used to check not got a mixed set
			if(self.data_from_lims_file["PLATES"][s_plt_idx]["STANDARDS_PARAMS"] in standards_types):
				standards_types[self.data_from_lims_file["PLATES"][s_plt_idx]["STANDARDS_PARAMS"]] 	+= 1
			else:
				standards_types[self.data_from_lims_file["PLATES"][s_plt_idx]["STANDARDS_PARAMS"]]	= 1

			curr_wells 					= self.data_from_lims_file["PLATES"][s_plt_idx]["WELLS"]

			# count numbers of samples and controls in this plate
			count_sample_wells 			= 0
			count_control_wells 		= 0
			well_index 					= 0

			while well_index < len(curr_wells):
				if(curr_wells[well_index]["ROLE"] == "SAMPLE"):
					count_sample_wells 	+= 1

				if(curr_wells[well_index]["ROLE"] == "CONTROL"):
					count_control_wells += 1

				well_index 				+= 1

			self.data_summary["plates"][s_plt_idx]["count_sample_wells"] 	= count_sample_wells
			self.data_summary["plates"][s_plt_idx]["count_control_wells"] 	= count_control_wells

			plate_index 				+= 1

		# validation check: there is only one type of standards plate required, if more than one it's an error
		if((len(standards_types) == 0) or (len(standards_types) > 1)):
			self.display_message(True, "ERROR: this group of plates contains %s types of standards plate requirements, which is not currently supported. Cannot continue." % len(standards_types))
			return False
		else:
			# fetch first and only key and store as the standards type for the group of plates
			self.data_summary["standards_type"] = standards_types.keys()[0]

		self.display_message(False, "LIMS file successfully validated. Plates found: %s" % str(self.data_summary["num_plates"]))

		if(args.debug == True):
			print(pprint(self.data_summary))

		return True

	def display_summary_of_plates(self):
		'''Displays a summary of the data extracted from the LIMS json file.

		Headings plus one summary line with key information per plate.
		'''

		if(args.debug == True):
			print("In QuantSetupGUI.%s" % inspect.currentframe().f_code.co_name)

		# define tags for formatting text
		self.txt_summary.tag_configure('tag_title', font='arial 18 bold', foreground=fg_ok_colour, relief='raised',justify='center', underline='True')
		self.txt_summary.tag_configure('tag_highlight', font='arial 16 bold', foreground=fg_ok_colour)

		# insert text for summary
		self.txt_summary.insert('1.0', "Summary of File\n", ('tag_title'))
		
		self.txt_summary.insert(END, "Run ID \t\t: ")
		self.txt_summary.insert(END, "%s\n" % self.data_summary.get('run_id'), ('tag_highlight'))
		
		self.txt_summary.insert(END, "Standards set \t\t: ")
		std_version_num = str(quant_standards[self.data_summary.get('standards_type')]['Version']['version_number'])
		self.txt_summary.insert(END, "%s  (version: %s)\n" % (self.data_summary.get('standards_type'), std_version_num), ('tag_highlight'))
		
		self.txt_summary.insert(END, "Library prep set \t\t: ")
		# lib_prep_version_num = str(lib_prep_standards[self.data_summary.get('lib_prep_type')]['Version']['version_number'])
		self.txt_summary.insert(END, "%s  (version: %s)\n" % ("????", "???"), ('tag_highlight'))
		# self.txt_summary.insert(END, "%s  (version: %s)\n" % (self.data_summary.get('lib_prep_type'), lib_prep_version_num), ('tag_highlight'))
		
		self.txt_summary.insert(END, "Number of plates \t\t: ")
		self.txt_summary.insert(END, "%s\n\n" % self.data_summary.get('num_plates'), ('tag_highlight'))

		# insert a row for each plate
		# TODO: add stripes; make every even numbered row light grey background to make it easier on the eye
		plate_index = 1
		for plate in self.data_summary.get('plates'):
			s_plt_idx = str(plate_index)
			self.txt_summary.insert(END, "[ %s ] \tBarcode \t: " % s_plt_idx)
			self.txt_summary.insert(END, "%s" % self.data_summary.get('plates').get(s_plt_idx).get('barcode'), ('tag_highlight'))
			self.txt_summary.insert(END, "\t\tNum Samples\t: ")
			self.txt_summary.insert(END, "%s" % self.data_summary.get('plates').get(s_plt_idx).get('count_sample_wells'), ('tag_highlight'))
			self.txt_summary.insert(END, "\t\tNum Controls\t: ")
			self.txt_summary.insert(END, "%s\n" % self.data_summary.get('plates').get(s_plt_idx).get('count_control_wells'), ('tag_highlight'))
			plate_index += 1

		# set required plates variable to num plates + 1 (for standards int plate)
		self.var_num_blk_plts_reqd.set(str(self.data_summary.get('num_plates') + 1))

		return

	def clear_screen(self):
		'''Clears the various GUI widgets'''

		if(args.debug == True):
			print("In QuantSetupGUI.%s" % inspect.currentframe().f_code.co_name)

		# clear the summary text and message panel
		self.txt_summary.delete('1.0', END)
		self.txt_msg_panel.delete('1.0', END)

		self.txt_summary.insert('1.0', "")
		self.txt_msg_panel.insert('1.0', "")

		# clear black plate reqd field and set combo back to zero
		self.var_num_blk_plts_reqd.set(0)
		self.num_blk_plates_deck.set(0)

		# disable the create files button
		self.btn_create_files.config(state = DISABLED)

		return

# -----------------------------------------------------------------------------
# Common GUI Elements
# -----------------------------------------------------------------------------
def setup_styles_and_themes():
	'''Create any common styles and themes for the GUI screens'''

	if(args.debug == True):
		print("In %s" % inspect.currentframe().f_code.co_name)

	# to get the combobox to display properly without a grey border we need to apply a style theme
	combobox_style 			= ttk.Style()
	combobox_style.theme_create('combobox_style', parent = 'alt',
				settings 	= {'TCombobox':
                                {'configure':
                                    {	'font':font_arial_normal,
                                    	'foreground':fg_colour,
                                     	'selectforeground': fg_colour,
                                     	'fieldbackground': bg_colour,
                                     	'background': bg_colour,
                                     	'selectbackground': bg_colour
                                  	}
                                }
                            }
                         )
	combobox_style.theme_use('combobox_style')

	return


def create_widget_button(widg_frame, widg_params):
	'''Create a button widget and add to a grid position within a frame

	Args:
		widg_frame 	- the frame in which the widget will be set.
		widg_params - dictionary of parameters for creating the widget and for
					setting the widget's position within the frame.
	'''

	# create the button
	button 	    			= Button(widg_frame,
				text 		= widg_params.get('widg_text', ""),
				textvariable= widg_params.get('widg_txt_var', None),
				width   	= widg_params.get('widg_width', 16),
				bg 			= widg_params.get('widg_bg', bg_colour),
				fg 			= widg_params.get('widg_fg', fg_colour),
				font 		= widg_params.get('widg_font', font_arial_normal),
				image   	= widg_params.get('widg_image', None),
				justify 	= widg_params.get('widg_justify', CENTER),
				anchor  	= widg_params.get('widg_anchor', CENTER),				
				command 	= widg_params.get('widg_command', None),
				state 		= widg_params.get('widg_state', NORMAL))

	# add the button to a grid
	add_widget_to_grid(button, widg_frame, widg_params)

	return button

def create_widget_combobox(widg_frame, widg_params):
	'''Create a combobox widget and add to a grid position within a frame

	Args:
		widg_frame 	- the frame in which the widget will be set.
		widg_params - dictionary of parameters for creating the widget and for
					setting the widget's position within the frame.
	'''

	if(args.debug == True):
		print("In %s" % inspect.currentframe().f_code.co_name)
		pprint(widg_params)

	# create the combobox (NB. style and theme set beforehand to make this display well)
	combobox 	    		= ttk.Combobox(widg_frame, 
				textvariable= widg_params.get('widg_txt_var', None),
				style       = 'SCG.TCombobox')
	combobox.config(values  = widg_params.get('widg_values', None),
				state   	= widg_params.get('widg_state', 'readonly'),
				font 		= widg_params.get('widg_font', font_arial_normal),
				width 		= widg_params.get('widg_width', 3))

	# combox doesn't recognise border so ensure disabled
	widg_params['grid_has_border'] = False

	# add the combobox to a grid (combobox is not a tkinter widget so is treated differently)
	add_widget_to_grid(combobox, widg_frame, widg_params)

	return combobox

def create_widget_label(widg_frame, widg_params):
	'''Create a label widget and add to a grid position within a frame

	Args:
		widg_frame 	- the frame in which the widget will be set.
		widg_params - dictionary of parameters for creating the widget and for
					setting the widget's position within the frame.
	'''

	if(args.debug == True):
		print("In %s" % inspect.currentframe().f_code.co_name)
		pprint(widg_params)

	# create the label
	label 	    			= Label(widg_frame,
				text 		= widg_params.get('widg_text', ""),
				textvariable= widg_params.get('widg_txt_var', None),
				width  		= widg_params.get('widg_width', None),
				bg 			= widg_params.get('widg_bg', bg_colour),
				fg 			= widg_params.get('widg_fg', fg_colour),
				font 		= widg_params.get('widg_font', font_arial_normal),
				image   	= widg_params.get('widg_image', None),
				justify 	= widg_params.get('widg_justify', LEFT),
				anchor  	= widg_params.get('widg_anchor', W),
				state 		= widg_params.get('widg_state', NORMAL))

	# add the label to a grid
	add_widget_to_grid(label, widg_frame, widg_params)

	return label

def create_widget_text(widg_frame, widg_params):
	'''Create a text widget and add to a grid position within a frame

	Args:
		widg_frame 	- the frame in which the widget will be set.
		widg_params - dictionary of parameters for creating the widget and for
					setting the widget's position within the frame.
	'''

	if(args.debug == True):
		print("In %s" % inspect.currentframe().f_code.co_name)
		pprint(widg_params)

	# create the text widget
	text 	    			= Text(widg_frame,
				wrap 		= widg_params.get("widg_wrap", WORD),
				height   	= widg_params.get('widg_height', 1),
				bg 			= widg_params.get('widg_bg', bg_colour),
				fg 			= widg_params.get('widg_fg', fg_colour),
				font 		= widg_params.get('widg_font', font_arial_normal))

	# write the text into the widget
	text.insert("1.0", widg_params.get("widg_text", ""))

	# have to set state last or a disabled state prevents text being written
	text.config(state = widg_params.get('widg_state', NORMAL))

	# add the text to a grid
	add_widget_to_grid(text, widg_frame, widg_params)

	return text

def add_widget_to_grid(widg, widg_frame, widg_params):
	'''Position a widget's grid position within a frame

	Called from the various create widget methods.
	Args:
		widg        - the widget
		widg_frame 	- the frame in which the widget will be set.
		widg_params - dictionary of parameters for creating the widget and for
					setting the widget's position within the frame.
	'''

	if(args.debug == True):
		print("In %s" % inspect.currentframe().f_code.co_name)
		pprint(widg_params)

	# add the widget to the grid
	widg.grid( 	row 		= widg_params.get('grid_row', None),
				column 		= widg_params.get('grid_col', None), 
				columns 	= widg_params.get('grid_cols', None),
				ipadx 		= widg_params.get('grid_ipadx', 2),
				ipady 		= widg_params.get('grid_ipady', 2),
				padx 		= widg_params.get('grid_padx', 2),
				pady 		= widg_params.get('grid_pady', 2),
				sticky 		= widg_params.get('grid_sticky', None))

	if(widg_params.get('grid_has_border', False)):
		widg.config(borderwidth=1, relief="solid")

	return

# -----------------------------------------------------------------------------
# Create GUI Screen Methods
# -----------------------------------------------------------------------------
def display_gui_quantsetup():
	'''Display the quantification setup GUI'''

	if(args.debug == True):
		print("In %s" % inspect.currentframe().f_code.co_name)
		print pprint(settings)

	root 	= Tk()
	root.geometry('%dx%d+%d+%d' % (settings.get('Common').get('gui_width'), settings.get('Common').get('gui_height'), settings.get('Common').get('gui_x_posn'), settings.get('Common').get('gui_y_posn')))
	root.title("Access System Script")
	app 	= QuantSetupGUI(root)

	root.mainloop()

	return

# -----------------------------------------------------------------------------
# Startup
# -----------------------------------------------------------------------------
def main():
	'''Entry point for the program

	Reads command line arguments and configuration file(s) and decides which 
	processing route to take.
	'''
	
	# read in the command line arguments
	parse_command_line_arguments()

	# read in the configuration file
	parse_access_system_config_file()

	# call a function based on the mode name (mode type validated already)
	method_name = 'process_' + str(args.mode)
	return getattr(sys.modules[__name__], method_name)()
    
if __name__ == '__main__':
    main()