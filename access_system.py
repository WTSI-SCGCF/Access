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

2. Standards file(s) that hold the layout of the relevant standards plate (which 
   file depends on the source plate sample type e.g. cDNA, gDNA etc.)

3. Library Prep Parameters file(s) that hold the parameters for the varying types
   of library preparation (which file depends on the source plate library prep 
   parameters type e.g. ??) 

4. Template RunDef file into which to insert run-specific parameters, some of
   which depend upon the numbers of plates.

Outputs:

1. Quantification RunDef file for the Labcyte Tempo software.

2. Echo transfer csv files specific to source DNA plates and standards plate.

3. Experiment run log to detail what was done.

4. Input LIMS file moved to experiment run directory to record what was 
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


TODO:
* add logging with a generic append to log function to the expt dir,
  this should log key information only, e.g. user choices, refs to files and where. 
  they have been written. Timestamp per entry row.

* can we replace the text file message in "Your plates are quantified.txt" with something else?


'''

import argparse 	# to parse command line arguments
import sys 			# for sys.exit
import os 			# for file directory selection
import json 		# for reading/writing json files
import csv 			# for reading/writing csv files
import inspect 		# for getting method name
import re 			# for regex expressions
import shutil 		# for file copying

from configobj 		import ConfigObj, ConfigObjError # for reading config files
from Tkinter        import * 				# for the GUI interfaces
import ttk 									# for the GUI interface widgets
from tkMessageBox 	import askyesno 		# for pop-up message boxes
from tkFileDialog 	import askopenfilename 	# for file selection


from pprint import pprint # for pretty printing e.g. lists and dictionaries

# -----------------------------------------------------------------------------
# Variables
# -----------------------------------------------------------------------------
script_version 				= "1.0"

config_filepath 			= '../Access_Configs/config/access_system.cfg' # configuration filename
valid_modes     			= ['quantsetup'] # valid program modes
args     					= {} # stores parsed command line arguments
settings 					= {} # stores parsed configuration settings
valid_quant_standards  		= ['SS2'] # list of the valid standards types
quant_standards				= {} # stores parsed quantification standards

# GUI colours and fonts
colour_white 				= 'white'
colour_black 				= 'black'
colour_yellow 	 			= 'yellow'
colour_blue  				= 'blue'
colour_green 				= '#32701e'
colour_red 					= 'red'
colour_lt_green 			= '#18631e'
colour_light_grey 			= '#e5e5e5'

font_arial_normal 			= 'Arial -16'
font_courier_normal			= 'Courier -12'
font_arial_huge_bold 		= 'Arial -24 bold'
font_arial_large_bold 		= 'Arial -18 bold'
font_arial_medium_bold 		= 'Arial -16 bold'
font_verdana_normal_italic 	= 'Verdana -12 italic'

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
		print_debug_message("DEBUG mode is active")		
		if(args.verbose == True):
			print_debug_message("Verbose mode is active")

	# check if mode is valid
	if args.mode in valid_modes:
		if(args.debug == True):
			print_debug_message("Mode is valid: " + args.mode)
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
		print_debug_message("In %s" % inspect.currentframe().f_code.co_name)

	# Read the configuration file
	config = read_configuration_file(config_filepath)

	# Verify that the options imported from the file match with an expected list maintained here
	options_list = {
		'Version':{'version_number':'flt'},
		'Common':{'dir_tempo_inbox':'str', 
				'filepath_ecp_384armadillo':'str',
				'filepath_ecp_384corningblack':'str',
				'filepath_ecp_384dest':'str',
				'dir_expt_root':'str',
				'dir_rundef_templates':'str',
				'src_plts_initial_stk_posn':'int',
				'gui_width':'int',
				'gui_height':'int',
				'gui_x_posn':'int',
				'gui_y_posn':'int'},
		'Quantification':{'quant_max_source_plates_allowed':'int', 
				'quant_dir_lims_file_network':'str', 
				'quant_dir_standards':'str',
				'quant_filename_standards_ss2':'str',
				'quant_filename_setup_rundef_template':'str',
				'quant_filename_sources_to_standards_csv':'str',
				'quant_filename_sources_to_black_plts_csv':'str',
				'quant_filename_standards_to_black_csv':'str',
				'quant_filename_setup_log':'str',
				'quant_cdna_min_conc':'flt',
				'quant_cdna_max_conc':'flt',
				'quant_cdna_min_percent_wells_ok':'flt',
				'quant_gdna_min_conc':'flt',
				'quant_gdna_max_conc':'flt',
				'quant_gdna_min_percent_wells_ok':'flt'
		}
	}

	global settings

	# Copy the options into the global settings list so that they are available throughout the program
	try:
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
						else:
							sys.exit("Access System Script: ERROR: Configuration file field type not understood for option group <%s> and option <%s>, cannot continue" % (str(opt_group), str(opt)))
					else:
						sys.exit("Access System Script: ERROR: Configuration file format option group <%s> is missing option <%s>, cannot continue" % (str(opt_group), str(opt)))
				else:
					sys.exit("Access System Script: ERROR: Configuration file format missing option group <%s>, cannot continue" % str(opt_group))
	except ValueError as ve:
		sys.exit("Access System Script: ValueError parsing config file into settings from filepath <%s>. Cannot continue. Message: %s" % (config_filepath, ve.message))
	except Exception as e:
		sys.exit("Access System Script: Exception parsing config file into settings from filepath <%s>. Cannot continue. Message: %s" % (config_filepath, e.message))

	# Print out settings if in debug mode
	if(args.debug == True):
		print_debug_message('-' * 80)
		print_debug_message("Settings from config file are as follows:")
		print_debug_message('-' * 80)
		pprint(settings)
		print_debug_message('-' * 80)

	return


def parse_quantification_standards_config_file(stnd_type):
	'''Parse the settings from the selected standards configuration file.

	There will be a version-controlled standards file for each type of standards plate.
	These standards files are located in a directory specified in the main config file,
	and have filenames specified in the main config file.
	'''

	if(args.debug == True):
		print_debug_message("In %s" % inspect.currentframe().f_code.co_name)
		print_debug_message("Input standards type: %s" % stnd_type)

	# Read the relevant configuration file
	standards_config_filename = ""
	if(stnd_type == 'SS2'):
		standards_config_filename = settings.get('Quantification').get('quant_filename_standards_ss2')

	if(args.debug == True):
		print_debug_message("Standards config filename : %s" % standards_config_filename)	

	standards_config_filepath = os.path.join(settings.get('Quantification').get('quant_dir_standards'), standards_config_filename)

	if(args.debug == True):
		print_debug_message("Standards config filepath : %s" % standards_config_filepath)

	config = read_configuration_file(standards_config_filepath)

	# Verify that the options imported from the file match with an expected list maintained here
	options_list = {
		'Version':{'version_number':'flt'},
		'Information':{'stnd_plt_stk_posn':'int'},
		'Kit':{'kit_name':'str',
				'manufacturer':'str',
				'order_number':'str'},
		'Ladder':{'num_of_ladder_wells':'int',
				'sheared_size_of_dna_kb':'flt',
				'num_of_ladder_reps_black_plate':'int'},
		'Pools':{'num_of_pool_reps_black_plate':'int',
				'vol_src_to_pool_nl':'int',
				'vol_src_to_black_plate_nl':'int',
				'vol_pool_to_black_plate_nl':'int',
				'max_num_sources':'int'}
	}

	global quant_standards
	quant_standards[stnd_type] = {}

	# copy the standards into the global standards list so that they are available throughout the program
	try:
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
						else:
							sys.exit("Access System Script: ERROR: Standards file field type not understood for option group <%s> and option <%s>, cannot continue" % (str(opt_group), str(opt)))
					else:
						sys.exit("Access System Script: ERROR: Standards file format option group <%s> is missing option <%s>, cannot continue" % (str(opt_group), str(opt)))
				else:
					sys.exit("Access System Script: ERROR: Standards file format missing option group <%s>, cannot continue" % str(opt_group))

		# retrieve sub-section information for the ladder
		lad_idx 			= 1
		num_ladder_wells 	= int(quant_standards[stnd_type]['Ladder']['num_of_ladder_wells'])

		quant_standards[stnd_type]['Ladder']['wells'] = {}

		while lad_idx <= num_ladder_wells:
			s_lad_idx 		= str(lad_idx)
			quant_standards[stnd_type]['Ladder']['wells'][s_lad_idx] 							= {}
			quant_standards[stnd_type]['Ladder']['wells'][s_lad_idx]['well_posn'] 				= str(config['Ladder'][s_lad_idx]['well_posn'])
			quant_standards[stnd_type]['Ladder']['wells'][s_lad_idx]['concentration_ng_ul']		= float(config['Ladder'][s_lad_idx]['concentration_ng_ul'])
			quant_standards[stnd_type]['Ladder']['wells'][s_lad_idx]['vol_to_dispense_nl'] 		= float(config['Ladder'][s_lad_idx]['vol_to_dispense_nl'])
			quant_standards[stnd_type]['Ladder']['wells'][s_lad_idx]['black_plt_well_locns'] 	= config['Ladder'][s_lad_idx]['black_plt_well_locns'] # creates list of well posn strings
			lad_idx 		+= 1

		# retrieve sub-section information for the pools
		src_idx 			= 1
		max_num_sources 	= int(quant_standards[stnd_type]['Pools']['max_num_sources'])

		quant_standards[stnd_type]['Pools']['sources'] = {}

		while src_idx <= max_num_sources:
			s_src_idx 		= str(src_idx)
			quant_standards[stnd_type]['Pools']['sources'][s_src_idx] 							= {}
			quant_standards[stnd_type]['Pools']['sources'][s_src_idx]['standards_plt_pool_locn']= str(config['Pools'][s_src_idx]['standards_plt_pool_locn'])
			quant_standards[stnd_type]['Pools']['sources'][s_src_idx]['black_plt_well_locns'] 	= config['Pools'][s_src_idx]['black_plt_well_locns'] # creates list of well posn strings
			src_idx 		+= 1
	except ValueError as ve:
		sys.exit("Access System Script: ValueError parsing standards config file from filepath <%s>. Cannot continue. Message: %s" % (standards_config_filepath, ve.message))
	except Exception as e:
		sys.exit("Access System Script: Exception parsing standards config file from filepath <%s>. Cannot continue. Message: %s" % (standards_config_filepath, e.message))


	# print out quant_standards if in debug mode
	if(args.debug == True):
		print_debug_message('-' * 80)
		print_debug_message("Standards are as follows:")
		print_debug_message('-' * 80)
		pprint(quant_standards)
		print_debug_message('-' * 80)

	return

def parse_library_prep_config_file():
	'''Parse the settings from the library preparation configuration file.


	'''

	if(args.debug == True):
		print_debug_message("In %s" % inspect.currentframe().f_code.co_name)

	return

# -----------------------------------------------------------------------------
# Processing Methods
# -----------------------------------------------------------------------------
def process_quantsetup():
	'''Entry method for processing the quantsetup mode'''

	if(args.debug == True):
		print_debug_message("In %s" % inspect.currentframe().f_code.co_name)

	# build the user interface to ask the user to the select files required 
	display_gui_quantsetup()

	return

def process_quantcalc():
	'''Entry method for processing the quantcalc mode'''

	if(args.debug == True):
		print_debug_message("In %s" % inspect.currentframe().f_code.co_name)

	return

# -----------------------------------------------------------------------------
# GUI Screen Classes
# -----------------------------------------------------------------------------
class QuantSetupGUI:
	'''GUI class for quantification setup'''

	# -----------------------------------------------------------------------------
	# Variables
	# -----------------------------------------------------------------------------
	lims_src_plt_grp_filepath		= "" # holds selected lims plate grouping file filepath
	data_lims_src_plt_grp 			= {} # holds selected lims plate grouping data once read from file
	data_summary 					= {} # holds summarised data

	def __init__(self, master):
		'''Initialise a new frame for the GUI'''

		if(args.debug == True):
			print_debug_message("In QuantSetupGUI.%s" % inspect.currentframe().f_code.co_name)

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
		s_script_info = "Single Cell Genomics Core Facility\nScript Version: %s" % script_version
		frame.grid_rowconfigure(0, weight=1)
		lbl_scgcf_params 		= {'widg_text':s_script_info, 
									'widg_fg':colour_lt_green,
									'widg_font':font_verdana_normal_italic,
									'grid_row':0,
									'grid_col':0,
									'grid_cols':2,
									'grid_sticky':NW}
		self.lbl_scgcf 			= create_widget_label(frame, lbl_scgcf_params)

		self.logo_image			= PhotoImage(file='images/sanger_logo.ppm')
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
									'widg_fg':colour_blue,
									'widg_font':font_arial_huge_bold,
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
		self.instr_text  		= "Click 'Select File' to choose the Quantification plate grouping file downloaded "\
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
		lbl_lims_fp_params 		= {'widg_bg':colour_yellow,
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
			print_debug_message("In QuantSetupGUI.%s" % inspect.currentframe().f_code.co_name)

		# clear error msg line in screen, summary text and disable create files button
		self.clear_screen()

		# normpath fixes path for different OSs
		limsdir = os.path.normpath(settings.get('Quantification').get('quant_dir_lims_file_network'))

		if(args.debug == True):
			print_debug_message("LIMS file directory: %s" % limsdir)

		# open file dialog and save selection N.B. if path not recognised it opens anyway
		self.lims_src_plt_grp_filepath = askopenfilename(title  = "Select LIMS plate grouping file", initialdir = limsdir)

		# update the lable to display the filepath
		self.lbl_lims_fp.configure(text = self.lims_src_plt_grp_filepath)

		if((len(self.lims_src_plt_grp_filepath) > 0) and (os.path.isfile(self.lims_src_plt_grp_filepath))):
			if(args.debug == True):
				print_debug_message("File chosen: %s" % self.lims_src_plt_grp_filepath)

			# read file into memory
			self.read_lims_plate_grouping_json_file()

			# validate data and create summary data
			if(self.validate_data_lims_src_plt_grp() == True):

				# at this point we need to know which standards we are using e.g. SS2 and load the relevant file
				parse_quantification_standards_config_file(self.data_summary['standards_type'])

				# TODO: need validation check to make sure standards layout can handle number of plates in this run
				# self.settings[Quantification][quant_max_source_plates_allowed]

				# display a summary of the data in the file (one line per plate) and ask user to confirm generation of rundef and csv files
				self.display_summary_of_plates()

				# TODO: we should also be able to calculate and summarise amount of reagents reqd
				# maybe have links to display what the standards plate looks like, incl where pools will go

				# Enable the Create Files Button
				self.btn_create_files.config(state = NORMAL)
		else:
			if(args.debug == True):
				print_debug_message("No file was selected!")

			# red error msg line in screen
			self.display_message(True, "No file was selected, please try again.")		

		return

	def create_access_files_button_callback(self):
		'''Triggered when user has pressed the Check File button'''

		if(args.debug == True):
			print_debug_message("In QuantSetupGUI.%s" % inspect.currentframe().f_code.co_name)

		# validate that there are enough black plates in the stack
		if(self.validate_number_of_black_plates_in_stack()):
			self.display_message(False, "Validated number of black plates in stack, now creating experiment directory, please wait...")
		else:
			if(args.debug == True):
				print_debug_message("Failed to validate number of black plates in stack")
			return

		# create expt dir using lims_plate_group_id
		if(self.create_experiment_directory()):
			self.display_message(False, "Created experiment directory, now generating Echo csv files, please wait...")
		else:
			if(args.debug == True):
				print_debug_message("Failed to create experiment directory")
			return

		# generate Echo csvs in expt dir
		if(self.generate_csv_files()):
			self.display_message(False, "ECHO csv files created, now generating RunDef file, please wait...")
		else:
			if(args.debug == True):
				print_debug_message("Failed to generate csv files")
			return

		# generate Access RunDef file
		if(self.generate_rundef_file()):
			self.display_message(False, "RunDef file created in experiment directory, now tidying up, please wait...")
		else:
			if(args.debug == True):
				print_debug_message("Failed to generate RunDef file")
			return

  		# create a new copy of the LIMS file and place it in the experiment directory
		try:
			lims_src_plt_grp_filename 	= os.path.basename(self.lims_src_plt_grp_filepath)
			new_expt_filepath 			= os.path.join(self.expt_directory, lims_src_plt_grp_filename)

			shutil.copyfile(self.lims_src_plt_grp_filepath, new_expt_filepath)
			self.display_message(False, "The RunDef file should now be in the Tempo Inbox and ready to run from Tempo.\nYou can now Quit or open another LIMS file.")
		except Exception as e:
			if(args.debug == True):
				print_debug_message("Failed to copy LIMS file into experiment directory")

			self.display_message(True, "ERROR: Exception copying the LIMS plate grouping file into the experiment directory.\nError Message: <%s>" % str(e))
			return False

		return

	def display_message(self, is_error, message):
		'''Display message on the screen'''

		if(args.debug == True):
			print_debug_message("In QuantSetupGUI.%s" % inspect.currentframe().f_code.co_name)
			if(is_error):
				print_debug_message("Input is Error message")
			else:
				print_debug_message("Input is Standard message")
			print_debug_message("Input message: %s" % message)			

		# delete current message content
		self.txt_msg_panel.delete('1.0', END)

		# set up tags
		self.txt_msg_panel.tag_configure('msg_error', 	 font=font_arial_normal, foreground=colour_red)
		self.txt_msg_panel.tag_configure('msg_standard', font=font_arial_normal, foreground=colour_green)

		# display message
		if(is_error):
			self.txt_msg_panel.insert('1.0', message, ('msg_error'))
		else:
			self.txt_msg_panel.insert('1.0', message, ('msg_standard'))

		return

	def quit_button_callback(self):
		'''Triggered when user has pressed the Quit button'''

		if(args.debug == True):
			print_debug_message("In QuantSetupGUI.%s" % inspect.currentframe().f_code.co_name)

		if askyesno('Verify', 'Are you sure you want to Quit?'):
			sys.exit()

		return

	def read_lims_plate_grouping_json_file(self):
		'''Reads the data from a json file'''

		if(args.debug == True):
			print_debug_message("In QuantSetupGUI.%s" % inspect.currentframe().f_code.co_name)

		# read the file N.B. the with clause closes the file as soon as we are out of the with context.
		# don't print or do anything in the with open clause as that keeps the file open longer.
		with open(self.lims_src_plt_grp_filepath) as data_file:    
			self.data_lims_src_plt_grp = json.load(data_file)
		
		if(args.debug == True):
			print_debug_message("LIMS file data:")
			print_debug_message('-' * 80)
			pprint(self.data_lims_src_plt_grp)
			print_debug_message('-' * 80)

		return

	def validate_data_lims_src_plt_grp(self):
		'''Validates the data extracted from the LIMS json file

		This file is created either by the LIMS or manually. It holds the details about the grouping of plates
		about to be processed on the Labcyte Access system. This information is used to generate the Access RunDef
		and Echo transfer csv files required for the quantification process.

		Method returns a boolean indicating validation success or failure.
		'''

		if(args.debug == True):
			print_debug_message("In QuantSetupGUI.%s" % inspect.currentframe().f_code.co_name)

		# validation check: there should be certain fields in the json file
		expected_fields = ['LIMS_PLATE_GROUP_ID', 'PLATES']
		for expected_field in expected_fields:
			if self.data_lims_src_plt_grp.has_key(expected_field):
				if(args.debug == True):
					print_debug_message("Field located : %s" % expected_field)
			else:
				self.display_message(True, "ERROR: Key field <%s> missing from this file. Cannot continue." % expected_field)
				return False

		self.data_summary['lims_plate_group_id'] 	= self.data_lims_src_plt_grp['LIMS_PLATE_GROUP_ID']
		self.data_summary['num_src_plts'] 			= len(self.data_lims_src_plt_grp['PLATES'])
		self.data_summary['plates'] 				= {}

		if(args.debug == True):
			print_debug_message("Run ID           = %s" % self.data_summary['lims_plate_group_id'])
			print_debug_message("Num plates found = %s" % self.data_summary['num_src_plts'])

		# validation check: there should be at least one plate
		if(self.data_summary['num_src_plts'] == 0):
			self.display_message(True, "ERROR: Unable to identify any plates in this file. Cannot continue.")
			return False

		# set black plates required (n + 1 for standards intermediate plate)
		self.num_black_plates_reqd = self.data_summary['num_src_plts'] + 1

		# cycle through the plates in the file and extract key information
		plate_index 					= 1
		num_standards_types 			= 0
		standards_types                 = {}

		while plate_index <= self.data_summary['num_src_plts']:
			s_plt_idx 					= str(plate_index)
			
			# store plate summary details
			self.data_summary['plates'][s_plt_idx] 							= {}
			self.data_summary['plates'][s_plt_idx]['barcode'] 				= self.data_lims_src_plt_grp['PLATES'][s_plt_idx]['BARCODE']
			self.data_summary['plates'][s_plt_idx]['library_prep_params'] 	= self.data_lims_src_plt_grp['PLATES'][s_plt_idx]['LIBRARY_PREP_PARAMS']
			self.data_summary['plates'][s_plt_idx]['standards_params'] 		= self.data_lims_src_plt_grp['PLATES'][s_plt_idx]['STANDARDS_PARAMS']

			if self.data_lims_src_plt_grp['PLATES'][s_plt_idx]['STANDARDS_PARAMS'] not in valid_quant_standards:
				self.display_message(True, "ERROR: plate with barcode <%s> has a standards type of <%s> which is not currently supported. Cannot continue." % (self.data_lims_src_plt_grp['PLATES'][s_plt_idx]['BARCODE'], self.data_lims_src_plt_grp['PLATES'][s_plt_idx]['STANDARDS_PARAMS']))
				return False

			# count number of types of standards plates, used to check not got a mixed set
			if(self.data_lims_src_plt_grp['PLATES'][s_plt_idx]['STANDARDS_PARAMS'] in standards_types):
				standards_types[self.data_lims_src_plt_grp['PLATES'][s_plt_idx]['STANDARDS_PARAMS']] 	+= 1
			else:
				standards_types[self.data_lims_src_plt_grp['PLATES'][s_plt_idx]['STANDARDS_PARAMS']]	= 1

			curr_wells 					= self.data_lims_src_plt_grp['PLATES'][s_plt_idx]['WELLS']

			# count numbers of samples and controls in this plate
			count_sample_wells 			= 0
			count_control_wells 		= 0
			well_index 					= 0

			while well_index < len(curr_wells):
				if(curr_wells[well_index]['ROLE'] == "SAMPLE"):
					count_sample_wells 	+= 1

				if(curr_wells[well_index]['ROLE'] == "CONTROL"):
					count_control_wells += 1

				well_index 				+= 1

			self.data_summary['plates'][s_plt_idx]['count_sample_wells'] 	= count_sample_wells
			self.data_summary['plates'][s_plt_idx]['count_control_wells'] 	= count_control_wells

			plate_index 				+= 1

		# validation check: there is only one type of standards plate required, if more than one it's an error
		if((len(standards_types) == 0) or (len(standards_types) > 1)):
			self.display_message(True, "ERROR: this group of plates contains %s types of standards plate requirements, which is not currently supported. Cannot continue." % len(standards_types))
			return False
		else:
			# fetch first and only key and store as the standards type for the group of plates
			self.data_summary['standards_type'] = standards_types.keys()[0]

		self.display_message(False, "LIMS file successfully validated. Plates found: %s\nPlease check and indicate how many black plates are in stack 4 and press 'Create Access Files'" % str(self.data_summary['num_src_plts']))

		if(args.debug == True):
			print_debug_message(pprint(self.data_summary))

		return True

	def display_summary_of_plates(self):
		'''Displays a summary of the data extracted from the LIMS json file.

		Headings plus one summary line with key information per plate.
		'''

		if(args.debug == True):
			print_debug_message("In QuantSetupGUI.%s" % inspect.currentframe().f_code.co_name)

		# define tags for formatting text
		self.txt_summary.tag_configure('tag_title', font=font_arial_large_bold, foreground=colour_green, relief='raised',justify='center', underline='True')
		self.txt_summary.tag_configure('tag_hl', font=font_arial_medium_bold, foreground=colour_green)
		self.txt_summary.tag_configure('tag_bg_grey', background=colour_light_grey)

		# insert text for summary
		self.txt_summary.insert('1.0', "Summary of File\n", ('tag_title'))
		
		self.txt_summary.insert(END, "LIMS plate grouping ID \t\t\t: ")
		self.txt_summary.insert(END, "%s\n" % self.data_summary.get('lims_plate_group_id'), ('tag_hl'))
		
		self.txt_summary.insert(END, "Standards set \t\t\t: ")
		std_version_num = str(quant_standards[self.data_summary.get('standards_type')]['Version']['version_number'])
		self.txt_summary.insert(END, "%s  (version: %s)\n" % (self.data_summary.get('standards_type'), std_version_num), ('tag_hl'))
		
		# TODO: add library prep set details 
		# self.txt_summary.insert(END, "Library prep set \t\t: ")
		# lib_prep_version_num = str(lib_prep_standards[self.data_summary.get('lib_prep_type')]['Version']['version_number'])
		# self.txt_summary.insert(END, "%s  (version: %s)\n" % (self.data_summary.get('lib_prep_type'), lib_prep_version_num), ('tag_hl'))
		
		self.txt_summary.insert(END, "Num DNA source plates \t\t\t: ")
		self.txt_summary.insert(END, "%s\n\n" % self.data_summary.get('num_src_plts'), ('tag_hl'))

		# insert a row for each plate
		# every odd numbered row has light grey background to make it easier on the eye
		plate_index = 1
		for plate in self.data_summary.get('plates'):
			s_plt_idx = str(plate_index)

			bg_tags = ()
			hl_tags = ('tag_hl')
			if(plate_index % 2 == 1):
				bg_tags = ('tag_bg_grey')
				hl_tags = ('tag_bg_grey', 'tag_hl')

			self.txt_summary.insert(END, "[ %s ] \tBarcode \t: " % s_plt_idx, bg_tags)
			self.txt_summary.insert(END, "%s" % self.data_summary.get('plates').get(s_plt_idx).get('barcode'), hl_tags)
			self.txt_summary.insert(END, "\t\tNum Samples\t: ", bg_tags)
			self.txt_summary.insert(END, "%s" % self.data_summary.get('plates').get(s_plt_idx).get('count_sample_wells'), hl_tags)
			self.txt_summary.insert(END, "\t\tNum Controls\t: ", bg_tags)
			self.txt_summary.insert(END, "%s\n" % self.data_summary.get('plates').get(s_plt_idx).get('count_control_wells'), hl_tags)
			plate_index += 1

		# set required plates variable to num plates + 1 (for standards int plate)
		self.var_num_blk_plts_reqd.set(self.data_summary.get('num_src_plts') + 1)

		return

	def clear_screen(self):
		'''Clears the various GUI widgets'''

		if(args.debug == True):
			print_debug_message("In QuantSetupGUI.%s" % inspect.currentframe().f_code.co_name)

		# clear the filepath 
		self.lbl_lims_fp.configure(text = "")

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

	def validate_number_of_black_plates_in_stack(self):
		'''Validate that the user has indicated that there are enough black plates in the stack'''

		if(args.debug == True):
			print_debug_message("In QuantSetupGUI.%s" % inspect.currentframe().f_code.co_name)

		# Clear the message panel
		self.txt_msg_panel.delete('1.0', END)

		num_blk_deck = self.num_blk_plates_deck.get()
		num_blk_reqd = self.var_num_blk_plts_reqd.get()

		if(args.debug == True):
			print_debug_message("Num plates required = %s" % str(num_blk_reqd))
			print_debug_message("Num plates in stack = %s" % str(num_blk_deck))

		if(num_blk_reqd == 0):
			# error should be > 0
			self.display_message(True, "ERROR: Number of black plates required is zero, should be 1 or more.")
			return False
		elif(num_blk_deck == 0):
			# error user should set to num in stack
			self.display_message(True, "ERROR: Please load sufficient black plates (%s or more required) into stack 4 "\
				"and use the dropdown entry field to set how many black plates there are now in the stack." % str(num_blk_reqd))
			return False
		elif(num_blk_deck < num_blk_reqd):
			# error user needs to load at least reqd plates to stack
			self.display_message(True, "ERROR: Insufficient black plates (%s or more required), please add more and use the "\
				"dropdown entry field to set how many black plates there are now in stack 4." % str(num_blk_reqd))
			return False
		if askyesno('Verify', 'This will generate the Access System RunDef and Echo files. Are you sure?'):
			self.display_message(False, "Starting file creation process, please wait...")

		if(args.debug == True):
			print_debug_message("Successfully validated number of black plates in stack")

		return True

	def create_experiment_directory(self):
		'''Create the experiment directory based on the lims plate grouping id'''

		if(args.debug == True):
			print_debug_message("In QuantSetupGUI.%s" % inspect.currentframe().f_code.co_name)

		self.expt_directory = os.path.join(settings.get('Common').get('dir_expt_root'), self.data_summary['lims_plate_group_id'])

		if(args.debug == True):
			print_debug_message("Attempting to create experiment dir = %s" % self.expt_directory)

		try:
			check_and_create_directory(self.expt_directory)
		except Exception as ex:
			self.display_message(True, "ERROR: Exception creating the experiment directory. Error Message: <%s>" % str(ex))
			return False

		if(args.debug == True):
			print_debug_message("Successfully created experiment directory")

		return True

	def generate_csv_files(self):
		'''Generate the csv files required for the quantification setup

		Creates the 3 csv files required for the Echo transfers:
			1. DNA Source plate wells to pool wells on Standards plate
			2. DNA Source plate wells to Corning Black plates
			3. Standards plate ladder and pool wells to Corning Black plate 
		'''

		if(args.debug == True):
			print_debug_message("In QuantSetupGUI.%s" % inspect.currentframe().f_code.co_name)

		if(not self.generate_sources_to_standards_csv_file()):
			return False

		if(not self.generate_sources_to_corning_black_csv_file()):
			return False

		if(not self.generate_standards_to_corning_black_csv_file()):
			return False

		if(args.debug == True):
			print_debug_message("Successfully generated csv files")

		return True


	def generate_sources_to_standards_csv_file(self):
		'''Generate the csv file for transferring DNA Source plate wells to pool wells on Standards plate

		Generation varies depending on number of pools required.
		'''

		if(args.debug == True):
			print_debug_message("In QuantSetupGUI.%s" % inspect.currentframe().f_code.co_name)

		stnd_type 							= self.data_summary['standards_type']
		
		# open csv file for writing, wb = write in binary format, overwrites the file if the file exists or creates a new file
		try:
			# csv filepath set in confguration file
			csv_filepath_sources_to_standards 	= os.path.join(self.expt_directory, settings.get('Quantification').get('quant_filename_sources_to_standards_csv'))

			with open(csv_filepath_sources_to_standards, 'wb') as csvfile:
				# for csv file open in binary format.  delimiter is defaulted to comma, quote character is defaulted to doublequote, quoting defaults to quote minimal
				csv_writer = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
				
				# write header row
				csv_writer.writerow(['Source Plate Name',
									'Source Plate Barcode',
									'Source Plate Type',
									'Source Well',
									'Transfer Volume',
									'Destination Plate Name',
									'Destination Plate Barcode',
									'Destination Plate Type',
									'Destination well'])

				# create rows for each src plate sample and control for each pool on the standards plate
				src_plt_type 				= '384PP_AQ_SP_High'
				src_transfer_vol 			= quant_standards[stnd_type]['Pools']['vol_src_to_pool_nl'] # volume from standards config file
				dest_plate_name 			= 'DNAQ_STD'
				dest_plate_barcode 			= None
				dest_plate_type 			= '384PP_Dest'

				i_plt_idx 					= 1
				while i_plt_idx <= self.data_summary['num_src_plts']:
					s_plt_idx 					= str(i_plt_idx)

					if(args.debug == True):
						print_debug_message("Plate indx = %s" % s_plt_idx)

					src_plt_name 				= "DNAQ_source_%s" % s_plt_idx
					src_plt_barcode 			= self.data_summary['plates'][s_plt_idx]['barcode']
					
					# fetch the wells for this plate
					curr_wells 					= self.data_lims_src_plt_grp['PLATES'][s_plt_idx]['WELLS']

					# fetch the pool position on the standards plate (N.B. limit of 1 intermediate pool currently )
					standards_plt_pool_locn 	= quant_standards[stnd_type]['Pools']['sources'][s_plt_idx]['standards_plt_pool_locn']

					if(args.debug == True):
						print_debug_message("standards_plt_pool_locn = %s" % str(standards_plt_pool_locn))

					dest_well 					= standards_plt_pool_locn

					# TODO: add validation for all fields before writing the csv row
					if(dest_well == None):
						# error, unable to continue
						self.display_message(True, "ERROR: Unable to determine valid pooling position for sources to standards plate csv creation. "\
							"Source plate barcode <%s>. Check Standards configuration file. Cannot continue." % src_plt_barcode)
						return False
					
					# create csv rows for each sample or control on the source plate
					i_src_well_idx 					= 0
					while i_src_well_idx < len(curr_wells):

						if(args.debug == True):
							print_debug_message("Well indx = %s" % str(i_src_well_idx))

						if(curr_wells[i_src_well_idx]['ROLE'] == "SAMPLE" or curr_wells[i_src_well_idx]['ROLE'] == "CONTROL"):

							# append line to csv
							csv_writer.writerow([src_plt_name,
								src_plt_barcode,
								src_plt_type,
								curr_wells[i_src_well_idx]['POSITION'],
								src_transfer_vol,
								dest_plate_name,
								dest_plate_barcode,
								dest_plate_type,
								dest_well])

						i_src_well_idx 				+= 1

					i_plt_idx 					+= 1

		except IOError as e:
			self.display_message(True, "ERROR: IOException writing Echo csv for the Standards intermediate plate to it's Black plate.\nError Code: <%s> Message: <%s>" % (str(e.errno), e.strerror))
			return False
		except Exception as ex:
			self.display_message(True, "ERROR: Exception writing Echo csv for the Standards intermediate plate to it's Black plate.\nError Message: <%s>" % str(ex))
			return False

		return True

	def generate_sources_to_corning_black_csv_file(self):
		'''Generate the csv file for transferring DNA Source plate wells to corning black plates

		This transfer is a straight map of any samples and controls to their equivalent wells on the corresponding black plate.
		'''

		if(args.debug == True):
			print_debug_message("In QuantSetupGUI.%s" % inspect.currentframe().f_code.co_name)

		stnd_type 							= self.data_summary['standards_type']

		# open csv file for writing, wb = write in binary format, overwrites the file if the file exists or creates a new file
		try:
			# csv filepath set in confguration file
			csv_filepath_sources_to_black_plts 	= os.path.join(self.expt_directory, settings.get('Quantification').get('quant_filename_sources_to_black_plts_csv'))

			with open(csv_filepath_sources_to_black_plts, 'wb') as csvfile:
				# for csv file open in binary format.  delimiter is defaulted to comma, quote character is defaulted to doublequote, quoting defaults to quote minimal
				csv_writer = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
				
				# write header row
				csv_writer.writerow(['Source Plate Name',
									'Source Plate Barcode',
									'Source Plate Type',
									'Source Well',
									'Transfer Volume',
									'Destination Plate Name',
									'Destination Plate Barcode',
									'Destination Plate Type',
									'Destination well'])

				# create rows for each src plate sample and control for each pool on the standards plate
				src_plt_type 				= '384PP_AQ_SP_High'
				src_transfer_vol 			= quant_standards[stnd_type]['Pools']['vol_src_to_black_plate_nl'] # volume from standards config file
				
				dest_plate_barcode 			= None
				dest_plate_type 			= 'Corning_384PS_Black'

				i_plt_idx 					= 1
				while i_plt_idx <= self.data_summary['num_src_plts']:
					s_plt_idx 					= str(i_plt_idx)

					if(args.debug == True):
						print_debug_message("s_plt_idx - %s" % s_plt_idx)

					src_plt_name 				= "DNAQ_source_%s" % s_plt_idx
					src_plt_barcode 			= self.data_summary['plates'][s_plt_idx]['barcode']

					dest_plate_name 			= "DNAQ_Black_%s" % s_plt_idx
					
					# fetch the wells for this plate
					curr_wells 					= self.data_lims_src_plt_grp['PLATES'][s_plt_idx]['WELLS']

					# create csv rows for each sample or control on the source plate
					i_src_well_idx 				= 0
					while i_src_well_idx < len(curr_wells):
						s_src_well_idx = str(i_src_well_idx)

						if(args.debug == True):
							print_debug_message("s_src_well_idx - %s" % s_src_well_idx)

						if(curr_wells[i_src_well_idx]['ROLE'] == "SAMPLE" or curr_wells[i_src_well_idx]['ROLE'] == "CONTROL"):

	 						# append line to csv
							csv_writer.writerow([src_plt_name,
								src_plt_barcode,
								src_plt_type,
								curr_wells[i_src_well_idx]['POSITION'],
								src_transfer_vol,
								dest_plate_name,
								dest_plate_barcode,
								dest_plate_type,
								curr_wells[i_src_well_idx]['POSITION']])

						i_src_well_idx 			+= 1

					i_plt_idx 				+= 1

		except IOError as e:
			self.display_message(True, "ERROR: IOException writing Echo csv for the Standards intermediate plate to it's Black plate.\nError Code: <%s> Message: <%s>" % (str(e.errno), e.strerror))
			return False
		except Exception as ex:
			self.display_message(True, "ERROR: Exception writing Echo csv for the Standards intermediate plate to it's Black plate.\nError Message: <%s>" % str(ex))
			return False

		return True


	def generate_standards_to_corning_black_csv_file(self):
		'''Generate the csv file for transferring the Standards plate wells to a Corning black plate

		Generation varies depending on number and location of source pools and the number and location of ladder wells.
		'''

		if(args.debug == True):
			print_debug_message("In QuantSetupGUI.%s" % inspect.currentframe().f_code.co_name)

		stnd_type 							= self.data_summary['standards_type']

		# open csv file for writing, wb = write in binary format, overwrites the file if the file exists or creates a new file
		try:
			# csv filepath set in confguration file
			csv_filepath_standards_to_black 	= os.path.join(self.expt_directory, settings.get('Quantification').get('quant_filename_standards_to_black_csv'))

			if(args.debug == True):
				print_debug_message("csv_filepath_standards_to_black = %s" % csv_filepath_standards_to_black)

			with open(csv_filepath_standards_to_black, 'wb') as csvfile:
				# for csv file open in binary format.  delimiter is defaulted to comma, quote character is defaulted to doublequote, quoting defaults to quote minimal
				csv_writer = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

				if(args.debug == True):
					print_debug_message("csv writer open")
				
				# write header row
				csv_writer.writerow(['Source Plate Name',
									'Source Plate Barcode',
									'Source Plate Type',
									'Source Well',
									'Transfer Volume',
									'Destination Plate Name',
									'Destination Plate Barcode',
									'Destination Plate Type',
									'Destination well'])

				if(args.debug == True):
					print_debug_message("written header row")

				# create rows for each src plate sample and control for each pool on the standards plate
				src_plt_name 				= "DNAQ_STD"
				src_plt_barcode 			= None
				src_plt_type 				= '384PP_AQ_SP_High'
				dest_plate_name 			= 'DNAQ_STD'
				dest_plate_barcode 			= None
				dest_plate_type 			= 'Corning_384PS_Black'

				
				# number of ladder wells and number of replicate ladders to make comes from standards configuration
				num_ladder_wells 			= quant_standards[stnd_type]['Ladder']['num_of_ladder_wells']
				num_ladder_reps 			= quant_standards[stnd_type]['Ladder']['num_of_ladder_reps_black_plate']

				if(args.debug == True):
					print_debug_message("num_ladder_wells = %s" % str(num_ladder_wells))
					print_debug_message("num_ladder_reps = %s" % str(num_ladder_reps))

				# perform ladder transfers num_of_ladder_reps_black_plate times
				i_ladder_rep_idx 			= 0
				while i_ladder_rep_idx < num_ladder_reps:

					# s_ladder_rep_idx = str(i_ladder_rep_idx)

					if(args.debug == True):
						print_debug_message("Ladder rep idx = %s" % str(i_ladder_rep_idx))
					
					# loop by ladder well, fetching source and destination well information and transfer volume from standards config file	
					i_ladder_idx 				= 1
					while i_ladder_idx <= num_ladder_wells:

						s_ladder_idx 				= str(i_ladder_idx)

						if(args.debug == True):
							print_debug_message("Ladder idx = %s" % s_ladder_idx)
										
						src_well 					= quant_standards[stnd_type]['Ladder']['wells'][s_ladder_idx]['well_posn']

						if(args.debug == True):
							print_debug_message("src_well = %s" % src_well)

						src_transfer_vol 			= quant_standards[stnd_type]['Ladder']['wells'][s_ladder_idx]['vol_to_dispense_nl']

						if(args.debug == True):
							print_debug_message("src_transfer_vol = %s" % src_transfer_vol)

						dest_well 					= quant_standards[stnd_type]['Ladder']['wells'][s_ladder_idx]['black_plt_well_locns'][i_ladder_rep_idx]

						if(args.debug == True):
							# print("src_well = %s" % src_well)
							# print("src_transfer_vol = %s" % src_transfer_vol)
							print_debug_message("dest_well = %s" % dest_well)

						# append line to csv
						csv_writer.writerow([src_plt_name,
							src_plt_barcode,
							src_plt_type,
							src_well,
							src_transfer_vol,
							dest_plate_name,
							dest_plate_barcode,
							dest_plate_type,
							dest_well])

						i_ladder_idx 				+= 1 # increment ladder position

					i_ladder_rep_idx 			+= 1 # increment ladder replicate


				# perform DNA source pool transfers, which depends on number of plates and number of pool replicates to make
				num_pool_replicates			= quant_standards[stnd_type]['Pools']['num_of_pool_reps_black_plate']

				# source transfer volumes are different depending on whether we are transferring ladder wells or DNA source pools
				src_transfer_vol 			= quant_standards[stnd_type]['Pools']['vol_pool_to_black_plate_nl'] # volume from standards config file

				if(args.debug == True):
					print_debug_message("num_pool_replicates = %s" % num_pool_replicates)
					print_debug_message("src_transfer_vol = %s" % src_transfer_vol)

				i_plt_idx 					= 1
				while i_plt_idx <= self.data_summary['num_src_plts']:
					s_plt_idx 				= str(i_plt_idx)

					if(args.debug == True):
						print_debug_message("s_plt_idx = %s" % s_plt_idx)

					# fetch the pool position on the standards plate 
					src_well 				= quant_standards[stnd_type]['Pools']['sources'][s_plt_idx]['standards_plt_pool_locn']

					if(args.debug == True):
						print_debug_message("src_well = %s" % src_well)

					# create csv row for each pool replicate
					i_pool_rep_idx 			= 0
					while i_pool_rep_idx < num_pool_replicates:

						if(args.debug == True):
							print_debug_message("i_pool_rep_idx = %s" % str(i_pool_rep_idx))
						
						# destination wells depend on number of replicates of pool required and the locations specified in the standards config file
						dest_well 				= quant_standards[stnd_type]['Pools']['sources'][s_plt_idx]['black_plt_well_locns'][i_pool_rep_idx]

						if(args.debug == True):
							print_debug_message("dest_well = %s" % dest_well)

						# append line to csv
						csv_writer.writerow([src_plt_name,
							src_plt_barcode,
							src_plt_type,
							src_well,
							src_transfer_vol,
							dest_plate_name,
							dest_plate_barcode,
							dest_plate_type,
							dest_well])

						i_pool_rep_idx 			+= 1

					i_plt_idx 				+= 1

		except IOError as e:
			self.display_message(True, "ERROR: IOException writing Echo csv for the Standards intermediate plate to it's Black plate.\nError Code: <%s> Message: <%s>" % (str(e.errno), e.strerror))
			return False
		except Exception as ex:
			self.display_message(True, "ERROR: Exception writing Echo csv for the Standards intermediate plate to it's Black plate.\nError Message: <%s>" % str(ex))
			return False

		return True

	def generate_rundef_file(self):
		'''Generate the Tempo RunDef file for operating the Access System'''

		if(args.debug == True):
			print_debug_message("In QuantSetupGUI.%s" % inspect.currentframe().f_code.co_name)

		# directories and filepaths
		try:
			rundef_template_filepath 	= os.path.join(settings.get('Common').get('dir_rundef_templates'), settings.get('Quantification').get('quant_filename_setup_rundef_template'))
			rundef_expt_filename 		= "dna_quantification_group_%s.rundef" % self.data_summary['lims_plate_group_id']
			rundef_expt_filepath 		= os.path.join(self.expt_directory, rundef_expt_filename)
			rundef_tempo_inbox_filepath = os.path.join(settings.get('Common').get('dir_tempo_inbox'), rundef_expt_filename)

			if(args.debug == True):
				print_debug_message("rundef_template_filepath 	= %s" % rundef_template_filepath)
				print_debug_message("rundef_expt_filename     	= %s" % rundef_expt_filename)
				print_debug_message("rundef_expt_filepath     	= %s" % rundef_expt_filepath)
				print_debug_message("rundef_tempo_inbox_filepath 	= %s" % rundef_tempo_inbox_filepath)

			# create a dictionary of search_string : value
			search_term_dict 	= self.generate_rundef_search_dictionary()
			search_list 		= search_term_dict.keys()

			if(args.debug == True):
				print_debug_message("-" * 80)
				print_debug_message("Rundef search term dictionary:")
				pprint(search_term_dict)
				print_debug_message("-" * 80)
				print_debug_message("Rundef search term list:")
				pprint(search_list)
				print_debug_message("-" * 80)

		except Exception as e:
			self.display_message(True, "ERROR: Exception creating the RunDef file when determining directory paths.\nError Message: <%s>" % str(e))
			return False

		# open the template rundef file, and the new output rubdef file, and copy lines across
		# replacing any terms matching those in the search_term_dict
		# N.B. the rundef file is a UTF-8 encoded XML file, so each line needs to be decoded then re-encoded to write
		i_line_num = 1
		try:
			with open(rundef_expt_filepath, 'w') as file_out:
				with open(rundef_template_filepath, 'r') as file_in:
					# check each line in the template
					for line in file_in:
						# decode the line from utf-8 so we can work with it
						u_line = line.decode('utf8')
						# see if the line contains ANY of the search terms
						for search_term in search_list:
							# modify the line to replace the search term with the new lines generated in the dictionary
							u_line = u_line.replace(search_term, search_term_dict[search_term])
						# encode the line to utf-8
						line = u_line.encode('utf8')
						# write the (potentially modified) line out to the new copy of the rundef file in the experiment directory
						file_out.write(line)
						i_line_num += 1
		except Exception as e:
			self.display_message(True, "ERROR: Exception creating the RunDef file at line <%s>.\nError Message: <%s>" % (str(i_line_num), str(e)))
			return False

		# create a new copy of the modified file and place it in the Tempo inbox directory
		try:
			shutil.copyfile(rundef_expt_filepath, rundef_tempo_inbox_filepath)
		except Exception as e:
			self.display_message(True, "ERROR: Exception copying the newly created RunDef file to the experiment directory.\nError Message: <%s>" % str(e))
			return False

		return True

	def generate_rundef_search_dictionary(self):
		'''Generate the dictionary of fields to be replaced in the quantification rundef file'''

		if(args.debug == True):
			print_debug_message("In QuantSetupGUI.%s" % inspect.currentframe().f_code.co_name)

		search_term_dict = {}

		# reference the LIMS plate group id
		search_term_dict['SSS_RUNSET_REFERENCE_ID_SSS'] 					= self.data_summary['lims_plate_group_id']

		# experiment directory root
		search_term_dict['SSS_EXPT_ROOT_DIR_SSS'] 							= self.expt_directory

		# ecp filepath for 384dest type plates
		search_term_dict['SSS_CHERRY_PICK_384DEST_ECP_FP_SSS'] 				= settings.get('Common').get('filepath_ecp_384dest')

		# ecp filepath for corning black type plates
		search_term_dict['SSS_CHERRY_PICK_CORNINGBLACK_ECP_FP_SSS'] 		= settings.get('Common').get('filepath_ecp_384corningblack')

		# csv filepath for sources being pooled into the standards plate
		search_term_dict['SSS_SOURCES_POOL_TO_STANDARDS_CSV_FP_SSS'] 		= os.path.join(self.expt_directory, settings.get('Quantification').get('quant_filename_sources_to_standards_csv'))

		# csv filepath for sources to corning black plates
		search_term_dict['SSS_SOURCES_TO_CORNINGBLACK_CSV_FP_SSS'] 			= os.path.join(self.expt_directory, settings.get('Quantification').get('quant_filename_sources_to_black_plts_csv'))

		# csv filepath for standards plate to corning black plate
		search_term_dict['SSS_STANDARDS_TO_CORNINGBLACK_CSV_FP_SSS'] 		= os.path.join(self.expt_directory, settings.get('Quantification').get('quant_filename_standards_to_black_csv'))

		i_srcs_init_stk_posn 				= settings.get('Common').get('src_plts_initial_stk_posn') # from main config
		stnd_type 							= self.data_summary['standards_type']
		s_stnd_plt_stk_posn 				= quant_standards[stnd_type]['Information']['stnd_plt_stk_posn'] # from standards config
		
		# create the source and destination plate rows for the sources pooling to the standards plate transfers
		# add the first source (DNA source) plate row
		s_src_stk_posn 						= str(i_srcs_init_stk_posn)
		s_pooling_platemap_rows 			= 	(u'\t\t<Plate EchoPlateID="1" PlateName="DNAQ_source_1" PlateType="384PP" Barcode="%s" LidType="" PlateCategory="Source" '
												'LocationURL="deck://Deck/1/%s/" FinalLocation="deck://Deck/1/%s/" PlateAccess="Sequential" PreRunActionSetName="Source" RunActionSetName="" '
												'PostRunActionSetName="" StorageDeviceSetName="" EchoTemplate="DNAQ_source_1" />\n' 
												% (self.data_summary['plates']['1']['barcode'], s_src_stk_posn, s_src_stk_posn))
		
		# next append the destination (Standards) plate row which is in a fixed position (will then stay in the Echo destination drawer while the transfers run)
		s_stnd_plt_id_num 					= str(self.data_summary['num_src_plts'] + 1) # standards destination set to 1 more than total source plates
		s_pooling_platemap_rows 			+= 	(u'\t\t<Plate EchoPlateID="%s" PlateName="DNAQ_STD" PlateType="384PP_Dest" Barcode="" LidType="" PlateCategory="Destination" '
												'LocationURL="deck://Deck/1/%s/" FinalLocation="deck://Deck/1/%s/" PlateAccess="Sequential" PreRunActionSetName="DNAQ_STD" RunActionSetName="" '
												'PostRunActionSetName="DNAQ_STD" StorageDeviceSetName="" EchoTemplate="DNAQ_STD" />\n' 
												% (s_stnd_plt_id_num, s_stnd_plt_stk_posn, s_stnd_plt_stk_posn))

		# append the remaining source (DNA source) plate rows
		i_src_stk_posn 						= i_srcs_init_stk_posn
		i_src_idx      						= 2 # start at second plate

		while i_src_idx <= self.data_summary['num_src_plts']:
			s_src_idx 						= str(i_src_idx)
			s_src_stk_posn 					= str(i_src_stk_posn)
			s_pooling_platemap_rows 		+= 	(u'\t\t<Plate EchoPlateID="%s" PlateName="DNAQ_source_%s" PlateType="384PP" Barcode="%s" LidType="" PlateCategory="Source" '
												'LocationURL="deck://Deck/1/%s/" FinalLocation="deck://Deck/1/%s/" PlateAccess="Sequential" PreRunActionSetName="Source" RunActionSetName="" '
												'PostRunActionSetName="" StorageDeviceSetName="" EchoTemplate="DNAQ_source_%s" />\n' 
												% (s_src_idx, s_src_idx, self.data_summary['plates'][s_src_idx]['barcode'], s_src_stk_posn, s_src_stk_posn, s_src_idx))
			i_src_stk_posn 					+= 1 # increment stack position
			i_src_idx 						+= 1 # increment source plate index

		search_term_dict['SSS_POOLING_PLATEMAP_ROWS_SSS'] 					= s_pooling_platemap_rows.encode('utf-8')

		# create the source and destination rows for the sources to black plates transfers
		# for each source add a destination row
		s_src_to_blks_rows 					= ""
		i_src_stk_posn 						= i_srcs_init_stk_posn
		i_src_idx      						= 1
		
		while i_src_idx <= self.data_summary['num_src_plts']:
			# add a source plate row
			s_src_idx 						= str(i_src_idx)
			s_src_stk_posn 					= str(i_src_stk_posn)
			s_plate_barcode 				= self.data_summary['plates'][s_src_idx]['barcode']
			s_src_to_blks_rows 				+= 	(u'\t\t<Plate EchoPlateID="%s" PlateName="DNAQ_source_%s" PlateType="384PP" Barcode="%s" LidType="" PlateCategory="Source" '
												'LocationURL="deck://Deck/1/%s/" FinalLocation="deck://Deck/1/%s/" PlateAccess="Sequential" PreRunActionSetName="" RunActionSetName="" '
												'PostRunActionSetName="" StorageDeviceSetName="" EchoTemplate="DNAQ_source_%s" />\n' 
												% (s_src_idx, s_src_idx, s_plate_barcode, s_src_stk_posn, s_src_stk_posn, s_src_idx))
			s_black_plate_id_num 			= str(self.data_summary['num_src_plts'] + i_src_idx)

			# add a corresponding destination plate row (these plates come from stack 4 and end up in stack 3)
			s_src_to_blks_rows 				+= 	(u'\t\t<Plate EchoPlateID="%s" PlateName="DNAQ_Black_%s" PlateType="Corning_384PS_Black" Barcode="" LidType="" PlateCategory="Destination" '
												'LocationURL="deck://Deck/4/*/" FinalLocation="deck://Deck/3/*/" PlateAccess="Sequential" PreRunActionSetName="" RunActionSetName="" '
												'PostRunActionSetName="Destination" StorageDeviceSetName="" EchoTemplate="DNAQ_Black_%s" />\n' 
												% (s_black_plate_id_num, s_src_idx, s_src_idx))
			i_src_stk_posn 					+= 1 # increment stack position
			i_src_idx 						+= 1 # increment source plate index

		search_term_dict['SSS_SOURCES_PLATEMAP_ROWS_SSS'] 					= s_src_to_blks_rows.encode('utf-8')

  		# TODO: What does the DataFileName here do?! Can we set this to fix the BMG read filepath?:
		# <PHERAstar_Read_Plate>
		# <ReadPlateParameters>Labcyte.POD.Devices.Pherastar.ReadPlateParameters</ReadPlateParameters>
		# <ProtocolName>AccuClear UHS</ProtocolName>
		# <DataFileName />
		# <LoadUnloadOnly>False</LoadUnloadOnly>
		# <Schedule>Immediately</Schedule>
		# <UserPriority>Normal</UserPriority>
		# <DeviceName>FLUOstar</DeviceName>
		# </PHERAstar_Read_Plate>

		# create the source and destination rows for the standards plate to its black plate transfer
		# add the standards source plate row (location is fixed)
		s_stnd_to_blk_rows 					= 	(u'\t\t<Plate EchoPlateID="1" PlateName="DNAQ_STD" PlateType="384PP" Barcode="" LidType="" PlateCategory="Source" '
												'LocationURL="deck://Deck/1/%s/" FinalLocation="deck://Deck/1/%s/" PlateAccess="Sequential" PreRunActionSetName="Source" RunActionSetName="" '
												'PostRunActionSetName="" StorageDeviceSetName="" EchoTemplate="DNAQ_STD" />\n' 
												% (s_stnd_plt_stk_posn, s_stnd_plt_stk_posn))

		# append the destination black plate row (these plates come from stack 4 and end up in stack 3)
		s_stnd_to_blk_rows 					+= 	(u'\t\t<Plate EchoPlateID="2" PlateName="DNAQ_Black" PlateType="Corning_384PS_Black" Barcode="" LidType="" PlateCategory="Destination" '
												'LocationURL="deck://Deck/4/*/" FinalLocation="deck://Deck/3/*/" PlateAccess="Sequential" PreRunActionSetName="" RunActionSetName="" '
												'PostRunActionSetName="DNAQ_Black" StorageDeviceSetName="" EchoTemplate="DNAQ_Black" />\n')

		search_term_dict['SSS_STANDARDS_PLATEMAP_ROWS_SSS'] 				= s_stnd_to_blk_rows.encode('utf-8')


		# PlateID,EchoPlateID,PlateName,PlateBarcode,PlateType,PlateCategory,LidType,Sealed,PlateStatus,OriginalLocation,FinalLocation,
		# CurrentLocation,RunLocation,CurrentRotation,PlateRename,PlateState,LidLocation

		# create the standards plate row for the plate storage map
		s_plate_storage_std_row  			= (u'1844,1,DNAQ_STD,,384PP_Dest,Destination,,False,Unknown,deck://Deck/1/1/,deck://Deck/1/1/,deck://Deck/1/1/,,0,,Unknown,,\n')

		search_term_dict['SSS_PLATE_STORAGE_STD_PLATE_SSS'] 					= s_plate_storage_std_row.encode('utf-8')

		# create the source plate rows for the plate storage map
		# e.g. 1845,0,DNAQ_source_1,,384PP,Source,,False,Unknown,deck://Deck/1/5/,deck://Deck/1/5/,deck://Deck/1/5/,,0,,Unknown,,
		# 1845 -> 1860 so 16 max
		s_plate_storage_src_rows			= ""
		i_src_stk_posn 						= i_srcs_init_stk_posn
		i_src_idx      						= 1
		i_plate_id 							= 1845
		
		while i_src_idx <= self.data_summary['num_src_plts']:
			s_src_idx 						= str(i_src_idx)
			s_src_stk_posn 					= str(i_src_stk_posn)
			s_plate_id 						= str(i_plate_id)
			s_plate_barcode 				= self.data_summary['plates'][s_src_idx]['barcode']

			s_plate_storage_src_rows 		+= (u'%s,%s,DNAQ_source_%s,%s,384PP,Source,,False,Unknown,deck://Deck/1/%s/,deck://Deck/1/%s/,deck://Deck/1/%s/,,0,,Unknown,,\n'
											% (s_plate_id, s_src_stk_posn, s_src_idx, s_plate_barcode, s_src_stk_posn, s_src_stk_posn, s_src_stk_posn))

			i_src_stk_posn 					+= 1 # increment stack position
			i_src_idx 						+= 1 # increment source plate index
			i_plate_id 						+= 1 # increment source plate id

		search_term_dict['SSS_PLATE_STORAGE_SRC_PLATES_SSS'] 					= s_plate_storage_src_rows.encode('utf-8')

		# create the destination plate rows for the plate storage map
		# e.g. 1861,0,DNAQ_Black,,Corning_384PS_Black,Destination,,False,Unknown,deck://Deck/4/1/,deck://Deck/3/*/,deck://Deck/4/1/,,0,,Unknown,,
		# 1861 -> 1880  NEED MAX 40!!
		s_plate_storage_dest_rows			= ""
		i_dest_idx      					= 1
		i_plate_id 							= 1861
		
		while i_dest_idx <= self.num_blk_plates_deck.get():
			s_plate_id 						= str(i_plate_id)
			s_dest_stk_posn 				= str(i_dest_idx)

			s_plate_storage_dest_rows 		+= (u'%s,0,DNAQ_Black,,Corning_384PS_Black,Destination,,False,Unknown,deck://Deck/4/%s/,deck://Deck/3/*/,deck://Deck/4/%s/,,0,,Unknown,,\n'
											% (s_plate_id, s_dest_stk_posn, s_dest_stk_posn))

			i_dest_idx 						+= 1 # increment source plate index
			i_plate_id 						+= 1 # increment source plate id

		search_term_dict['SSS_PLATE_STORAGE_DEST_PLATES_SSS'] 					= s_plate_storage_dest_rows.encode('utf-8')

		return search_term_dict

# -----------------------------------------------------------------------------
# Utility Methods
# -----------------------------------------------------------------------------
def check_and_create_directory(directory):
	'''Check whether a directory exists and create it if not'''

	if(args.debug == True):
		print_debug_message("In %s" % inspect.currentframe().f_code.co_name)

	if(args.debug == True):
		print_debug_message("Attempting to create directory = %s" % directory)

	# check if the directory exists and create it if not
	if not os.path.exists(directory):
		os.makedirs(directory)

	return

def print_debug_message(message):
	'''Print a debug message'''

	print("DEBUG: " + str(message))

	return

def read_configuration_file(filepath):
	'''Reads a configuration file and checks for errors'''

	if(args.debug == True):
		print_debug_message("In %s" % inspect.currentframe().f_code.co_name)

	try:
		print_debug_message("Attempting to read config")
		config = ConfigObj(filepath)
	except ValueError as ve:
		print_debug_message("Caught ValueError")
		if(args.debug == True):
			print_debug_message("-" * 80)
			print_debug_message("ValueError raised:")
			print_debug_message("Line: <%s> %s" % (ve.line_number, ve.line))
			print_debug_message("Msg : %s" % ve.message)
			print_debug_message("-" * 80)

		sys.exit("Access System Script: ValueError parsing config file from filepath  <%s>. Cannot continue. Message: %s" % (filepath, ve.message))
	except ConfigObjError as ce:
		print_debug_message("Caught ConfigObjError")
		if(args.debug == True):
			print_debug_message("-" * 80)
			print_debug_message("Part of Config that parsed successfully:")
			pprint(ce.config)
			print_debug_message("-" * 80)
			print_debug_message("Errors raised:")
			for er in ce.errors:
				print_debug_message("- " * 40)
				print_debug_message("Line: <%s> %s" % (er.line_number, er.line))
				print_debug_message("Msg : %s" % er.message)
			print_debug_message("-" * 80)

		sys.exit("Access System Script: ConfigObjError parsing config file from filepath  <%s>. Cannot continue. Message: %s" % (filepath, er.message))

	return config

def regex_replace_field_in_xml(xml_string, re_string, replacement_string):
	'''Use a Regular Expression to replace one value in an XML string with another'''

	if(args.debug == True):
		print_debug_message("In %s" % inspect.currentframe().f_code.co_name)
		print_debug_message("XML string: %s" % xml_string)
		print_debug_message("RegEx string: %s" % re_string)
		print_debug_message("Replacement string: %s" % replacement_string)


	# Substitute in a string to replace the regex string currently there
	# substitutions		= {re_string: replacement_string, ...}
	substitutions		= {re_string: replacement_string}
	pattern 			= re.compile(r'%([^%]+)%')
	xml_string_modified = re.sub(pattern, lambda m: substitutions[m.group(1)], xml_string)

	if(args.debug == True):
		print_debug_message("Modified XML: %s" % xml_string_modified)

	# return modified XML string
	return xml_string_modified

def append_to_log(log_filepath, message):
	'''Append a message line to a specified log file'''
	
	# quant_filename_setup_log

	return

# -----------------------------------------------------------------------------
# Common GUI Elements
# -----------------------------------------------------------------------------
def setup_styles_and_themes():
	'''Create any common styles and themes for the GUI screens'''

	if(args.debug == True):
		print_debug_message("In %s" % inspect.currentframe().f_code.co_name)

	# to get the combobox to display properly without a grey border we need to apply a style theme
	combobox_style 			= ttk.Style()
	combobox_style.theme_create('combobox_style', 
								parent = 'alt',
								settings 	= {'TCombobox':
				                                {'configure':
				                                    {	'font':font_arial_normal,
				                                    	'foreground':colour_black,
				                                     	'selectforeground': colour_black,
				                                     	'fieldbackground': colour_white,
				                                     	'background': colour_white,
				                                     	'selectbackground': colour_white
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
				bg 			= widg_params.get('widg_bg', colour_white),
				fg 			= widg_params.get('widg_fg', colour_black),
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
		print_debug_message("In %s" % inspect.currentframe().f_code.co_name)
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
		print_debug_message("In %s" % inspect.currentframe().f_code.co_name)
		pprint(widg_params)

	# create the label
	label 	    			= Label(widg_frame,
				text 		= widg_params.get('widg_text', ""),
				textvariable= widg_params.get('widg_txt_var', None),
				width  		= widg_params.get('widg_width', None),
				bg 			= widg_params.get('widg_bg', colour_white),
				fg 			= widg_params.get('widg_fg', colour_black),
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
		print_debug_message("In %s" % inspect.currentframe().f_code.co_name)
		pprint(widg_params)

	# create the text widget
	text 	    			= Text(widg_frame,
				wrap 		= widg_params.get("widg_wrap", WORD),
				height   	= widg_params.get('widg_height', 1),
				bg 			= widg_params.get('widg_bg', colour_white),
				fg 			= widg_params.get('widg_fg', colour_black),
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
		print_debug_message("In %s" % inspect.currentframe().f_code.co_name)
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
		print_debug_message("In %s" % inspect.currentframe().f_code.co_name)

	root 		= Tk()
	root.geometry('%dx%d+%d+%d' % (settings.get('Common').get('gui_width'), settings.get('Common').get('gui_height'), settings.get('Common').get('gui_x_posn'), settings.get('Common').get('gui_y_posn')))
	root.title("Access System Script")
	app 		= QuantSetupGUI(root)

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