
animal_name: eg'JB_181017_4'
patching_date: eg '181127LGN'
experimentator: eg 'SW'
cellname: eg '0002'
eye_inj_ord: 1 (right eye red) or 0
slice_nr: 1 
brain_contra_ipsi: 1 (contra red, ipsi green) or 0
ocular_category: [0:6] refering to cat {non responseive, contra, ipsi, bino, contra silent ipsi, ...
		ipsi silent contra, one eye only has AMPA but no NMDA}
MD: 0 or 1
clear_sl_nr: eg 2
hemisphere: eg 'LH'
photodiode_flag: 0 (one means problem)
step_red: data on redpulse (first pulse) of both red step protocol and double pulse step ...
		protocol (first column is always red only step protocol). neg referes to -70mV ...
		measurements and _pos refers to +40mV. other info about second step protocol and ...
		analyis is also contained here. 
step_red.neg_peak1: [2×11 double] (protocol,step) peak amp during first pulse of -70mV protocols 
step_red.neg_mean1: [2×11 double] 
step_red.neg_fail1: [2×11 logical] test of response based on 4 * baseline std (not used)
step_red.neg_integ1: [2×11 double] area under curve (not used)
step_red.neg_PD: [2×11 double] photodiode signal 
step_red.neg_irr_red: [2×11 double] double] irradiance of red pulse
step_red.neg_laser_amp: [2×11 double]  ?
step_red.neg_Rs: [2×11 double] series reistance
step_red.neg_Cm: [2×11 double] membrane capasitance
step_red.neg_base_mean1: [2×11 double] baseline mean
step_red.pos_peak1: [2×11 double] 
step_red.pos_mean1: [2×11 double]
step_red.pos_fail1: [2×11 logical]
step_red.vpos_PD: [2×11 double]
step_red.pos_irr_red: [2×11 double]
step_red.pos_laser_amp: [2×11 double]
step_red.pos_Rs: [2×11 double] 
step_red.pos_Cm: [2×11 double] 
step_red.pos_base_mean1: [2×11 double] 
step_red.ephys_traces_70: [1000×11×2 double] 
step_red.fit_traces_70: [1000×11×2 double] 
step_red.diff_traces_70: [1000×11×2 double] 
step_red.ephys_traces_40: [1000×11×2 double] 
step_red.fit_traces_40: [1000×11×2 double] 
step_red.diff_traces_40: [1000×11×2 double] 
step_red.go_fit: [11×4 double]
step_red.fit_param: {11×4 cell}
step_red.resp_tests: [1×1 struct]
step_red.red_resp_AMPA: eg. 429.2371, mean of response
step_red.blue_resp_AMPA: eg. 0, mean of response (is set to 0 if noparametric test fails)
step_red.red_resp_AMPA_peak: 834.7846, peak of response
step_red.blue_resp_AMPA_peak: 0, peak of response (is set to 0 if noparametric test fails)
step_red.steps_use_AMPA: 4, indicates the step(s) that is used to caluclate the peak/mean response (for red this is alway just 1 value but for green this can be ether 1 or 6)
step_red.red_resp_NMDA: 122.3569
step_red.blue_resp_NMDA: 7.6363
step_red.red_resp_NMDA_peak: 174.3744
step_red.blue_resp_NMDA_peak: 20.7162
step_red.steps_use_NMDA: 4
step_blue: data on blue pulse (second pulse) of both red step protocol and double pulse step ...
		protocol (first column is always red only step protocol). neg referes to -70mV ...
		measurements and _pos refers to +40mV. See comments for step_red subfields for explanations 
step_blue.neg_peak2: [2×11 double] 
step_blue.neg_mean2: [2×11 double]
step_blue.neg_fail2: [2×11 logical]
step_blue.neg_PD: [2×11 double] 
step_blue.neg_irr_blue: [2×11 double]
step_blue.neg_laser_amp: [2×11 double] 
step_blue.neg_Rs: [2×11 double] 
step_blue.neg_Cm: [2×11 double] 
step_blue.neg_base_mean2: [2×11 double] 
step_blue.pos_peak2: [2×11 double]
step_blue.pos_mean2: [2×11 double]
step_blue.pos_fail2: [2×11 double]
step_blue.pos_integ2: [2×11 double] 
step_blue.pos_PD: [2×11 double] 
step_blue.pos_irr_blue: [2×11 double]
step_blue.pos_laser_amp: [2×11 double]
step_blue.pos_Rs: [2×11 double]
step_blue.pos_Cm: [2×11 double]
step_blue.pos_base_mean2: [2×11 double]
ODI_AMPA_step: 1, AMPA ODI based on mean during response period
ODI_NMDA_step: 0.8825
ODI_AMPA_step_peak: 1, AMPA ODI based on peak during response period
ODI_NMDA_step_peak: 0.8825
AMPA_NMDA_r_red: 4.7873, AMPA response peak / NMDA response peak to red light
AMPA_NMDA_r_blue: 0
red_failure_AMPA: [1×1 struct] field containing all info on failure recording at -70mV with red light stim
red_failure_AMPA.traces_all: [1000×50 double]
red_failure_AMPA.peaks: [1×50 double]
red_failure_AMPA.test: 1
red_failure_AMPA.resp_thresh: [1×50 double]
red_failure_AMPA.resp_idx: [1x50 logical] 
red_failure_AMPA.steady_state: 6
red_failure_AMPA.resp_prop: 0.7778
red_failure_AMPA.avg_amp: 152.2819
red_failure_AMPA.photodiode_mean_std: [2×1000 double]
red_failure_AMPA.IR_pulse1: [50×1 double]
red_failure_AMPA.IR_pulse2: [50×1 double]
red_failure_AMPA.Rs_cell: [50×1 double]
red_failure_AMPA.Rm_cell: [50×1 double]
red_failure_AMPA.Cm_cell: [50×1 double]
blue_failure_AMPA: [1×1 struct] field containing all info on failure recording at -70mV with blue light stim (see red_failure_AMPA for subfield explanations)
red_failure_NMDA: [1×1 struct] field containing all info on failure recording at +40mV with red light stim (see red_failure_AMPA for subfield explanations)
blue_failure_NMDA: [1×1 struct] field containing all info on failure recording at +40mV with blue light stim (see red_failure_AMPA for subfield explanations)
scracm_blue: [] 
scracm_red: []
schematic_file_loc: 'I:\Martin Fernholz\LGN _Project_Common\cell_loc_info2\LGN_schematic_stack.tif'
schematic_loc_xyz: [64 54 146], xyz coordinates of cell in dLGN Schematic 
schematic_imagej_zipfile_loc: 'I:\Martin Fernholz\LGN _Project_Common\cell_loc_info2\JB_181017_4\ephys Simon\slice 1\Schematic\cell_loc_schematic.zip'
schematic_loc_zone: 3, zone number in which cell lies within schematic 
Overview2p: [512×512×3 double] xyz coordinates of cell in 2 photon stack
FlourLogRatio_mean: -0.301, red/green axon flourescence ratio around the cell 
FlourLogRatio_data: [1x1 struct], stack of axon flourescence arround cell
FlourLogRatio_data.mouse: 'JB_181126_1'
FlourLogRatio_data.slice: '190124MF1'
FlourLogRatio_data.cell_name: 'MF0003'
FlourLogRatio_data.pos_xyz: [1012 541 29]
FlourLogRatio_data.LogRatioF: -0.3010, average of Stack_LogRatioF
FlourLogRatio_data.Stack_red: [187x187x67 double], red stack around cell
FlourLogRatio_data.Stack_green: [187x187x67 double], blue stack around cell
FlourLogRatio_data.Stack_LogRatioF: [187x187x67 double], flourescence log ratio stack around cell
FlourLogRatio_data.arbour_radius: 150
FlourLogRatio_data.pixel_res: [1.6100 1.6100 4]
Confocal_data_loc: {1x3 cell}, location of data for FlourLogRatio_mean & FlourLogRatio_data





