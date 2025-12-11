import numpy as np
import h5py as h5

"""
This function applies the reverse calibration to A-STEP simulations.

Steps
    1: Define file names and parameters
    2: Read .sim file
    3: Transform hits in (x,y,z) to layer, chip, row, column
    4: Determine clock timestamps
    5: Transform from energy to ToT, smear, and apply thresholds
    6: Define output array
    7: Write to output file

"""

class ASTEP_RevCal:
    def __init__(self,sim_name,calib_name,is_h5,seed):
        self.sim_name   = sim_name
        self.calib_name = calib_name
        self.is_h5      = is_h5
        self.seed       = seed
        
        if is_h5:
            self.out_name = sim_name + '.ASTEP.h5'
        else:
            self.out_name = sim_name + '.ASTEP.csv'
        
        self.Chip_Offset_X = (0.1475 + .018 + .06)/2
        self.Chip_Offset_Z = (.069 + 2*.06)/2
        
        self.Chip_Width_X = 1.75
        self.Chip_Width_Z = 1.75
        
        self.Chip_Spacing_X = self.Chip_Width_X + 2*self.Chip_Offset_X
        self.Chip_Spacing_Z = self.Chip_Width_Z + 2*self.Chip_Offset_Z
        
        self.N_Chip_Xs = 2
        self.N_Chip_Zs = 2
        
        self.N_pixel_Xs = 35
        self.N_pixel_Zs = 35
        
        self.pixel_size_X = self.Chip_Width_X/self.N_pixel_Xs
        self.pixel_size_Z = self.Chip_Width_Z/self.N_pixel_Zs
        
        self.Layer_Offset_Y = 0
        self.Layer_Spacing_Y = 0.965
        
        self.x0 = -(self.N_Chip_Xs - 1)*self.Chip_Offset_X - (self.N_Chip_Xs/2)*self.Chip_Width_X + self.pixel_size_X/2
        self.z0 = -(self.N_Chip_Zs - 1)*self.Chip_Offset_Z - (self.N_Chip_Zs/2)*self.Chip_Width_Z + self.pixel_size_Z/2
        
        self.AstroPix_Clock_Freq = 1000000 # 1 MHz
        self.AstroPix_Max_Clock = 2**8
        
        self.AstroPix_ToT_Clock_Freq = 100000000 # 100 MHz
        self.AstroPix_ToT_Max_Clock = 2**12
        
        self.FPGA_Clock_Freq = 1000000 # 1 MHz
        self.FPGA_Max_Clock = 2**32
        self.FPGA_Clock_Offset = 1e-6 # s
        
        self.FPGA_Rollover_Buffer = 5e6
        self.FPGA_diff_cutoff = 1e6
        
        self.FPGA_readout_cycles = 42
    
    def clean_FPGA_times(self,array):
    
        # This function will search for extra hits that might be related to having a way to high event rate
        # These extra hits have jumbled FPGA timestamps and will be removed
        
        edit_array = array.copy()
        complete_drops = False
    
        while not complete_drops:
    
            FPGA_times = edit_array[:,-1].copy()
            FPGA_time_diffs = np.diff(FPGA_times)
            FPGA_time_diffs_rollovers = np.where(FPGA_time_diffs < (-self.FPGA_Max_Clock + self.FPGA_Rollover_Buffer))[0]+1
    
            for rollover_idx in FPGA_time_diffs_rollovers:
                FPGA_times[rollover_idx:] += self.FPGA_Max_Clock
    
            # Look for outliers
            FPGA_times_diff_up = (np.where(np.diff(FPGA_times) > self.FPGA_diff_cutoff)[0]) # this is just before shooting up
            FPGA_times_diff_down = (np.where(np.diff(FPGA_times) < -self.FPGA_diff_cutoff)[0]) # this is just before shooting down
    
            # Drop outliers
            droppable_idx = []
            for up_idx in FPGA_times_diff_up:
                if up_idx + 1 in FPGA_times_diff_down:
                    droppable_idx.append(up_idx + 1)
                elif up_idx - 1 in FPGA_times_diff_down:
                    droppable_idx.append(up_idx)
            droppable_idx = np.array(droppable_idx)
    
            if len(droppable_idx) == 0:
                complete_drops = True
            else:
                mask = np.ones(len(edit_array),dtype = bool)
                mask[droppable_idx] = False
                edit_array = edit_array[mask,:]
                
        return edit_array
    
    def read_sim(self):
        TKR_hits = []
        for line in open(self.sim_name,'r'):
            if 'ID' == line[:2]:
                eid = float(line.split()[1])
            if 'TI' == line[:2]:
                time = float(line.split()[1])
            if 'HTsim 1' == line[:7]:
                split_line = line.split(';')
                x = float(split_line[1])
                y = float(split_line[2])
                z = float(split_line[3])
                e = float(split_line[4])
                TKR_hits.append([eid,time,x,y,z,e])
        self.TKR_hits = np.array(TKR_hits)
    
    def pinpoint(self):
        # A-STEP is flat in the zx plane
        # -x, -z = Chip 0
        # -x, +z = Chip 1
        # +x, -z = Chip 2
        # +x, +z = Chip 3
        # ^similar setup for the rows/columns
        
        # Assume (x,z) = (0,0) is in the middle of a layer
        # Look for the smallest x,z and call those x0,z0
        # 
        
        xs = self.TKR_hits[:,2]
        ys = self.TKR_hits[:,3]
        zs = self.TKR_hits[:,4]
        
        Chip_Xs = (xs-self.x0)//self.Chip_Spacing_X
        Chip_Zs = (zs-self.z0)//self.Chip_Spacing_Z
        
        self.Chip_IDs = (self.N_Chip_Zs*Chip_Xs + Chip_Zs).astype(int)
        
        # Now look for row & column number within the pixel
        self.rows = np.round((xs-self.x0)%self.Chip_Spacing_X/self.pixel_size_X,5).astype(int)
        self.cols = np.round((zs-self.z0)%self.Chip_Spacing_Z/self.pixel_size_Z,5).astype(int)
        
        # Now look for layer number
        self.Layer_IDs = np.round((ys - self.Layer_Offset_Y)/self.Layer_Spacing_Y,5).astype(int)
    
    def RevCal(self):
        ToT_us = np.zeros_like(self.TKR_hits[:,5])
        ToT_us_row_smear = np.zeros_like(ToT_us)
        ToT_us_col_smear = np.zeros_like(ToT_us)
        
        subthresh_row = np.ones_like(self.TKR_hits[:,5])
        subthresh_col = np.ones_like(self.TKR_hits[:,5])

        
        calib_f = h5.File(self.calib_name,'r')
        calib_layers = calib_f.keys()
        
        if self.seed >= 0:
            rng = np.random.default_rng(seed = self.seed)
        else:
            rng = np.random.default_rng()
        
        for layer_str in calib_layers:
            layer_group = calib_f[layer_str]
            calib_chips = layer_group.keys()
            
            for chip_str in calib_chips:
                chip_group = layer_group[chip_str]
                
                layer = int(layer_str.split('Layer')[1][0])
                chip  = int(chip_str.split('Chip')[1][0])
        
                calib_array = chip_group['Calibration'][...] # row, col, slope, intercept, should be calibrated to row ToT
                res_array = chip_group['Resolution'][...] # row, col, slope, intercept, minimum in ToT_us
                try:
                    offset_array = chip_group['Offset'][...] # row, col, <col minus row in us>
                except:
                    offset_array = np.zeros((len(calib_array),3))
                    offset_array[:,:2] = calib_array[:,:2]
                    
                try:
                    thresh_array = chip_group['Threshold'][...] # row, col, thresh in us
                except:
                    thresh_array = np.zeros((len(calib_array),3))
                    thresh_array[:,:2] = calib_array[:,:2]
        
                cond1 = (self.Layer_IDs == layer) & (self.Chip_IDs == chip)
                
                for i in range(len(calib_array)):
                    i_row = calib_array[i,0]
                    i_col = calib_array[i,1]
                    
                    cond2 = (self.rows == i_row) & (self.cols == i_col)
                    
                    cond = cond1 & cond2
                    
                    ToT_us[cond] = (self.TKR_hits[cond,5] - calib_array[i,3])/calib_array[i,2] # should be in us now
                    ToT_sigma = res_array[i,2]*ToT_us[cond] + res_array[i,3]
                    ToT_sigma[ToT_sigma < res_array[i,4]] = res_array[i,4]
                    ToT_us_row_smear[cond] = rng.normal(loc = ToT_us[cond], scale = ToT_sigma)
                    ToT_us_col_smear[cond] = ToT_us_row_smear[cond] + offset_array[i,2]
                    
                    # check for threshold
                    subthresh_row[cond] = ToT_us_row_smear[cond] > thresh_array[i,2]
                    subthresh_col[cond] = ToT_us_col_smear[cond] > thresh_array[i,2]
                    
        calib_f.close()
        
        ToT_us_row_smear = ((ToT_us_row_smear*100)//1)/100
        ToT_tot_row = (ToT_us_row_smear / 1e6 * self.AstroPix_ToT_Clock_Freq).astype(int) # should be in AstroPix ToT clock units
        
        ToT_us_col_smear = ((ToT_us_col_smear*100)//1)/100
        ToT_tot_col = (ToT_us_col_smear / 1e6 * self.AstroPix_ToT_Clock_Freq).astype(int) # should be in AstroPix ToT clock units
        
        self.ToT_tot_row = ToT_tot_row % self.AstroPix_ToT_Max_Clock
        self.ToT_us_row = self.ToT_tot_row * 1e6 / self.AstroPix_ToT_Clock_Freq
        self.ToT_lsb_row = self.ToT_tot_row%2**8
        self.ToT_msb_row = (self.ToT_tot_row - self.ToT_lsb_row) >> 8
        
        self.ToT_tot_col = ToT_tot_col % self.AstroPix_ToT_Max_Clock
        self.ToT_us_col = self.ToT_tot_col * 1e6 / self.AstroPix_ToT_Clock_Freq
        self.ToT_lsb_col = self.ToT_tot_col%2**8
        self.ToT_msb_col = (self.ToT_tot_col - self.ToT_lsb_col) >> 8
        
        self.subthresh_row = subthresh_row
        self.subthresh_col = subthresh_col
        
    def get_clock_times(self):
        self.AstroPix_times = ((self.TKR_hits[:,1]*self.AstroPix_Clock_Freq)%self.AstroPix_Max_Clock).astype(int)
        
        self.FPGA_row_times = (((self.TKR_hits[:,1] + 1e-6*self.ToT_us_row + self.FPGA_Clock_Offset)*self.FPGA_Clock_Freq)%self.FPGA_Max_Clock).astype(int)
        self.FPGA_col_times = (((self.TKR_hits[:,1] + 1e-6*self.ToT_us_col + self.FPGA_Clock_Offset)*self.FPGA_Clock_Freq)%self.FPGA_Max_Clock).astype(int)
        
    def make_out_array(self):
        self.out_header = 'dec_ord,readout,layer,chipID,payload,location,isCol,timestamp,tot_msb,tot_lsb,tot_total,tot_us,fpga_ts'
        out_array = np.zeros([2*len(self.TKR_hits),13])
        out_subthresh = np.zeros(2*len(self.TKR_hits))
        
        # dec_ord & payload are constant
        out_array[:,0] = np.zeros(2*len(self.TKR_hits)) # dec_ord
        out_array[:,4] = 4*np.ones(2*len(self.TKR_hits)) # payload
        
        for i in range(len(self.TKR_hits)):
            
            # dec_ord already done
            
            # readout
            out_array[2*i  ,1] = 0
            out_array[2*i+1,1] = 0
            
            # layer
            out_array[2*i  ,2] = self.Layer_IDs[i]
            out_array[2*i+1,2] = self.Layer_IDs[i]
            
            # chip ID
            out_array[2*i  ,3] = self.Chip_IDs[i]
            out_array[2*i+1,3] = self.Chip_IDs[i]
            
            # payload already done
            
            # location
            out_array[2*i  ,5] = self.rows[i]
            out_array[2*i+1,5] = self.cols[i]
            
            # location
            out_array[2*i  ,6] = 0 # row
            out_array[2*i+1,6] = 1 # col
            
            # timestamp
            out_array[2*i  ,7] = self.AstroPix_times[i]
            out_array[2*i+1,7] = self.AstroPix_times[i] 
            
            # tot_msb
            out_array[2*i  ,8] = self.ToT_msb_row[i]
            out_array[2*i+1,8] = self.ToT_msb_col[i]
            
            # tot_lsb
            out_array[2*i  ,9] = self.ToT_lsb_row[i]
            out_array[2*i+1,9] = self.ToT_lsb_col[i]
            
            # tot_tot
            out_array[2*i  ,10] = self.ToT_tot_row[i]
            out_array[2*i+1,10] = self.ToT_tot_col[i]
            
            # tot_us
            out_array[2*i  ,11] = self.ToT_us_row[i]
            out_array[2*i+1,11] = self.ToT_us_col[i]
            
            # FPGA timestamp
            out_array[2*i  ,12] = self.FPGA_row_times[i]
            out_array[2*i+1,12] = self.FPGA_col_times[i]
            
            # sub-threshold hits
            out_subthresh[2*i  ] = self.subthresh_row[i]
            out_subthresh[2*i+1] = self.subthresh_col[i]
            
        # Remove sub-threshold hits
        out_array = out_array[out_subthresh == 1]
        
        self.out_array = out_array.copy()
        
    def write_output(self):
        if self.is_h5:
            out_file = h5.File(self.out_name,'w')
            out_file.create_dataset('Data',data = self.out_array)
            out_file.create_dataset('Column_Names',data = self.out_header.split(','))
            out_file.close()
        
        else:
            out_file = open(self.out_name,'w')
            out_file.write(self.out_header + '\n')
            for row in self.out_array:
                string = ""
                for i,el in enumerate(row):
                    if i == 11:
                        string = string + f'{el:.2f},'
                    elif i == 12:
                        string = string + f'{int(el)}\n'
                    else:
                        string = string + f'{int(el)},'
                out_file.write(string)
            out_file.close()
    
    def process(self):
        self.read_sim()
        self.pinpoint()
        self.RevCal()
        self.get_clock_times()
        self.make_out_array()
        self.write_output()