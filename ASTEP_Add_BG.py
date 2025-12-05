import numpy as np
import h5py as h5

"""
This function applies the reverse calibration to A-STEP simulations.

Steps
    1: Define file names and parameters
    2: Read .csv background file
    3: Clean FPGA timestamps
    4: Combine background array & simulated source array
    5: Sort by FPGA timestamp
    6: Write to output file

"""

class ASTEP_Add_BG:
    def __init__(self,ARC,BG_name,is_h5):
        self.ARC        = ARC
        self.BG_name    = BG_name
        self.is_h5      = is_h5
        
        if is_h5:
            self.out_name = ARC.sim_name + '.ASTEP_wBG.h5'
        else:
            self.out_name = ARC.sim_name + '.ASTEP_wBG.csv'
                
    def read_BG(self):
        # Read in BG file
        output = []
        first_line = True
        for line in open(self.BG_name, 'r'):
            if not first_line:
                split_line = line.split(',')
                if (split_line[6].strip() == 'True') or (split_line[6].strip() == '1'):
                    split_line[6] = '1'
                elif (split_line[6].strip() == 'False') or (split_line[6].strip() == '0'):
                    split_line[6] = '0'
                else:
                    print('Bad is_col entry')
                
                output.append(split_line)
            first_line = False
        output = np.array(output,dtype = float)
        output = output[output[:,4] == 4] # filter on payload == 4
        self.BG_array = output
    
    def sort_FPGA_times(self):
        # Sort BG FGPA times
        self.BG_array = self.ARC.clean_FPGA_times(self.BG_array)
        
        BG_FPGA_times = self.BG_array[:,-1]
        BG_FPGA_time_diffs = np.diff(BG_FPGA_times)
        BG_FPGA_time_diffs_rollovers = np.where(BG_FPGA_time_diffs < (-self.ARC.FPGA_Max_Clock + self.ARC.FPGA_Rollover_Buffer))[0]+1
        
        BG_FPGA_times_corrected = BG_FPGA_times.copy()
        for rollover_idx in BG_FPGA_time_diffs_rollovers:
            BG_FPGA_times_corrected[rollover_idx:] += self.ARC.FPGA_Max_Clock
        
        # Set clock to start at 0
        self.BG_FPGA_times_corrected = BG_FPGA_times_corrected - BG_FPGA_times_corrected[0]
        
        drop_BG_mask = (self.BG_FPGA_times_corrected > 0)
        self.BG_FPGA_times_corrected = self.BG_FPGA_times_corrected[drop_BG_mask]
        self.BG_array = self.BG_array[drop_BG_mask]
        
        self.BG_array[:,12] = (self.BG_FPGA_times_corrected%self.ARC.FPGA_Max_Clock)
        
        self.ARC.out_array = self.ARC.clean_FPGA_times(self.ARC.out_array)
        
        out_FPGA_times = self.ARC.out_array[:,-1]
        out_FPGA_time_diffs = np.diff(out_FPGA_times)
        out_FPGA_time_diffs_rollovers = np.where(out_FPGA_time_diffs < (-self.ARC.FPGA_Max_Clock + self.ARC.FPGA_Rollover_Buffer))[0]+1
        
        out_FPGA_times_corrected = out_FPGA_times.copy()
        for rollover_idx in out_FPGA_time_diffs_rollovers:
            out_FPGA_times_corrected[rollover_idx:] += self.ARC.FPGA_Max_Clock
        
        self.out_FPGA_times_corrected = out_FPGA_times_corrected
        self.ARC.out_array[:,12] = (out_FPGA_times_corrected%self.ARC.FPGA_Max_Clock)
        
        # Get maximum FPGA timestamp time
        self.max_FPGA_time = min(max(self.BG_FPGA_times_corrected),max(self.out_FPGA_times_corrected).astype(int))
        
    def combine_arrays(self):

        combined_array = np.zeros((len(self.ARC.out_array) + len(self.BG_array),np.shape(self.ARC.out_array)[1]))
        combined_array[:len(self.ARC.out_array)] = self.ARC.out_array
        combined_array[len(self.ARC.out_array):] = self.BG_array
        
        combined_times = np.zeros((len(self.ARC.out_array) + len(self.BG_array)))
        combined_times[:len(self.ARC.out_array)] = self.out_FPGA_times_corrected
        combined_times[len(self.ARC.out_array):] = self.BG_FPGA_times_corrected
        
        combined_array = combined_array[combined_times < self.max_FPGA_time,:]
        combined_times = combined_times[combined_times < self.max_FPGA_time]
        
        sort_order = np.argsort(combined_times)
        self.combined_array_sorted = combined_array[sort_order,:]
        
        combined_time_sorted = combined_times[sort_order]

        
    def write_output(self):
        if self.is_h5:
            out_file = h5.File(self.out_name,'w')
            out_file.create_dataset('Data',data = self.combined_array_sorted)
            out_file.create_dataset('Column_Names',data = self.ARC.out_header.split(','))
            out_file.close()
        
        else:
            out_file = open(self.out_name,'w')
            out_file.write(self.ARC.out_header + '\n')
            for row in self.combined_array_sorted:
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
        self.read_BG()
        self.sort_FPGA_times()
        self.combine_arrays()
        self.write_output()