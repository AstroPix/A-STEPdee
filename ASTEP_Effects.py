import numpy as np
import h5py as h5

"""
This function applies the reverse calibration to A-STEP simulations.

Steps
    1: Define file names and parameters
    2: Check for coincident entries
    3: Separate coincident entries by ARC.FPGA_readout_cycles
    4: Define output array
    5: Write to output file

"""

class ASTEP_Effects:
    def __init__(self,sim_name,in_array,is_h5,with_BG,ARC):
        self.sim_name   = sim_name
        self.is_h5      = is_h5
        self.in_array   = in_array
        
        if is_h5:
            if with_BG:
                self.out_name = sim_name + '.ASTEP_wBG_wEff.h5'
            else:
                self.out_name = sim_name + '.ASTEP_wEff.h5'
        else:
            if with_BG:
                self.out_name = sim_name + '.ASTEP_wBG_wEff.csv'
            else:
                self.out_name = sim_name + '.ASTEP_wEff.csv'
                
        self.ARC = ARC
    
    def handle_coinc(self,array,times):
    
        new_array = np.zeros_like(array)
        new_times = np.zeros_like(times)
        
        is_col = array[:,6]
        ToT_us = array[:,11]
        FPGA_time = array[:,12]
        start_times = times/self.ARC.FPGA_Clock_Freq - ToT_us*1e-6 # transform to seconds
        
        # do rows first
        rows = np.where(is_col == 0)[0]
        sort_rows = np.argsort(start_times[rows])
        for i,idx in enumerate(sort_rows):
            new_array[i] = array[rows[idx]]
            new_times[i] = times[rows[idx]]
        
        # now do cols
        cols = np.where(is_col == 1)[0]
        sort_cols = np.argsort(start_times[cols])
        sort_cols = sort_cols[::-1]
        for i,idx in enumerate(sort_cols):
            new_array[i+len(sort_rows)] = array[cols[idx]]
            new_times[i+len(sort_rows)] = times[cols[idx]]
        
        for i in range(1,len(times)):
            new_times[i] = new_times[0] + i*42
            new_array[i,-1] = (new_array[0,-1] + i*42)%self.ARC.FPGA_Max_Clock
            
        return new_array, new_times
        
    def sort_FPGA_timestamps(self):
        # Sort FPGA times & handle rollover
        self.in_array = self.ARC.clean_FPGA_times(self.in_array)
        
        in_FPGA_times = self.in_array[:,-1]
        in_FPGA_time_diffs = np.diff(in_FPGA_times)
        in_FPGA_time_diffs_rollovers = np.where(in_FPGA_time_diffs < (-self.ARC.FPGA_Max_Clock + self.ARC.FPGA_Rollover_Buffer))[0]+1
        
        in_FPGA_times_corrected = in_FPGA_times.copy()
        for rollover_idx in in_FPGA_time_diffs_rollovers:
            in_FPGA_times_corrected[rollover_idx:] += self.ARC.FPGA_Max_Clock
        
        self.in_FPGA_times = in_FPGA_times_corrected
        self.in_array[:,12] = (in_FPGA_times_corrected%self.ARC.FPGA_Max_Clock)
        
    def coincidence_hits(self):
        # Look for hits within 42 FPGA clock cycles
        self.out_time = self.in_FPGA_times
        self.out_array = self.in_array.copy()
        
        n_coin = sum(np.diff(self.out_time) < self.ARC.FPGA_readout_cycles)
        
        while n_coin != 0:
        
            print(f'Starting Coin Handling with N = {n_coin}')
        
            start_idx = 0
            end_idx = 1
            while end_idx < len(self.out_time):
                # check that there is a coincidence
                coinc = (self.out_time[end_idx] - self.out_time[start_idx]) < self.ARC.FPGA_readout_cycles
        
                # if there is a coincidence, check if there are multiple coincidences
                if coinc:
                    end_of_cluster = False
                    while not end_of_cluster:
                        if end_idx == len(self.out_time)-1:
                            end_of_cluster = True
                        elif (self.out_time[end_idx+1] - self.out_time[start_idx]) <= self.ARC.FPGA_readout_cycles*(end_idx+1 - start_idx):
                            end_idx += 1
                        else:
                            end_of_cluster = True
                    # Handle
                    #print(f'Indices:   {start_idx}\t{end_idx}')
                    self.out_array[start_idx:end_idx+1,:], self.out_time[start_idx:end_idx+1] = self.handle_coinc(self.out_array[start_idx:end_idx+1,:], self.out_time[start_idx:end_idx+1])
                    #print()
                    
                # move on to the next start_idx
                start_idx = end_idx
                end_idx = start_idx + 1
                
            # Check again
            n_coin = sum(np.diff(self.out_time) < self.ARC.FPGA_readout_cycles)
        
    def write_output(self):
        if self.is_h5:
            out_file = h5.File(self.out_name,'w')
            out_file.create_dataset('Data',data = self.out_array)
            out_file.create_dataset('Column_Names',data = self.ARC.out_header.split(','))
            out_file.close()
        
        else:
            out_file = open(self.out_name,'w')
            out_file.write(self.ARC.out_header + '\n')
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
        self.sort_FPGA_timestamps()
        self.coincidence_hits()
        self.write_output()