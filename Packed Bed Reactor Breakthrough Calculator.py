import numpy as np
import pandas as pd
from datetime import datetime
import os, glob, pathlib
'''------------------------------------------------------------------------'''
'''------------------------------------------------------------------------'''
def PBR_Import(file, c_max, stage, Mass_Loss):
    '''--------------------------------------------------------------------'''
    '''Imports data from PBR .csv files.-----------------------------------'''
    '''Returns time data in min/g and concentration data in cc/g.----------'''
    '''file = .csv file from PBR instrument.-------------------------------'''
    '''c_max = value of highest concentration of target gas in ppm.--------'''
    '''stage = value in 'STAGE' column of .csv file; usually 'Adsorption'.-'''
    '''Mass_Loss = % loss of weight from activation; estimate from TGA.----'''
    '''--------------------------------------------------------------------'''
    #Determine index values where STAGE column = stage
    s_ = pd.read_csv(file).filter(like='STAGE',axis=1)
    index = s_.index[s_['STAGE'] == stage].tolist()
    index_length = len(index)
    #Extract ppm data from within stage region and normalize
    d_ = pd.read_csv(file).filter(like='PLC.CO2_ANALYZER_ppm',axis=1)
    data = d_.iloc[index,:]
    data = data.reset_index().iloc[:,1:]
    data = data.dropna()
    data.columns = [0]
    if c_max == None:
        c_c0 = (data-data.min())/(data.max()-data.min())
        c_max = data.max()
    else:
        c_c0 = (data-data.min())/(c_max-data.min()) 
        c_max = c_max
    #Extract the mass of the loaded sample, approximate activated weight
    m_ = pd.read_csv(file).filter(like='Sorbent_Quantity',axis=1)
    mass = m_.iloc[index,:]
    mass = mass.reset_index().iloc[:,1:]
    mass = mass.dropna()
    mass.columns = [0]
    mass = mass.iloc[0,0] - (mass.iloc[0,0]*(Mass_Loss/100))
    #Convert datatime to timestamp and normalize to experiment
    #Reduce experiment time (s) to minutes per g
    t_ = pd.read_csv(file,parse_dates=['Time'],
                     date_parser=lambda x: pd.to_datetime(
                         x, format='%m/%d/%Y %H:%M:%S')
                    ).filter(like='Time',axis=1)
    time = t_.iloc[index,:]
    time = time.reset_index().iloc[:,1:]
    time = time.dropna()
    time = time.astype(np.int64) // 10**9
    time.columns = [0]
    reduced_time = []
    for i in range(index_length):
        j = time.iloc[i,:] - time.iloc[0,:]
        j = (j/60)/mass
        reduced_time.append(j)
    m_g = pd.DataFrame(reduced_time)
    m_g = m_g.reset_index().iloc[:,1:]
    return(m_g, c_c0, c_max)
'''------------------------------------------------------------------------'''
'''------------------------------------------------------------------------'''
def Axial_Dispersion(m_g, c_c0):
    '''--------------------------------------------------------------------'''
    '''Finds stoichiometric time i.e. where C/C0 = 0.5 and subtracts this--'''
    '''from all other time points to compare changes in mass transfer and--'''
    '''evaluate axial dispersion effects.----------------------------------'''
    '''Takes m_g and c_c0 from PBR import as data.-------------------------'''
    '''--------------------------------------------------------------------'''
    #Find time where c_c0=0.5, subtract this time from all other time points
    c_c0_round = c_c0.round(1)
    index = c_c0_round.index[c_c0_round[0] == 0.5].values
    AxD = m_g - m_g.iloc[index[0],:]
    return(AxD)
'''------------------------------------------------------------------------'''
'''------------------------------------------------------------------------'''
def Tau_Break(m_g, c_c0):
    '''--------------------------------------------------------------------'''
    '''--------------------------------------------------------------------'''
    b = c_c0[0].where(c_c0[0].between(0.4,0.6))
    b = b.dropna()
    index = b.index
    a = m_g.iloc[index,:].values.flatten().tolist()
    b = b.values.tolist()
    slope, intercept = np.polyfit(a, b, 1)
    tau_break = (-intercept)/slope
    return(tau_break)
'''--------------------------------------------------------------------'''
def Breakthrough_Capacity(total_flow_rate, c_max, tau_break):
    gas_flow = total_flow_rate*(c_max/1000000)
    breakthrough_capacity = gas_flow * tau_break
    return(gas_flow, breakthrough_capacity)
'''--------------------------------------------------------------------'''
def Dimensionless_Breakthrough(m_g, c_c0, gas_flow, c_max):
    total_area = m_g.iloc[-1,0]
    cumulative_trapezoid = []   
    for i in range(len(c_c0)-1):
        trapezoid = (
            (c_c0.iloc[i+1,0]+c_c0.iloc[i,0])/2
        )*(
            (m_g.iloc[i+1,0]-m_g.iloc[i,0])
        )
        cumulative_trapezoid.append(trapezoid)       
    peak_area = sum(cumulative_trapezoid)
    dimensionless_breakthrough = total_area - peak_area
    equilibrium_capacity = dimensionless_breakthrough*gas_flow
    return(dimensionless_breakthrough, equilibrium_capacity)
'''--------------------------------------------------------------------'''   

def Breakthrough_Analysis(path, c_max=None, stage='Adsorption', Mass_Loss=0, total_flow_rate=400):
    '''--------------------------------------------------------------------'''
    
    
    '''--------------------------------------------------------------------'''
    outdir = path+'results\\'
    if not os.path.exists(outdir):
        os.mkdir(outdir)
     
    
    file = glob.glob(path + '*csv')
    File_List = []
    for i in file:
        File_List.append(i)
        
    Time = pd.DataFrame()
    Concentration = pd.DataFrame()
    Axial_Disp = pd.DataFrame()
    
    Breakthrough_Time = []
    Dimensionless_Breakthrough_Time = []
    Breakthrough_Loading = []
    Equilibrium_Loading = []
    
    
    for i in File_List:
        m_g, c_c0, c_max = PBR_Import(file=i, c_max=c_max, stage=stage, Mass_Loss=Mass_Loss)
        AxD = Axial_Dispersion(m_g,c_c0)
        tau_break = Tau_Break(m_g, c_c0)
        gas_flow, breakthrough_capacity = Breakthrough_Capacity(total_flow_rate, c_max, tau_break)
        dimensionless_breakthrough, equilibrium_capacity = Dimensionless_Breakthrough(m_g, c_c0, gas_flow, c_max)
        
        Time = pd.concat([Time, m_g], ignore_index=True, axis=1)
        Concentration = pd.concat([Concentration, c_c0], ignore_index=True, axis=1)
        Axial_Disp = pd.concat([Axial_Disp, AxD], ignore_index=True, axis=1)
        
        Breakthrough_Time.append(tau_break)
        Dimensionless_Breakthrough_Time.append(dimensionless_breakthrough)
        Breakthrough_Loading.append(breakthrough_capacity)
        Equilibrium_Loading.append(equilibrium_capacity)
        print('Analysis Complete on:', i)
    
    Time.columns = np.arange(0,len(Time.T))
    Concentration.columns = np.arange(0,len(Concentration.T))
    Axial_Disp.columns = np.arange(0,len(Axial_Disp.T))
    
    Results = {'Breakthrough Time / min/g' : Breakthrough_Time, 
               'Dimensionless Breakthrough Time / a.u.': Dimensionless_Breakthrough_Time, 
               'Breakthrough Capacity / cc/g' : Breakthrough_Loading, 
               'Equilibrium Capacity / cc/g' : Equilibrium_Loading}
    
    Results = pd.DataFrame(Results)
    

    Time.to_csv(outdir+'Time_min-per-g.csv')
    Concentration.to_csv(outdir+'Concentration_C-over-C0.csv')
    Axial_Disp.to_csv(outdir+'Axial_Dispersion_Time.csv')
    Results.to_csv(outdir+'Results.csv')

    
    return(Time, Concentration, Axial_Dispersion, Results)
