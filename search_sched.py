import urllib.request
import numpy as np
import json
import pandas as pd
import astropy.units as u
from astropy.coordinates import Angle, EarthLocation, SkyCoord, AltAz
from astropy.time import Time

from astroplan import Observer

ATCA = EarthLocation(lat=-30.31277778*u.deg, lon=149.55000000*u.deg, height=236.87*u.m)
ATCA_Obs = Observer(location=ATCA, name='ATCA')

def get_schedule(semester, filepath=None):
    sched_url = 'https://www.narrabri.atnf.csiro.au/observing/schedules/%sSem/atca_maint.json'%(semester)
    
    if filepath is None:
        filepath = '%s_atca.json'%(semester)
        
    urllib.request.urlretrieve(sched_url,filepath)
    
    return filepath
    

def load_schedule(filepath):
    df = pd.read_json(filepath)
    
    for time in ['start','end']:
        lst = time+'LST'
        
        df[lst] = Angle(df[lst], unit=u.hourangle)
        df[time] = pd.to_datetime(df[time])
        
        t = Time(df[time], location=ATCA)
        
        df[lst] = t.sidereal_time('apparent')
        
    return df
  
def apply_obs_constraints(df,good_arrays=['6A','6B','6C','6D','1.5A','1.5B','1.5C','1.5D'], greentime=True):
    array_constraint = df['array'].isin(good_arrays)
    greentime_constraint = df['className']=='Green Time'

    good_df = df[array_constraint & greentime_constraint]
    
    return good_df


def get_overlaps(df, lst_range):
    source_rise = Angle(lst_range[0], unit=u.hourangle)
    source_set = Angle(lst_range[1], unit=u.hourangle)
    
    nrows = df.shape[0]
    ones = np.ones(nrows)
    
    obs_start = np.maximum(df['startLST'].values, source_rise.value)
    obs_end = np.minimum(df['endLST'].values, source_set.value)
    obs_length = obs_end-obs_start
    
    sched_times = df.loc[:, ['startLST','endLST']]
    for i in range(nrows):
        length = obs_length[i]
        if length < 3:
            continue
        print(source_rise.value, source_set.value)
        print(sched_times.iloc[i,:].values, obs_start[i],obs_end[i])
        print((obs_end-obs_start)[i])
        print()
    
    
def get_AltAz(df,coord, horizon=12*u.deg):
    obs_frame_start = AltAz(obstime=df['start'], location=ATCA)
    obs_frame_end = AltAz(obstime=df['end'], location=ATCA)
    
    altaz_start = coord.transform_to(obs_frame_start)
    altaz_end = coord.transform_to(obs_frame_end)
    
    earliest_start = df['start'].copy()
    below_horizon = altaz_start.alt < horizon
    
    rise_times = ATCA_Obs.target_rise_time(earliest_start[below_horizon], coord, which='next', horizon=horizon).utc.isot.astype('datetime64')
    
    earliest_start[below_horizon] = rise_times
    
    latest_end = df['end'].copy()
    below_horizon = altaz_end.alt < horizon
    set_times = ATCA_Obs.target_set_time(latest_end[below_horizon], coord, which='previous', horizon=horizon).utc.isot.astype('datetime64')
    
    latest_end[below_horizon] = set_times
    
    df['obs_start'] = earliest_start
    df['obs_end'] = latest_end
    df['obs_length'] = df['obs_end']-df['obs_start']
    
    lengthy = df['obs_length'] >= pd.Timedelta(6,'h')
    obs_blocks = df[lengthy][['obs_start','obs_end', 'obs_length','array']]
    num_obs = len(obs_blocks)
    
    if num_obs > 0:
        print(obs_blocks)
        
        return True
        
    return False


def analyse_sources():
    sources = pd.read_csv('sources.csv')#, encoding = "ISO-8859-1")
    
    print(sources)
    sc = SkyCoord(sources['ra'],sources['dec'], unit=u.deg)
    
    sched_filepath = get_schedule('2019Oct')
    full_schedule = load_schedule(sched_filepath)
    good_schedule = apply_obs_constraints(full_schedule)
    
    num_observable = 0
    for i,coord in enumerate(sc):
        print(sources['name'].iloc[i])
        observable =get_AltAz(good_schedule,coord)
        print()
        
        num_observable += observable
    print(num_observable)
        
        
analyse_sources()
