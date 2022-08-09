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
  
def apply_obs_constraints(df,
                          good_arrays,
                          greentime=True,
                          start_dt=None,
                          end_dt=None
                          ):
 
    array_constraint = df['array'].isin(good_arrays)
    greentime_constraint = df['className']=='Green Time'
    
    overall_constraint = array_constraint & greentime_constraint
    if start_dt is not None:
        overall_constraint = overall_constraint & (df['end'] > start_dt)
    if end_dt is not None:
        overall_constraint = overall_constraint & (df['start'] < end_dt)
    
    good_df = df[overall_constraint]
    
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
    
    
def get_AltAz(df, coord, min_obs_length, horizon=12*u.deg):
    obs_frame_start = AltAz(obstime=df['start'], location=ATCA)
    obs_frame_end = AltAz(obstime=df['end'], location=ATCA)
    
    altaz_start = coord.transform_to(obs_frame_start)
    altaz_end = coord.transform_to(obs_frame_end)
    
    earliest_start = df['start'].copy()
    below_horizon = altaz_start.alt < horizon
    
    rise_times = ATCA_Obs.target_rise_time(earliest_start[below_horizon], coord, which='next', horizon=horizon).utc.isot.astype('datetime64')
    
    earliest_start.iloc[np.where(below_horizon)] = pd.to_datetime(rise_times, utc=True)
    
    
    df['obs_start'] = earliest_start
    
    latest_end = df['end'].copy()
    below_horizon = altaz_end.alt < horizon
    set_times = ATCA_Obs.target_set_time(latest_end[below_horizon], coord, which='previous', horizon=horizon).utc.isot.astype('datetime64')
    
    latest_end.iloc[np.where(below_horizon)] = pd.to_datetime(set_times, utc=True)
    
    df['obs_end'] = latest_end
    df['obs_length'] = df['obs_end']-df['obs_start']
    
    lengthy = df['obs_length'] >= min_obs_length
    obs_blocks = df[lengthy][['obs_start','obs_end', 'obs_length','array']]
    num_obs = len(obs_blocks)
    
    if num_obs > 0:
        print(obs_blocks[['obs_start','obs_length']])
        
        return True
        
    return False


def analyse_sources(sc=None,
                    filename=None,
                    sched_filepath=None,
                    semester=None,
                    good_arrays=['6A','6B','6C','6D','1.5A','1.5B','1.5C','1.5D'],
                    start_dt=pd.Timestamp("today", tz='UTC'),
                    end_dt=None,
                    min_obs_length=pd.Timedelta(1,'h')
                    ):
    if sc is not None and filename is not None:
        print("Received skycoords and catalogue - pick one!")
    elif sc is None and filename is None:
        print("Must supply a skycoord or catalogue!")
        
    if filename is not None:
        sources = pd.read_csv(filename)
        units = (u.deg, u.deg)
        if ':' in sources['ra'].iloc[0]:
            units = (u.hourangle, u.deg)
        sc = SkyCoord(sources['ra'],sources['dec'], unit=units)
    
    if sched_filepath is None:
        sched_filepath = get_schedule(semester)
    full_schedule = load_schedule(sched_filepath)
    good_schedule = apply_obs_constraints(full_schedule,
                                          good_arrays,
                                          start_dt,
                                          end_dt)
    
    num_observable = 0
    for i,coord in enumerate(sc):
        if filename is not None:
            print(sources['name'].iloc[i])
        else:
            print(f"Source {i}")
        observable = get_AltAz(good_schedule,coord, min_obs_length)
        print()
        
        num_observable += observable
    print('{} sources are observable'.format(num_observable))


#analyse_sources(filename='sn2012dy.csv', semester='2022Apr')

sc = SkyCoord(['06:23:09.28'], ['-04:56:22.79'], unit=(u.hourangle, u.deg))
analyse_sources(sc, semester='2022Apr', good_arrays=['H168'])
