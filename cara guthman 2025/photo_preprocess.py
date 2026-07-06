import glob, os, scipy, math, pickle, ast, re, sys
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
import statistics as stat
import statsmodels.api as sm
import tmac.models as tm
from datetime import datetime
import seaborn as sns 
from matplotlib.patches import Rectangle
import matplotlib as mpl

def locate_data(dir): 
    # this function locates the data files in the directory. each experiment has multiple files and file types associated with them. 
    file_pattern = r'(.*)_(ts|pd)\.csv|(.*)\.csv|(.*)\.txt'
        # Scan the date directory for files
    for filename in os.listdir(dir):
        match = re.match(file_pattern, filename)
        if match:
            if match.group(2) == 'ts':
                time_file = filename
            elif match.group(2) == 'pd':
                power_file = filename
            elif match.group(3):
                data_file = filename
            elif match.group(4):
                metadata_file = filename
    try:
        file_dict = {'time_file': time_file, 'power_file': power_file, 'data_file': data_file, 'metadata_file': metadata_file}
    except:
        pass
    finally: 
        file_dict = {'time_file': time_file, 'data_file': data_file, 'metadata_file': metadata_file}
    return file_dict

def check_output(metadata_file):
    animal_id =  metadata_file.split('_')[1]
    condition = metadata_file.split('_')[0]
    date = metadata_file.split('_')[2].split('.')[0]
    file_name = '{}_{}_{}_tmacNeural.npy'.format(condition, animal_id, date)S
    out_dir = r'\\cup\falkner\Sae\seeER\outputs\preprocessed'
    out_dir = os.path.join(out_dir, animal_id, date,file_name)
    return out_dir

#extract metadata, date, rig, animal_id, condition, rois
def extract_metadata(metadata_file):
    print('extracting metadata from ', os.getcwd())
    txt = open(metadata_file, 'r')
    lines = txt.readlines()
    metadata = {}
    for line in lines: 
        key, value = line.strip().split(':')
        metadata[key] = value
    metadata['roi_ap_labels'] = ast.literal_eval(metadata['roi_ap_labels'])
    metadata['roi_perm'] = ast.literal_eval(metadata['roi_perm'])
    metadata['zeitgeber_time'] = float(metadata['zeitgeber_time'])
    metadata['rig'] = ast.literal_eval(metadata['rig'])
    rois = [None] * len(metadata['roi_perm'])
    x = 0
    for number in metadata['roi_perm']: 
        rois[number] = metadata['roi_ap_labels'][x]
        x = x + 1
    #extract other metadata
    animal_id =  metadata_file.split('_')[1]
    condition = metadata_file.split('_')[0]
    date = metadata_file.split('_')[2].split('.')[0]
    #organize into two dictionaries, data_dict carries data and identifirs, metadata carries all metadata minus the data itself
    if 'estrous' in metadata: 
        metadata['estrous'] = ast.literal_eval(metadata['estrous'])
        data_dict = {'animal_id': animal_id, 'condition': condition, 'date': date, 'estrous': metadata['estrous']}
    elif 'inj_time' in metadata: 
        metadata['inj_time'] = float(metadata['inj_time'])
        data_dict = {'animal_id': animal_id, 'condition': condition, 'date': date, 'inj_time': metadata['inj_time']}
    else:
        data_dict = {'animal_id': animal_id, 'condition': condition, 'date': date}
    for key in data_dict.keys():
        metadata[key] = data_dict[key]
    data_dict['rois'] = rois
    return metadata, data_dict

def extract_indices(data_file):
    ''' load data file and adjust for any dropped frame '''
    #read the data
    data = pd.read_csv(data_file)
    headers = data.columns.tolist()
    frames = data['FrameCounter'].to_numpy()
    led = data['LedState'].to_numpy()
    #pull out the index of the frames which correspond to LED-blue (2) or LED-green (4)
    #the order is 4 0 0 0 2
    led_green = np.array(frames[led == 4])
    led_blue = np.array(frames[led == 2])
    print('blue frames: {}, green frames: {}'.format(len(led_blue), len(led_green)))
    #edit the frames we pull from just in case there are any dropped frames
    #here, we check if there's an extra frame taken at the beginning (2 comes before the first 4)
    if led_blue[0] < led_green[0]: 
        led_blue = led_blue[1:]
        print('the first blue frame was cut')
    #here, we check if there's an extra frame at the end, and if there are any dropped frames
    missing_frames = np.setdiff1d(led_green, led_blue -4)
    othermissing_frames = np.setdiff1d(led_blue, led_green +4)
    if len(missing_frames) > 0:
        print('led_blue frame is shorter, and led_green has extra frame(s)', missing_frames)
        for i in missing_frames:
            try:
                led_green = led_green.tolist()
            except:
                pass
            led_green.remove(i)
    if len(othermissing_frames) > 0:
        print('there are missing blue frames,', othermissing_frames)
        for i in othermissing_frames: 
            try:
                led_blue = led_blue.tolist()
            except:
                pass
            led_blue.remove(i)
    print('led_blue is now:', len(led_blue), 'led_green is now:', len(led_green))
    return data, led_blue, led_green, headers

def extract_time(time_file, data_file, metadata, data_dict, led_blue):
    rig = metadata['rig']
    print('rig is ursula:', rig == 'Ursula')
    # timefiles are saved slightly differently on different bonsai versions
    if rig == 'Ursula': 
        ts = pd.read_csv(time_file, header=None)
        times = ts[1].to_numpy() /1000 #convert timestamps to seconds from ms
    else: 
        data = pd.read_csv(data_file)
        times = data['ComputerTimestamp']/1000  
    if isinstance(times, np.ndarray):
        signal_time = times[led_blue]
    else:
        signal_time = times.to_numpy()[led_blue]
    #we convert the times into zeigtgeiber time (where the 24 hour period is shifted to whenever light and dark is)
    start_time = metadata['zeitgeber_time']
    seconds = start_time * 3600
    adj_timestamp = (signal_time - seconds)/3600
    for i in range(len(adj_timestamp)):
        if adj_timestamp[i] < 0:
            adj_timestamp[i] = adj_timestamp[i] + 24 
        else: 
            adj_timestamp[i] = adj_timestamp[i]
    #extract the indices for when the light goes on and off for the mouse
    light_off = scipy.signal.find_peaks(adj_timestamp)[0] #active time
    if len(light_off)>= 2: 
        light_on = (light_off - (light_off[1] - light_off[0])/2).astype(int)
    else: 
        light_on = light_off - 12*12
        print(light_on, ":less than a full cycle")
    #inj_hours = how many hours after your start time did you inject your hormone?
    if 'inj_time' in metadata:
        inj_idx = float(metadata['inj_time'])*12 #hours injected after start time
        inj_label = '{0:.1f}'.format(float(metadata['inj_time']))
        data_dict['inj_idx'] = inj_idx
        data_dict['inj_label'] = inj_label
    else: 
        pass
    data_dict['timedata'] = adj_timestamp
    data_dict['raw_time'] = signal_time
    data_dict['light_off'] = np.sort(light_off)
    metadata['start_time'] = signal_time[0]
    try:
        if light_on[0] < 0: 
            light_on[0] = 0
        else:
            light_on[0] = min(range(len(adj_timestamp)), key=lambda i: abs(adj_timestamp[i]-12))
        data_dict['light_on'] = np.sort(light_on)
    except:
        print('this data set is shroter than a full 24 hours...')
    return data_dict, metadata

def extract_ledpower(power_file, metadata):
    import pandas as pd
    #extract the power_file
    data = pd.read_csv(power_file)
    headers = data.columns.tolist()
    power = [data['PD415'].to_numpy()]
    power = np.append(power, [data['PD470'].to_numpy()], axis = 0)
    power = np.append(power, [data['PD560'].to_numpy()], axis = 0)
    #save onto metadata for now
    metadata['raw_led_power'] = power
    #extract the peak data
    return metadata

def led_correlation(metadata, data_dict):
    from scipy.stats import pearsonr
    rois = data_dict['rois']
    data = data_dict['corrected_signal']
    led_power = metadata['raw_led_power']
    peak_array = []
    properties_array = []
    blue_led_power = led_power[1,:]
    green_led_power = led_power[2,:]
    print('testing correlation with 470nm')
    peaks, properties = scipy.signal.find_peaks(green_led_power, distance = 100)
    peak_array.append(peaks)
    properties_array.append(properties)
    peak_values = green_led_power[peaks]
    peaks_dff = (peak_values - np.mean(peak_values))/np.mean(peak_values)
    for i, x  in zip(range(len(rois)), rois):
        corr, p_value = pearsonr(data[i,:], peaks_dff[1:])
    if p_value < 0.05: 
        print('{} is significantly correlated with a p_value of {}'.format(x, p_value))
        print('r value = {:.3f} '.format(corr))
    else: 
        print('{} is not correlated'.format(x))
    print('testing correlation with 560nm')
    peaks, properties = scipy.signal.find_peaks(blue_led_power, distance = 100)
    peak_array.append(peaks)
    properties_array.append(properties)
    peak_values = blue_led_power[peaks]
    peaks_dff = (peak_values - np.mean(peak_values))/np.mean(peak_values)
    for i, x  in zip(range(len(rois)), rois):
        corr, p_value = pearsonr(data[i,:], peaks_dff[1:])
    if p_value < 0.05: 
        print('{} is significantly correlated with a p_value of {}'.format(x, p_value))
        print('r value = {:.3f} '.format(corr))
    else: 
        print('{} is not correlated'.format(x))
    return None


def split_channels(metadata, data_dict, data, headers, led_blue, led_green):
    #extract photometry data 
    rig = metadata['rig']
    rois = data_dict['rois']
    GFP_dict = {}
    RFP_dict = {}
    x = 0
    if rig == 'Ursula': 
        print('using ursula index')
        for j in rois: 
            GFP_dict[j] = data[headers[8+2*x]].to_numpy()[led_blue]
            RFP_dict[j] = data[headers[9+2*x]].to_numpy()[led_green]
            x = x+1
    else: 
        print('using oogway index')
        for j in rois:
            #pilot data headers are in a different order from newest data 
            GFP_dict[j] = data[headers[16 + x]].to_numpy()[led_blue]
            RFP_dict[j] = data[headers[4 + x]].to_numpy()[led_green]
            x = x + 1
    #convert this dictionary into a 12 x whatever length array
    raw_data = []
    rfp_data = []
    for i in metadata['roi_ap_labels']:
        raw_data.append(GFP_dict[i])
        rfp_data.append(RFP_dict[i])
    raw_data.extend(rfp_data)
    raw_data = np.vstack(raw_data)
    data_dict['raw_data'] = raw_data
    #extract dark regions 
    if 'G25' in headers:
        dark_gfp = []
        dark_rfp = []
        for i in range(2):
            dark_gfp.append(data[headers[29+2*i]].to_numpy()[led_blue])
            dark_rfp.append(data[headers[28+2*i]].to_numpy()[led_green])
        dark_rois = np.vstack((dark_gfp, dark_rfp))
        means = np.array(np.mean(dark_rois, axis = 1))
        detrend = scipy.signal.detrend(dark_rois) + means[:,np.newaxis]
        data_dict['dark_channels'] = detrend
        data_dict['dark_rois'] = ['green_in', 'green_out', 'red_in', 'red_out']
    else: 
        pass
    return data_dict

def linear_detrend(data_dict): 
    raw_data = data_dict['raw_data']
    means = np.array(np.mean(raw_data, axis = 1))
    detrend = scipy.signal.detrend(raw_data) + means[:,np.newaxis]
    data_dict['detrended_data'] = detrend
    return data_dict

def run_tmac(data_dict, metadata):
    if 'inj_time' in metadata.keys():
        data = data_dict['raw_data'] 
    else:
        data = data_dict['detrended_data']
    rois = data_dict['rois']
    motion_artifact = []
    corrected_signal = []
    for i in range(len(rois)):
        tmac_variables = tm.tmac_ac(data[i + 12], data[i])
        motion_artifact.append(tmac_variables['m'].flatten())
        corrected_signal.append(tmac_variables['a'].flatten())
    corrected_signal.extend(motion_artifact)
    corrected_signal = np.vstack(corrected_signal)
    data_dict['corrected_signal'] = corrected_signal
    print('tmac has run')
    return data_dict

def save_data(data_dict, metadata, save_dir):
    condition = data_dict['condition']
    animal_id = data_dict['animal_id']
    date = data_dict['date']
    #save data
    filename = '/{}_{}_{}_datadict.npy'.format(condition, animal_id, date)
    #check that folder exists 
    if not os.path.exists(data_dir):
        print('data folder does not exist. creating directory')
        os.makedirs(save_dir)
    # if not os.path.exists(data_dir + '/{}'.format(date)):
    #     print('date folder does not exist... creating directory now...')
    #     os.makedirs(data_dir + '/{}'.format(date))
    data_dir = '//cup/falkner/Sae/seeER/outputs/preprocessed/{}'.format(animal_id)     
    np.save(os.path.join(save_dir, filename), data_dict, allow_pickle = True)
    #check that the file saved properly
    if os.path.isfile(data_dir + filename) == True: 
        print('{} has been saved onto {}'. format(filename, data_dir))
    else: 
        print('{} has not been saved'.format(filename))
    #save metadata
    filename = '/{}_{}_{}_metadata.npy'.format(condition, animal_id, date)
    np.save(data_dir + filename, metadata, allow_pickle = True)
    #check that things have been saved properly
    if os.path.isfile(data_dir + filename) == True: 
        print('{} has been saved onto {}'. format(filename, data_dir))
    else: 
        print('{} has not been saved'.format(filename))
    return None

def load_data(meta_dir, data_dir): 
    metadata = np.load(meta_dir, allow_pickle = True).item()
    data_dict = np.load(data_dir, allow_pickle = True).item()
    return metadata, data_dict

def dark_channel_tmac(data_dict, metadata):
    dark_regions = data_dict['dark_channels']
    dark_rois = data_dict['dark_rois']
    rois = data_dict['rois']
    corrected_signal = data_dict['corrected_signal']
    motion_artifact = []
    dark_corrected_signal = []
    for i in range(len(dark_rois)):
        for j in range(len(rois)):
            tmac_variables = tm.tmac_ac(dark_regions[i], corrected_signal[j])
            motion_artifact = np.append(motion_artifact, tmac_variables['m'].flatten())
            dark_corrected_signal = np.append(dark_corrected_signal, tmac_variables['a'].flatten())
            print('this is the dark corrected signals:', np.shape(dark_corrected_signal))
        dark_corrected_signal = np.append(dark_corrected_signal, motion_artifact)
        dark_corrected_signal = np.vstack(dark_corrected_signal)
        data_dict['{}_corrected'.format(dark_rois[i])] = dark_corrected_signal
        key = '{}_corrected'.format(dark_rois[i])
        trace = 'GFP'
        optional_save = input('would you like to save {} corrected graph and processed data?'.format(dark_rois[i]))
        if optional_save == 'y': 
            print('saving data and figure...')
            save_data(data_dict, metadata)
        else: 
            print('skipping data and figure save...')
    return data_dict, metadata

def linear_plus_tmac(data_dict):
    raw_data = data_dict['raw_data']
    rois = data_dict['rois']
    means = np.array(np.mean(raw_data, axis = 1))
    x = np.arange(0, len(raw_data[0]), 1)
    all_detrended = []
    green_detrended = []
    red_detrended = []
    for i in range(len(rois)): 
        slope, intercept = np.polyfit(x, raw_data[i+12], 1)
        trendline_green = slope * x + intercept
        detrended_green = raw_data[i] - trendline_green + means[i]
        detrended_red = raw_data[i + 12] - trendline_green + means[i + 12]
        green_detrended.append(detrended_green)
        red_detrended.append(detrended_red)
    green_detrended.extend(red_detrended)
    all_detrended = np.vstack(green_detrended)
    #run tmac
    data_dict['detrended_to_red'] = all_detrended
    motion_artifact = []
    corrected_signal = []
    for i in range(len(rois)):
        tmac_variables = tm.tmac_ac(all_detrended[i + 12], all_detrended[i])
        motion_artifact.append(tmac_variables['m'].flatten() - np.min(tmac_variables['m'].flatten()))
        corrected_signal.append(tmac_variables['a'].flatten() - np.min(tmac_variables['a'].flatten()))
    corrected_signal.extend(motion_artifact)
    corrected_signal = np.vstack(corrected_signal)
    data_dict['tmac_detrended_signal'] = corrected_signal
    return data_dict

def linear_separate_tmac(data_dict, metadata):
    data_dict = linear_detrend(data_dict)
    data_dict = run_tmac(data_dict, metadata)
    return data_dict


def only_tmac (data_dict):
    data = data_dict['raw_data'] 
    rois = data_dict['rois']
    motion_artifact = []
    corrected_signal = []
    for i in range(len(rois)):
        tmac_variables = tm.tmac_ac(data[i + 12], data[i])
        motion_artifact.append(tmac_variables['m'].flatten() - np.min(tmac_variables['m'].flatten()))
        corrected_signal.append(tmac_variables['a'].flatten() - np.min(tmac_variables['a'].flatten()))
    corrected_signal.extend(motion_artifact)
    corrected_signal = np.vstack(corrected_signal)
    data_dict['tmac_signal'] = corrected_signal
    return data_dict

def seeER_parallel_preprocess(file_dict): 
    metadata_file = file_dict['metadata_file']
    data_file = file_dict['data_file']
    time_file = file_dict['time_file']
    metadata, data_dict = extract_metadata(metadata_file)
    data, led_blue, led_green, headers =  extract_indices(data_file)
    data_dict, metadata = extract_time(time_file, data_file, metadata, data_dict, led_blue)
    data_dict = split_channels(metadata, data_dict, data, headers, led_blue, led_green)
    #parallel processing: do only tmac and do linear detrend then tmac 
    data_dict = linear_plus_tmac(data_dict)
    data_dict = only_tmac(data_dict)
    data_dict = linear_separate_tmac(data_dict, metadata)
    if 'power_file' in file_dict: 
        power_file = file_dict['power_file']
        metadata = extract_ledpower(power_file, metadata)
        led_correlation(metadata, data_dict)
    else:
        pass
    save_data(data_dict, metadata)
    return data_dict, metadata

#find raw files
def fetch_photometry_data(photo_raw_dir):
    #copy the path of the photometry raw data to a photometry data .txt file
    #SeeER/raw_data
    photometry_paths = []
    for root, dirs, files in os.walk(photo_raw_dir):
        for file in files:
            if file.endswith('.csv') and not file.endswith('_ts.csv') and not file.endswith('_pd.csv'):
            # if file.endswith('.txt'):
                photometry_paths.append(os.path.join(root, file))
    txt_file = os.path.join(photo_raw_dir, 'photo_raw_data.txt')
    with open(txt_file, 'w') as f:
        for p in photometry_paths:
            f.write(p + '\n')
    print('raw data paths have been saved')
    return None

#find and make a list of preprocessed directories
def find_pp_data(pp_data_dir):
    #SeeER/outputs/preprocessed
    pp_dirs = []
    for root, dirs, files in os.walk(pp_data_dir):
        for file in files: 
            if file.endswith('datadict.npy') and not file.endswith('_metadata.npy'):
                pp_dirs.append(file)
    txt_file = os.path.join(pp_data_dir, 'pp_data.txt')
    with open(txt_file, 'w') as f:
        for p in pp_dirs:
            f.write(p + '\n')
    print('preprocessed data have been saved')
    return None

#compare which files exists and make a third .txt file which lists the files that needs preprocessing
def compare_files(photo_raw_dir, pp_data_dir):
    #SeeER/raw_data/photo_raw_data.txt vs SeeER/outputs/preprocessed/pp_data.txt
    raw_files = open(os.path.join(photo_raw_dir, 'photo_raw_data.txt'), 'r')
    pp_files = open(os.path.join(pp_data_dir, 'pp_data.txt'), 'r')
    raw = raw_files.readlines()
    pp = pp_files.readlines()
    raw_files.close()
    pp_files.close()
    to_pp_files = []
    for r in raw:
        raw_file = os.path.split(r)[-1]
        raw_file = raw_file.split('.csv')[0] + '_datadict.npy\n'
        if raw_file not in pp:
            print('raw file:', raw_file)
            to_pp_files.append(r)
    txt_file = os.path.join(photo_raw_dir, 'files_to_preprocess.txt')
    with open(txt_file, 'w') as f:
        for t in to_pp_files:
            f.write(os.path.split(t)[0] + '\n')
    print('unpreprocessed list saved')
    return None
    
#only run preprocessing on those files or if you updated your preprocessing script run on all 
#feed the .txt file with the files to preprocess 
def run_preprocess(photo_raw_dir, run_all = True):
    failed_pp = []
    txt_file = os.path.join(photo_raw_dir, 'failed_preprocess.txt')
    if run_all is True:
        print('running all files in the raw data directory')
        raw_files = open(os.path.join(photo_raw_dir, 'photo_raw_data.txt'), 'r')
        raw = raw_files.readlines()
        raw_files.close()
        for r in raw:
            try: 
                print('now preprocessing {}'.format(r))
                file_dict = locate_data(os.path.split(r)[0])
                os.chdir(os.path.split(r)[0])
                seeER_parallel_preprocess(file_dict)
            except: 
                current_dict = r
                print('preprocessing failed for {}'.format(current_dict))
                failed_pp.append(current_dict)
        with open(txt_file, 'w') as f:
            for t in failed_pp:
                f.write(t)
    else: 
        to_pp_files = open(os.path.join(photo_raw_dir, 'files_to_preprocess.txt'), 'r')
        to_pp = to_pp_files.readlines()
        to_pp_files.close()
        for t in to_pp:
            try: 
                print('now preprocessing {}'.format(t.split('\n')[0]))
                file_dict = locate_data(t.split('\n')[0])
                os.chdir(t.split('\n')[0])
                seeER_parallel_preprocess(file_dict)
            except Exception as e:
                print('preprocessing failed for {}'.format(t))
                print('error:', e)
                failed_pp.append(t)
        with open(txt_file, 'w') as f:
            for t in failed_pp:
                print(t)
                f.write(t)
    return None

def run_preprocess_all(raw_data_dir, preprocessed_dir, run_all = True):
    fetch_photometry_data(raw_data_dir)
    find_pp_data(preprocessed_dir)
    compare_files(raw_data_dir, preprocessed_dir)
    current_dict = run_preprocess(raw_data_dir, run_all)
    return None

def run_preprocess_single(file_dir, file_dict): 
    file_dict = locate_data(file_dir)
    os.chdir(file_dir)
    datadict, metadata = seeER_parallel_preprocess(file_dict)
    return None

