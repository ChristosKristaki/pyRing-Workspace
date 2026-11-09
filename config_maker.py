from matplotlib import text
import numpy as np
import configparser 
from gwosc import datasets
from gwosc.datasets import query_events
from gwosc.datasets import event_detectors
import glob
from pathlib import Path
import os
from scipy import constants
from astropy import constants as const
import astropy.units as u
import lal
from collections import defaultdict


config = configparser.ConfigParser()
config.optionxform = str

def AQ_220_config(gw_event, trigtime, detectors, mf_median, fix_t, fix_ra, fix_dec, channel_H1, channel_L1, channel_V1, start_time_flag, gw_catalogue_flag, t_min):
    
    input_section = {

        'run-type' : 'full',
        'pesummary' : '0',
        'screen-output' : '0',
        'trigtime' : trigtime,
        'template' : 'Kerr',
        'reference-amplitude' : '1E-21',
        'amp-non-prec-sym' : '1',
        'channel-H1' : channel_H1,
        'channel-L1' : channel_L1,
        'channel-V1' : channel_V1,
        'sampling-rate' : '4096.0',
        'detectors' : detectors,
        'sky-frame' : 'equatorial',
        'download-data' : '1',
        'datalen-download' : '64.0',
        'noise-chunksize' : '4.0',
        'signal-chunksize' : '4.0',
        'f-min-bp' : '20.0',
        'f-max-bp' : '2043.0',
        'window-onsource' : '0',
        'analysis-duration' : '0.2',
        'qnm-fit' : '0',
        'ignore-data-filename' : '0',
        'ref-det' : 'H1',
        'truncate' : '1',
        'gw-data-find' : '0',
        'spheroidal' : '0',
        'run-tag' : '', 
        'kerr-modes' : '[(2,2,2,0)]',
        'area-quantization' : '1',
        'output' : f'/groups/astro/ctj304/catalogue_events_prod/catalogue_events_prod_{start_time_flag}/GW{gw_event}_220_AQ_{start_time_flag}',
        #'data-H1' : data_H1,
        #'data-L1' : data_L1,
        #'data-V1' : data_V1,
    }

    config['input'] = input_section

    config['Sampler settings'] = {
    
        'nlive' : '4096',
        'poolsize' : '100',
        'maxmcmc' : '4096',
        'nthreads' : '8',
        'verbose' : '2',
        'seed' : '1234',
    }

    config['Priors'] = {

        'mf-time-prior' : mf_median,
        'gr-time-prior' : '1',
        't-min' : t_min,
        't-max' : '60',
        'fix-t' : fix_t,
        'fix-ra' : fix_ra,
        'fix-dec' : fix_dec,
        'fix-phi' : '0',
        'alpha-min' : '0',
        'alpha-max' : '50',
    }

    config['Injection'] = {}

    config['Plot'] = {

        'imr-samples': '',
    }

    filename = f'/groups/astro/ctj304/pyring/pyRing/config_files/catalogue_{start_time_flag}/{gw_catalogue_flag}/config_GW_{gw_event}_220_AQ_prod_download_{start_time_flag}.ini'
    with open(filename, 'w') as file:
        config.write(file)
    print(f"Written {filename}")


def GR_220_config(gw_event, trigtime, detectors, mf_median, fix_t, fix_ra, fix_dec, channel_H1, channel_L1, channel_V1, start_time_flag, gw_catalogue_flag, t_min):

    input_section = {

        'run-type' : 'full',
        'pesummary' : '0',
        'screen-output' : '0',
        'trigtime' : trigtime,
        'template' : 'Kerr',
        'reference-amplitude' : '1E-21',
        'amp-non-prec-sym' : '1',
        'channel-H1' : channel_H1,
        'channel-L1' : channel_L1,
        'channel-V1' : channel_V1,
        'sampling-rate' : '4096.0',
        'detectors' : detectors,
        'sky-frame' : 'equatorial',
        'download-data' : '1',
        'datalen-download' : '64.0',
        'noise-chunksize' : '4.0',
        'signal-chunksize' : '4.0',
        'f-min-bp' : '20.0',
        'f-max-bp' : '2043.0',
        'window-onsource' : '0',
        'analysis-duration' : '0.2',
        'qnm-fit' : '0',
        'ignore-data-filename' : '0',
        'ref-det' : 'H1',
        'truncate' : '1',
        'gw-data-find' : '0',
        'spheroidal' : '0',
        'run-tag' : '', 
        'kerr-modes' : '[(2,2,2,0)]',
        'area-quantization' : '0',
        'output' : f'/groups/astro/ctj304/catalogue_events_prod/catalogue_events_prod_{start_time_flag}/GW{gw_event}_220_GR_{start_time_flag}',
        #'data-H1' : data_H1,
        #'data-L1' : data_L1,
        #'data-V1' : data_V1,
    }

    config['input'] = input_section

    config['Sampler settings'] = {

        'nlive' : '4096',
        'poolsize' : '100',
        'maxmcmc' : '4096',
        'nthreads' : '8',
        'verbose' : '2',
        'seed' : '1234',
    }

    config['Priors'] = {
    
        'mf-time-prior' : mf_median,
        'gr-time-prior' : '1',
        't-min' : t_min,
        't-max' : '60',
        'fix-t' : fix_t,
        'fix-ra' : fix_ra,
        'fix-dec' : fix_dec,
        'fix-phi' : '0',
        }

    config['Injection'] = {}

    config['Plot'] = {
    
        'imr-samples': '',
    }

    filename = f'/groups/astro/ctj304/pyring/pyRing/config_files/catalogue_{start_time_flag}/{gw_catalogue_flag}/config_GW_{gw_event}_220_GR_prod_download_{start_time_flag}.ini'
    with open(filename, 'w') as file:
        config.write(file)
    print(f"Written {filename}")


def AQ_221_config(gw_event, trigtime, detectors, mf_median, fix_t, fix_ra, fix_dec, channel_H1, channel_L1, channel_V1, start_time_flag, gw_catalogue_flag, t_min):

    input_section = {

        'run-type' : 'full',
        'pesummary' : '0',
        'screen-output' : '0',
        'trigtime' : trigtime,
        'template' : 'Kerr',
        'reference-amplitude' : '1E-21',
        'amp-non-prec-sym' : '1',
        'channel-H1' : channel_H1,
        'channel-L1' : channel_L1,
        'channel-V1' : channel_V1,
        'sampling-rate' : '4096.0',
        'detectors' : detectors,
        'sky-frame' : 'equatorial',
        'download-data' : '1',
        'datalen-download' : '64.0',
        'noise-chunksize' : '4.0',
        'signal-chunksize' : '4.0',
        'f-min-bp' : '20.0',
        'f-max-bp' : '2043.0',
        'window-onsource' : '0',
        'analysis-duration' : '0.2',
        'qnm-fit' : '0',
        'ignore-data-filename' : '0',
        'ref-det' : 'H1',
        'truncate' : '1',
        'gw-data-find' : '0',
        'spheroidal' : '0',
        'run-tag' : '', 
        'kerr-modes' : '[(2,2,2,0), (2,2,2,1)]',
        'area-quantization' : '1',
        'output' : f'/groups/astro/ctj304/catalogue_events_prod/catalogue_events_prod_{start_time_flag}/GW{gw_event}_221_AQ_{start_time_flag}',
        #'data-H1' : data_H1,
        #'data-L1' : data_L1,
        #'data-V1' : data_V1,
    }

    config['input'] = input_section    

    config['Sampler settings'] = {

        'nlive' : '4096',
        'poolsize' : '100',
        'maxmcmc' : '4096',
        'nthreads' : '8',
        'verbose' : '2',
        'seed' : '1234',
    }

    config['Priors'] = {
    
        'mf-time-prior' : mf_median,
        'gr-time-prior' : '1',
        't-min' : t_min,
        't-max' : '60',
        'fix-t' : fix_t,
        'fix-ra' : fix_ra,
        'fix-dec' : fix_dec,
        'fix-phi' : '0',
        'alpha-min' : '0',
        'alpha-max' : '50',
        }

    config['Injection'] = {}

    config['Plot'] = {
    
        'imr-samples': '',
    }

    filename = f'/groups/astro/ctj304/pyring/pyRing/config_files/catalogue_{start_time_flag}/{gw_catalogue_flag}/config_GW_{gw_event}_221_AQ_prod_download_{start_time_flag}.ini'
    with open(filename, 'w') as file:
        config.write(file)
    print(f"Written {filename}")

def GR_221_config(gw_event, trigtime, detectors, mf_median, fix_t, fix_ra, fix_dec, channel_H1, channel_L1, channel_V1, start_time_flag, gw_catalogue_flag, t_min):

    input_section = {

        'run-type' : 'full',
        'pesummary' : '0',
        'screen-output' : '0',
        'trigtime' : trigtime,
        'template' : 'Kerr',
        'reference-amplitude' : '1E-21',
        'amp-non-prec-sym' : '1',
        'channel-H1' : channel_H1,
        'channel-L1' : channel_L1,
        'channel-V1' : channel_V1,
        'sampling-rate' : '4096.0',
        'detectors' : detectors,
        'sky-frame' : 'equatorial',
        'download-data' : '1',
        'datalen-download' : '64.0',
        'noise-chunksize' : '4.0',
        'signal-chunksize' : '4.0',
        'f-min-bp' : '20.0',
        'f-max-bp' : '2043.0',
        'window-onsource' : '0',
        'analysis-duration' : '0.2',
        'qnm-fit' : '0',
        'ignore-data-filename' : '0',
        'ref-det' : 'H1',
        'truncate' : '1',
        'gw-data-find' : '0',
        'spheroidal' : '0',
        'run-tag' : '', 
        'kerr-modes' : '[(2,2,2,0), (2,2,2,1)]',
        'area-quantization' : '0',
        'output' : f'/groups/astro/ctj304/catalogue_events_prod/catalogue_events_prod_{start_time_flag}/GW{gw_event}_221_GR_{start_time_flag}',
        #'data-H1' : data_H1,
        #'data-L1' : data_L1,
        #'data-V1' : data_V1,
    }

    #if data_H1 : 
    #        input_section['data-H1'] = data_H1
    #if data_L1 : 
    #        input_section['data-L1'] = data_L1
    #if data_V1 : 
    #        input_section['data-V1'] = data_V1
    
    
    config['input'] = input_section

    config['Sampler settings'] = {

        'nlive' : '4096',
        'poolsize' : '100',
        'maxmcmc' : '4096',
        'nthreads' : '8',
        'verbose' : '2',
        'seed' : '1234',
    }

    config['Priors'] = {

        'mf-time-prior' : mf_median,
        'gr-time-prior' : '1',
        't-min' : t_min,
        't-max' : '60',
        'fix-t' : fix_t,
        'fix-ra' : fix_ra,
        'fix-dec' : fix_dec,
        'fix-phi' : '0',
    }

    config['Injection'] = {}

    config['Plot'] = {
    
        'imr-samples': '',
    }

    filename = f'/groups/astro/ctj304/pyring/pyRing/config_files/catalogue_{start_time_flag}/{gw_catalogue_flag}/config_GW_{gw_event}_221_GR_prod_download_{start_time_flag}.ini'
    with open(filename, 'w') as file:
        config.write(file)
    print(f"Written {filename}")    

gw_catalogue_GWTC4 = [
    '230601_224134',
    '230609_064958',
    '230628_231200',
    '230811_032116',
    '230814_061920',
    '230824_033047',
    '230914_111401',
    '230922_020344',
    '230922_040658',
    '230924_124453',
    '230927_043729',
    '230927_153832',
    '230928_215827',
    '231001_140220',
    '231028_153006',
    '231102_071736',
    '231108_125142',
    '231213_111417',
    '231223_032836',
    '231226_101520',
]

gw_catalogue_GWTC3 = [
    '150914',
    '170104',
    '170814',
    '170823',
    '190408_181802',
    '190512_180714',
    '190513_205428',
    '190519_153544',
    '190521_030229',
    '190521_074359',
    '190602_175927',
    '190706_222641',
    '190708_232457',
    '190727_060333',
    '190828_063405',
    '190910_112807',
    '190915_235702',
    '191109_010717',
    '191222_033537',
    '200129_065458',
    '200224_222234',
    '200311_115853',
]

#CHOOSE CATALOGUE AND START TIME. ALWAYS CHECK BEFORE RUNNING
gw_catalogue = gw_catalogue_GWTC4
gw_catalogue_flag = 'GWTC4'
sp_character = '_'
start_time = 10 #in units of M
start_time_flag = str(start_time)+'M'

same_day = defaultdict(list)
for ev in gw_catalogue:
    same_day[ev.split('_')[0]].append(ev) #DICT WITH EACH EVENT CORRESPONDING TO A 6DIGIT VALUE, SAME DAY EVENTS IN ONE VALUE
print(same_day)
suffix_map = {} #will contain only events which happened on the same day
for day, evs in same_day.items():
    if len(evs) > 1:
        for letter, ev in zip('AB', sorted(evs)):   # sorted by time string -> earliest = A
            suffix_map[ev] = letter
print(suffix_map)
for gw_event in gw_catalogue:

    config.clear()
    detectors = []

    if gw_catalogue is gw_catalogue_GWTC3:
        if '_' in gw_event and gw_event in suffix_map:
            name = gw_event.split('_')[0]
            temp = name + suffix_map[gw_event]
        elif '_' in gw_event and gw_event not in suffix_map:
            temp = gw_event.split('_')[0]
        else:
            temp = gw_event

        base = '/groups/astro/ctj304/GWTC3/GWTC3_configs'
        file = (config.read(f'{base}/config_GW{temp}A_GWTC3-RR1_TEOB_22_0M_SEED_1234.ini')
                or config.read(f'{base}/config_GW{temp}A_GWTC3-RR2_TEOB_22_0M_SEED_1234.ini')
                or config.read(f'{base}/config_GW{temp}B_GWTC3-RR1_TEOB_22_0M_SEED_1234.ini')
                or config.read(f'{base}/config_GW{temp}B_GWTC3-RR2_TEOB_22_0M_SEED_1234.ini')
                or config.read(f'{base}/config_GW{temp}_GWTC3-RR2_TEOB_22_0M_SEED_1234.ini')
                or config.read(f'{base}/config_GW{temp}_GWTC3-RR1_TEOB_22_0M_SEED_1234.ini'))
        if not file:
            raise FileNotFoundError(f'No GWTC3 config for {gw_event} (tried {temp})')
    elif gw_catalogue is gw_catalogue_GWTC4:
            file_read = config.read(f'/groups/astro/ctj304/GWTC4/pyring/samples/GW{gw_event}/GW{gw_event}_Kerr_220_np_10M/config_GW{gw_event}_Kerr_220_np_10M.ini')
            print(file_read)
    
    trigtime = config.getfloat("input", "trigtime")
    mf_median = config.getfloat("Priors", "mf-time-prior")
    
    if gw_catalogue is gw_catalogue_GWTC4:
        if start_time in (0, 5, 15, 20):
            temp_t = start_time * lal.MTSUN_SI * mf_median
            fix_t = f"{temp_t:.6f}"
            t_min = start_time
        elif start_time == 10:
            fix_t = config.getfloat("Priors", "fix-t")
            t_min = 10
    elif gw_catalogue is gw_catalogue_GWTC3:
        if start_time in (0, 5, 15, 20):
            temp_t = start_time * lal.MTSUN_SI * mf_median
            fix_t = f"{temp_t:.6f}"
            t_min = start_time
        elif start_time == 10:
            fix_t = config.getfloat("Priors", "fix-t")
            t_min = 10

    fix_ra = config.getfloat("Priors", "fix-ra")
    fix_dec = config.getfloat("Priors", "fix-dec")
    if config.has_option("input", "channel-H1"):
        channel_H1 = config.get("input", "channel-H1")
        detectors.append('H1')
    else:
        channel_H1 = ""
        
    
    if config.has_option("input", "channel-L1"):
        channel_L1 = config.get("input", "channel-L1")
        detectors.append('L1')
    else:
        channel_L1 = ""
    
    if config.has_option("input", "channel-V1"):
        channel_V1 = config.get("input", "channel-V1")
        detectors.append('V1')
    else:
        channel_V1 = ""

    ref_dets = [d.strip() for d in config.get("input", "detectors", fallback="").split(",") if d.strip()]

    detectors = []
    if "H1" in ref_dets: detectors.append("H1")
    if "L1" in ref_dets: detectors.append("L1")
    if "V1" in ref_dets: detectors.append("V1")

    AQ_220_config(gw_event = gw_event, trigtime = trigtime, detectors = ','.join(detectors), mf_median = mf_median, fix_t=fix_t, fix_dec=fix_dec, fix_ra=fix_ra, channel_H1 = channel_H1, channel_L1 = channel_L1, channel_V1 = channel_V1, start_time_flag=start_time_flag, gw_catalogue_flag=gw_catalogue_flag, t_min=t_min)
    GR_220_config(gw_event = gw_event, trigtime = trigtime, detectors = ','.join(detectors), mf_median = mf_median, fix_t=fix_t, fix_dec=fix_dec, fix_ra=fix_ra, channel_H1 = channel_H1, channel_L1 = channel_L1, channel_V1 = channel_V1, start_time_flag=start_time_flag, gw_catalogue_flag=gw_catalogue_flag, t_min=t_min)
    AQ_221_config(gw_event = gw_event, trigtime = trigtime, detectors = ','.join(detectors), mf_median = mf_median, fix_t=fix_t, fix_dec=fix_dec, fix_ra=fix_ra, channel_H1 = channel_H1, channel_L1 = channel_L1, channel_V1 = channel_V1, start_time_flag=start_time_flag, gw_catalogue_flag=gw_catalogue_flag, t_min=t_min)
    GR_221_config(gw_event = gw_event, trigtime = trigtime, detectors = ','.join(detectors), mf_median = mf_median, fix_t=fix_t, fix_dec=fix_dec, fix_ra=fix_ra, channel_H1 = channel_H1, channel_L1 = channel_L1, channel_V1 = channel_V1, start_time_flag=start_time_flag, gw_catalogue_flag=gw_catalogue_flag, t_min=t_min)

