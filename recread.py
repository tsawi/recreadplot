''' Record Reading Plots
Program to extract and plot global waveforms for a specified seismic event and 
maps regional focal mechanisms. This program creates an interactive plot to 
view global waveforms for a user-specified seismic event. The program also has 
options to plot regional maps of focal mechanisms from the CMT catalog 
globalcmt.org).
This is a python adaptation to the original record reading waveform plot 
program written by Ge Jin (jinwar@gmail.com) in Matlab.
Contributors: Janine Birnbaum, Theresa Sawi, Christopher Carchedi, and Michelle 
Lee
Instructions:
1. Event Selection: Select event to plot global waveforms for by inputting 
location and magnitude parameters of desired event
2. (User steps to download the data?)
3. Plotting the waveforms: The global waveforms will be be displayed on an 
interactive plot allowing users to zoom in/out, change component, change 
frequency band, ... (and other features)
4. Plotting focal mechanisms: (info on what user needs to do for this)
Use the ``bokeh serve`` command to run the example by executing:
    bokeh serve --show recread.py
at your command prompt. Then navigate to the URL
    http://localhost:5006/recread
in your browser.
'''

###############################################################################
# Imports
###############################################################################
import sys
import os
import warnings
from io import BytesIO
import base64

import numpy as np
import pandas as pd
from pandas.core.common import SettingWithCopyWarning
from matplotlib import pyplot as plt

warnings.simplefilter(action="ignore", category=SettingWithCopyWarning) # suppresses a warning in binning and plotting
warnings.filterwarnings('ignore')

from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource, CustomJS, Div, OpenURL, Slider, ColorBar
from bokeh.models import DatetimeTickFormatter
from bokeh.models import TextInput, HoverTool, TapTool, Select, RangeSlider, MultiChoice, RadioButtonGroup, CheckboxGroup, FreehandDrawTool, PolyDrawTool, BooleanFilter, CDSView, Range1d
from bokeh.models.widgets import Button, Panel, Tabs
from bokeh.plotting import figure, curdoc
from bokeh.tile_providers import get_provider
from bokeh import events
from bokeh.palettes import Viridis256, grey
from bokeh.transform import linear_cmap
from bokeh.io.export import get_screenshot_as_png

from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from obspy.core import AttribDict
from obspy.io.sac import SACTrace
from obspy.imaging.mopad_wrapper import beach
from obspy.clients.fdsn.mass_downloader import CircularDomain, Restrictions, MassDownloader
from obspy.geodetics.base import gps2dist_azimuth
from obspy.taup import TauPyModel

from pyproj import Transformer

###############################################################################
# Tab 1: Select event
###############################################################################

# load initial inputs from file
try:
    meta = pd.read_hdf('eventdat',key='meta')

    # Selecting an event
    lat = meta.lat[0];
    lon = meta.lon[0];
    ID = meta.ID[0];
    depth = meta.depth[0];
    t11 = UTCDateTime(meta.time[0])
    mag = meta.mag[0]
except: 
    lat = 0
    lon = 0
    ID = None
    depth = 0
    t11 = UTCDateTime(0)
    mag = 5
    
search_time_range = 72 # hours
t22 = t11 + 3600*search_time_range

# Setup search parameter inputs
input_lat = TextInput(value=str(lat), title='Latitude (-90 to 90):')
input_lon = TextInput(value=str(lon), title='Longitude (-180 to 180):')
input_time = TextInput(value=str(t11),title='Start time UTC (YYYY-MM-DD HH:MM:SS)')
search_time_range = Slider(start=0,end=72,value=72, title='Search time range (hrs)')
input_mag = RangeSlider(start=0,end=10,value=(mag-2,mag+2),title='Magnitude') # Can probably add an if condition to change the slider bar boundaries for other magnitude scales
mag_type = Select(title='Magnitude type', value='Mw',
               options=['Mw'])
search_rad = Slider(start=0,end=20,value=10,title='Search radius (deg)')
input_webservice = Select(title='Catalog', value='IRIS',
               options=['IRIS'])
input_velmodel = Select(title='Velocity Model', value='iasp91',
                options=['iasp91','ak135','ak135f','prem'])
input_ID = TextInput(value=str(ID))
input_depth = TextInput(value=str(depth))

# define coordinate transformations from lat/lon to Web Mercator
latlon2webmercator = Transformer.from_crs('EPSG:4326', 'EPSG:3857')
webmercator2latlon = Transformer.from_crs('EPSG:3857', 'EPSG:4326')

[x0, y0] = latlon2webmercator.transform(lat, lon) # transform search coordinate for plotting  
xcoord = TextInput(value=str(x0))
ycoord = TextInput(value=str(y0))

# Fetch initital data
evlon = np.array([])
evlat = np.array([])
evdepth = np.array([])
evmag = np.array([])
evtime = np.array([])
ID = np.array([])
    
try:
    client = Client(input_webservice.value)
    eventlist = (client.get_events(starttime=t11-1,endtime=t22,latitude=lat,
                               longitude=lon,minradius=0,maxradius=10,
                               minmagnitude=mag-0.5,
                               maxmagnitude=mag+0.5,
                               magnitudetype=mag_type.value,
                               catalog='GCMT',orderby='time-asc'))    

    for iev, event in enumerate(eventlist):
        evlon = np.append(evlon,eventlist[iev].origins[0].longitude)
        evlat = np.append(evlat,eventlist[iev].origins[0].latitude)
        evdepth = np.append(evdepth,eventlist[iev].origins[0].depth)
        evmag = np.append(evmag,eventlist[iev].magnitudes[0].mag)
        evtime = np.append(evtime,eventlist[iev].origins[0].time.strftime('%y/%m/%d %H:%M'))
        ID = np.append(ID,str(eventlist[iev].resource_id))

    # plot initial data
    [x1, y1] = latlon2webmercator.transform(evlat,evlon) # transform event coordinates  
    source = ColumnDataSource(data={'x': x1,'y': y1,'lat':evlat,'lon':evlon,'depth':evdepth/1000,
                                    'mag':evmag,'size':4*evmag,'time':evtime, 'id':ID})
except Exception as e:
    print(e)
    source = ColumnDataSource(data={'x':[],'y': [],'lat':[],'lon':[],'depth':[],
				   'mag':[],'size':[],'time':[], 'id':[]})

# callback for updating search parameters and plot
def update_search_params(attrname, old, new):

    eventlist = []
    
    [x_cntr, y_cntr] = latlon2webmercator.transform(float(input_lat.value),float(input_lon.value))
    
    p.x_range.start = x_cntr-padding
    p.x_range.end   = x_cntr+padding
    p.y_range.start = y_cntr-padding
    p.y_range.end   = y_cntr+padding
    # Get the current slider values
    t11 = UTCDateTime(str(input_time.value).replace('-','').replace(' ','').replace(':',''))
    t22 = t11 + 3600*search_time_range.value

    # Fetch data
    client = Client(input_webservice.value)
    
    try:
        eventlist = (client.get_events(starttime=t11,endtime=t22,latitude=float(input_lat.value),longitude=float(input_lon.value),minradius=0,maxradius=float(search_rad.value),
                minmagnitude=float(input_mag.value[0])-0.5,maxmagnitude=float(input_mag.value[1])+0.5,magnitudetype=mag_type.value,catalog='GCMT',
                orderby='time-asc'))
    
        evlon = np.array([])
        evlat = np.array([])
        evdepth = np.array([])
        evmag = np.array([])
        evtime = np.array([])
        ID = np.array([])

        for iev, event in enumerate(eventlist) :
            evlon = np.append(evlon,eventlist[iev].origins[0].longitude)
            evlat = np.append(evlat,eventlist[iev].origins[0].latitude)
            evdepth = np.append(evdepth,eventlist[iev].origins[0].depth)
            evmag = np.append(evmag,eventlist[iev].magnitudes[0].mag)
            evtime = np.append(evtime,eventlist[iev].origins[0].time.strftime('%y/%m/%d %H:%M'))
            ID = np.append(ID,str(eventlist[iev].resource_id))

        [x1, y1] = latlon2webmercator.transform(evlat,evlon) # transform event coordinates  
        source.data = {'x': x1,'y': y1,'lat':evlat,'lon':evlon,'depth':evdepth/1000,
                       'mag':evmag,'size':4*evmag,'time':evtime, 'id':ID}
    
    except Exception as e:
        print(e)
        source.data = {'x':[],'y': [],'lat':[],'lon':[],'depth':[],
				'mag':[],'size':[],'time':[], 'id':[]}

for w in [input_lat, input_lon, input_time, input_mag, mag_type, search_rad, input_webservice, input_ID]:
    w.on_change('value', update_search_params)

[x0, y0] = latlon2webmercator.transform(lat, lon) # transform search coordinate for plotting  
    
# draw a map
# figure bounds supplied in web mercator coordinates
padding = 0.5*10**6 # about 5 degrees
p = figure(x_range=(x0-padding, x0+padding), y_range=(y0-padding, y0+padding),
                x_axis_type='mercator', y_axis_type='mercator',
		tools='tap,pan,wheel_zoom,box_zoom,reset',
           	active_drag='pan',active_scroll='wheel_zoom')
p.add_tile(get_provider('ESRI_IMAGERY'))
   

# hover tool
# Add a hover tool, that selects the circle
TOOLTIPS = [
    ('Lat', '@lat{0.00}'),
    ('Lon', '@lon{0.00}'),
    ('Depth (km)','@depth{0.0}'),
    ('Mw','@mag{0.0}'),
    ('Time (UTC)','@time')
]

cr = p.circle('x', 'y', size='size',source=source,color='red',line_color='black',hover_color='yellow',hover_line_color='black')
callbackhover = CustomJS(args={'circle': cr.data_source})
p.add_tools(HoverTool(tooltips=TOOLTIPS, callback=callbackhover, renderers=[cr], mode='mouse'))

# Tap tool
def callbacktap(xcoord,ycoord, attributes=[]):
    """
    Function to build a suitable CustomJS to display the current event
    in the div model.
    """
    
    return CustomJS(args={'xcoord':xcoord,'ycoord':ycoord}, code="""
        const attrs = %s;
        const args = [];
        const val = JSON.stringify(cb_obj[attrs[0]], function(key, val) {
            return val.toFixed ? Number(val.toFixed(2)) : val;
        })
        args.push(val)
        const text = "<span>" + cb_obj.event_name + args.join(", ") + "</span>";
        const values = val.split('}')[0].split(',')
        const x_select = values[3].split(':')[1]
        const y_select = values[4].split(':')[1]
        xcoord.value = x_select;
        ycoord.value = y_select;
        
    """ % (attributes))
        
def callbacktap2(attrname, old, new):
    [lat,lon] = webmercator2latlon.transform(float(xcoord.value),float(ycoord.value))
    input_lat.value = '%.2f' % lat
    input_lon.value = '%.2f' % lon
    t11 = UTCDateTime(str(input_time.value).replace('-','').replace(' ','').replace(':',''))
    t22 = t11 + 3600*search_time_range.value

    # Fetch data
    client = Client(input_webservice.value)
    
    try:
        eventlist = (client.get_events(starttime=t11,endtime=t22,latitude=float(input_lat.value),longitude=float(input_lon.value),minradius=0,maxradius=0.5,
                minmagnitude=float(input_mag.value[0]),maxmagnitude=float(input_mag.value[1]),magnitudetype=mag_type.value,catalog='GCMT',
                orderby='time-asc'))
        input_time.value = str(eventlist[0].origins[0].time).split('.')[0].replace('T', ' ')
        input_mag.value  = (eventlist[0].magnitudes[0].mag - 0.1,eventlist[0].magnitudes[0].mag + 0.1)
        input_ID.value   = str(eventlist[0].resource_id)
        evtime = eventlist[0].origins[0].time.strftime('%y/%m/%d %H:%M')
        input_depth.value = str(eventlist[0].origins[0].depth/1000)
        
        source.data = {'x': [float(xcoord.value)],'y': [float(ycoord.value)],'lat':[lat],'lon':[lon],'depth':[eventlist[0].origins[0].depth/1000],
                       'mag':[eventlist[0].magnitudes[0].mag],'size':[4*eventlist[0].magnitudes[0].mag],'time':[evtime], 'id':[input_ID.value]}
    except Exception as e: 
        print(e)     

p.js_on_event(events.SelectionGeometry, callbacktap(xcoord,ycoord, attributes=['geometry']))
ycoord.on_change('value', callbacktap2)
    
# Save event info
def save_button_callback():
    eventdat = pd.DataFrame(data={'ID':str(input_ID.value),'lat':float(input_lat.value), 
                                  'lon':float(input_lon.value), 'time':str(input_time.value), 
                                  'mag':float(input_mag.value[0]), 'depth':str(input_depth.value)}, index=[0])
    eventdat.to_hdf('eventdat',key='meta')

save_button = Button(label='Save', button_type='success')
save_button.on_click(save_button_callback)

# Button to stop the server
def quit_button_callback():
    sys.exit()  # Stop the server
    
quit_button = Button(label='Quit', button_type='success')
quit_button.on_click(quit_button_callback)

# figure caption
figure1 = Div(text='Figure 1: Hover over points to see event details. Click desired event to highlight. Click save to write event details to file.')

# Display plots
layout1 = column(Div(text='<h1> Select event <h1>'),
                 (row(column([input_lat,
                              input_lon,
                              search_rad,
                              input_time,
                              search_time_range,
                              input_mag,mag_type,
                              input_webservice,
                              input_velmodel,
                              save_button,
                              figure1,
                              quit_button]),p)))
panel1 = Panel(child=layout1, title='Select event')

###############################################################################
# Tab 2: Download data
###############################################################################
network_options = [('0', '_US-ALL'), ('1', '_GSN'), ('2', 'TA')]
network = MultiChoice(title='Station network', value=['0', '1'],
               options=network_options)
add_network_text = TextInput(value='')
add_network_button = Button(label='Add network',button_type='success')

epicentral_distance = RangeSlider(start=0,end=180,value=(0,180),title='Epicentral distance (deg)')

min_before = TextInput(value='10',title='Minutes before:')
min_after = TextInput(value='120',title='Minutes after:')

unit_options = ['Frequncy (Hz)', 'Period (s)']
radio_button_group = RadioButtonGroup(labels=unit_options, active=1)
lowfilter_left = TextInput(value='200')
lowfilter_right = TextInput(value='30')
midfilter_left = TextInput(value='25')
midfilter_right = TextInput(value='5')
highfilter_left = TextInput(value='5')
highfilter_right = TextInput(value='2')
resample_delta = TextInput(value='0.5',title='Resample rate:')
phase_shift = Select(title="Align:", value="None", options=['None'])

load_stations = Button(label='Load station map', button_type='success')
download_data = Button(label='Download data', button_type='success')

try: 
    df
except: 
    try:
        df = pd.read_hdf('eventdat',key='data')
    except: 
        df = pd.DataFrame(data={'Network':['None'],
                                'Station':['None'], 
                                'Lat':[0],
                                'Lon':[0],
                                'Azimuth':[0], 
                                'Distance':[1], 
                                'Frequency':['Raw'],
                                'Channel':['BHZ'],
                                'Time':[np.array([
                                        str(UTCDateTime(0)).split('.')[0].replace('T', ' '), 
                                        str(UTCDateTime(0)+1).split('.')[0].replace('T', ' ')])], 
                                'Data':[np.array([0,0])], 
                                'url':['None']})

stat_lat_from_file = df['Lat'].values
stat_lon_from_file = df['Lon'].values
[stat_x_from_file, stat_y_from_file] = latlon2webmercator.transform(stat_lat_from_file, stat_lon_from_file)
stations_source = ColumnDataSource(data={'stat_x': stat_x_from_file,'stat_y': stat_y_from_file,
                                             'stat_lat':stat_lat_from_file,'stat_lon':stat_lon_from_file,
                                             'stat_name':df['Station'].values})


# download some station info based on networks
# Stations map (could maybe pull station info from the networks and allow a select tool on a map)
def add_network_option():
    network_options.append((str(len(network_options)),str(add_network_text.value)))
    network.options = network_options
    network.value.append(str(len(network_options)-1))

add_network_button.on_click(add_network_option)

#def networkscallback(attrs,new,old):
#    print(','.join([network_options[i][1] for i in np.array(network.value).astype(int)]))
#network.on_change('value', networkscallback)

def select_by_dist(attrs,new,old):
    [x_cntr,y_cntr] = latlon2webmercator.transform(input_lat.value,input_lon.value)
    xpadding = latlon2webmercator.transform(0,180)[0]/180*float(epicentral_distance.value[1])
    ypadding = latlon2webmercator.transform(85,0)[1]/90*float(epicentral_distance.value[1])

    p2.x_range.start=x_cntr-xpadding
    p2.x_range.end  =x_cntr+xpadding
    p2.y_range.start=np.max([y_cntr-ypadding,latlon2webmercator.transform(-85,0)[1]])
    p2.y_range.end  =np.min([y_cntr+ypadding,latlon2webmercator.transform(85,0)[1]])

[x_cntr,y_cntr] = latlon2webmercator.transform(input_lat.value,input_lon.value)
xpadding = latlon2webmercator.transform(0,180)[0]/180*float(epicentral_distance.value[1])
ypadding = latlon2webmercator.transform(85,0)[1]/90*float(epicentral_distance.value[1])

p2 = figure(x_range=(x_cntr-xpadding, x_cntr+xpadding), 
           y_range=(np.max([y_cntr-ypadding,latlon2webmercator.transform(-85,0)[1]]), np.min([y_cntr+ypadding,latlon2webmercator.transform(85,0)[1]])),
           x_axis_type='mercator', y_axis_type='mercator',
           tools='pan,wheel_zoom,box_zoom,reset',
           active_drag='pan',active_scroll='wheel_zoom')

p2.circle('x', 'y', size='size',source=source,color='red',line_color='black')
p2.add_tile(get_provider('ESRI_IMAGERY'))

epicentral_distance.on_change('value',select_by_dist)

def toggle_units(attrs,new,old):
    lowfilter_left.value   = str(1/float(lowfilter_left.value))
    lowfilter_right.value  = str(1/float(lowfilter_right.value))
    midfilter_left.value   = str(1/float(midfilter_left.value))
    midfilter_right.value  = str(1/float(midfilter_right.value))
    highfilter_left.value  = str(1/float(highfilter_left.value))
    highfilter_right.value = str(1/float(highfilter_right.value))
    resample_delta.value   = str(1/float(resample_delta.value))
    
radio_button_group.on_change('active',toggle_units)

def load_stations_callback():
    stat_lat = np.array([])
    stat_lon = np.array([])
    stat_name = np.array([])
    t11 = UTCDateTime(str(input_time.value).replace('-','').replace(' ','').replace(':',''))
    try:
        inventory = client.get_stations(network=','.join([network_options[i][1] for i in np.array(network.value).astype(int)]),
                                    station="*",
                                    latitude=float(input_lat.value), 
                                    longitude=float(input_lon.value),
                                    minradius=float(epicentral_distance.value[0]), 
                                    maxradius=float(epicentral_distance.value[1]),
                                    starttime=t11-float(min_before.value)*60,
                                    endtime=t11+float(min_after.value)*60)
        for net in inventory.networks:
            for station in net.stations:
                stat_lat = np.append(stat_lat,station.latitude)
                stat_lon = np.append(stat_lon,station.longitude)
                stat_name = np.append(stat_name,station.code)
                
    except Exception as e: 
        print(e)   
    try:
        [stat_x,stat_y] = latlon2webmercator.transform(stat_lat,stat_lon)
        stat_x[stat_x<(x_cntr-latlon2webmercator.transform(0,180)[0])] += 2*latlon2webmercator.transform(0,180)[0]
        stat_x[stat_x>(x_cntr+latlon2webmercator.transform(0,180)[0])] -= 2*latlon2webmercator.transform(0,180)[0]
    except: 
        stat_x = np.array([])
        stat_y = np.array([])
    stations_source.data = {'stat_x': stat_x,'stat_y': stat_y,
                            'stat_lat':stat_lat,'stat_lon':stat_lon,
                            'stat_name':stat_name}
            
load_stations.on_click(load_stations_callback)
p2.triangle('stat_x', 'stat_y', size=10,source=stations_source,color='red',line_color='black')

def download_data_callback():
    t11 = UTCDateTime(str(input_time.value).replace('-','').replace(' ','').replace(':',''))
    inventory = client.get_stations(network=','.join([network_options[i][1] for i in np.array(network.value).astype(int)]), 
                                    station="*",
                                    latitude=float(input_lat.value), 
                                    longitude=float(input_lon.value),
                                    minradius=float(epicentral_distance.value[0]), 
                                    maxradius=float(epicentral_distance.value[1]),
                                    starttime=t11-float(min_before.value)*60,
                                    endtime=t11+float(min_after.value)*60,
                                    level='response')    
    bulk = []
    bulk_stat = []
    for net in inventory.networks:
        for stat in net.stations:
            chan = np.array([i.code for i in stat.channels])
            location = np.array([i.location_code for i in stat.channels])
            stored = False
            for loc_i in np.unique(location): # pull only one location from each station with preference for BH* and then HH*
                if ('BHZ' in chan[location==loc_i]) & \
                (('BHN' in chan[location==loc_i]) | ('BH1' in chan[location==loc_i])) & \
                (('BHE' in chan[location==loc_i]) | ('BH2' in chan[location==loc_i])) & (~stored):
                    bulk.append((net.code,stat.code,loc_i,'BH*',t11-float(min_before.value)*60,t11+float(min_after.value)*60))
                    bulk_stat.append(stat.code)
                    stored = True
            if ~stored:
                for loc_i in np.unique(location): # pull only one location from each station with preference for BH* and then HH*
                    if ('HHZ' in chan[location==loc_i]) & \
                    (('HHN' in chan[location==loc_i]) | ('HH1' in chan[location==loc_i])) & \
                    (('HHE' in chan[location==loc_i]) | ('HH2' in chan[location==loc_i])) & (~stored):
                        bulk.append((net.code,stat.code,loc_i,'HH*',t11-float(min_before.value)*60,t11+float(min_after.value)*60))
                        bulk_stat.append(stat.code)
                        stored = True
    # Download data
    print('Begin data download')
    st = client.get_waveforms_bulk(bulk)
    print(st)
    print('Download complete; Pre-processing data')
    
    if radio_button_group.active == 0: 
        freq_resample = float(resample_delta.value)
        freqmin_high  = np.min([float(highfilter_left.value),float(highfilter_right.value)])
        freqmax_high  = np.max([float(highfilter_left.value),float(highfilter_right.value)])
        freqmin_mid   = np.min([float(midfilter_left.value), float(midfilter_right.value)])
        freqmax_mid   = np.max([float(midfilter_left.value), float(midfilter_right.value)])
        freqmin_low   = np.min([float(lowfilter_left.value), float(lowfilter_right.value)])
        freqmax_low   = np.max([float(lowfilter_left.value), float(lowfilter_right.value)])
    else:
        freq_resample = 1/float(resample_delta.value)
        freqmin_high  = np.min([1/float(highfilter_left.value),1/float(highfilter_right.value)])
        freqmax_high  = np.max([1/float(highfilter_left.value),1/float(highfilter_right.value)])
        freqmin_mid   = np.min([1/float(midfilter_left.value), 1/float(midfilter_right.value)])
        freqmax_mid   = np.max([1/float(midfilter_left.value), 1/float(midfilter_right.value)])
        freqmin_low   = np.min([1/float(lowfilter_left.value), 1/float(lowfilter_right.value)])
        freqmax_low   = np.max([1/float(lowfilter_left.value), 1/float(lowfilter_right.value)])
        
    if freq_resample<(2*np.max([freqmin_high,freqmax_high,freqmin_mid,freqmax_mid,freqmin_low,freqmax_low])):
        print('Filter frequency exceeds Nyquist frequency')
    
    st.resample(freq_resample) # downsample
    st.detrend() # detrend
    #st.rotate('->ZNE',inventory=inventory) # rotates channels **1 and **2 to **N and **E
    # divide between raw, high, mid, and low frequency bands
    st_raw = st.copy()
    # Initialize data storage structures
    meta = np.array([['Network','Station','Lat','Lon','Azimuth','Distance',
                      'Frequency','Channel']]) # channel metadata
    time = np.array([[pd.to_datetime(str(st_raw[0].stats.starttime)) + pd.Timedelta(st_raw[0].stats.delta*i,unit='s') for i in np.arange(st_raw[0].stats.npts)]]) # time data
    data = np.array([np.arange(st_raw[0].stats.npts)*st_raw[0].stats.delta]) # time data
    #worked = 0
    per_done = 0
    for k,stat in enumerate(bulk_stat): # loop over each station
        # get metadata
        try:
            st_raw_i = st_raw.select(station=stat).rotate('->ZNE',inventory=inventory)
            if len(st_raw_i)==3:
                net = st_raw_i[0].stats.network
                inv = inventory.select(network=net, station=stat)
                station_info = inv.networks[0].stations[0]
                #station_info = inventory.select(network=net,station=stat).networks[0].stations[0]
                stlat = station_info.latitude 
                stlon = station_info.longitude
        
                # calculate distance, azimuth, abd back-azimuth
                [dist, az, back_az] = gps2dist_azimuth(float(input_lat.value),float(input_lon.value),stlat,stlon)
    
                # rotate channels to ZRT
                st_raw_i = st_raw_i.rotate('NE->RT',back_azimuth=back_az)
                st_high_i = st_raw_i.copy().filter('bandpass',freqmin=freqmin_high, freqmax=freqmax_high)
                st_mid_i = st_raw_i.copy().filter('bandpass',freqmin=freqmin_mid, freqmax=freqmax_mid)
                st_low_i = st_raw_i.copy().filter('bandpass',freqmin=freqmin_low, freqmax=freqmax_low)
    
                # loop over channels
                for channel in ['BHZ','BHR','BHT']:
                    for freq,freq_name in zip([st_raw_i,st_high_i,st_mid_i,st_low_i],['Raw','High','Mid','Low']):
                        meta = np.append(meta,[[net,stat,stlat,stlon,az,dist/111139,freq_name,channel]],axis=0)   
                        time = np.append(time,
                                         [np.array([pd.to_datetime(str(st_raw_i[0].stats.starttime)) + pd.Timedelta(st_raw_i[0].stats.delta*i,unit='s') for i in np.arange(st_raw_i[0].stats.npts)]
                                         ,dtype='datetime64')],
                                         axis=0)
                        data = np.append(data,[freq.select(channel=channel.replace('B','*'))[0].data],axis=0)
        except Exception as e:
            bulk_stat.remove(stat)
            print(e)

        if (k/len(bulk_stat)*100) > (per_done + 10 ):
            print("%2.0f %% done"% (per_done))
            per_done += 10

    # store in pandas dataframe
    df = pd.DataFrame(meta[1:,:], columns = np.array(meta[0,:],
                      dtype='str')).astype({'Lat':'float32','Lon':'float32',
                      'Azimuth':'float32',
                      'Distance':'float32'}).join(pd.Series(list(time[1:,:]),
                      name="Time")).join(pd.Series(list(data[1:,:]),name="Data"))
    
    print('100%% done')
    print('Gathering station data')
                                                    
    df['url'] = ''
    for index, stat in df.groupby('Station').first().iterrows():
        x,y = latlon2webmercator.transform(stat['Lat'],stat['Lon'])
        station_data.data = {'x':[x],'y':[y]}
        p3b.x_range.start = x-padding
        p3b.x_range.end = x+padding
        p3b.y_range.start = y-padding
        p3b.y_range.end = y+padding
    
        img = get_screenshot_as_png(p3b)
        im_file = BytesIO()
        img.save(im_file, format='png')
        im_bytes = im_file.getvalue()  # im_bytes: image in binary format.
        url = 'data:image/png;base64,' + base64.b64encode(im_bytes).decode('utf-8')
        df.loc[df['Station']==index,'url'] = url

    df = df.dropna()
    df.to_hdf('eventdat',key='data') # save to 'eventdat.h5'
    eventdat = pd.DataFrame(data={'ID':str(input_ID.value),
                                  'lat':float(input_lat.value), 
                                  'lon':float(input_lon.value), 
                                  'time':str(input_time.value), 
                                  'mag':float(input_mag.value[0]),
                                  'depth':float(input_depth.value)}, index=[0])
    eventdat.to_hdf('eventdat',key='meta')
    print('Loaded ' + str(int((df.shape[0])/12)) + ' stations')
    print('Saving to file')
    
    stat_lat = df['Lat'].values
    stat_lon = df['Lon'].values
    [stat_x, stat_y] = latlon2webmercator.transform(stat_lat, stat_lon)
    stations_source.data = {'stat_x': stat_x,'stat_y': stat_y,
                                             'stat_lat':stat_lat,'stat_lon':stat_lon,
                                             'stat_name':df['Station'].values}
    
    df['SNR'] = df['Data'].apply(lambda x: np.mean(x**2)/np.mean(x[:100]**2))
    
    full = df.copy()
    full = full.groupby(['Frequency','Channel'])
    
    bins = np.linspace(df['Distance'].min()-0.1,df['Distance'].max()+0.1,60)
    if np.isnan(bins).any():
        bins = np.linspace(0,180,60)
    df['binned'] = pd.cut(df['Distance'], bins) # bin

    binned_dist = df.loc[df.groupby(['Frequency','Channel','binned'])['SNR'].agg(
                    lambda x : np.nan if x.count() == 0 else x.idxmax()
                ).dropna().sort_values().values].drop(columns=['binned']).groupby(
                    ['Frequency','Channel'])
    
    bins = np.linspace(df['Azimuth'].min()-0.1,df['Azimuth'].max()+0.1,60)
    if np.isnan(bins).any():
        bins = np.linspace(0,360,60)
    df['binned'] = pd.cut(df['Azimuth'], bins) # bin

    binned_az = df.loc[df.groupby(['Frequency','Channel','binned'])['SNR'].agg(
                    lambda x : np.nan if x.count() == 0 else x.idxmax()
                ).dropna().sort_values().values].drop(columns=['binned']).groupby(
                    ['Frequency','Channel'])
    
    [group,mapper] = update_plotting_data(freq_select,filled_select,amplitude_slider,
                         normalize_select,sort_opts,sort_select,color_by_az,
                         azimuth_range,binning_select,phase_cheatsheet,color_bar,
                         full, binned_dist, binned_az)
    
    source_records.data = group
    
    model = TauPyModel(model=input_velmodel.value)
    arr_dist = {}
    arr_time = {}

    branch_dist = {}
    branch_time = {}

    [dists, dd]  = np.linspace(epicentral_distance.value[0],
                               epicentral_distance.value[1],100,retstep=True)
    
    tol = 1e-5
    for dist in dists:
        arrivals = model.get_travel_times(source_depth_in_km=float(input_depth.value[0]),
                                  distance_in_degree=dist)
        for arrival in arrivals: 
            try:
                arr_time[arrival.name] = np.append(arr_time[arrival.name],arrival.time)
                arr_dist[arrival.name] = np.append(arr_dist[arrival.name],arrival.distance)
            except: 
                arr_time[arrival.name] = arrival.time 
                arr_dist[arrival.name] = arrival.distance 
    
    for item in arr_dist: 
        dict_t = {}
        ind = 0
        if type(arr_dist[item])==np.ndarray:
            unique, counts = np.unique(arr_dist[item], return_counts=True)
            for d, n in zip(unique, counts): 
                dict_t[d] = arr_time[item][ind:ind+n]
                ind += n
        else:
            dict_t[arr_dist[item]] = np.array([arr_time[item]])
            unique = [arr_dist[item]]
            counts = [1]
        
        N = -1
        no = 0
        do = 0
        map_branch = {0:0}
        for (d,n) in zip(unique,counts):
            if N == -1: # initial
                N = 0
                for i in np.arange(n):
                    map_branch[i] = i
            elif (d - do - dd) > tol: # gap
                more = True
                while more: 
                    try: 
                        branch_dist[item + '_' + str(N)]
                        N += 1
                    except: 
                        more = False
                map_branch = {0:0}
                for i in np.arange(n):
                    map_branch[i] = i
            elif (no - n) > (1 - tol): # decrease in branches
                tf = np.zeros((no,n))
                map_brancho = map_branch
                map_branch = {}
                for i in np.arange(no):
                    try:
                        v = (branch_time[item + '_' + str(N+i)][-1] - branch_time[item + '_' + str(N+i)][0])/(branch_dist[item + '_' + str(N+i)][-1] - branch_dist[item + '_' + str(N+i)][0])
                        if np.isnan(v):
                            v = 0
                    except: 
                        v = 0
                    tf[i,:] = v*dd + branch_time[item + '_' + str(N+i)][-1] - dict_t[d] 
                for i in np.arange(n):
                    (old,new) = np.unravel_index(np.nanargmin(np.abs(tf), axis=None), tf.shape)
                    map_branch[new] = map_brancho[old] # add to the branch with the closest slope
                    tf[old,:] = np.nan
                    tf[:,new] = np.nan

            elif (no - n) < (-1 + tol): # increase in branches
                tf = np.zeros((no,n))
                map_brancho = map_branch
                map_branch = {}
                for i in np.arange(no):
                    try:
                        v = (branch_time[item + '_' + str(N+i)][-1] - branch_time[item + '_' + str(N+i)][0])/(branch_dist[item + '_' + str(N+i)][-1] - branch_dist[item + '_' + str(N+i)][0])
                        if np.isnan(v):
                            try:
                                v = dict_t[d+dd] - dict_t[d]/(-dd)
                            except:
                                v = 0
                    except: 
                        v = 0
                    tf[i,:] = v*dd + branch_time[item + '_' + str(N+i)][-1] - dict_t[d]
                for i in np.arange(no):
                    (old,new) = np.unravel_index(np.nanargmin(np.abs(tf), axis=None), tf.shape)
                    map_branch[new] = map_brancho[old] # add to the branch with the closest slope
                    tf[old,:] = np.nan
                    tf[:,new] = np.nan
                k = 0
                for j in np.arange(n):
                    if j not in map_branch:
                        more = True
                        i = 1
                        while more: 
                            try: 
                                branch_dist[item + '_' + str(N + i)]
                                i += 1
                            except: 
                                map_branch[j] = i + k
                                more = False
                        k += 1
                        
            for i in np.arange(n): 
                try:
                    branch_dist[item + '_' + str(N+map_branch[i])] = np.append(branch_dist[item + '_' + str(N+map_branch[i])],[d])
                    branch_time[item + '_' + str(N+map_branch[i])] = np.append(branch_time[item + '_' + str(N+map_branch[i])],[dict_t[d][i]])
                except: 
                    branch_dist[item + '_' + str(N+map_branch[i])] = np.array([d]) # starts new branch
                    branch_time[item + '_' + str(N+map_branch[i])] = np.array([dict_t[d][i]])   
            no = n
            do = d
         
    df_arr = pd.concat([pd.DataFrame.from_dict({k: [v] for k, v in branch_dist.items()},orient='index',columns=['Distance']),
           pd.DataFrame.from_dict({k: [v] for k, v in branch_time.items()},orient='index',columns=['Time'])], axis=1)
    df_arr.reset_index(inplace=True)
    df_arr['Phase'] = df_arr['index'].apply(lambda x: x.split('_')[0])
    df_arr['Time'] = df_arr['Time'].apply(lambda x:np.array([pd.to_datetime(str(t11)) + pd.Timedelta(i,unit='sec') for i in x]))
    df_arr.to_hdf('eventdat',key='arrivals')       
    df = df_arr
    update_plotting_data(freq_select,filled_select,amplitude_slider,
                         normalize_select,sort_opts,sort_select,color_by_az,
                         azimuth_range,binning_select,phase_cheatsheet,color_bar,
                         full, binned_dist, binned_az)
    
    print('Downloading regional seismicity')
    lat_bot = float(input_lat.value)-5 #meta.lat[0] - r
    lat_top = float(input_lat.value)+5 #meta.lat[0] + r
    lon_left = float(input_lon.value)-5 #meta.lon[0] - r
    lon_right = float(input_lon.value)+5 #meta.lon[0] + r

    df_reg = 0
    for CMTurl in CMTurls:
        df_reg_i = focal_mech(CMTurl,lat_bot,lat_top,lon_left,lon_right,Mw_min,latlon2webmercator,ax)
        if type(df_reg) == pd.core.frame.DataFrame:
            df_reg = df_reg.append(df_reg_i)
        else: 
            df_reg = df_reg_i.copy()
    df_regq = focal_mech(CMTurl_quick,lat_bot,lat_top,lon_left,lon_right,Mw_min,latlon2webmercator,ax)
    
    df_mt = 0
    for CMTurl in np.append(CMTurls,CMTurl_quick):
        df_mt_i = focal_mech(CMTurl,float(input_lat.value)-0.05,float(input_lat.value)+0.05,
                             float(input_lon.value)-0.05,float(input_lon.value)+0.05,Mw_min,latlon2webmercator,ax)
        if type(df_mt) == pd.core.frame.DataFrame:
            df_mt = df_mt.append(df_mt_i)
        else: 
            df_mt = df_mt_i.copy()
    
    if len(df_reg)>0:
        df_reg['P_x'] = df_reg.apply(lambda x: [-sind(x['P'][1])*cosd(x['P'][0])*stax+x['x'], sind(x['P'][1])*cosd(x['P'][0])*stax+x['x']],axis=1)
        df_reg['P_y'] = df_reg.apply(lambda x: [-cosd(x['P'][1])*cosd(x['P'][0])*stax+x['y'], cosd(x['P'][1])*cosd(x['P'][0])*stax+x['y']],axis=1)
        df_reg['T_x'] = df_reg.apply(lambda x: [-sind(x['T'][1])*cosd(x['T'][0])*stax+x['x'], sind(x['T'][1])*cosd(x['T'][0])*stax+x['x']],axis=1)
        df_reg['T_y'] = df_reg.apply(lambda x: [-cosd(x['T'][1])*cosd(x['T'][0])*stax+x['y'], cosd(x['T'][1])*cosd(x['T'][0])*stax+x['y']],axis=1)

    if len(df_regq)>0:
        df_regq['P_x'] = df_regq.apply(lambda x: [-sind(x['P'][1])*cosd(x['P'][0])*stax+x['x'], sind(x['P'][1])*cosd(x['P'][0])*stax+x['x']],axis=1)
        df_regq['P_y'] = df_regq.apply(lambda x: [-cosd(x['P'][1])*cosd(x['P'][0])*stax+x['y'], cosd(x['P'][1])*cosd(x['P'][0])*stax+x['y']],axis=1)
        df_regq['T_x'] = df_regq.apply(lambda x: [-sind(x['T'][1])*cosd(x['T'][0])*stax+x['x'], sind(x['T'][1])*cosd(x['T'][0])*stax+x['x']],axis=1)
        df_regq['T_y'] = df_regq.apply(lambda x: [-cosd(x['T'][1])*cosd(x['T'][0])*stax+x['y'], cosd(x['T'][1])*cosd(x['T'][0])*stax+x['y']],axis=1)

    if len(df_mt)>0:
        df_mt['P_x'] = df_mt.apply(lambda x: [-sind(x['P'][1])*cosd(x['P'][0])*stax+x['x'], sind(x['P'][1])*cosd(x['P'][0])*stax+x['x']],axis=1)
        df_mt['P_y'] = df_mt.apply(lambda x: [-cosd(x['P'][1])*cosd(x['P'][0])*stax+x['y'], cosd(x['P'][1])*cosd(x['P'][0])*stax+x['y']],axis=1)
        df_mt['T_x'] = df_mt.apply(lambda x: [-sind(x['T'][1])*cosd(x['T'][0])*stax+x['x'], sind(x['T'][1])*cosd(x['T'][0])*stax+x['x']],axis=1)
        df_mt['T_y'] = df_mt.apply(lambda x: [-cosd(x['T'][1])*cosd(x['T'][0])*stax+x['y'], cosd(x['T'][1])*cosd(x['T'][0])*stax+x['y']],axis=1)

    df_reg['radius'] = df_reg['Mw']*np.min([lat_top-lat_bot,lon_right-lon_left])*10**tensor_size_exp.value
    df_regq['radius'] = df_regq['Mw']*np.min([lat_top-lat_bot,lon_right-lon_left])*10**tensor_size_exp.value
    df_mt['radius'] = df_mt['Mw']*np.min([lat_top-lat_bot,lon_right-lon_left])*10**tensor_size_exp.value
    
    df_reg.to_hdf('eventdat',key='regional')
    df_regq.to_hdf('eventdat',key='regional_quick')
    df_mt.to_hdf('eventdat',key='moment_tensor')
    
    source_reg.data = df_reg
    source_regq.data = df_regq
    source_regq.data = df_mt
    
    [xl, yb] = latlon2webmercator.transform(lat_bot, lon_left)
    [xr, yt] = latlon2webmercator.transform(lat_top, lon_right)
    
    p_reg.x_range=Range1d(xl, xr)
    p_reg.y_range=Range1d(yb, yt)

    print('Complete')
           
download_data.on_click(download_data_callback)

layout2 = row(column(Div(text='<h1>Downloading the data parameters<h1>'),
                     network,
                     row(add_network_text,add_network_button),
                     epicentral_distance,
                     Div(text='Download waveform length:'),
                     row(children=[min_before, min_after,phase_shift],sizing_mode='scale_width'),
                     Div(text='Waveform processing parameters:'),
                     radio_button_group,
                     Div(text='Low filter:'),
                     row(lowfilter_left,lowfilter_right),
                     Div(text='Mid filter:'),
                     row(midfilter_left,midfilter_right),
                     Div(text='High filter:'),
                     row(highfilter_left,highfilter_right),
                     resample_delta),
                                                             column(p2,Div(text='Figure 2: Map of selected stations and event. Move the slider to change the radius from the event.'),load_stations,download_data))
panel2 = Panel(child=layout2,title='Download data')

###############################################################################
# Display records
###############################################################################
#Load data if reading from saved file
channel_select = Select(value="Vertical (BHZ)", 
                        options=['Vertical (BHZ)', 'Radial (BHR)',
                                 'Transverse (BHT)'],
                        width_policy='min', min_width=150,
                        margin=(5, 5, 5, 100))

if radio_button_group.active == 0: 
    freq_resample = float(resample_delta.value)
    freqmin_high  = np.min([float(highfilter_left.value),float(highfilter_right.value)])
    freqmax_high  = np.max([float(highfilter_left.value),float(highfilter_right.value)])
    freqmin_mid   = np.min([float(midfilter_left.value), float(midfilter_right.value)])
    freqmax_mid   = np.max([float(midfilter_left.value), float(midfilter_right.value)])
    freqmin_low   = np.min([float(lowfilter_left.value), float(lowfilter_right.value)])
    freqmax_low   = np.max([float(lowfilter_left.value), float(lowfilter_right.value)])
    freq_options = ['Raw', 
                    'High (' + str(freqmin_high) + '-' + str(freqmax_high) + 'Hz)', 
                    'Mid (' + str(freqmin_mid) + '-' + str(freqmax_mid) + 'Hz)',
                    'Low (' + str(freqmin_low) + '-' + str(freqmax_low) + 'Hz)']
else:
    freq_resample = 1/float(resample_delta.value)
    freqmin_high  = np.min([1/float(highfilter_left.value),1/float(highfilter_right.value)])
    freqmax_high  = np.max([1/float(highfilter_left.value),1/float(highfilter_right.value)])
    freqmin_mid   = np.min([1/float(midfilter_left.value), 1/float(midfilter_right.value)])
    freqmax_mid   = np.max([1/float(midfilter_left.value), 1/float(midfilter_right.value)])
    freqmin_low   = np.min([1/float(lowfilter_left.value), 1/float(lowfilter_right.value)])
    freqmax_low   = np.max([1/float(lowfilter_left.value), 1/float(lowfilter_right.value)])
    freq_options = ['Raw', 
                    'High (' + str(np.min([float(highfilter_left.value),float(highfilter_right.value)])) + '-' + str(np.max([float(highfilter_left.value),float(highfilter_right.value)])) + 's)', 
                    'Mid (' + str(np.min([float(midfilter_left.value),float(midfilter_right.value)])) + '-' + str(np.max([float(midfilter_left.value),float(midfilter_right.value)])) + 's)',
                    'Low (' + str(np.min([float(lowfilter_left.value),float(lowfilter_right.value)])) + '-' + str(np.max([float(lowfilter_left.value),float(lowfilter_right.value)])) + 's)']

freq_select = Select(value=freq_options[0], options=freq_options, width_policy='min',
                     min_width=150)
reduction_velocity = TextInput(value='0', width_policy='min', min_width=100,margin=(5, 10, 5, 10))
filled_select = CheckboxGroup(labels=['Fill above zero'], active=[],margin=(15, 10, 5, 10))
amplitude_slider = Slider(start=0.1, end=10, value=1, step=1, title='Relative amplitude',
                          orientation='vertical', direction='rtl',margin=(50, 5, 5, 10))
normalize_select = Select(options=['Individual', 'Global'], value='Individual',
                          width_policy='min',min_width=100,margin=(15, 5, 5, 5))
sort_opts = ['Distance', 'Azimuth']
sort_select = Select(options=sort_opts, value=sort_opts[0], width_policy='min',
                     min_width=100,margin=(300, 5, 5, 5))
color_by_az = CheckboxGroup(labels=['Color by azimuth'], active=[],margin=(15, 10, 5, 10))
azimuth_range = RangeSlider(start=0,end=360,value=(0,360),title='Azimuth range (deg)',margin=(15, 10, 5, 10))
binning_select = CheckboxGroup(labels=['Binned'], active=[0],margin=(15, 10, 5, 10)) # binned data plots one trace for each bin that has the lowest SNR

phase_cheatsheet = CheckboxGroup(labels=['Phase cheatsheet'], active=[],margin=(15, 10, 5, 10))
phases_short = ['P','Pdiff','S','Sdiff','SP','PP','PPP','SS','SSS','SSP','PSP',
                'SKS','SKP','SKKS','SKKKS','ScSScS','ScSScSScS','ScSScSScSScS',
                'ScS','PcP','PKIKP','PKKP','PKIKKIKP','PKP','PKPPKP','PKIKP']
    
df['SNR'] = df['Data'].apply(lambda x: np.mean(x**2)/np.mean(x[:30]**2))
    
full = df.copy()
full = full.groupby(['Frequency','Channel'])
    
bins = np.linspace(df['Distance'].min()-0.1,df['Distance'].max()+0.1,60)
if np.isnan(bins).any():
    bins = np.linspace(0,180,60)
df['binned'] = pd.cut(df['Distance'], bins) # bin
try:
    binned_dist = df.loc[df.groupby(['Frequency','Channel','binned'])['SNR'].agg(
                lambda x : np.nan if x.count() == 0 else x.idxmax()
            ).dropna().sort_values().values].drop(columns=['binned']).groupby(
                ['Frequency','Channel'])
except: 
    binned_dist = df.copy().groupby(['Frequency','Channel'])
    
bins = np.linspace(df['Azimuth'].min()-0.1,df['Azimuth'].max()+0.1,60)
if np.isnan(bins).any():
    bins = np.linspace(0,360,60)
df['binned'] = pd.cut(df['Azimuth'], bins) # bin

try:
    binned_az = df.loc[df.groupby(['Frequency','Channel','binned'])['SNR'].agg(
              lambda x : np.nan if x.count() == 0 else x.idxmax()
            ).dropna().sort_values().values].drop(columns=['binned']).groupby(
                    ['Frequency','Channel'])
except: 
    binned_az = df.copy().groupby(['Frequency','Channel'])
    
binned = binned_dist
group = binned.get_group((freq_select.value.split(' ')[0],channel_select.value[-4:-1]))

station_data = ColumnDataSource(data={'x':[],'y':[]})

try:
    arrival_data = pd.read_hdf('eventdat',key='arrivals')
except: 
    arrival_data = {'index':[],'Distance':[], 'Time':[], 'Phase':[]} 
source_arrivals = ColumnDataSource({'index':[],'Distance':[], 'Time':[], 'Phase':[]})

p3b = figure(plot_height=200,plot_width=200,x_axis_type='mercator',y_axis_type= 'mercator',
                 x_range=(-padding, padding),y_range=(-padding, padding),tools='')
p3b.add_tile(get_provider('ESRI_IMAGERY'))
p3b.triangle('x','y',source=station_data,size=20,color='red',line_color='black')
    
mapper = linear_cmap(field_name='Azimuth', palette=grey(1) ,low=df['Azimuth'].min() ,high=df['Azimuth'].max())
color_bar = ColorBar(color_mapper=mapper['transform'],visible=False)

def update_plotting_data(freq_select,filled_select,amplitude_slider,
                         normalize_select,sort_opts,sort_select,color_by_az,
                         azimuth_range,binning_select,phase_cheatsheet,color_bar,
                         full, binned_dist, binned_az):
    
    sort_by = sort_select.value
    
    if 0 in binning_select.active:
        if 'Distance' in sort_by:
            bins = np.linspace(df['Distance'].min()-0.1,df['Distance'].max()+0.1,60)
            if np.isnan(bins).any():
                bins = np.linspace(0,180,60)
            binned = binned_dist
        else:
            bins = np.linspace(df['Azimuth'].min()-0.1,df['Azimuth'].max()+0.1,60)
            if np.isnan(bins).any():
                bins = np.linspace(0,360,60)
            binned = binned_az
    else: 
        binned = full
        bins = np.linspace(0,180,60)
    
    try:
        group = binned.get_group((freq_select.value.split(' ')[0],channel_select.value[-4:-1]))
    except:
        group = full.get_group((freq_select.value.split(' ')[0],channel_select.value[-4:-1]))
    group = group[group['Azimuth'].between(azimuth_range.value[0], azimuth_range.value[1], inclusive=True)]
    
    arrivals_temp = arrival_data.copy()
    
    norm_stretch = 0.02*(bins[-1]-bins[0])*amplitude_slider.value

    if 'Individual' in normalize_select.value:
        group['plot_trace'] = group[sort_by] + group['Data'].apply(lambda x: -x*norm_stretch/np.max(np.abs(x))) # normalize by max value of each trace
    else: 
        global_max = group['Data'].agg(lambda x: np.max(np.abs(x))).max()
        group['plot_trace'] = group[sort_by] + group['Data'].apply(lambda x: -x*norm_stretch/global_max)
    
         
    vel = float(reduction_velocity.value)/111.139
    if ('Distance' in sort_by) & (vel>1e-5):
        group['Time'] = group.apply(lambda x: x['Time'] - pd.Timedelta(x['Distance']/vel,unit='sec'),axis=1)
        source_drawline.data['xs'] = [[source_drawline.data['xs'][0][0],source_drawline.data['xs'][0][0]]]
        arrivals_temp['Time'] = arrivals_temp.apply(lambda x: [t - pd.Timedelta(d/vel,unit='sec') for (t,d) in zip(x['Time'],x['Distance'])],axis=1)
        
    if 0 in filled_select.active:
        group['Fills'] = (group['plot_trace']-group[sort_by]).apply(lambda x: np.max([x,np.zeros_like(x)],axis=0)) + group[sort_by]
    else:
        group['Fills'] = 0*group['plot_trace'] + group[sort_by]
        
    if 0 in color_by_az.active:
        mapper['transform'].palette=Viridis256
        mapper['transform'].low=group['Azimuth'].min()
        mapper['transform'].high=group['Azimuth'].max()
        color_bar.visible = True
    else: 
        mapper['transform'].palette=grey(1)
        mapper['transform'].low=group['Azimuth'].min()
        mapper['transform'].high=group['Azimuth'].max()
        color_bar.visible = False
        
    return group, mapper

def update_cheatsheet_plotting(phase_cheatsheet):
    vel = float(reduction_velocity.value)/111.139
    arrivals_temp = arrival_data.copy()
    #if ('Distance' in sort_by) & (vel>1e-5):
    #    arrivals_temp['Time'] = arrivals_temp.apply(lambda x: [t - pd.Timedelta(d/vel,unit='sec') for (t,d) in zip(x['Time'],x['Distance'])],axis=1)
        
    if ('Distance' in sort_by) & (0 in phase_cheatsheet.active):
        arrivals_group = arrivals_temp.groupby(['Phase'])
        phases = []
        for phase in phases_short: 
            try:
                phases = phases.append(arrivals_group.get_group(phase))
            except: 
                try: 
                    phases = arrivals_group.get_group(phase)
                except:
                    phase
        source_arrivals.data = phases
    else: 
        source_arrivals.data = {'index':[],'Distance':[], 'Time':[], 'Phase':[]}
    return group, mapper


[group,mapper] = update_plotting_data(freq_select,filled_select,amplitude_slider,
                         normalize_select,sort_opts,sort_select,color_by_az,
                         azimuth_range,binning_select,phase_cheatsheet,color_bar,
                         full, binned_dist, binned_az)
source_records = ColumnDataSource(group)
    
sort_by = sort_select.value
p3 = figure(plot_height=600,plot_width=1200,
            x_range=[pd.to_datetime(str(t11))-pd.Timedelta(float(min_before.value)*60,unit='sec'), 
                     pd.to_datetime(str(t11))+pd.Timedelta(float(min_after.value)*60,unit='sec')],         
            x_axis_type='datetime',y_axis_label= 'degrees',
            x_axis_label='Time (min)',
            tools='box_zoom,undo,redo,reset,save', active_drag='box_zoom')
p3.y_range.flipped = True
p3.xaxis[0].formatter = DatetimeTickFormatter(minutes=['%H:%M'])
p3.patches(xs='Time', ys='Fills', source=source_records, fill_color='blue',line_color='white',line_alpha=1)
tr = p3.multi_line(xs='Time', ys='plot_trace', line_width=1, line_color=mapper,
             source=source_records)
p3.add_layout(color_bar, 'right')
ar = p3.multi_line('Time', 'Distance',line_color='red', line_width=1,
             source=source_arrivals)

TOOLTIPS3 = """
    <div>
        <div>
            <img
                src="@url" height="100" alt="@url" width="100"
                style="float: left; margin: 0px 5px 5px 0px;"
            ></img>
        </div>
        <div>
            <span style="font-size: 12px"> @Network.@Station</span>
        </div>
        <div>
            <span style="font-size: 12px"> Lat: @Lat</span>
        </div>
        <div>
            <span style="font-size: 12px"> Lon: @Lon</span>
        </div>
        <div>
            <span style="font-size: 12px"> Az: @Azimuth</span>
        </div>
        <div>
            <span style="font-size: 12px"> Dis: @Distance</span>
        </div>
    </div>
"""

TOOLTIPS3c = [
    ('Phase', '@Phase')
]

p3.add_tools(HoverTool(tooltips=TOOLTIPS3, renderers=[tr], mode='mouse',description='Station information'))
p3.toolbar.active_inspect = None

source_drawfree = ColumnDataSource({'xs':[],'ys':[]})
draw_rfree = p3.multi_line('xs', 'ys',line_color='red', line_width=3,
             source=source_drawfree)
source_drawline = ColumnDataSource({'xs':[],'ys':[]})
draw_rline = p3.multi_line('xs', 'ys',line_color='red', line_width=3,
             source=source_drawline)
freehanddraw = FreehandDrawTool(renderers=[draw_rfree])
polydraw = PolyDrawTool(renderers=[draw_rline],num_objects=1,description='Velocity reduction')
p3.add_tools(freehanddraw,polydraw)

p3.add_tools(HoverTool(tooltips=TOOLTIPS3c, renderers=[ar], mode='mouse',description='Phase cheatsheet'))

def velocity_reduce(attrs,new,old):
    try: 
        source_drawline.data['xs'][0]
        if len(source_drawline.data['xs'][0])>2:
            source_drawline.data['xs'] = [source_drawline.data['xs'][0][-2:]]
            source_drawline.data['ys'] = [source_drawline.data['ys'][0][-2:]]
        try: 
            vel = (source_drawline.data['ys'][0][1]-source_drawline.data['ys'][0][0])*111.139/((source_drawline.data['xs'][0][1]-source_drawline.data['xs'][0][0])/1000)
        except:
            vel = 0
        if (vel/111.139)>1e-5:
            reduction_velocity.value = '%.2f' % vel
    except: 
        reduction_velocity.value = '0'

source_drawline.on_change('data',velocity_reduce)

def update_display(attrname, old, new):
    
    [group,mapper] = update_plotting_data(freq_select,filled_select,amplitude_slider,
                         normalize_select,sort_opts,sort_select,color_by_az,
                         azimuth_range,binning_select,phase_cheatsheet,color_bar,
                         full, binned_dist, binned_az)
    source_records.data = group
    
def update_cheatsheet(attrname, old, new):
    
    [group,mapper] = update_cheatsheet_plotting(phase_cheatsheet)
    source_records.data = group
        
for u in [reduction_velocity,channel_select,freq_select,sort_select,
          amplitude_slider,normalize_select,azimuth_range]:
    u.on_change('value', update_display)
for v in [filled_select,color_by_az,binning_select]:
    v.on_change('active',update_display)
    
reduction_velocity.on_change('value',update_cheatsheet)
phase_cheatsheet.on_change('active',update_cheatsheet)
    
#change_y_axis_label = CustomJS(args=dict(plot=p3, source=source_records, sort_by=sort_select, sort_opts = sort_opts, ax=p3.yaxis), code="""
#    ax[0].axis_label = sort_opts[sort_by.active] + ' (deg)';
#    source.change.emit();
#""")

#sort_select.js_on_change('active', change_y_axis_label)
    
layout3 = row(column(sort_select),
    column(row(channel_select,freq_select),
           p3,
           Div(text='<hr style="width:1200px">'),
           row(column(Div(text='Velocity reduction'),
               row(reduction_velocity,Div(text='km/s'))),
               azimuth_range,
               column(color_by_az,
                      filled_select),
               column(binning_select,
                      phase_cheatsheet))
                ), 
    column(amplitude_slider,
           Div(text='Normalize:'),
                 normalize_select
           ))

panel3 = Panel(child=layout3,title='Display records')

#layout4 = column(Div(text='Velocity reduction'),
#                 row(reduction_velocity,Div(text='km/s')),
#                 Div(text='Virtualization'),
#                 filled_select,
#                 amplitude_slider,
#                 Div(text='Normalize:'),
#                 normalize_select,
#                 Div(text='Sort by:'),
#                 sort_select,
#                 color_by_az,
#                 azimuth_range,
#                 binning_select,
#                 Div(text='Phase cheatsheet'),
#                 phase_cheatsheet,
#                 Div(text='Sonification'))
#panel4 = Panel(child=layout4, title='Display controls')

###############################################################################
# Display regional maps
###############################################################################
def cosd(theta):
    return np.cos(np.deg2rad(theta))
    
def sind(theta):
    return np.sin(np.deg2rad(theta))

lat_bot = float(input_lat.value)-5 #meta.lat[0] - r
lat_top = float(input_lat.value)+5 #meta.lat[0] + r
lon_left = float(input_lon.value)-5 #meta.lon[0] - r
lon_right = float(input_lon.value)+5 #meta.lon[0] + r

# Setup search parameter inputs
min_year = TextInput(value='2010', title='Starting year:')
max_year = TextInput(value='2022', title='Ending year:')
min_depth = TextInput(value='0', title='Minimum depth (km):')
max_depth = TextInput(value='300', title='Maximum depth (km):')
min_mag = TextInput(value='5', title='Minimum Mw:')
max_mag = TextInput(value='9', title='Maximum Mw:')
Mw_min = 4

tensor_size_exp = Slider(start=1,end=5,value=3.3, title='Tensor size')
tensor_size = np.min([lat_top-lat_bot,lon_right-lon_left])*10**tensor_size_exp.value
stax = np.min([lat_top-lat_bot,lon_right-lon_left])*10**(tensor_size_exp.value+0.5)

CMTurls = ['http://www.ldeo.columbia.edu/~gcmt/projects/CMT/catalog/jan76_dec20.ndk']

months = ['jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec']
UTCDateTime.utcnow()
for yr in np.arange(2021,UTCDateTime.utcnow().year+1):
    for mon in np.arange(0,UTCDateTime.utcnow().month):
        cat = 'https://www.ldeo.columbia.edu/~gcmt/projects/CMT/catalog/NEW_MONTHLY/' + str(
            UTCDateTime.utcnow().year) + '/' + months[mon] + str(UTCDateTime.utcnow().year)[2:] + '.ndk'
        CMTurls = np.append(CMTurls,cat)
        
CMTurl_quick = 'https://www.ldeo.columbia.edu/~gcmt/projects/CMT/catalog/NEW_QUICK/qcmt.ndk'

def focal_mech(CMTurl,lat_bot,lat_top,lon_left,lon_right,Mw_min,latlon2webmercator,ax):
    try: 
        events = pd.read_table(CMTurl).values
        while np.remainder(len(events),5) > 0: #there's a missing line in the 1976-2017 catalog
            events = np.append(['0000/0'],events)
        events = events.reshape(-1,5)
        def in_bounds(lat,lon,lat_bot,lat_top,lon_left,lon_right): 
            if (lon_left<-180) | (lon_right>180):
                in_x = (lon_i>(lon_left+360))|(lon_i<(lon_right-360))
            else:
                in_x = (lon_i>lon_left)&(lon_i<lon_right)
                in_y = (lat_i>lat_bot)&(lat_i<lat_top)
            return (in_x & in_y)
    
        df_reg_i = pd.DataFrame(data={'year':[],'lat':[],'lon':[],'depth':[],'Mw':[],'mt':[],'P':[],'T':[]})
        for ev_row in events:
            try:
                year_i = np.array(ev_row[0].split('/')[0][-4:],dtype=float)
                lat_i, lon_i, depth_i = np.array(ev_row[2].split()[3:8:2],dtype=float)
                x_i, y_i = latlon2webmercator.transform(lat_i,lon_i)
                Mw_i = 2/3*np.log10(float(ev_row[4].split()[10]) * 10 ** int(ev_row[3].split()[0]))-10.7
                if in_bounds(lat_i,lon_i,lat_bot,lat_top,lon_left,lon_right) & (Mw_i>=Mw_min):
                    df_reg_i = df_reg_i.append({'year':year_i,'lat':lat_i,'lon':lon_i,'depth':depth_i,'Mw':Mw_i,
                                'mt':np.array(ev_row[3].split()[1::2],dtype=float),
                                'np1':np.array(ev_row[4].split()[11:14],dtype=int),
                                'P':np.array(ev_row[4].split()[8:10],dtype=int),
                                'T':np.array(ev_row[4].split()[2:4],dtype=int)},ignore_index=True)
            except Exception as e: 
                print(e)
    
        url = []
        if len(df_reg_i)>0:
            for index, ev_row in df_reg_i.iterrows():
                ax.clear()
                if 'QUICK' in CMTurl:
                    ball_color = 'r'
                else:
                    ball_color = 'b'
                if np.abs(lat_bot - lat_top)<=0.2:
                    ball_color = 'g'
                ball = beach(ev_row['mt'], xy=(0,0),width=2,size=50,facecolor=ball_color,zorder=0)
                ball2 = beach(ev_row['np1'], xy=(0,0),width=2,size=50,nofill=True,linewidth=6,zorder=1)
                ball2.set_linewidth=20
                ax.add_collection(ball,autolim=True)
                ax.add_collection(ball2, autolim=True)
                ax.set_aspect("equal")
                ax.set_xlim(-1,1)  
                ax.set_ylim(-1,1)
    
                ax.scatter(sind(ev_row['P'][1])*cosd(ev_row['P'][0]),cosd(ev_row['P'][1])*cosd(ev_row['P'][0]),1000,'k',zorder=2)
                ax.scatter(sind(ev_row['T'][1])*cosd(ev_row['T'][0]),cosd(ev_row['T'][1])*cosd(ev_row['T'][0]),1000,'k',zorder=2)

                fig.patch.set_facecolor('w') # instead of fig.patch.set_facecolor
                fig.patch.set_alpha(0)
                plt.axis('off')
                im_file = BytesIO()
                plt.savefig(im_file, facecolor=([1,1,1,0]),transparent=True)
                im_bytes = im_file.getvalue()
                url = url + ['data:image/png;base64,' + base64.b64encode(im_bytes).decode('utf-8')]
                ax.clear()
            df_reg_i['x'],df_reg_i['y'] = latlon2webmercator.transform(df_reg_i['lat'],df_reg_i['lon'])
        df_reg_i['url'] = url
    except Exception as e: 
        print(e)
        df_reg_i = pd.DataFrame(data={'year':[],'lat':[],'lon':[],'depth':[],'Mw':[],
                                'mt':[],'np1':[],
                                'P':[],'T':[],
                                'x':[],'y':[],
                                'url':[]})
    return df_reg_i

fig,ax = plt.subplots(figsize=(10,10))
try: 
    df_reg = pd.read_hdf('eventdat',key='regional')
    df_regq = pd.read_hdf('eventdat',key='regional_quick')
    df_mt = pd.read_hdf('eventdat',key='moment_tensor')
except: 
    df_reg = pd.DataFrame(data={'year':[],'lat':[],'lon':[],'depth':[],'Mw':[],
                                'mt':[],'np1':[],
                                'P':[],'T':[],
                                'x':[],'y':[],
                                'url':[],
                                'P_x':[],'P_y':[],
                                'T_x':[],'T_y':[],
                                'radius':[]})
    df_regq = pd.DataFrame(data={'year':[],'lat':[],'lon':[],'depth':[],'Mw':[],
                                'mt':[],'np1':[],
                                'P':[],'T':[],
                                'x':[],'y':[],
                                'url':[],
                                'P_x':[],'P_y':[],
                                'T_x':[],'T_y':[],
                                'radius':[]})
    df_mt = pd.DataFrame(data={'year':[],'lat':[],'lon':[],'depth':[],'Mw':[],
                                'mt':[],'np1':[],
                                'P':[],'T':[],
                                'x':[],'y':[],
                                'url':[],
                                'P_x':[],'P_y':[],
                                'T_x':[],'T_y':[],
                                'radius':[]})

source_reg = ColumnDataSource(df_reg)
source_regq = ColumnDataSource(df_regq)
source_mt = ColumnDataSource(df_mt)

bool_year_reg = [True if (np.float(year)>=np.float(min_year.value))&(np.float(year)<=np.float(max_year.value)) else False for year in source_reg.data['year']]
bool_depth_reg = [True if (np.float(depth)>=np.float(min_depth.value))&(np.float(depth)<=np.float(max_depth.value)) else False for depth in source_reg.data['depth']]
bool_mag_reg = [True if (np.float(mag)>=np.float(min_mag.value))&(np.float(mag)<=np.float(max_mag.value)) else False for mag in source_reg.data['Mw']]
view_reg = CDSView(source=source_reg,filters=[BooleanFilter(bool_year_reg),BooleanFilter(bool_depth_reg),BooleanFilter(bool_mag_reg)])

bool_year_regq = [True if (np.float(year)>=np.float(min_year.value))&(np.float(year)<=np.float(max_year.value)) else False for year in source_regq.data['year']]
bool_depth_regq = [True if (np.float(depth)>=np.float(min_depth.value))&(np.float(depth)<=np.float(max_depth.value)) else False for depth in source_regq.data['depth']]
bool_mag_regq = [True if (np.float(mag)>=np.float(min_mag.value))&(np.float(mag)<=np.float(max_mag.value)) else False for mag in source_regq.data['Mw']]
view_regq = CDSView(source=source_regq,filters=[BooleanFilter(bool_year_regq),BooleanFilter(bool_depth_regq),BooleanFilter(bool_mag_regq)])

bool_year_mt = [True if (np.float(year)>=np.float(min_year.value))&(np.float(year)<=np.float(max_year.value)) else False for year in source_mt.data['year']]
bool_depth_mt = [True if (np.float(depth)>=np.float(min_depth.value))&(np.float(depth)<=np.float(max_depth.value)) else False for depth in source_mt.data['depth']]
bool_mag_mt = [True if (np.float(mag)>=np.float(min_mag.value))&(np.float(mag)<=np.float(max_mag.value)) else False for mag in source_mt.data['Mw']]
view_mt = CDSView(source=source_mt,filters=[BooleanFilter(bool_year_mt),BooleanFilter(bool_depth_mt),BooleanFilter(bool_mag_mt)])

[xl, yb] = latlon2webmercator.transform(lat_bot, lon_left)
[xr, yt] = latlon2webmercator.transform(lat_top, lon_right)


def update_regional_seismicity(attr,old,new):
    bool_year_reg = [True if (np.float(year)>=np.float(min_year.value))&(np.float(year)<=np.float(max_year.value)) else False for year in source_reg.data['year']]
    bool_depth_reg = [True if (np.float(depth)>=np.float(min_depth.value))&(np.float(depth)<=np.float(max_depth.value)) else False for depth in source_reg.data['depth']]
    bool_mag_reg = [True if (np.float(mag)>=np.float(min_mag.value))&(np.float(mag)<=np.float(max_mag.value)) else False for mag in source_reg.data['Mw']]
    view_reg.filters=[BooleanFilter(bool_year_reg),BooleanFilter(bool_depth_reg),BooleanFilter(bool_mag_reg)]

    bool_year_regq = [True if (np.float(year)>=np.float(min_year.value))&(np.float(year)<=np.float(max_year.value)) else False for year in source_regq.data['year']]
    bool_depth_regq = [True if (np.float(depth)>=np.float(min_depth.value))&(np.float(depth)<=np.float(max_depth.value)) else False for depth in source_regq.data['depth']]
    bool_mag_regq = [True if (np.float(mag)>=np.float(min_mag.value))&(np.float(mag)<=np.float(max_mag.value)) else False for mag in source_regq.data['Mw']]
    view_regq.filters=[BooleanFilter(bool_year_regq),BooleanFilter(bool_depth_regq),BooleanFilter(bool_mag_regq)]

    bool_year_mt = [True if (np.float(year)>=np.float(min_year.value))&(np.float(year)<=np.float(max_year.value)) else False for year in source_mt.data['year']]
    bool_depth_mt = [True if (np.float(depth)>=np.float(min_depth.value))&(np.float(depth)<=np.float(max_depth.value)) else False for depth in source_mt.data['depth']]
    bool_mag_mt = [True if (np.float(mag)>=np.float(min_mag.value))&(np.float(mag)<=np.float(max_mag.value)) else False for mag in source_mt.data['Mw']]
    view_mt.filters=[BooleanFilter(bool_year_mt),BooleanFilter(bool_depth_mt),BooleanFilter(bool_mag_mt)]

for w in [min_year,max_year,min_depth,max_depth,min_mag,max_mag]:
    w.on_change('value', update_regional_seismicity)

p_reg = figure(x_range=(xl, xr), y_range=(yb, yt),
                x_axis_type='mercator', y_axis_type='mercator',
                tools='pan,wheel_zoom,box_zoom,reset',
           active_drag='pan',active_scroll='wheel_zoom')
p_reg.add_tile(get_provider('ESRI_IMAGERY'))

p_reg.image_url(url='url', x='x', y='y', w='radius', h='radius', anchor='center', source=source_reg, view=view_reg)
p_reg.image_url(url='url', x='x', y='y', w='radius', h='radius', anchor='center', source=source_regq, view=view_regq)
p_reg.image_url(url='url', x='x', y='y', w='radius', h='radius', anchor='center', source=source_mt, view=view_mt)

p_regP = figure(x_range=p_reg.x_range, y_range=p_reg.y_range,
                x_axis_type='mercator', y_axis_type='mercator',
                title="P-axes")
p_regP.add_tile(get_provider('ESRI_IMAGERY'))

p_regP.multi_line(xs='P_x',ys='P_y',source=source_reg,color='yellow',width=2,view=view_reg)
p_regP.multi_line(xs='P_x',ys='P_y',source=source_regq,color='yellow',width=2,view=view_regq)
p_regP.multi_line(xs='P_x',ys='P_y',source=source_mt,color='red',width=2,view=view_mt)
p_regP.circle(x='x',y='y',source=source_reg,radius=tensor_size/5,line_color='yellow',fill_color='black',line_width=2,view=view_reg)
p_regP.circle(x='x',y='y',source=source_regq,radius=tensor_size/5,line_color='yellow',fill_color='black',line_width=2,view=view_regq)
p_regP.circle(x='x',y='y',source=source_mt,radius=tensor_size/5,line_color='red',fill_color='black',line_width=2,view=view_mt)

p_regT = figure(x_range=p_reg.x_range, y_range=p_reg.y_range,
                x_axis_type='mercator', y_axis_type='mercator',
                title="T-axes")
p_regT.add_tile(get_provider('ESRI_IMAGERY'))
p_regT.multi_line(xs='T_x',ys='T_y',source=source_reg,color='yellow',width=2,view=view_reg)
p_regT.multi_line(xs='T_x',ys='T_y',source=source_regq,color='yellow',width=2,view=view_regq)
p_regT.multi_line(xs='T_x',ys='T_y',source=source_mt,color='red',width=2,view=view_mt)
p_regT.circle(x='x',y='y',source=source_reg,radius=tensor_size/5,line_color='yellow',fill_color='black',line_width=2,view=view_reg)
p_regT.circle(x='x',y='y',source=source_regq,radius=tensor_size/5,line_color='yellow',fill_color='black',line_width=2,view=view_regq)
p_regT.circle(x='x',y='y',source=source_mt,radius=tensor_size/5,line_color='red',fill_color='black',line_width=2,view=view_mt)

layout4 = column(Div(text='<h1> Regional seismicity <h1>'),
                               row(min_year, max_year, min_depth, max_depth, min_mag, max_mag),
                               row(p_reg, p_regP, p_regT))

panel4 = Panel(child=layout4,title='Regional seismicity')

###############################################################################
# Compile tabs
###############################################################################
tabs = Tabs(tabs=[panel1, panel2, panel3, panel4])

bokeh_doc = curdoc()
bokeh_doc.add_root(tabs)
bokeh_doc.title = 'Record Reading Server App'
