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

warnings.simplefilter(action="ignore", category=SettingWithCopyWarning) # suppresses a warning in binning and plotting

from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource, CustomJS, Div, OpenURL, Slider, ColorBar
from bokeh.models import TextInput, HoverTool, TapTool, Select, RangeSlider, MultiChoice, RadioButtonGroup, CheckboxGroup, FreehandDrawTool
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
from obspy.clients.fdsn.mass_downloader import CircularDomain, Restrictions, MassDownloader
from obspy.geodetics.base import gps2dist_azimuth

from pyproj import Transformer

###############################################################################
# Tab 1: Select event
###############################################################################
# load initial inputs from file
meta = pd.read_hdf('eventdat',key='meta')

# Selecting an event
lat = meta.lat[0];
lon = meta.lon[0];
ID = meta.ID[0];
t11 = UTCDateTime(meta.time[0])
search_time_range = 72 # hours
t22 = t11 + 3600*search_time_range

# Setup search parameter inputs
input_lat = TextInput(value=str(lat), title='Latitude (-90 to 90):')
input_lon = TextInput(value=str(lon), title='Longitude (-180 to 180):')
input_time = TextInput(value=meta.time[0],title='Start time UTC (YYYY-MM-DD HH:MM:SS)')
search_time_range = Slider(start=0,end=72,value=72, title='Search time range (hrs)')
input_mag = RangeSlider(start=0,end=9,value=(meta.mag[0]-0.5,meta.mag[0]+0.5),title='Magnitude') # Can probably add an if condition to change the slider bar boundaries for other magnitude scales
mag_type = Select(title='Magnitude type', value='Mw',
               options=['Mw'])
search_rad = Slider(start=0,end=20,value=10,title='Search radius (deg)')
input_webservice = Select(title='Catalog', value='IRIS',
               options=['IRIS'])
input_ID = TextInput(value=str(ID[0]))

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
                               minmagnitude=meta.mag[0]-0.5,
                               maxmagnitude=meta.mag[0]+0.5,
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
                minmagnitude=float(input_mag.value[0]),maxmagnitude=float(input_mag.value[1]),magnitudetype=mag_type.value,catalog='GCMT',
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
                x_axis_type='mercator', y_axis_type='mercator',tools='tap')
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
        evtime = eventlist[iev].origins[0].time.strftime('%y/%m/%d %H:%M')
        
        source.data = {'x': [float(xcoord.value)],'y': [float(ycoord.value)],'lat':[lat],'lon':[lon],'depth':[eventlist[0].origins[0].depth/1000],
                       'mag':[eventlist[0].magnitudes[0].mag],'size':[4*eventlist[0].magnitudes[0].mag],'time':[evtime], 'id':[input_ID.value]}
    except Exception as e: 
        print(e)     

p.js_on_event(events.SelectionGeometry, callbacktap(xcoord,ycoord, attributes=['geometry']))
ycoord.on_change('value', callbacktap2)
    
# Save event info
def save_button_callback():
    eventdat = pd.DataFrame(data={'ID':str(input_ID.value),'lat':float(input_lat.value), 'lon':float(input_lon.value), 'time':str(input_time.value), 'mag':float(input_mag.value[0])}, index=[0])
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
phase_shift = Select(title="Align:", value="None", options=['None', 'P-wave'])

load_stations = Button(label='Load station map', button_type='success')
download_data = Button(label='Download data', button_type='success')

try: 
    df
except: 
    df = pd.read_hdf('eventdat',key='data')

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
    st.rotate('->ZNE',inventory=inventory) # rotates channels **1 and **2 to **N and **E
    # divide between raw, high, mid, and low frequency bands
    st_raw = st.copy()
    st_high = st.copy().filter('bandpass',freqmin=freqmin_high, freqmax=freqmax_high)
    st_mid  = st.copy().filter('bandpass',freqmin=freqmin_mid,  freqmax=freqmax_mid)
    st_low  = st.copy().filter('bandpass',freqmin=freqmin_low,  freqmax=freqmax_low)
    # Initialize data storage structures
    meta = np.array([['Network','Station','Lat','Lon','Azimuth','Distance',
                      'Frequency','Channel']]) # channel metadata
    time = np.array([[pd.to_datetime(str(st_raw[0].stats.starttime)) + pd.Timedelta(st_raw[0].stats.delta*i,unit='s') for i in np.arange(st_raw[0].stats.npts)]]) # time data
    data = np.array([np.arange(st_raw[0].stats.npts)*st_raw[0].stats.delta]) # time data
    for i,stat in enumerate(bulk_stat): # loop over each station
        # get metadata
        st_raw_i = st_raw.select(station=stat)
        try:
            net = st_raw_i[0].stats.network
            station_info = inventory.select(network=net,station=stat).networks[0].stations[0]
            stlat = station_info.latitude 
            stlon = station_info.longitude
        
            # calculate distance, azimuth, abd back-azimuth
            [dist, az, back_az] = gps2dist_azimuth(float(input_lat.value),float(input_lon.value),stlat,stlon)
    
            # rotate channels to ZRT
            st_raw_i.rotate('NE->RT',back_azimuth=back_az)
            st_high_i = st_high.select(station=stat).rotate('NE->RT',back_azimuth=back_az)
            st_mid_i = st_mid.select(station=stat).rotate('NE->RT',back_azimuth=back_az)
            st_low_i = st_low.select(station=stat).rotate('NE->RT',back_azimuth=back_az)
    
            # loop over channels
            for channel in ['BHZ','BHR','BHT']:
                for freq,freq_name in zip([st_raw_i,st_high_i,st_mid_i,st_low_i],['Raw','High','Mid','Low']):
                    meta = np.append(meta,[[net,stat,stlat,stlon,az,dist/111139,freq_name,channel]],axis=0)   
                    time = np.append(time,
                    [np.array([pd.to_datetime(str(st_raw_i[0].stats.starttime)) + pd.Timedelta(st_raw_i[0].stats.delta*i,unit='s') for i in np.arange(st_raw_i[0].stats.npts)]
                     ,dtype='datetime64')],
                                 axis=0)
                    data = np.append(data,[freq.select(channel=channel.replace('B','*'))[0].data],axis=0)
        except:
            bulk_stat.remove(stat)
    
    # store in pandas dataframe
    df = pd.DataFrame(meta[1:,:], columns = np.array(meta[0,:],
                      dtype='str')).astype({'Lat':'float32','Lon':'float32',
                      'Azimuth':'float32',
                      'Distance':'float32'}).join(pd.Series(list(time[1:,:]),
                      name="Time")).join(pd.Series(list(data[1:,:]),name="Data"))
    df = df.dropna()
    df.to_hdf('eventdat',key='data') # save to 'eventdat.h5'
    eventdat = pd.DataFrame(data={'ID':str(input_ID.value),
                                  'lat':float(input_lat.value), 
                                  'lon':float(input_lon.value), 
                                  'time':str(input_time.value), 
                                  'mag':float(input_mag.value[0])}, index=[0])
    eventdat.to_hdf('eventdat',key='meta')
    print('Loaded ' + str((df.shape[0])/12) + ' stations')
    
    stat_lat = df['Lat'].values
    stat_lon = df['Lon'].values
    [stat_x, stat_y] = latlon2webmercator.transform(stat_lat, stat_lon)
    stations_source.data = {'stat_x': stat_x,'stat_y': stat_y,
                                             'stat_lat':stat_lat,'stat_lon':stat_lon,
                                             'stat_name':df['Station'].values}
    
    df['SNR'] = df['Data'].apply(lambda x: np.mean(x**2)/np.mean(x[:100]**2))
    
    [binned,mapper] = update_plotting_data(freq_select,filled_select,amplitude_slider,
                         normalize_select,sort_opts,sort_select,color_by_az,
                         azimuth_range,binning_select,phase_cheatsheet,color_bar,
                         df)
    
    source_records.data = binned
                    
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
channel_select = Select(value="Vertical (BHZ)", options=['Vertical (BHZ)', 'Radial (BHR)','Transverse (BHT)'])

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

freq_select = Select(value=freq_options[0], options=freq_options)
filled_select = CheckboxGroup(labels=['Fill above zero'], active=[])
amplitude_slider = Slider(start=0.1, end=10, value=1, step=.1, title='Relative amplitude')
normalize_select = RadioButtonGroup(labels=['Individual', 'Global'], active=0)
sort_opts = ['Distance', 'Azimuth']
sort_select = RadioButtonGroup(labels=sort_opts, active=0)
color_by_az = CheckboxGroup(labels=['Color by azimuth'], active=[])
azimuth_range = RangeSlider(start=0,end=360,value=(0,360),title='Azimuth range (deg)')
binning_select = CheckboxGroup(labels=['Binned'], active=[0]) # binned data plots one trace for each bin that has the lowest SNR

phase_cheatsheet = CheckboxGroup(labels=['Phase cheatsheet'], active=[])
    
df['SNR'] = df['Data'].apply(lambda x: np.mean(x**2)/np.mean(x[:100]**2))

station_data = ColumnDataSource(data={'x':[],'y':[]})
p3b = figure(plot_height=200,plot_width=200,x_axis_type='mercator',y_axis_type= 'mercator',
                 x_range=(-padding, padding),y_range=(-padding, padding),tools='')
p3b.add_tile(get_provider('ESRI_IMAGERY'))
p3b.triangle('x','y',source=station_data,size=20,color='red',line_color='black')
    
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
    
mapper = linear_cmap(field_name='Azimuth', palette=grey(1) ,low=df['Azimuth'].min() ,high=df['Azimuth'].max())
color_bar = ColorBar(color_mapper=mapper['transform'],visible=False)

def update_plotting_data(freq_select,filled_select,amplitude_slider,
                         normalize_select,sort_opts,sort_select,color_by_az,
                         azimuth_range,binning_select,phase_cheatsheet,color_bar,df):
    
    group = df.groupby(['Frequency','Channel']).get_group((freq_select.value.split(' ')[0],channel_select.value[-4:-1]))
    group = group[group['Azimuth'].between(azimuth_range.value[0], azimuth_range.value[1], inclusive=True)]
    
    sort_by = sort_opts[sort_select.active]
    
    bins = np.linspace(group[sort_by].min()-0.1,group[sort_by].max()+0.1,179)
    if np.isnan(bins).any():
        bins = np.linspace(0,180*(1 + sort_select.active),179)
    norm_stretch = 20*(bins[1]-bins[0])*amplitude_slider.value

    if normalize_select.active == 0:
        group['plot_trace'] = group[sort_by] + group['Data'].apply(lambda x: -x*norm_stretch/np.max(np.abs(x))) # normalize by max value of each trace
    else: 
        global_max = group['Data'].agg(lambda x: np.max(np.abs(x))).max()
        group['plot_trace'] = group[sort_by] + group['Data'].apply(lambda x: -x*norm_stretch/global_max)
    
    if 0 in binning_select.active:
        group['binned'] = pd.cut(group[sort_by], bins) # bin
        binned = group.loc[group.groupby('binned')['SNR'].agg(
                    lambda x : np.nan if x.count() == 0 else x.idxmax()
                ).dropna().sort_values().values].drop(columns=['binned'])
    else:
        binned = group
        
    if 0 in filled_select.active:
        binned['Fills'] = (binned['plot_trace']-binned[sort_by]).apply(lambda x: np.max([x,np.zeros_like(x)],axis=0)) + binned[sort_by]
    else:
        binned['Fills'] = 0*binned['plot_trace'] + binned[sort_by]
        
    if 0 in color_by_az.active:
        mapper['transform'].palette=Viridis256
        mapper['transform'].low=binned['Azimuth'].min()
        mapper['transform'].high=binned['Azimuth'].max()
        color_bar.visible = True
    else: 
        mapper['transform'].palette=grey(1)
        mapper['transform'].low=binned['Azimuth'].min()
        mapper['transform'].high=binned['Azimuth'].max()
        color_bar.visible = False
    return binned, mapper

[binned,mapper] = update_plotting_data(freq_select,filled_select,amplitude_slider,
                         normalize_select,sort_opts,sort_select,color_by_az,
                         azimuth_range,binning_select,phase_cheatsheet,color_bar,df)
source_records = ColumnDataSource(binned)
    
sort_by = sort_opts[sort_select.active]
p3 = figure(plot_height=600,plot_width=1200,x_axis_type='datetime',y_axis_label= sort_by + ' (deg)',
           x_axis_label='Time (min)',tools='box_zoom,undo,redo,reset,save', active_drag='box_zoom')
p3.y_range.flipped = True
p3.patches(xs='Time', ys='Fills', source=source_records, fill_color='blue',line_color='white',line_alpha=1)
tr = p3.multi_line(xs='Time', ys='plot_trace', line_width=1, line_color=mapper,
             source=source_records)
p3.add_layout(color_bar, 'right')

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

p3.add_tools(HoverTool(tooltips=TOOLTIPS3, mode='mouse'))
p3.toolbar.active_inspect = None

draw_r = p.multi_line('Time', 'plot_trace', source=source_records)
freehand_draw = FreehandDrawTool(renderers=[draw_r])

def update_display(attrname, old, new):
    
    [binned,mapper] = update_plotting_data(freq_select,filled_select,amplitude_slider,
                         normalize_select,sort_opts,sort_select,color_by_az,
                         azimuth_range,binning_select,phase_cheatsheet,color_bar,df)
    source_records.data = binned

for u in [channel_select,freq_select,amplitude_slider,azimuth_range]:
    u.on_change('value', update_display)
for v in [filled_select,normalize_select,sort_select,color_by_az,binning_select]:
    v.on_change('active',update_display)
    
change_y_axis_label = CustomJS(args=dict(plot=p3, source=source_records, sort_by=sort_select, sort_opts = sort_opts, ax=p3.yaxis), code="""
    ax[0].axis_label = sort_opts[sort_by.active] + ' (deg)';
    source.change.emit();
""")


sort_select.js_on_change('active', change_y_axis_label)
    
layout3 = column(row(channel_select,freq_select),
                 p3,
                 Div(text='Figure 3: Seismic traces'))
panel3 = Panel(child=layout3,title='Display records')

layout4 = column(Div(text='Velocity reduction'),
                 Div(text='Virtualization'),
                 filled_select,
                 amplitude_slider,
                 Div(text='Normalize:'),
                 normalize_select,
                 Div(text='Sort by:'),
                 sort_select,
                 color_by_az,
                 azimuth_range,
                 binning_select,
                 Div(text='Phase cheatsheet'),
                 phase_cheatsheet,
                 Div(text='Sonification'))
panel4 = Panel(child=layout4, title='Display controls')

###############################################################################
# Compile tabs
###############################################################################
tabs = Tabs(tabs=[panel1, panel2, panel3, panel4])

bokeh_doc = curdoc()
bokeh_doc.add_root(tabs)
bokeh_doc.title = 'Record Reading Server App'