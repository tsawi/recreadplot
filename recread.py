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

import numpy as np
import pandas as pd

from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource, CustomJS, Div, OpenURL, Slider
from bokeh.models import TextInput, HoverTool, TapTool, Select, RangeSlider, MultiChoice, RadioButtonGroup
from bokeh.models.widgets import Button, Panel, Tabs
from bokeh.plotting import figure, curdoc
from bokeh.tile_providers import get_provider
from bokeh import events

from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from obspy.core import AttribDict
from obspy.io.sac import SACTrace

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
t11 = UTCDateTime(str(meta.time[0]).replace('-','').replace(' ','').replace(':',''))
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
    ('Depth','@depth{0.0}'),
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
network_options = [('1', '_US-ALL'), ('2', '_GSN'), ('3', 'TA')]
network = MultiChoice(title='Station network', value=['1', '2'],
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
resample_delta = TextInput(value='0.1',title='Resample rate:')

load_stations = Button(label='Load station map', button_type='success')
download_data = Button(label='Download data', button_type='success')

stations_source = ColumnDataSource(data={'stat_x': [],'stat_y': [],
                                             'stat_lat':[],'stat_lon':[],
                                             'stat_name':[]})


# download some station info based on networks
# Stations map (could maybe pull station info from the networks and allow a select tool on a map)
def add_network_option():
    network_options.append((str(len(network_options)+1),str(add_network_text.value)))
    network.options = network_options
    network.value.append(str(len(network_options)))

add_network_button.on_click(add_network_option)

def networkscallback(attrs,new,old):
    networks = network.value
network.on_change('value', networkscallback)

def select_by_dist(attrs,new,old):
    xpadding = latlon2webmercator.transform(0,180)[0]/180*float(epicentral_distance.value[1])
    ypadding = latlon2webmercator.transform(85,0)[1]/90*float(epicentral_distance.value[1])
    p2.x_range.start=x_cntr-xpadding
    p2.x_range.end  =x_cntr+xpadding
    p2.y_range.start=np.max([y_cntr-ypadding,latlon2webmercator.transform(-85,0)[1]])
    p2.y_range.end  =np.min([y_cntr+ypadding,latlon2webmercator.transform(85,0)[1]])

x_cntr = float(xcoord.value)
y_cntr = float(ycoord.value)
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

def load_stations_callback():
    stat_lat = np.array([])
    stat_lon = np.array([])
    stat_name = np.array([])
    t11 = UTCDateTime(str(input_time.value).replace('-','').replace(' ','').replace(':',''))
    for opt in network.value:
        try:
            inventory = client.get_stations(network=str(network_options[int(opt)-1][1]), 
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

layout2 = row(column(Div(text='<h1>Downloading the data parameters<h1>'),
                     network,
                     row(add_network_text,add_network_button),
                     epicentral_distance,
                     Div(text='Download waveform length:'),
                     row(min_before, min_after),
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
# Compile tabs
###############################################################################
tabs = Tabs(tabs=[panel1, panel2])

bokeh_doc = curdoc()
bokeh_doc.add_root(tabs)
bokeh_doc.title = 'Record Reading Server App'