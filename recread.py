''' Record Reading Plots
Program to extract and plot global waveforms for a specified seismic event and maps regional focal mechanisms. This program creates an interactive plot to view global waveforms for a user-specified seismic event. The program also has options to plot regional maps of focal mechanisms from the CMT catalog (globalcmt.org).
This is a python adaptation to the original record reading waveform plot program written by Ge Jin (jinwar@gmail.com) in Matlab.
Contributors: Janine Birnbaum, Theresa Sawi, Christopher Carchedi, and Michelle Lee

Instructions:
1. Event Selection: Select event to plot global waveforms for by inputting location and magnitude parameters of desired event
2. (User steps to download the data?)
3. Plotting the waveforms: The global waveforms will be be displayed on an interactive plot allowing users to zoom in/out, change component, change frequency band, ... (and other features)
4. Plotting focal mechanisms: (info on what user needs to do for this)

Use the ``bokeh serve`` command to run the example by executing:
    bokeh serve recread.py
at your command prompt. Then navigate to the URL
    http://localhost:5006/recread
in your browser.
'''
import sys

import numpy as np
import pandas as pd

from bokeh.layouts import column, row, widgetbox
from bokeh.models import ColumnDataSource, CustomJS, Div, Slider, TextInput, Circle, HoverTool, TapTool, Select, RangeSlider
from bokeh.models.widgets import Button, Panel, Tabs
from bokeh.plotting import figure, curdoc
from bokeh.tile_providers import get_provider, Vendors

import obspy
from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from obspy.core import AttribDict
from obspy.io.sac import SACTrace

from pyproj import Transformer

# load initial inputs from file

meta = pd.read_hdf('eventdat',key='meta')

# Selecting an event
lat = meta.lat;
lon = meta.lon;

ID = meta.ID;

#tstart = "2020-07-22T00:00:00"
t11 = UTCDateTime(meta.time)
search_time_range = 72 # hours
t22 = t11 + 3600*search_time_range

def change_input_lat(attr,old,new):
    lat = attr

input_lat = TextInput(value=str(lat), title='Latitude')
input_lon = TextInput(value=str(lon), title='Longitude')
input_time = TextInput(value=meta.time, title='Time')
search_time_range = Slider(start=0,end=72,value=72, title='Search time (hrs)')
input_mag = RangeSlider(start=0,end=10,value=(meta.mag-0.5,meta.mag+0.5),title='Magnitude')
mag_type = Select(title="Magnitude type", value="Mw",
               options=['Mw'])
search_rad = Slider(start=0,end=20,value=10,title='Search radius (deg)')
input_webservice = Select(title="Catalog", value='IRIS',
               options=['IRIS'])
input_ID = TextInput(value=str(ID))

# Fetch initital data
client = Client(input_webservice.value)
eventlist = (client.get_events(starttime=t11,endtime=t22,latitude=lat,
                               longitude=lon,minradius=0,maxradius=10,
                               minmagnitude=meta.mag-0.5,
                               maxmagnitude=meta.mag+0.5,
                               magnitudetype=mag_type.value,
                               catalog='GCMT',orderby='time-asc'))
    
    # HELPFUL THINGS TO PRINT?
    # print(eventlist.__str__(print_all=True))
    # eventlist.plot(projection="global")
    # print(eventlist)
    # eventlist[0].origins has time,long,lat,depth,creation_info
    # eventlist[0].event_descriptions has location description text
    # eventlist[0].magnitudes lists magnitude and magnitude type
    
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
    evtime = np.append(evtime,eventlist[iev].origins[0].time.strftime("%y/%m/%d %H:%M"))
    ID = np.append(ID,str(eventlist[iev].resource_id))
    
# define coordinate transformations from lat/lon to Web Mercator
latlon2webmercator = Transformer.from_crs('EPSG:4326', 'EPSG:3857')
webmercator2latlon = Transformer.from_crs('EPSG:3857', 'EPSG:4326')

# plot initial data
[x1, y1] = latlon2webmercator.transform(evlat,evlon) # transform event coordinates  
source = source = ColumnDataSource(data={'x': x1,'y': y1,'lat':evlat,'lon':evlon,'depth':evdepth/1000,
				'mag':evmag,'size':4*evmag,'time':evtime, 'id':ID})


def update_search_params(attrname, old, new):

    eventlist = []
    
    [x_cntr, y_cntr] = latlon2webmercator.transform(float(input_lat.value),float(input_lon.value))
    
    p.x_range.start = x_cntr-padding
    p.x_range.end   = x_cntr+padding
    p.y_range.start = y_cntr-padding
    p.y_range.end   = y_cntr+padding
    # Get the current slider values
    t11 = UTCDateTime(input_time.value)
    t22 = t11 + 3600*search_time_range.value

    # Fetch data
    client = Client(input_webservice.value)
    
    try:
        eventlist = (client.get_events(starttime=t11,endtime=t22,latitude=float(input_lat.value),longitude=float(input_lon.value),minradius=0,maxradius=float(search_rad.value),
                minmagnitude=float(input_mag.value[0]),maxmagnitude=float(input_mag.value[1]),magnitudetype=mag_type.value,catalog='GCMT',
                orderby='time-asc'))
     
    # HELPFUL THINGS TO PRINT?
    # print(eventlist.__str__(print_all=True))
    # eventlist.plot(projection="global")
    # print(eventlist)
    # eventlist[0].origins has time,long,lat,depth,creation_info
    # eventlist[0].event_descriptions has location description text
    # eventlist[0].magnitudes lists magnitude and magnitude type
    
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
            evtime = np.append(evtime,eventlist[iev].origins[0].time.strftime("%y/%m/%d %H:%M"))
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
                x_axis_type="mercator", y_axis_type="mercator")
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

# Hover tool
cr = p.circle('x', 'y', size='size',source=source,color='red',line_color='black',hover_color='yellow',hover_line_color='black')
callbackhover = CustomJS(args={'circle': cr.data_source})

# Add tools
p.add_tools(HoverTool(tooltips=TOOLTIPS, callback=callbackhover, renderers=[cr],mode='mouse'))

# Save event info
def save_button_callback():
    eventdat = pd.Series(data={'ID':input_ID.value,'lat':input_lat.value, 'lon':input_lon.value, 'time':input_time.value, 'mag':input_mag.value[0]})
    eventdat.to_hdf('eventdat',key='meta')

save_button = Button(label='Save', button_type='success')
save_button.on_click(save_button_callback)

# Button to stop the server
def quit_button_callback():
    sys.exit()  # Stop the server
    
quit_button = Button(label="Quit", button_type="success")
quit_button.on_click(quit_button_callback)

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
                              quit_button]),p)))
panel1 = Panel(child=layout1, title='Select event')
tabs = Tabs(tabs=[panel1])

bokeh_doc = curdoc()
bokeh_doc.add_root(tabs)
bokeh_doc.title = "Record Reading Server App"