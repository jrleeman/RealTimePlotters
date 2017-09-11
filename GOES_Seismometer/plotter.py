from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from obspy.clients.fdsn import Client
from obspy import UTCDateTime
import numpy as np
from matplotlib.dates import DateFormatter, DayLocator, HourLocator
from siphon.catalog import TDSCatalog
from netCDF4 import Dataset
from matplotlib import patheffects
import cartopy.feature as cfeat
import cartopy.crs as ccrs
from metpy.plots import add_logo


def get_station_location(t0, net, st0, loc):
    """
    Get the station latitude, longitude, and elevation.

    Given a time, duration, loc code, and station network/name, get
    station information from IRIS.  Return a list containing the
    lat, lon, and elevation.
    """
    client = Client('IRIS')
    st0 = client.get_stations(starttime=t0, endtime=t0+timedelta(seconds=60),
                              network=net, station=st0, level='station')
    slat = st0[0][0].latitude
    slon = st0[0][0].longitude
    selev = st0[0][0].elevation
    return [slat, slon, selev]


def get_seismometer_data(t0, net, st0, loc, ch, duration):
    """
    Download data from the IRIS datacenter and output
    with the instrument response removed and calibrated.
    Return a station object.
    """
    client = Client("IRIS")
    st = client.get_waveforms(net, st0, loc, ch, t0,
                              t0+timedelta(minutes=duration), attach_response=True)
    st.detrend(type='linear')
    st.detrend(type='constant')
    st.taper(max_percentage=0.01)
    st.remove_response(output='DISP')
    return st


def open_GOES_dataset(date, channel, idx):
    """
    Open and return a netCDF Dataset object for a given date, channel, and image index
    of GOES-16 data from THREDDS test server.
    """
    cat = TDSCatalog('http://thredds-test.unidata.ucar.edu/thredds/catalog/satellite/goes16/GOES16/'
             'Mesoscale-1/Channel{:02d}/{:%Y%m%d}/catalog.xml'.format(channel, date)) #Mesoscale-1
    dataset = cat.datasets[idx]
    ds = Dataset(dataset.access_urls['OPENDAP'])
    return ds

#
# Only change these parameters
#

history = 48 # How much seismic history to plot in hours
n_goes_images = 24 * 60 * 2 # Number of GOES-16 images to plot
goes_skip = 10 # Only plot every n images from goes
start_number = 288

# Station to be used in plot
seismic_station = {'station': 'AGPR', 'network': 'PR', 'location': '--', 'channel': 'HNZ'}
seismic_station = {'station': 'SJG', 'network': 'IU', 'location': '00', 'channel': 'BHZ'}

#
# Change nothing below here!
#

# Based on the current time, calculate the start and end time to the nearest
# hour to the present.
now = datetime.utcnow()
end_time = datetime(now.year, now.month, now.day, now.hour)
start_time = end_time - timedelta(hours=history)

# Get the seismic data from IRIS
station = seismic_station['station']
network = seismic_station['network']
location = seismic_station['location']
channel = seismic_station['channel']

st = get_seismometer_data(start_time, network, station, location, channel, history * 60)
st.filter('highpass', freq=0.02)
tr = st[0]
tr = tr.decimate(10)
trace_data = tr.data

trace_time = []
for i in range(len(trace_data)):
    trace_time.append(start_time + timedelta(seconds=i * 1./st[0].stats.sampling_rate))

# Set up a feature for the state/province lines. Tell cartopy not to fill in the polygons
state_boundaries = cfeat.NaturalEarthFeature(category='cultural',
                                             name='admin_1_states_provinces_lakes',
                                             scale='50m', facecolor='none')

station_latitude, station_longitude, _ = get_station_location(start_time, network, station, location)

for i, goes_idx in enumerate(range(-1*n_goes_images, 0)[::goes_skip]):
    print('Making image: {}'.format(i))
    goes_ds = open_GOES_dataset(datetime(2017,9,7), 14, goes_idx)
    proj_var = goes_ds.variables[goes_ds.variables['Sectorized_CMI'].grid_mapping]
    #print(proj_var)
    x = goes_ds.variables['x'][:]
    y = goes_ds.variables['y'][:]
    goes_image = goes_ds.variables['Sectorized_CMI'][:]
    time = datetime.strptime(goes_ds.start_date_time, '%Y%j%H%M%S')
    # Create a Globe specifying a spherical earth with the correct radius
    globe = ccrs.Globe(ellipse='sphere', semimajor_axis=proj_var.semi_major,
                       semiminor_axis=proj_var.semi_minor)


    if proj_var.grid_mapping_name == 'lambert_conformal_conic':
        # If the projection is LCC
        proj = ccrs.LambertConformal(central_longitude=proj_var.longitude_of_central_meridian,
                                     central_latitude=proj_var.latitude_of_projection_origin,
                                     standard_parallels=[proj_var.standard_parallel],
                                     globe=globe)
    else:
        # If the projection if Mercator
        proj = ccrs.Mercator(central_longitude=proj_var.longitude_of_projection_origin,
                             globe=globe, latitude_true_scale=proj_var.standard_parallel)

    # Make the plot
    fig = plt.figure(figsize=(8, 8.5))
    ax1 = plt.subplot2grid((4, 1), (0, 0), projection=proj, rowspan=3)
    ax2 = plt.subplot2grid((4, 1), (3, 0))
    #ax1 = fig.add_subplot(4, 1, 1, projection=proj)
    #ax2 = fig.add_subplot(4, 1, 4)

    plt.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.99, wspace=0, hspace=0.05)


    # Plot the map fiducials and GOES image
    ax1.coastlines(resolution='50m', color='black')
    ax1.add_feature(state_boundaries, linestyle=':', edgecolor='black')
    ax1.add_feature(cfeat.BORDERS, linewidth='2', edgecolor='black')
    im = ax1.imshow(goes_image, extent=(x.min(), x.max(), y.min(), y.max()), origin='upper')
    im.set_cmap('Greys')
    im.set_norm(plt.Normalize(200, 330))

    # Plot the seismometer on the map
    ax1.scatter(station_longitude, station_latitude, s=80, color='r', marker='*', transform=ccrs.PlateCarree())

    # Plot text on the satellite image
    text_time = ax1.text(0.99, 0.01, time.strftime('%d %B %Y %H%MZ'),
                   horizontalalignment='right', transform=ax1.transAxes,
                   color='white', fontsize='x-large', weight='bold')

    text_channel = ax1.text(0.5, 0.97, 'Experimental GOES-16 Channel 14',
                   horizontalalignment='center', transform=ax1.transAxes,
                   color='white', fontsize='large', weight='bold')

    # Make the text stand out even better using matplotlib's path effects
    outline_effect = [patheffects.withStroke(linewidth=2, foreground='black')]
    text_time.set_path_effects(outline_effect)
    text_channel.set_path_effects(outline_effect)

    # Plot the seismic data
    ax2.plot(trace_time, trace_data, color='k')
    ax2.axvline(x=time, color='red', linewidth=2)

    ax2.xaxis.set_major_formatter(DateFormatter('%m/%d'))
    ax2.xaxis.set_minor_formatter(DateFormatter('%HZ'))
    ax2.xaxis.set_major_locator(DayLocator())
    ax2.xaxis.set_minor_locator(HourLocator([6, 12, 18]))

    add_logo(fig, size='small', x=45, y=270)
    ax2.get_yaxis().set_ticks([])

    plt.savefig('plots/{:04d}.png'.format(i+start_number))
    plt.close()
