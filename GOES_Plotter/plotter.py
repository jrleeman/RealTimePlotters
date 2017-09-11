from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.dates import DateFormatter, DayLocator, HourLocator
from siphon.catalog import TDSCatalog
from netCDF4 import Dataset
from matplotlib import patheffects
import cartopy.feature as cfeat
import cartopy.crs as ccrs
from metpy.plots import add_logo


def open_GOES_dataset(date, region, abschannel, idx):
    """
    Open and return a netCDF Dataset object for a given date, channel, and image index
    of GOES-16 data from THREDDS test server.
    """
    cat = TDSCatalog('http://thredds-test.unidata.ucar.edu/thredds/catalog/satellite/goes16/GOES16/'
             '{}/Channel{:02d}/{:%Y%m%d}/catalog.xml'.format(region, channel, date)) #Mesoscale-1
    dataset = cat.datasets[idx]
    ds = Dataset(dataset.access_urls['OPENDAP'])
    return ds

def get_number_datasets(date, region, channel):
    """
    Return the number of datasets available for a given channel, day, region.
    """
    cat = TDSCatalog('http://thredds-test.unidata.ucar.edu/thredds/catalog/satellite/goes16/GOES16/'
             '{}/Channel{:02d}/{:%Y%m%d}/catalog.xml'.format(region, channel, date))
    return len(cat.datasets)


#
# Only change these parameters
#
start_number = 0  # Frame numbering starts from here
start_day = datetime(2017, 8, 27)  # Inclusive
end_day = datetime(2017, 9, 11)  # Inclusive
channel = 13  # ABI Channel Number
region = 'Mesoscale-1'  # Region (CONUS, Mesoscale-1, Mesoscale-2, FullDisk)
skipframes = 10  # Number of frames to skip, 1 for plotting all frames

#
# Change nothing below here!
#

# Make range of days to get data from
num_days = (end_day - start_day).days + 1
dates = [start_day + timedelta(days=x) for x in range(0, num_days)]

# Set up a feature for the state/province lines. Tell cartopy not to fill in the polygons
state_boundaries = cfeat.NaturalEarthFeature(category='cultural',
                                             name='admin_1_states_provinces_lakes',
                                             scale='50m', facecolor='none')

# Loop over each day and plot all of the data available for that day
frame_number = 0
for day in dates:
        num_datasets = get_number_datasets(day, region, channel)
        for idx in range(num_datasets)[::skipframes]:
            print('{:%Y-%m-%d} - {}/{}'.format(day, idx, num_datasets))
            goes_ds = open_GOES_dataset(day, region, channel, idx)
            proj_var = goes_ds.variables[goes_ds.variables['Sectorized_CMI'].grid_mapping]
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
            fig = plt.figure(figsize=(18, 15.7))
            ax1 = fig.add_subplot(1, 1, 1, projection=proj)

            plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)


            # Plot the map fiducials and GOES image
            ax1.coastlines(resolution='50m', color='black')
            ax1.add_feature(state_boundaries, linestyle=':', edgecolor='black')
            ax1.add_feature(cfeat.BORDERS, linewidth='2', edgecolor='black')
            im = ax1.imshow(goes_image, extent=(x.min(), x.max(), y.min(), y.max()), origin='upper')
            im.set_cmap('Greys')
            im.set_norm(plt.Normalize(200, 330))

            # Plot text on the satellite image
            text_time = ax1.text(0.99, 0.01, time.strftime('%d %B %Y %H%MZ'),
                           horizontalalignment='right', transform=ax1.transAxes,
                           color='white', fontsize='x-large', weight='bold')

            text_channel = ax1.text(0.5, 0.97, 'Experimental GOES-16 Channel 13',
                           horizontalalignment='center', transform=ax1.transAxes,
                           color='white', fontsize='large', weight='bold')

            # Make the text stand out even better using matplotlib's path effects
            outline_effect = [patheffects.withStroke(linewidth=2, foreground='black')]
            text_time.set_path_effects(outline_effect)
            text_channel.set_path_effects(outline_effect)


            add_logo(fig, size='large')

            plt.savefig('plots/{:06d}.png'.format(frame_number+start_number))
            plt.close()
            frame_number += 1
