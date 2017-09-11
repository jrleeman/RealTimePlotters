from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from obspy.clients.fdsn import Client
from obspy import UTCDateTime
import numpy as np
from matplotlib.dates import DateFormatter, DayLocator, HourLocator

def GetData(t0, net, st0, loc, ch, duration):
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
    st.taper(max_percentage=0.1)
    st.remove_response(output='DISP')
    return st

def calculate_amplitude(st):
    """
    Calculate a noise proxy for the first trace in that station.
    """
    res = plt.psd(st[0].data,Fs=st[0].stats.sampling_rate,NFFT=4096,noverlap=4096/2)
    amp_sum = np.sum(res[0])
    return amp_sum

def plot_amplitude_timeseries(amplitudes, station, channel):
    """
    Make a plot of the amplitude time series and save it with a unique
    filename that will be overwritten if the script is run again later.
    (i.e. not having time in the name)
    """
    fig = plt.figure(figsize=(10,7.5))
    ax1 = plt.subplot(111)

    # Set labels and tick sizes
    ax1.set_xlabel(r'Time', fontsize=16)
    ax1.set_ylabel(r'Relative Seismic Noise', fontsize=16)

    # Turn off top and right tick marks
    ax1.get_xaxis().tick_bottom()
    ax1.get_yaxis().tick_left()

    # Make ticks readable
    ax1.tick_params(axis='both', which='major', labelsize=14)

    # Turn off top and right splines
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # Plotting
    ax1.plot(amplitudes[:,0], amplitudes[:,1]/np.min(amplitudes[:,1]))

    # Setup the x axis time ticking
    ax1.xaxis.set_major_formatter(DateFormatter('%m/%d'))
    ax1.xaxis.set_minor_formatter(DateFormatter('%HZ'))
    ax1.xaxis.set_major_locator(DayLocator())
    ax1.xaxis.set_minor_locator(HourLocator(range(2, 24, 2)))

    # Add a title
    plt.title('Station {} {}'.format(station, channel), fontsize=18)

    # Save the figure
    plt.savefig('plots/seismic_amplitude_{}_{}.png'.format(station, channel), bbox_inches='tight');

#
# Only change these parameters
#

duration = 30 # How long each analysis window is in minutes
history = 24 # How much history to plot in hours

# List of station property dictionaries - a plot will be made for each
stations = [{'station': 'SJ07', 'network': 'PR', 'location': '--', 'channel': 'HHN'},
            {'station': 'SJ07', 'network': 'PR', 'location': '--', 'channel': 'HHE'},
            {'station': 'SJ07', 'network': 'PR', 'location': '--', 'channel': 'HHZ'},
            {'station': 'SC01', 'network': 'DR', 'location': '--', 'channel': 'BHE'},
            {'station': 'SC01', 'network': 'DR', 'location': '--', 'channel': 'BHN'},
            {'station': 'SC01', 'network': 'DR', 'location': '--', 'channel': 'BHZ'},
            {'station': 'CAMR', 'network': 'CW', 'location': '00', 'channel': 'HHE'},
            {'station': 'CAMR', 'network': 'CW', 'location': '00', 'channel': 'HHN'},
            {'station': 'CAMR', 'network': 'CW', 'location': '00', 'channel': 'HHZ'},
            {'station': 'GRTK', 'network': 'CU', 'location': '00', 'channel': 'BH1'},
            {'station': 'GRTK', 'network': 'CU', 'location': '00', 'channel': 'BH2'},
            {'station': 'GRTK', 'network': 'CU', 'location': '00', 'channel': 'BHZ'},
            {'station': '061Z', 'network': 'N4', 'location': '--', 'channel': 'BHE'},
            {'station': '061Z', 'network': 'N4', 'location': '--', 'channel': 'BHN'},
            {'station': '061Z', 'network': 'N4', 'location': '--', 'channel': 'BHZ'},
            {'station': 'SMRT', 'network': 'NA', 'location': '--', 'channel': 'BHE'},
            {'station': 'SMRT', 'network': 'NA', 'location': '--', 'channel': 'BHN'},
            {'station': 'SMRT', 'network': 'NA', 'location': '--', 'channel': 'BHZ'}]

stations = [{'station': 'AGPR', 'network': 'PR', 'location': '--', 'channel': 'HNZ'}]

#
# Change nothing below here!
#

# Based on the current time, calculate the start and end time to the nearest
# hour to the present.
now = datetime.utcnow()
end_time = datetime(now.year, now.month, now.day, now.hour)
start_time = end_time - timedelta(hours=history)

# Loop over each item in the stations list
for item in stations:

    station = item['station']
    network = item['network']
    location = item['location']
    channel = item['channel']

    print('Analyzing {} - {}'.format(station, channel))

    data_start_time = start_time
    amplitudes = []
    while data_start_time < end_time:
        try:
            st = GetData(data_start_time, network, station, location, channel, 30)
            amplitude_sum = calculate_amplitude(st)
            amplitudes.append([data_start_time + timedelta(minutes=duration/2), amplitude_sum])
        except:
            # Can't get the data for some reason, just skip over it
            pass
        # Increment
        data_start_time += timedelta(minutes=duration)

    amplitudes = np.array(amplitudes)

    if len(amplitudes) > 0:
        plot_amplitude_timeseries(amplitudes, station, channel)
