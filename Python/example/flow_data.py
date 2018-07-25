"""class to read flow data from text"""
import numpy as np
import re

class data:
    """
    Read and store flow cytometry data into a workable format
    each instance of this object will have the data, and list of
    channels
    """
    def __init__(self, filename):
        self.filename = filename
        self.data = self._Altload()

    def _load(self):
        '''
        Reads tab delimited file of flow cytometry data

        Returns :
        FACs data from text file, note all values <= 0 ignored
        '''
        # read txt file
        fid = open(self.filename, 'r')
        # split based upon delimiter, all
        data = re.split('\t|\n', fid.read())
        fid.close()

        # find header
        k = 0
        while re.findall('[A-Za-z]+', data[k]) != []:
            k += 1
        # set channels
        self.channels = data[:k]

        # eliminate empty lines
        j = -1
        while data[j] == '':
            j -= 1
        # change data to numpy array and reshape
        data = np.array(data[k:j+1], dtype='f')
        data = np.reshape(data,
                          (data.size/k, k))
        data[data <= 0] = np.nan
        return data[np.isnan(np.sum(data, 1)) == 0, :]

    def _Altload(self):
        '''
        Reads tab delimited file of flow cytometry data
        Returns :
        FACs data from text file, note all values <= 0 ignored
        '''
        import csv
        # read txt file
        fid = open(self.filename, 'r')
        
                # split based upon delimiter, all
        lol = list(csv.reader(fid, delimiter='\t'))
        self.channels = lol[0]

        # change data to numpy array and reshape
        data = np.array(lol[1:], dtype='f')
        data[data <= 0] = np.nan
        return data[np.isnan(np.sum(data, 1)) == 0, :]



    def _get_channel(self, marker):
        '''
        Retrieve channel idex from marker names.
        '''
        marker = marker.lower()
        chan_count = 0
        # split text by non-character thing
        while re.split('\s', self.channels[chan_count])[-1].lower() != marker:
            chan_count += 1
        return chan_count

    def get_data(self, marker):
        '''
        Retrieve data by marker name.
        '''
        return self.data[:, self._get_channel(marker)]

    def print_channels(self):
        count = 0
        for wchan in self.channels:
            print('{0}) {1}'.format(count, wchan))
            count += 1
