"""Module to track information during computations."""

import multiprocessing
import Queue
import collections
import pandas
import tempfile
import webbrowser
from urllib import pathname2url
import Parameters


class InfoTracker(object):
    """Tracks information during computations and displays it as table on a html page.

    Attributes:
        name (str): Name/title.

        parameter_names_list (list(str)): List with the parameter and property names.

        data_dict (dict(tuple -> tuple/list)): Dictionary (parameters -> properties) containing the data.

        queue (Queue): Queue to push new data.

        p (multiprocessing.Process): Process which runs the information trakcer.

        f (tempfile.NamedTemporaryFile): Temporary html file.

        html_path (path): Path to the temporary html file.

        url (url): Url.
    """
    def __init__(self, name):
        self.name = name
        self.parameter_names_list = []
        self.data_dict = collections.OrderedDict()
        self.queue = multiprocessing.Queue()
        self.p = None
        self.f = tempfile.NamedTemporaryFile(delete=False)
        self.html_path = self.f.name + '.html'
        self.url ='file:{}'.format(pathname2url(self.html_path))

    def set_parameter_names_list(self, parameter_names_list):
        """Sets the name of parameters and properties.

        :param parameter_names_list: list(str): List with the parameter and property names.
        """
        self.parameter_names_list = parameter_names_list

    def get_parameter_names_list(self):
        """Returns a list with the names of parameters and properties.

        :return: parameter_list: list(str): List with the parameter and property names.
        """
        return self.parameter_names_list

    def get_data_dict(self):
        """Returns the data dictionary.

        :return: dict(tuple -> tuple/list): Dictionary (parameters -> properties) containing the data.
        """
        return self.data_dict

    def get_data_list(self):
        """Returns a list of data entries.

        :return: list(list(int)): Data in form of a list of lists.
        """
        data_list = []
        for (params, properties) in self.data_dict.items():
            data_list.append(list(params) + list(properties))
        return data_list

    def update_data(self, data_dict):
        """Updates the data dictionary with data_dict.

        :param data_dict: dict(tuple -> tuple/list): Dictionary (parameters -> properties) with new data to update the data dictionary.
        """
        self.data_dict.update(data_dict)

    def get_pandas_df(self):
        """Returns a pandas data frame with the information table.

        :return: pandas.DataFrame: Data Frame with the stored information.
        """
        return pandas.DataFrame(data=self.get_data_list(), columns=self.get_parameter_names_list())

    def update_html_file(self):
        """Updates the html file."""
        self.get_pandas_df().to_html(self.html_path)

    def update(self, data_dict):
        """Updates the data dictionary with data_dict and updates the html file.

        :param data_dict: dict: New data to update the data dictionary.
        """
        self.update_data(data_dict)
        self.update_html_file()

    def show_in_browser(self):
        """Shows the information as table in a html file. Opens the file in the webbrowser."""
        webbrowser.open(self.url)

    def get_url(self):
        """Returns the url for the webpage.

        :return: url: Webpage with the information table.
        """
        return self.url

    def get_queue(self):
        """Returns the queue for the information input.

        Push updated information to the que in the form of a dict(tuple -> tuple/list) dictionary (parameters -> properties).

        :return: Queue: Queue for information input.
        """
        return self.queue

    def start(self):
        """Starts tracking information.

        Opens the html file with the table of information in the webbrowser.
        """
        self.update_html_file()
        self.show_in_browser()
        self.p = multiprocessing.Process(target=self.track)
        self.p.start()

    def stop(self):
        """Stops tracking information."""
        self.queue.put('stop')
        self.p.join()

    def track(self, timeout=Parameters.data_tracker_timeout):
        """Tracks information.

        Waits for new data from the queue and updates the html file.

        :param timeout: positive float: Timeout to get new data from queue.
        """
        while True:
            try:
                message = self.queue.get(timeout=timeout)
                if str(message) == 'stop':
                    self.update_html_file()
                    break
                else:
                    self.update(message)
            except Queue.Empty:
                continue
