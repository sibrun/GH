"""Module to track information during computations."""

import multiprocessing
import queue
import collections
import pandas
import tempfile
import webbrowser
from urllib.request import pathname2url
import Parameters
import StoreLoad
import os


def plot_info(data_list, header_list, path, to_html=True, to_csv=False):
    """Write a data list to a html or csv file.

    :param data_list:
    :type data_list: list(list)
    :param header_list: List with the parameter and property names.
    :type header_list: list(str)
    :param path: Path to the file without suffix.
    :type path: path
    :param to_html: Option to write data to a html file (Default: True).
    :type to_html: bool
    :param to_csv: Option to write data to a csv file (Default: False).
    :type to_csv: bool
    :return:
    """
    StoreLoad.generate_path(path)
    data_frame = pandas.DataFrame(data=data_list, columns=header_list)
    if to_html:
        html_path = path + '.html'
        data_frame.to_html(html_path)
    if to_csv:
        csv_path = path + '.csv'
        data_frame.to_csv(csv_path)


def display_html(shtml):
    """Display the given html code in a browser"""
    s = 'temp.html'
    temp_path = os.path.join(Parameters.temp_folder, s)
    with open(temp_path, "w") as html_file:
        html_file.write(shtml)
    url = 'file:{}'.format(pathname2url(os.path.abspath(temp_path)))
    webbrowser.open_new_tab(url)

def display_html_body(s):
    """Displays a html page with the string s as body"""
    pre='''<!DOCTYPE html>
    <html>
    <head>
    <meta
    http - equiv = "refresh"
    content = "10">
    <title> graph
    list </title>
    <style>
    p
    {
        text - align: center;
    }
    div
    {
        display: inline - block;
    }
    img
    {
        height: 200px;
    max - width:200
    px;
    width: expression(this.width > 200 ? 200: true);
    }
    svg
    {
        height: 200px;
    max - width:200
    px;
    width: expression(this.width > 200 ? 200: true);
    }
    </style>
    </head>
    <body>'''

    post='</body></html>'
    display_html(pre+s+post)


class InfoTracker(object):
    """Track information during computations and displays it as table on a html page.

    Attributes:
        - name (str): Name/title.
        - parameter_names_list (list(str)): List with the parameter and property names.
        - data_dict (dict(tuple -> tuple/list)): Dictionary (parameters -> properties) containing the data.
        - queue (Queue): Queue to push new data.
        - p (multiprocessing.Process): Process which runs the information trakcer.
        - f (tempfile.NamedTemporaryFile): Temporary html file.
        - html_path (path): Path to the temporary html file.
        - url (url): Url of the html file to open in the webbrowser.
    """
    def __init__(self, name):
        """Initialize the info tracker.

        :param name: Name describing the info tracker.
        :type name: str
        """
        self.name = name
        self.header_list = []
        self.data_dict = collections.OrderedDict()
        self.queue = multiprocessing.Queue()
        self.p = None
        self.f = tempfile.NamedTemporaryFile(delete=False)
        self.html_path = self.f.name + '.html'
        self.url ='file:{}'.format(pathname2url(self.html_path))

    def set_header_list(self, header_list):
        """Set the names of parameters and properties.

        :param header_list: List with the parameter and property names.
        :type header_list: list(str)
        """
        self.header_list = header_list

    def get_header_list(self):
        """Return a list with the names of parameters and properties.

        :return: parameter_list: List with the parameter and property names.
        :rtype list(str)
        """
        return self.header_list

    def get_data_dict(self):
        """Return the data dictionary.

        :return: Dictionary (parameters -> properties) containing the data.
        :rtype: dict(tuple -> tuple/list)
        """
        return self.data_dict

    def get_data_list(self):
        """Return a list of data entries.

        :return: Data in form of a list of lists.
        :rtype: list(list)
        """
        data_list = []
        for (params, properties) in self.data_dict.items():
            data_list.append(list(params) + list(properties))
        return data_list

    def update_data(self, data_dict):
        """Update the data dictionary with data_dict.

        :param data_dict: Dictionary (parameters -> properties) with new data to update the data dictionary.
        :type data_dict: dict(tuple -> tuple/list)
        """
        self.data_dict.update(data_dict)

    def get_pandas_df(self):
        """Return a pandas data frame with the information table.

        :return: Data Frame with the stored information.
        :rtype: pandas.DataFrame
        """
        return pandas.DataFrame(data=self.get_data_list(), columns=self.get_header_list())

    def update_html_file(self):
        """Update the html file."""
        self.get_pandas_df().to_html(self.html_path)

    def update(self, data_dict):
        """Update the data dictionary with data_dict and updates the html file.

        :param data_dict: New data to update the data dictionary.
        :type data_dict: dict
        """
        self.update_data(data_dict)
        self.update_html_file()

    def show_in_browser(self):
        """Show the information as table in a html file. Opens the file in the webbrowser."""
        webbrowser.open(self.url)

    def get_url(self):
        """Return the url for the webpage.

        :return: Webpage with the information table.
        :rtype: url
        """
        return self.url

    def get_queue(self):
        """Return the queue for the information input.

        :Note: Push updated information to the que in the form of a dict(tuple -> tuple/list)
               dictionary (parameters -> properties).

        :return: Queue for information input.
        :rtype: Queue
        """
        return self.queue

    def start(self):
        """Start tracking information.

        Open the html file with the table of information in the webbrowser.
        """
        self.update_html_file()
        self.show_in_browser()
        self.p = multiprocessing.Process(target=self.track)
        self.p.start()

    def stop(self):
        """Stop tracking information."""
        self.queue.put('stop')
        self.p.join()

    def track(self, timeout=Parameters.data_tracker_timeout):
        """Track information.

        Wait for new data from the queue and updates the html file.

        :param timeout: Timeout to get new data from queue.
        :type timeout: float
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
