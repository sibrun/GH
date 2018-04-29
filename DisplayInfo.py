import time
import multiprocessing as mp
import pandas
import tempfile
import webbrowser
from urllib import pathname2url
import Parameters


class DataFrame(object):
    def __init__(self, name, parameters_dict):
        self.name = name
        self.parameters_dict = parameters_dict
        self.manager = mp.Manager()
        self.data_dict = self.manager.dict()
        self.f = tempfile.NamedTemporaryFile(delete=False)
        self.html_path = self.f.name + '.html'
        self.url ='file:{}'.format(pathname2url(self.html_path))

    def get_parameters_dict(self):
        return self.parameters_dict

    def get_parameters_list(self):
        return self.parameters_dict.items()

    def get_dict(self):
        return self.data_dict

    def get_list(self):
        return self.data_dict.items()

    def get_pandas_df(self):
        return pandas.DataFrame(data=self.get_list(), columns=self.get_parameters_list())

    def update_html_file(self):
        self.get_pandas_df().to_html(self.html_path)

    def show_in_browser(self):
        webbrowser.open(self.url)

    def track(self, delay_time = Parameters.data_tracker_delay):
        while True:
            self.update_html_file()
            time.sleep(delay_time)


