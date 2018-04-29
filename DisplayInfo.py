import multiprocessing as mp
import Queue
import pandas
import tempfile
import webbrowser
from urllib import pathname2url
import Parameters


class InfoTracker(object):
    def __init__(self, name):
        self.name = name
        self.parameter_list = []
        self.data_dict = dict()
        self.queue = mp.Queue()
        self.p = None
        self.f = tempfile.NamedTemporaryFile(delete=False)
        self.html_path = self.f.name + '.html'
        self.url ='file:{}'.format(pathname2url(self.html_path))

    def set_parameter_list(self, parameter_list):
        self.parameter_list = parameter_list

    def get_parameter_list(self):
        return self.parameter_list

    def get_dict(self):
        return self.data_dict

    def get_list(self):
        L = self.data_dict.items()

    def update_data(self, data_dict):
        self.data_dict.update(data_dict)

    def get_pandas_df(self):
        return pandas.DataFrame(data=self.get_list(), columns=self.get_parameter_list())

    def update_html_file(self):
        self.get_pandas_df().to_html(self.html_path)

    def update(self, data_dict):
        self.update_data(data_dict)
        self.update_html_file()

    def show_in_browser(self):
        webbrowser.open(self.url)

    def get_url(self):
        return self.url

    def get_queue(self):
        return self.queue

    def start(self):
        self.update_html_file()
        self.p = mp.Process(target=self.track)
        self.p.start()

    def stop(self):
        self.queue.put('stop')
        self.p.join()

    def track(self, timeout=Parameters.data_tracker_timeout):
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
