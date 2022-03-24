import yaml
from enum import Enum
from pathlib import Path
import logging
import sys
logging.getLogger('matplotlib').setLevel(logging.ERROR)

CONFIG_TYPE = Enum("CONFIG_TYPE", 'REPORT BENCH COMPARE')


class Config:
    def __init__(self, config_file):
        self.config_file = config_file
        self.experiment_dir = str(Path(config_file).parent.resolve())
        self.log_dir = self.experiment_dir + "/logs/"
        self.report_dir = self.experiment_dir + "/reports/"

        # setup the experiment directory structure
        Path(self.log_dir).mkdir(parents=True, exist_ok=True)
        Path(self.report_dir).mkdir(parents=True, exist_ok=True)

        # logging
        self.log_file = self.log_dir + 'main.log'
        # noinspection PyArgumentList
        logging.basicConfig(format='[%(levelname)s] %(message)s', level=logging.DEBUG,
                            handlers=[logging.FileHandler(self.log_file, mode='w'), logging.StreamHandler(sys.stdout)])

        # shared configs
        # ...
        logging.info(self)

    def __str__(self):
        s = " ===== Config ====="
        s += "\n\tYAML config file: " + self.config_file
        s += "\n\tExperiment directory: " + self.experiment_dir
        s += "\n\tMain LOG file: " + str(self.log_file) + "\n\t"
        return s


class ReportConfig(Config):
    def __init__(self, config_file, **entries):
        self.__dict__.update(entries)
        super().__init__(config_file)
        # additional initialization

    def __str__(self):
        s = super().__str__()
        s += '\n\t'.join("{}: {}".format(k, v) for k, v in self.__dict__.items())
        return s

def load_config(fname, config_type):
    """
    Load a YAML configuration file
    """
    with open(fname) as file:
        config = yaml.load(file, Loader=yaml.FullLoader)
    if config_type == CONFIG_TYPE.REPORT:
        return ReportConfig(fname, **config)
    else:
        return None
