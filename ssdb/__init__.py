import os
from ftplib import FTP
from tqdm import tqdm, tqdm_notebook

from typhon.arts.xml import load, save

if "SSDB_PATH" in os.environ:
    ssdb_path = os.environ["SSDB_PATH"]
else:
    ssdb_path = os.path.join(os.environ["HOME"], "ssdb")

    if not os.path.exists(ssdb_path):
        os.makedirs(ssdb_path)

ftp = FTP("ftp-projects.cen.uni-hamburg.de")
ftp.login()
ftp.cwd("/arts/ArtsScatDbase/v1.0.0")

def create_directories(ssdb_path):
    """
    Create SSDB folder structure.
    """
    if not os.path.exists(ssdb_path):
        os.makedirs(ssdb_path)

    paths = [os.path.join(ssdb_path, "StandardHabits"),
             os.path.join(ssdb_path, "StandardHabits", "FullSet")]
    for p in paths:
        if not os.path.exists(p):
            os.makedirs(p)

def get_standard_habits():
    """
    Returns list of standard habits.
    """
    ftp.cwd("/arts/ArtsScatDbase/v1.0.0/StandardHabits/FullSet")
    names = ftp.nlst()
    names = set([s.split(".")[0] for s in names])
    return names

class FileDownload:
    def __init__(self, path):
        self.path =path
        self.directory, self.filename = os.path.split(path)

    def __enter__(self):
        self.file = open(self.path, "wb")
        ftp.voidcmd('TYPE I')
        self.size = ftp.size(self.filename)
        self.bar  = tqdm(total = self.size)
        self.stored = 0
        return self

    def __exit__(self, type, value, traceback):
        self.file.close()

    def callback(self, b):
        self.file.write(b)
        self.stored += len(b) / 8
        self.bar.update(len(b))


class Habit:

    @property
    def path(self):
        return os.path.join(ssdb_path, "StandardHabits", "FullSet", self.name + ".xml")

    @property
    def meta_path(self):
        return os.path.join(ssdb_path, "StandardHabits", "FullSet", self.name + ".meta.xml")

    @property
    def local(self):
        return os.path.exists(self.path) and os.path.exists(self.path + ".bin") \
            and os.path.exists(self.meta_path)

    @property
    def data(self):
        if Habit.local:
            return load(self.path)
        else:
            raise Exception("Data for {0} is not available locally.".format(self.name))

    @property
    def meta_data(self):
        if Habit.local:
            return load(self.meta_path)
        else:
            raise Exception("Data for {0} is not available locally.".format(self.name))

    def download(self):
        ftp.cwd("/arts/ArtsScatDbase/v1.0.0/StandardHabits/FullSet")

        with FileDownload(self.path) as fp:
            ftp.retrbinary("RETR " + self.name + ".xml", fp.callback)

        with FileDownload(self.path + ".bin") as fp:
           ftp.retrbinary("RETR " + self.name + ".xml.bin", fp.callback)

        with FileDownload(self.meta_path) as fp:
            ftp.retrbinary("RETR " + self.name + ".meta.xml", fp.callback)

    def __init__(self, name):
        self.name = name

for n in get_standard_habits():
    globals()[n] = Habit(n)
