import numpy as np
import glob
import re
import netCDF4
import bisect
import xarray

def copy_variables(src, dst):
    """
    Helper function to copy all variables in a NetCDF group.
    Params:
        src: netCDF4 group or file handle whose content to copy.
        dest: netCDF4 group or file handle into which to copy
            the variables.
    """
    dst = dst.createGroup("SingleScatteringData")
    for name, dimension in src.dimensions.items():
        dst.createDimension(
            name, (len(dimension) if not dimension.isunlimited() else None))
        # copy all file data except for the excluded
    for name, variable in src.variables.items():
        x = dst.createVariable(name, variable.datatype, variable.dimensions)
        dst[name][:] = src[name][:]

class ScatteringData:
    """
    The ScatteringData class holds the scattering data of a particle for
    a specific frequency and temperature. The data thus depends only
    on the incoming and outgoing angles.
    """
    def __init__(self,
                 azimuth_angles_incoming,
                 zenith_angles_incoming,
                 azimuth_angles_scattering,
                 zenith_angles_scattering,
                 phase_matrix_data,
                 extinction_matrix_data,
                 absorption_vector_data):
        self._azimuth_angles_incoming = azimuth_angles_incoming
        self._zenith_angles_incoming = zenith_angles_incoming
        self._azimuth_angles_scattering = azimuth_angles_scattering
        self._zenith_angles_scattering = zenith_angles_scattering
        self._phase_matrix_data = phase_matrix_data
        self._extinction_matrix_data = extinction_matrix_data
        self._absorption_vector_data = absorption_vector_data

    @property
    def lon_inc(self):
        return self._azimuth_angles_incoming

    @property
    def lat_inc(self):
        return self._zenith_angles_incoming

    @property
    def lon_scat(self):
        return self._azimuth_angles_scattering

    @property
    def lat_scat(self):
        return self._zenith_angles_scattering

    @property
    def phase_matrix(self):
        return self._phase_matrix_data

    @property
    def extinction_matrix(self):
        return self._extinction_matrix_data

    @property
    def absorption_vector(self):
        return self._absorption_vector_data

class Particle:
    """
    The particle class represents scattering data for a specific particle. Its
    purpose is to give access to the scattering properties at the different frequencies
    and temperatures that are available.
    """
    def _parse_temperatures_and_frequencies(self):
        group_names = self.file_handle.groups.keys()
        temps = []
        freqs = []
        for gn in group_names:
            match = re.match('Freq([0-9\.]*)GHz_T([0-9\.]*)K', gn)
            freq = match.group(1)
            temp = match.group(2)
            temps.append(temp)
            freqs.append(freq)
        self.temperatures = np.array(temps, dtype=np.float)
        self.frequencies = np.array(freqs, dtype=np.float)

        key = lambda i: (freqs[i], temps[i])
        indices = list(range(self.frequencies.size))
        indices.sort(key=key)

        self.frequencies = self.frequencies[indices]
        self.temperatures = self.temperatures[indices]
        self.frequencies = list(set(self.frequencies))
        self.frequencies.sort()
        self.temperatures = list(set(self.temperatures))
        self.temperatures.sort()
        self.keys = [(f, t) for f in self.frequencies for t in self.temperatures]

    def __init__(self, filename):
        """
        Create a particle object by reading a NetCDF4 file from the ARTS SSDB.
        Args:
            filename: The path of the NetCDF4 file containign the data.
        """
        self.filename = filename
        self.file_handle = netCDF4.Dataset(filename)
        self._parse_temperatures_and_frequencies()

    def __len__(self):
        """
        Returns:
            The number of available unique frequency and temperature
            keys.
        """
        return len(self.frequencies) * len(self.temperatures)

    def __getitem__(self, *args):
        """
        Access scattering data for given frequency and temperature index.
        """
        if len(args) == 2:
            i, j = args
            return self.get_scattering_data(self.frequencies[i],
                                            self.temperatures[j])
        else:
            f, t = self.keys[args[0]]
            return self.get_scattering_data(f, t)

    def get_scattering_data(self, frequency, temperature):
        """
        Return scattering data for given frequency and temperature.


        """
        requested = (frequency, temperature)
        index = bisect.bisect_left(self.keys, requested)
        found = self.keys[index]
        print(index)
        print(found)
        if not (np.all(np.isclose(requested, found))):
            raise Exception("Could not find scattering data for given temperature and frequency")

        group_name = list(self.file_handle.groups.keys())[index]
        group = self.file_handle[group_name]["SingleScatteringData"]

        azimuth_angles_incoming = group["aa_inc"][:]
        zenith_angles_incoming = group["za_inc"][:]
        azimuth_angles_scattering = group["aa_scat"][:]
        zenith_angles_scattering = group["za_scat"][:]
        phase_matrix_data = group["phaMat_data"][:]
        extinction_matrix_data = group["extMat_data"][:]
        absorption_vector_data = group["absVec_data"][:]

        return ScatteringData(azimuth_angles_incoming,
                              zenith_angles_incoming,
                              azimuth_angles_scattering,
                              zenith_angles_scattering,
                              phase_matrix_data,
                              extinction_matrix_data,
                              absorption_vector_data)

class Habit:

    @staticmethod
    def get_particle_props(path):
        filename = os.path.basename(path)
        match = re.match('Dmax([0-9]*)um_Dveq([0-9]*)um_Mass([-0-9\.e]*)kg\.nc', filename)
        dmax = np.float32(match.group(1)) * 1e-6
        dveq = np.float32(match.group(2)) * 1e-6
        mass = np.float32(match.group(3))
        return (dmax, dveq, mass)

    def __init__(self,
                 name,
                 phase,
                 kind,
                 riming,
                 orientation,
                 path):
        self._name = name
        self._phase = phase
        self._kind = kind
        self._riming = riming
        self._orientation = orientation

        self.files = glob.glob(os.path.join(path, "Dmax*Dveq*Mass*.nc"))
        properties = np.array([Habit.get_particle_props(f) for f in self.files])

        self.d_eq = properties[:, 0]
        self.d_max = properties[:, 1]
        self.mass = properties[:, 2]

        indices = np.argsort(self.d_eq)
        self.d_eq = self.d_eq[indices]
        self.d_max = self.d_max[indices]
        self.mass = self.mass[indices]
        self.files = [self.files[i] for i in indices]

    @property
    def name(self):
        return self._name


    def __repr__(self):
        s = f"SSDB Particle: {self._name}\n"
        s += f"\t Phase: {self._phase}, Kind: {self._kind}, "
        s += f"Riming: {self._riming}, Orientation: {self._orientation}"
        return s

class SSDBReader:

    @staticmethod
    def particle_from_path(path):
        match = re.match('.*/([^/]*)/([^/]*)/([^/]*)/(.*)_Id[0-9]*/([^/]*).*', path)
        phase = match.group(1)
        kind = match.group(2)
        riming = match.group(3)
        name = match.group(4)
        orientation = match.group(5)
        return Habit(name, phase, kind, riming, orientation, path)

    def __init__(self, path):
        search_path = os.path.join(path, "**/Dmax*Dveq*Mass*.nc")
        files = glob.glob(search_path, recursive=True)
        self.particle_paths = list(set([os.path.dirname(f) for f in files]))
        self.habits = [SSDBReader.particle_from_path(p) for p in self.particle_paths]

    def get_available_habits(self):
        return [h.name for h in self.habits]
def get_habit(self, name):
        return self.habits[[h.name for h in self.habits].index(name)]


def extract_data(particle_path, filename):
    """
    Function to extract test data from particle data.
    Args:
        particle_path: Path to particle data.
        filename: Filename for output
    """
    p = Particle(particle_path)
    temps = list(set([k[1] for k in p.keys]))
    freqs = list(set(k[0] for k in p.keys))
    temps = [temps[0], temps[-1]]
    freqs = [freqs[0], freqs[-1]]
    new_file = netCDF4.Dataset(filename, "w")
    for f in freqs:
        for t in temps:
            group_name = f"Freq{f:08.3f}GHz_T{t:6.2f}K"
            src = p.file_handle[group_name + "/SingleScatteringData"]
            dst = new_file.createGroup(group_name)
            copy_variables(src, dst)
    new_file.close()

# particle_spherical_1 = "/home/simonpf/Dendrite/SSDB/SSDB_EUMETSAT/SSD_TRO/TotallyRandom/Ice/SingleCrystals/Pristine/IceSphere_Id24/TotallyRandom/Dmax38098um_Dveq38098um_Mass2.65425e-02kg.nc"
# particle_spherical_2 = "/home/simonpf/Dendrite/SSDB/SSDB_EUMETSAT/SSD_TRO/TotallyRandom/Ice/SingleCrystals/Pristine/IceSphere_Id24/TotallyRandom/Dmax00004um_Dveq00004um_Mass3.33452e-14kg.nc"

# particle_random_1 = "/home/simonpf/Dendrite/SSDB/SSDB_EUMETSAT/SSD_TRO/TotallyRandom/Ice/SingleCrystals/Pristine/SectorSnowflake_Id3/TotallyRandom/Dmax03369um_Dveq00771um_Mass2.19881e-07kg.nc"
# particle_random_2 = "/home/simonpf/Dendrite/SSDB/SSDB_EUMETSAT/SSD_TRO/TotallyRandom/Ice/SingleCrystals/Pristine/SectorSnowflake_Id3/TotallyRandom/Dmax03369um_Dveq00771um_Mass2.19881e-07kg.nc"

# particle_ar_1 = "/home/simonpf/Dendrite/SSDB/SSDB_EUMETSAT/SSD_ARO/AzimuthallyRandom/Ice/SingleCrystals/Pristine/PlateType1_Id9/AzimuthallyRandom_beta010.0deg/Dmax04151um_Dveq01257um_Mass9.53208e-07kg.nc"
# particle_ar_2 = "/home/simonpf/Dendrite/SSDB/SSDB_EUMETSAT/SSD_ARO/AzimuthallyRandom/Ice/SingleCrystals/Pristine/PlateType1_Id9/AzimuthallyRandom_beta010.0deg/Dmax00590um_Dveq00251um_Mass7.59425e-09kg.nc"

# extract_data(particle_spherical_1, "test_data_spherical_1.nc")
# extract_data(particle_spherical_2, "test_data_spherical_2.nc")
# extract_data(particle_random_1, "test_data_random_1.nc")
# extract_data(particle_random_2, "test_data_random_2.nc")
# extract_data(particle_ar_1, "test_data_azimuthally_random_1.nc")
# extract_data(particle_ar_2, "test_data_azimuthally_random_2.nc")




#d = p.get_scattering_data(886.4, 270.0)
#
#dims = ("azimuth_angles_incoming", "zenith_angles_incoming", "azimuth_angles_scattering", "zenith_angles_scattering")
#ds = xarray.Dataset({"phase_matrix": (("phase_matrix_elements",) + dims, d._phase_matrix_data.data),
#                     "extinction_matrix": (("absorption_vector_elements",) + dims[:2], d._extinction_matrix_data.data),
#                     "absorption_vector": (("absorption_vector_elements",) + dims[:2], d._absorption_vector_data.data)})
#ds["azimuth_angles_incoming"] = d._azimuth_angles_incoming
#ds["zenith_angles_incoming"] = d._zenith_angles_incoming
#ds["azimuth_angles_scattering"] = d._azimuth_angles_scattering
#ds["zenith_angles_scattering"] = d._zenith_angles_scattering

