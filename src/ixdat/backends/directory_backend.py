import json
import numpy as np
from .memory_backend import BackendBase
from ..config import CFG


char_substitutions = {
    "/": "_DIV_",  # slash (divided by)
    "\\": "_BKSL_",  # backslash
    ".": "_DOT_",  # decimal
    "^": "_CFLX_",  # circumflex accent (raised-to-the)
    "<": "_LTS_",  # less-than sign
    ">": "_GTS_",  # greater-than sign
}
# TODO: consider implementing some kind of general solution with a tmp dir
#    see: https://github.com/ixdat/ixdat/pull/5#discussion_r565075588


def fix_name_for_saving(name):
    """Replace problematic characters in name with the substitutions defined above"""
    for bad_char, substitution in char_substitutions.items():
        name = name.replace(bad_char, substitution)
    return name


def id_from_path(path):
    """Return the id (int) of the row represented by given path to an ixdat file"""
    try:
        return int(path.stem.split("_")[0])
    except ValueError:
        print(f"couldn't find id in {path}")  # debugging
        return None


def name_from_path(path):
    """Return the name (str) of the row represented by given path to an ixdat file"""
    return path.stem.split("_", 1)[1]


class DirBackend(BackendBase):
    """A database backend that loads and saves .ix files from a directory

    TODO: refactor i with the standard id_ inside methods
        see github: https://github.com/ixdat/ixdat/pull/1#discussion_r546400226
    """

    def __init__(
        self,
        directory=CFG.standard_data_directory,
        project_name=CFG.default_project_name,
        metadata_suffix=CFG.standard_metadata_suffix,
        data_suffix=CFG.standard_data_suffix,
    ):
        """Initialize a directory database backend with the directory as Path

        Args:
            directory (Path): the main ixdat directory
            project_name (str): the name of the project (ixdat subdirectory)
            metadata_suffix (str): The suffix to use for JSON-formatted metadata files
            data_suffix (str): The suffix to use for numpy-formatted data files
        """
        self.project_directory = directory / project_name
        try:
            self.project_directory.mkdir(parents=True, exist_ok=True)
        except Exception:
            raise  # TODO, figure out what gets raised, then except with line below
            # raise ConfigError(f"Cannot make dir '{self.standard_data_directory}'")

        self.metadata_suffix = metadata_suffix
        self.data_suffix = data_suffix
        super().__init__()

    @property
    def name(self):
        return f"DirBackend({self.project_directory})"

    def save(self, obj):
        """Save the Saveable object as a file corresponding to a row in a table"""
        if obj.data_objects:
            # save any data objects first as this may change the references
            for data_obj in obj.data_objects:
                self.save_data_obj(data_obj)
        table_name = obj.table_name
        obj_as_dict = obj.as_dict()
        i = self.add_row(obj_as_dict, table_name=table_name)
        obj.set_id(i)
        obj.set_backend(self)
        return i

    def get(self, cls, i):
        """Open a Saveable object represented as row i of table cls.table_name"""
        table_name = cls.table_name
        obj_as_dict = self.get_row_as_dict(table_name, i)
        i = obj_as_dict.pop("id", i)
        obj = cls.from_dict(obj_as_dict)
        obj.set_backend(self)
        obj.set_id(i)
        return obj

    def contains(self, table_name, i):
        """Check if id `i` is already a principle key in the table named `table_name`"""
        return i in self.get_id_list(table_name)

    def load_obj_data(self, obj):
        """Return the data for an object loaded from its .ixdata file"""
        path_to_row = self.get_path_to_row(obj.table_name, obj.id)
        try:
            return np.load(path_to_row.with_suffix(self.data_suffix))
        except FileNotFoundError:
            # there's no data to be got.
            print(f"could not find file {path_to_row}")
            return

    def save_data_obj(self, data_obj):
        """Save the object as a .ix for metadata and .ixdata for numerical data"""
        table_name = data_obj.table_name
        if data_obj.backend == self and self.contains(table_name, data_obj.id):
            return data_obj.id  # already saved!
        obj_as_dict = data_obj.as_dict()
        data = obj_as_dict["data"]
        obj_as_dict["data"] = None
        # first we save the metadata and set the object's id:
        i = self.add_row(obj_as_dict, table_name=table_name)
        data_obj.set_id(i)
        data_obj.set_backend(self)
        #  ... and now we save the data
        folder = self.project_directory / table_name
        fixed_name = fix_name_for_saving(data_obj.name)
        data_file_name = f"{data_obj.id}_{fixed_name}{self.data_suffix}"
        np.save(folder / data_file_name, data)
        return i

    def add_row(self, obj_as_dict, table_name):
        """Save object's serialization to the folder table_name (like adding a row)"""
        folder = self.project_directory / table_name
        if not folder.exists():
            folder.mkdir()
        i = self.get_next_available_id(table_name)
        obj_as_dict.update({"id": i})
        fixed_name = fix_name_for_saving(obj_as_dict["name"])
        file_name = f"{i}_{fixed_name}{self.metadata_suffix}"

        with open(folder / file_name, "w") as f:
            json.dump(obj_as_dict, f, indent=4)
        return i

    def get_row_as_dict(self, table_name, i):
        """Return the serialization of the object represented in row i of table_name"""
        path_to_row = self.get_path_to_row(table_name, i)
        with open(path_to_row, "r") as file:
            obj_as_dict = json.load(file)
        return obj_as_dict

    def get_path_to_row(self, table_name, i):
        """Return the Path to the file representing row i of the table `table_name`"""
        folder = self.project_directory / table_name
        for p in folder.iterdir():
            if id_from_path(p) == i and p.suffix == self.metadata_suffix:
                return p
        print(f"could not find row with id={i} in table '{table_name}'")
        print(f"looking in folder: {folder}")  # debugging
        return None  # if that row is not in the table.

    def get_id_list(self, table_name):
        """List the principle keys of the existing rows of a given table"""
        folder = self.project_directory / table_name
        id_list = []
        for file in folder.iterdir():
            if not file.is_dir():
                try:
                    i = int(file.name.split("_")[0])
                    id_list.append(i)
                except TypeError:
                    pass
        return id_list

    def get_next_available_id(self, table_name):
        """Return the next available id for a given table"""
        id_list = self.get_id_list(table_name)
        if not id_list:
            return 1
        return max(self.get_id_list(table_name)) + 1

    def __eq__(self, other):
        """Two DirBackends are equivalent if they refer to the same directory"""
        if other is self:
            return True
        if (
            hasattr(other, "directory")
            and other.project_directory.resolve() == self.project_directory.resolve()
            and other.project_directory.lstat() == self.project_directory.lstat()
            and other.__class__ is self.__class__
        ):
            return True
        return False
