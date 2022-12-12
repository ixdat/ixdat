"""This module implements a local json-file-based representation of a relational db

FIXME: Saving and/or loading get quite slow when the number of rows in a table (usually
  data_series) grows to hundreds. How to figure out why that happens?
  # see https://github.com/ixdat/ixdat/pull/11#discussion_r663468719
"""


import json
import numpy as np
from .backend_base import BackendBase
from ..config import config, prompt_for_permission


char_substitutions = {  # substitutions needed to name .json file with data series name
    "/": "_DIV_",  # slash (divided by)
    "\\": "_BKSL_",  # backslash
    ".": "_DOT_",  # decimal
    "^": "_CFLX_",  # circumflex accent (raised-to-the)
    "<": "_LTS_",  # less-than sign
    ">": "_GTS_",  # greater-than sign
}
# TODO: consider implementing some kind of general solution with a tmp dir
#   see: https://github.com/ixdat/ixdat/pull/5#discussion_r565075588


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

    backend_type = "directory"

    def __init__(
        self,
        directory=config.standard_data_directory,
        project_name=config.default_project_name,
        metadata_suffix=config.standard_metadata_suffix,
        data_suffix=config.standard_data_suffix,
    ):
        """Initialize a directory database backend with the directory as Path

        Args:
            directory (Path): the main ixdat directory
            project_name (str): the name of the project (ixdat subdirectory)
            metadata_suffix (str): The suffix to use for JSON-formatted metadata files
            data_suffix (str): The suffix to use for numpy-formatted data files
        """
        self.project_directory = directory / project_name
        self.project_directory.mkdir(parents=True, exist_ok=True)

        self.metadata_suffix = metadata_suffix
        self.data_suffix = data_suffix
        super().__init__()

    @property
    def address(self):
        """The directory containing the tables (folders) and rows (.ix files)"""
        return str(self.project_directory)

    def save(self, obj, force=False, no_updates=True):
        """Save the Savable object as a file corresponding to a row in a table

        Args:
            obj (Savable): an object
            force (bool): Whether to force updates if the object is already saved
            no_updates (bool): Whether to allow updates if the object is already saved.
                If both force and no_updates are False, the user will be prompted on
                whether to save.
        """
        # First, we save any objects referenced by this object that need to survive a
        # save-load cycle. These are listed in obj.child_attrs. They need to be saved
        # first, so that they get their id's in this backend for the main object to
        # correctly reference. This is done recursively.
        if obj.child_attrs:
            for child_list_name in obj.child_attrs:
                # save any data objects first as this may change the references
                child_list = getattr(obj, child_list_name) or []
                for child_obj in child_list:
                    self.save(child_obj, force=force, no_updates=True)
        # Now we're ready to save the main object.
        # The table_name is the table, the as_dict is the info for the row in the table.
        table_name = obj.table_name
        obj_as_dict = obj.as_dict()
        # check if it's already saved and decide what to do if so:
        if obj.backend is self and self.contains(table_name, obj.id):
            okay_to_update = not no_updates
            update_the_row = force or (
                okay_to_update
                and prompt_for_permission(
                    f"Are you sure you would like to overwrite "
                    f"{self} table={table_name} id={obj.id} with {obj}? "
                    f"(You can use save() with force=True to suppress this.)"
                )
            )
            if update_the_row:
                self.update_row(table_name, obj.id, obj_as_dict)
                return obj.id  # return the id of the updated row
            else:
                return  # return nothing since nothing was done
        else:
            i = self.add_row(table_name, obj_as_dict)
            obj.set_id(i)
            obj.set_backend(self)
            return i

    def save_data(self, data, table_name, i, fixed_name=None):
        """Save the data item of a given row, by default as .ix.npy

        Args:
            data (Array): the numerical data
            table_name (str): The name of the table to save in
            i (int): The id of the row to save in
            fixed_name (the name of the data, just used for the file name
        """
        folder = self.project_directory / table_name
        data_file_name = f"{i}_{fixed_name}{self.data_suffix}"
        np.save(folder / data_file_name, data)

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

    def add_row(self, table_name, obj_as_dict):
        """Save object's serialization to the folder table_name (like adding a row)"""
        folder = self.project_directory / table_name
        if not folder.exists():
            folder.mkdir(parents=True)
        i = self.get_next_available_id(table_name)
        obj_as_dict.update({"id": i})
        fixed_name = fix_name_for_saving(obj_as_dict["name"])
        if "data" in obj_as_dict:
            self.save_data(obj_as_dict["data"], table_name, i, fixed_name)
            obj_as_dict["data"] = None  # FIXME this could instead point to the data.
        file_name = f"{i}_{fixed_name}{self.metadata_suffix}"
        with open(folder / file_name, "w") as f:
            json.dump(obj_as_dict, f, indent=4)
        return i

    def update_row(self, table_name, i, obj_as_dict):
        """Update a file specified by `i` in the folder specified by `table_name`"""
        folder = self.project_directory / table_name
        if not folder.exists():
            folder.mkdir()
        obj_as_dict.update({"id": i})
        fixed_name = fix_name_for_saving(obj_as_dict["name"])
        if "data" in obj_as_dict:
            self.save_data(obj_as_dict["data"], table_name, i, fixed_name)
            obj_as_dict["data"] = None  # FIXME this could instead point to the data.
        file_name = f"{i}_{fixed_name}{self.metadata_suffix}"
        with open(folder / file_name, "w") as f:
            json.dump(obj_as_dict, f, indent=4)

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

    def get_next_available_id(self, table_name, obj=None):
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
