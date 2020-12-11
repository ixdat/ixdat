"""This module contains the classes which pass on database functionality"""

import json
import numpy as np
from .config import STANDARD_DATA_DIRECTORY


def id_from_path(path):
    try:
        return int(path.stem.split("_")[0])
    except ValueError:
        return None


def name_from_path(path):
    return "".join(path.stem.split("_")[1:])


class DirBackend:
    """A database backend that loads and saves .ix files from a directory"""

    def __init__(
        self, directory=STANDARD_DATA_DIRECTORY, suffix="ix", data_suffix="ixdata"
    ):
        """Initialize a directory database backend with the directory as Path"""
        self.directory = directory
        self.suffix = suffix
        self.data_suffix = data_suffix

    @property
    def name(self):
        return f"DirBackend({self.directory})"

    def save(self, obj):
        try:
            table_name = obj.table_name
        except AttributeError:
            table_name = str(type(obj))
        obj_as_dict = obj.as_dict()
        obj.id = self.add_row(obj_as_dict, table_name=table_name)
        obj.backend_name = self.name

        if hasattr(obj, "data_objects"):
            for data_obj in obj.data_objects:
                self.save_data_obj(data_obj)
        return obj.id

    def open(self, cls, i):
        table_name = cls.table_name
        obj_as_dict = self.open_serialization(table_name, i)
        obj = cls.from_dict(obj_as_dict)
        obj.backend = self
        return obj

    def contains(self, table_name, i):
        return i in self.get_id_list(table_name)

    def get_obj_data(self, obj):
        if (not hasattr(obj, "data")) or (obj.data is not None):
            # there's no data to be got or obj already has its data
            return
        path_to_row = self.get_path_to_row(obj.table_name, obj.id)
        try:
            obj.data = np.load(path_to_row.with_suffix(self.data_suffix))
            return obj
        except FileNotFoundError:
            # there's no data to be got.
            return

    def save_data_obj(self, data_obj):
        table_name = data_obj.table_name
        if data_obj.backend == self and self.contains(table_name, data_obj.id):
            return data_obj.id  # already saved!
        obj_as_dict = data_obj.as_dict()
        data = obj_as_dict["data"]
        obj_as_dict["data"] = None
        data_obj.id = self.add_row(obj_as_dict, table_name=table_name)
        folder = self.directory / table_name
        data_obj.backend = self
        data_file_name = f"{data_obj.id}_{data_obj.name}.{self.data_suffix}"
        np.save(folder / data_file_name, data)
        return data_obj.id

    def add_row(self, obj_as_dict, table_name):
        folder = self.directory / table_name
        if not folder.exists():
            folder.mkdir()
        i = self.get_next_available_id(table_name)

        name = obj_as_dict["name"]
        file_name = f"{id}_{name}.{self.suffix}"

        with open(folder / file_name, "w") as f:
            json.dump(obj_as_dict, f)
        return i

    def open_serialization(self, table_name, i):
        path_to_row = self.get_path_to_row(table_name, i)
        with open(path_to_row, "r") as file:
            obj_as_dict = json.load(file)
        return obj_as_dict

    def get_path_to_row(self, table_name, i):
        folder = self.directory / table_name
        try:
            return next(
                p
                for p in folder.iterdir()
                if id_from_path(p) == i and p.suffix == self.suffix
            )
        except StopIteration:
            return None

    def get_id_list(self, table_name):
        folder = self.directory / table_name
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
        return max(self.get_id_list(table_name)) + 1

    def __eq__(self, other):
        if other is self:
            return True
        if (
            hasattr(other, "directory")
            and other.directory.resolve() == self.directory.resolve()
            and other.directory.lstat() == self.directory.lstat()
        ):
            return True
        return False


class DataBase:
    def __init__(self, backend=None):
        """Initialize the database with its backend"""
        backend = backend or DirBackend()
        self.backend = backend

    def save(self, obj):
        return self.backend.save(obj)

    def open(self, cls, i):
        return self.backend.open(cls, i)

    def get_obj_data(self, obj):
        return self.backend.get_obj_data(obj)

    def set_backend(self, db_kind, **db_kwargs):
        if db_kind == "dir":
            BackendClass = DirBackend
        else:
            raise NotImplementedError(
                f"ixdat doresn't recognize db_kind = '{db_kind}'. If this is a new"
                "database backend, make sure it is explicitly handled in "
                "DataBase.set_backend() in the db.py module."
                "Or manually set it directly with DB.backend = <my_backend>"
            )
        self.backend = BackendClass(**db_kwargs)


DB = DataBase()


# THIS is proposed as the main mechanism for changing backend, to make
# the shared global nature of it explicit. And in any case, the user
# will never have to deal with the db, except when changing it away
# from the default. This function should probably be exposed in the
# top name space.


def change_database(db_kind, **db_kwargs):
    DB.set_backend(db_kind, **db_kwargs)


def get_database_kind():
    return DB.backend.__class__.__name__
