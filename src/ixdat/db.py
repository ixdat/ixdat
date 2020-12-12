"""This module contains the classes which pass on database functionality"""

import json
import numpy as np
from .config import STANDARD_DATA_DIRECTORY
from .exceptions import DataBaseError


def id_from_path(path):
    try:
        return int(path.stem.split("_")[0])
    except ValueError:
        return None


def name_from_path(path):
    return "".join(path.stem.split("_")[1:])


class LocalBackend:
    def __init__(self):
        self.next_available_ids = {}

    def get_next_available_id(self, table_name):
        if table_name in self.next_available_ids:
            i = self.next_available_ids[table_name]
            self.next_available_ids[table_name] += 1
        else:
            i = 1
            self.next_available_ids[table_name] = 2
        return i


local_backend = LocalBackend()


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
        if obj.data_objects:
            # save any data objects first as this may change the references
            for data_obj in obj.data_objects:
                self.save_data_obj(data_obj)
        try:
            table_name = obj.table_name
        except AttributeError:
            table_name = str(type(obj))
        obj_as_dict = obj.as_dict()
        obj.id = self.add_row(obj_as_dict, table_name=table_name)
        obj.backend_name = self.name
        return obj.id

    def open(self, cls, i):
        table_name = cls.table_name
        obj_as_dict = self.open_serialization(table_name, i)
        obj = cls.from_dict(obj_as_dict)
        obj.backend = self
        return obj

    def contains(self, table_name, i):
        return i in self.get_id_list(table_name)

    def load_obj_data(self, obj):
        if not hasattr(obj, "data"):
            # there's no data to be got or obj already has its data
            return
        path_to_row = self.get_path_to_row(obj.table_name, obj.id)
        try:
            return np.load(path_to_row.with_suffix(self.data_suffix))
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
        obj = self.backend.open(cls, i)
        obj.backend_name = self.backend.name
        return obj

    def load_obj_data(self, obj):
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


class Saveable:
    db = DB
    table_name = None  # THIS MUST BE OVERWRITTEN IN INHERITING CLASSES
    column_attrs = None  # THIS SHOULD BE OVERWRITTEN IN INHERITING CLASSES
    extra_column_attrs = None  # THIS CAN BE OVERWRITTEN IN INHERITING CLASSES
    extra_linkers = None  # THIS CAN BE OVERWRITTEN IN INHERITING CLASSES
    data_objects = None  # THIS SHOULD BE OVERWRITTEN IN CLASSES WITH DATA REFERENCES

    def __init__(self, **self_as_dict):
        for attr, value in self_as_dict:
            setattr(self, attr, value)
        if self_as_dict and not self.column_attrs:
            self.column_attrs = {attr: attr for attr in self_as_dict.keys()}
        self.backend_name = "local"  # SHOULD BE SET AFTER __INIT__ FOR LOADED OBJECT
        self._id = None  # SHOULD BE SET AFTER THE __INIT__ OF INHERITING CLASSES
        self.name = None  # MUST BE SET IN THE __INIT__ OF INHERITING CLASSES

    def __repr__(self):
        return f"{self.__class__}(id={self.id}, name={self.name})"

    @property
    def id(self):
        if not self._id:
            if self.backend_name == "local":
                self._id = local_backend.get_next_available_id(self.table_name)
            else:
                raise DataBaseError(
                    f"{self} comes from {self.backend_name} "
                    f"but did not get an id from its backend."
                )
        return self._id

    @id.setter
    def id(self, i):
        self._id = i

    def get_main_dict(self):
        if self.column_attrs is None:
            raise DataBaseError(
                f"{self} can't be seriealized because the class {self.__class__} "
                f"hasn't defined column_attrs"
            )
        self_as_dict = {
            column: getattr(self, attr) for column, attr in self.column_attrs.items()
        }
        return self_as_dict

    def as_dict(self):
        main = self.get_main_dict()
        self_as_dict = main
        if self.extra_column_attrs:
            aux_tables_dict = {
                table_name: {column: getattr(self, attr) for column, attr in extras}
                for table_name, extras in self.extra_column_attrs.items()
            }
            self_as_dict.update(**aux_dict for aux_dict in aux_tables_dict.values())
        if self.extra_linkers:
            linker_tables_dict = {
                table_name: {column: getattr(self, attr) for column, attr in linkers}
                for table_name, linkers in self.extra_linkers.items()
            }
            self_as_dict.update(**linker[1] for linker in linker_tables_dict.values())
        return self_as_dict

    def save(self):
        return self.db.save(self)

    @classmethod
    def from_dict(cls, obj_as_dict):
        return cls(**obj_as_dict)

    @classmethod
    def open(cls, i):
        return cls.db.open(cls, i)

    def load_data(self):
        return self.db.load_obj_data(self)


# THIS is proposed as the main mechanism for changing backend, to make
# the shared global nature of it explicit. And in any case, the user
# will never have to deal with the db, except when changing it away
# from the default. This function should probably be exposed in the
# top name space.

def change_database(db_kind, **db_kwargs):
    DB.set_backend(db_kind, **db_kwargs)


def get_database_kind():
    return DB.backend.__class__.__name__
