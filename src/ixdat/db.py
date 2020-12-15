"""This module contains the classes which pass on database functionality"""

from .exceptions import DataBaseError
from .backends import DATABASE_BACKENDS


MemoryBackend = DATABASE_BACKENDS["memory"]
memory_backend = MemoryBackend()
DirBackend = DATABASE_BACKENDS["directory"]


class BackendBase:
    def save(self, obj):
        raise NotImplementedError

    def open(self, cls, i):
        raise NotImplementedError

    def load_obj_data(self, obj):
        raise NotImplementedError

    def get_next_available_id(self, table_name):
        raise NotImplementedError


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
        self.backend_name = "memory"  # SHOULD BE SET AFTER __INIT__ FOR LOADED OBJECT
        self._id = None  # SHOULD BE SET AFTER THE __INIT__ OF INHERITING CLASSES
        self.name = None  # MUST BE SET IN THE __INIT__ OF INHERITING CLASSES

    def __repr__(self):
        return f"{self.__class__}(id={self.id}, name={self.name})"

    @property
    def id(self):
        if not self._id:
            if self.backend_name == "memory":
                self._id = memory_backend.get_next_available_id(self.table_name)
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


class PlaceHolderObject:
    def __init__(self, i, cls):
        self.id = i
        self.cls = cls

    def get_object(self):
        return self.cls.open(self.id)


# THIS is proposed as the main mechanism for changing backend, to make
# the shared global nature of it explicit. And in any case, the user
# will never have to deal with the db, except when changing it away
# from the default. This function should probably be exposed in the
# top name space.

def change_database(db_kind, **db_kwargs):
    DB.set_backend(db_kind, **db_kwargs)


def get_database_kind():
    return DB.backend.__class__.__name__
