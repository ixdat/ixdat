"""This module contains the classes which pass on database functionality"""

from .exceptions import DataBaseError
from .backends import DATABASE_BACKENDS


MemoryBackend = DATABASE_BACKENDS["memory"]  # The default for a local not-yet-saved obj
memory_backend = MemoryBackend()  # The backend to assign id's for in-memory objects
DirBackend = DATABASE_BACKENDS["directory"]  # The default backend for saving


class DataBase:
    """This class is a kind of middle-man between a Backend and a Saveabe class

    The reason for a middle man here is that it enables different databases (backends)
    to be switched between and kept track of in a single ixdat session.

    The DataBase should be initialized with a backend, and by default uses DirBackend,
    which saves to a folder.
    """

    def __init__(self, backend=None):
        """Initialize the database with its backend"""
        backend = backend or DirBackend()
        self.backend = backend

    def save(self, obj):
        """Save a Saveable object with the backend"""
        return self.backend.save(obj)

    def open(self, cls, i):
        """Open and return an object of a Saveable class from the backend"""
        obj = self.backend.open(cls, i)
        obj.backend = self.backend  # How we keep track with multiple backends
        return obj

    def load_obj_data(self, obj):
        """Load and return the numerical data (obj.data) for a Saveable object"""
        return self.backend.get_obj_data(obj)

    def set_backend(self, backend_name, **db_kwargs):
        """Change backend to the class given by backend_name initiated with db_kwargs"""
        if backend_name in DATABASE_BACKENDS:
            BackendClass = DATABASE_BACKENDS[backend_name]
        else:
            raise NotImplementedError(
                f"ixdat doresn't recognize db_kind = '{backend_name}'. If this is a new"
                "database backend, make sure it is added to the DATABASE_BACKENDS "
                "constant in ixdat.backends."
                "Or manually set it directly with DB.backend = <my_backend>"
            )
        self.backend = BackendClass(**db_kwargs)


DB = DataBase()  # initate the database. It functions as a global "constant"


# THIS is proposed as the main mechanism for changing backend, to make
# the shared global nature of it explicit. And in any case, the user
# will never have to deal with the db, except when changing it away
# from the default. This function should probably be exposed in the
# top name space.


def change_database(db_kind, **db_kwargs):
    """Change the backend specifying which database objects are saved to/loaded from"""
    DB.set_backend(db_kind, **db_kwargs)


def get_database_kind():
    """Return the name of the class of which the database backend is an instance"""
    return DB.backend.__class__.__name__


class Saveable:
    """Base class for table-representing classes implementing database functionality.

    This enables seamless interoperability between databse tables and ixdat classes.
    Classes inheriting from this need to provide just a bit of info to define the
    corresponding table, and then then saving and loading should just work.

    At a minimum, the `table_name` and `column_attrs` class attributes need to be
    overwritten in inheriting classes to define the name and columns of the main
    corresponding table. If an auxiliary table is needed to store lists of references
    as rows, this should be represented in `linkers`. Doubly-inheriting classes can use
    `extra_column_attrs` to add extra columns via an auxiliary table without changing
    the main table name.

    ixdat is lazy, only loading things when needed. Correspondingly, all of the columns
    of table mentioned above should refer to (lists of) id's and not actual objects of
    other ixdat classes.

    The class attributes (defined before __init__) and object attributes (defined in
    __init__) are described here. See the other methods and the relevant inheriting
    classes for more info.

    Class attributes:
        db (DataBase): the database, DB, which has the save, open, and load_data methods
        table_name (str): The name of the table or folder in which objects are saved
        column_attrs (dict): {column: attr} where column (str) is the name of the
            column in the table and attr (str) is the name of the attribute of the
            class. The two names column and attr are often but not always the same.
        extra_column_attrs (dict): {table_name: {column: attr}} for auxiliary tables
            to represent the "extra" attributes, for double-inheriting classes.
        linkers (dict): {table_name: (reference_table, {column: attr})} for defining
            the connections between objects.

    Object attributes:
        backend (Backend): the backend where the object is saved. For a
            new, un-saved, object, this is "memory".
        _id (int): the principle key of the object, also accessible as `id`. This should
            be set explicitly in the backend when loading an object. For objects
            initiated directly in the session, it will become the id provided by the
            memory backend, which just counts objects of each table starting with 1.
        name (str): The name of the object. `name` should be a column in ixdat tables.
    """

    db = DB
    table_name = None  # THIS MUST BE OVERWRITTEN IN INHERITING CLASSES
    column_attrs = None  # THIS SHOULD BE OVERWRITTEN IN INHERITING CLASSES
    extra_column_attrs = None  # THIS CAN BE OVERWRITTEN IN INHERITING CLASSES
    extra_linkers = None  # THIS CAN BE OVERWRITTEN IN INHERITING CLASSES
    data_objects = None  # THIS SHOULD BE OVERWRITTEN IN CLASSES WITH DATA REFERENCES

    def __init__(self, **self_as_dict):
        """Initialize a Saveable object from its dictionary serialization

        This the default behavior, and should be overwritten using an argument-free
        call to super().__init__() in inheriting classes.

        Args:
            self_as_dict: all key-word arguments are by default set to object attributes
        """
        for attr, value in self_as_dict:
            setattr(self, attr, value)
        if self_as_dict and not self.column_attrs:
            self.column_attrs = {attr: attr for attr in self_as_dict.keys()}
        self.backend = MemoryBackend  # SHOULD BE SET AFTER __INIT__ FOR LOADED OBJECT
        self._id = None  # SHOULD BE SET AFTER THE __INIT__ OF INHERITING CLASSES
        self.name = None  # MUST BE SET IN THE __INIT__ OF INHERITING CLASSES

    def __repr__(self):
        return f"{self.__class__}(id={self.id}, name={self.name})"

    @property
    def backend_name(self):
        """The name of the backend in which self has been saved to / loaded from"""
        return self.backend.name

    @property
    def id(self):
        """The principle-key identifier. Set by backend or counted in memory."""
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
        """Backends should explicitly set obj.id after loading a Saveable obj"""
        self._id = i

    def get_main_dict(self):
        """Return dict: serializition only of the row of the object's main table"""
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
        """Return dict: serialization of the object main and auxiliary tables"""
        main = self.get_main_dict()
        self_as_dict = main
        if self.extra_column_attrs:
            aux_tables_dict = {
                table_name: {column: getattr(self, attr) for column, attr in extras}
                for table_name, extras in self.extra_column_attrs.items()
            }
            for aux_dict in aux_tables_dict.values():
                self_as_dict.update(**aux_dict)
        if self.extra_linkers:
            linker_tables_dict = {
                table_name: {column: getattr(self, attr) for column, attr in linkers}
                for table_name, linkers in self.extra_linkers.items()
            }
            for linker in linker_tables_dict.values():
                self_as_dict.update(**linker[1])
        return self_as_dict

    def save(self):
        """Save self and return the id. This sets self.backend_name and self.id"""
        return self.db.save(self)

    @classmethod
    def from_dict(cls, obj_as_dict):
        """Return an object built from its serialization."""
        return cls(**obj_as_dict)

    @classmethod
    def open(cls, i):
        """Open an object given its id (the table is cls.table_name)"""
        return cls.db.open(cls, i)

    def load_data(self):
        """Load the data of the object, if ixdat in its laziness hasn't done so yet"""
        return self.db.load_obj_data(self)


class PlaceHolderObject:
    """A tool for ixdat's laziness, instances sit in for Saveable objects."""

    def __init__(self, i, cls):
        """Initiate a PlaceHolderObject with info for loading the real obj when needed

        Args:
            i (int): The id (principle key) of the object represented
            cls (class): Class inheriting from Saveabe and thus specifiying the table
        """
        self.id = i
        self.cls = cls

    def get_object(self):
        """Return the loaded real object represented by the PlaceHolderObject"""
        return self.cls.open(self.id)
