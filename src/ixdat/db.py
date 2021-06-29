"""This module contains the classes which pass on database functionality

Note on terminology:
    In ixdat, we seek to use the following naming conventions:
        `load` grabs *an object* from a database backend given its class or table name
            and the name of the specific object desired (see DataBase.load).
        `load_xxx` grabs `xxx` from a database backend given the object for which
            xxx is desired (see DataBase.load_object_data).
        `get` grabs *an object* from a database backend given its class or table name and
            the princple key (the id) of the row in the corresponding table
        `get_xxx` grabs `xxx` from a database backend given the principle key (the id) of
            the row in the table corresponding to xxx (see `Database.get`)

        So `load` works by `name` or existing object, while `get` works by `id`.
        `get_xx` is also used as the counterpart to `set_xx` to grab `xx`, typically a
        managed attribute, from an object in memory.

        `load` and `get` convention holds vertically - i.e. the Backend, the DataBase,
            up through the Saveable parent class for all ixdat classes corresponding to
            database tables have `load` and `get` methods which call downwards.
    see: https://github.com/ixdat/ixdat/pull/1#discussion_r546400793
"""

from .exceptions import DataBaseError
from .backends import BACKEND_CLASSES, database_backends

MemoryBackend = BACKEND_CLASSES["memory"]  # The default for a local not-yet-saved obj
memory_backend = MemoryBackend()  # The backend to assign id's for in-memory objects
DirBackend = BACKEND_CLASSES["directory"]  # The default backend for saving


class DataBase:
    """This class is a kind of middle-man between a Backend and a Savealbe class

    The reason for a middle man here is that it enables different databases (backends)
    to be switched between and kept track of in a single ixdat session.

    The DataBase should be initialized with a backend, and by default uses DirBackend,
    which saves to a folder.
    """

    def __init__(self, backend=None):
        """Initialize the database with its backend"""
        self.backend = backend or DirBackend()

    def save(self, obj):
        """Save a Saveable object with the backend"""
        return self.backend.save(obj)

    def get(self, cls, i):
        """Select and return object of Saveable class cls with id=i from the backend"""
        obj = self.backend.get(cls, i)
        # obj will already have obj.id = i and obj.backend = self.backend from backend
        return obj

    def load(self, cls, name):
        """Select and return object of Saveable class cls with name=name from backend"""

    def load_obj_data(self, obj):
        """Load and return the numerical data (obj.data) for a Saveable object"""
        return self.backend.load_obj_data(obj)

    def set_backend(self, backend_name, **db_kwargs):
        """Change backend to the class given by backend_name initiated with db_kwargs"""
        if backend_name in BACKEND_CLASSES:
            BackendClass = BACKEND_CLASSES[backend_name]
        else:
            raise NotImplementedError(
                f"ixdat doesn't recognize db_name = '{backend_name}'. If this is a new"
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


def change_database(db_name, **db_kwargs):
    """Change the backend specifying which database objects are saved to/loaded from"""
    DB.set_backend(db_name, **db_kwargs)


def get_database_name():
    """Return the name of the class of which the database backend is an instance"""
    return DB.backend.__class__.__name__


class Saveable:
    """Base class for table-representing classes implementing database functionality.

    This enables seamless interoperability between database tables and ixdat classes.
    Classes inheriting from this need to provide just a bit of info to define the
    corresponding table, and then saving and loading should just work.

    At a minimum, the `table_name` and `column_attrs` class attributes need to be
    overwritten in inheriting classes to define the name and columns of the main
    corresponding table. If an auxiliary table is needed to store lists of references
    as rows, this should be represented in `linkers`. Sub-sub classes can use
    `extra_column_attrs` to add extra columns via an auxiliary table without changing
    the main table name.

    ixdat is lazy, only loading things when needed. Correspondingly, all of the columns
    of table mentioned above should refer to (lists of) id's and not actual objects of
    other ixdat classes.

    The class attributes (defined before __init__) and object attributes (defined in
    __init__) are described here. See the other methods and the relevant inheriting
    classes for more info.

    Class attributes:
        db (DataBase): the database, DB, which has the save, get, and load_data methods
        table_name (str): The name of the table or folder in which objects are saved
        column_attrs (set of str): {attr} where attr is the name of the column in the
            table and also the name of the attribute of the class.
        extra_column_attrs (dict): {table_name: {attr}} for auxiliary tables
            to represent the "extra" attributes, for double-inheriting classes.
        linkers (dict): {table_name: (reference_table, attr)} for defining
            the connections between objects.

    Object attributes:
        backend (Backend): the backend where the object is saved. For a
            new, un-saved, object, this is "memory".
        _id (int): the principle key of the object, also accessible as `id`. This should
            be set explicitly in the backend when loading an object. For objects
            initiated directly in the session, it will become the id provided by the
            memory backend, which just counts objects of each table starting with 1.
            TODO: consider renaming.
                See: https://github.com/ixdat/ixdat/pull/1#discussion_r546434676
        name (str): The name of the object. `name` should be a column in ixdat tables.
    """

    db = DB
    table_name = None  # THIS MUST BE OVERWRITTEN IN INHERITING CLASSES
    column_attrs = None  # THIS SHOULD BE OVERWRITTEN IN INHERITING CLASSES
    extra_column_attrs = None  # THIS CAN BE OVERWRITTEN IN INHERITING CLASSES
    extra_linkers = None  # THIS CAN BE OVERWRITTEN IN INHERITING CLASSES
    child_attrs = None  # THIS SHOULD BE OVERWRITTEN IN CLASSES WITH DATA REFERENCES

    def __init__(self, backend="none", **self_as_dict):
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
        self._backend = backend
        self._id = None  # SHOULD BE SET AFTER THE __INIT__ FOR LOADED OBJECTS
        self.name = None  # MUST BE SET IN THE __INIT__ OF INHERITING CLASSES

    def __repr__(self):
        return f"{self.__class__.__name__}(id={self.id}, name='{self.name}')"

    @property
    def id(self):
        """The principle-key identifier. Set by backend or counted in memory."""
        if not self._id:
            if self.backend_type in ("none", "memory"):
                self._id = self.backend.get_next_available_id(self.table_name, obj=self)
                # TODO: Wouldn't it be better if the backend was always asked for the
                #   ID by Saveable.__init__ ?
            else:
                raise DataBaseError(
                    f"{self} comes from {self.backend_name} "
                    "but did not get an id from its backend."
                )
        return self._id

    @property
    def identity(self):
        return self.backend, self.id

    @property
    def full_identity(self):
        return self.backend_type, self.backend.address, self.table_name, self.id

    @property
    def backend(self):
        """The backend the Saveable object was loaded from or last saved to."""
        if not self._backend:
            self._backend = database_backends["none"]
        elif isinstance(self._backend, str):
            if self._backend in database_backends:
                self._backend = database_backends[self._backend]
            elif self._backend in BACKEND_CLASSES:
                self._backend = BACKEND_CLASSES[self._backend]()
            else:
                print(f"WARNING! {self} has unrecognized backend = {self._backend}")
        return self._backend

    @property
    def backend_name(self):
        """The name of the backend in which self has been saved to / loaded from"""
        return self.backend.name

    @property
    def backend_type(self):
        return self.backend.backend_type

    def set_id(self, i):
        """Backends set obj.id here after loading/saving a Saveable obj"""
        self._id = i

    def set_backend(self, backend):
        """Backends set obj.backend here after loading/saving a Saveable obj"""
        self._backend = backend

    def get_main_dict(self):
        """Return dict: serializition only of the row of the object's main table"""
        if self.column_attrs is None:
            raise DataBaseError(
                f"{self} can't be serialized because the class "
                f"{self.__class__.__name__} hasn't defined column_attrs"
            )
        self_as_dict = {attr: getattr(self, attr) for attr in self.column_attrs}
        return self_as_dict

    def as_dict(self):
        """Return dict: serialization of the object main and auxiliary tables"""
        self_as_dict = self.get_main_dict()
        if self.extra_column_attrs:
            aux_tables_dict = {
                table_name: {attr: getattr(self, attr) for attr in extras}
                for table_name, extras in self.extra_column_attrs.items()
            }
            for aux_dict in aux_tables_dict.values():
                self_as_dict.update(**aux_dict)
        if self.extra_linkers:
            linker_tables_dict = {
                (table_name, linked_table_name): {attr: getattr(self, attr)}
                for table_name, (linked_table_name, attr) in self.extra_linkers.items()
            }
            for linked_attrs in linker_tables_dict.values():
                self_as_dict.update(**linked_attrs)
        return self_as_dict

    def save(self, db=None):
        """Save self and return the id. This sets self.backend_name and self.id"""
        db = db or self.db
        return db.save(self)

    @classmethod
    def from_dict(cls, obj_as_dict):
        """Return an object built from its serialization."""
        return cls(**obj_as_dict)

    @classmethod
    def get(cls, i, db=None):
        """Open an object of cls given its id (the table is cls.table_name)"""
        db = db or cls.db
        return db.get(cls, i)

    def load_data(self, db=None):
        """Load the data of the object, if ixdat in its laziness hasn't done so yet"""
        db = db or self.db
        return db.load_obj_data(self)


class PlaceHolderObject:
    """A tool for ixdat's laziness, instances sit in for Saveable objects."""

    def __init__(self, i, cls):
        """Initiate a PlaceHolderObject with info for loading the real obj when needed

        Args:
            i (int): The id (principle key) of the object represented
            cls (class): Class inheriting from Saveable and thus specifiying the table
        """
        self.id = i
        self.cls = cls

    def get_object(self):
        """Return the loaded real object represented by the PlaceHolderObject"""
        return self.cls.get(self.id)


def fill_object_list(object_list, obj_ids, cls=None):
    """Add PlaceHolderObjects to object_list for any unrepresented obj_ids.

    Args:
        object_list (list of objects or None): The objects already known,
            in a list. This is the list to be appended to. If None, an empty
            list will be appended to.
        obj_ids (list of ints or None): The id's of objects to ensure are in
            the list. Any id in obj_ids not already represented in object_list
            is added to the list as a PlaceHolderObject
        cls (Saveable class): the class remembered by any PlaceHolderObjects
            added to the object_list, so that eventually the right object will
            be loaded.
    """
    cls = cls or object_list[0].__class__
    object_list = object_list or []
    provided_series_ids = [s.id for s in object_list]
    if not obj_ids:
        return object_list
    for i in obj_ids:
        if i not in provided_series_ids:
            object_list.append(PlaceHolderObject(i=i, cls=cls))
    return object_list
