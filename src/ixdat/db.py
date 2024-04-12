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
            database tables have `load` and `get` methods which call downwards. TODO.
    see: https://github.com/ixdat/ixdat/pull/1#discussion_r546400793
"""

from .exceptions import DataBaseError
from .backends import BACKEND_CLASSES, database_backends
from .tools import thing_is_close


class DataBase:
    """This class is a kind of middle-man between a Backend and a Saveable class

    The reason for a middle man here is that it enables different databases (backends)
    to be switched between and kept track of in a single ixdat session.

    The DataBase should be initialized with a backend, and by default uses DirBackend,
    which saves to a folder.
    """

    def __init__(self, backend=None):
        """Initialize the database with its backend"""
        self.backend = backend or database_backends["directory"]
        self.new_object_backend = "none"

    def save(self, obj):
        """Save a Saveable object with the backend"""
        return self.backend.save(obj)

    def get(self, cls, i, backend=None):
        """Select and return object of Saveable class cls with id=i from the backend"""
        backend = backend or self.backend
        obj = backend.get(cls, i)
        # obj will already have obj.id = i and obj.backend = self.backend from backend
        return obj

    def load(self, cls, name):
        """Select and return object of Saveable class cls with name=name from backend"""

    def load_obj_data(self, obj):
        """Load and return the numerical data (obj.data) for a Saveable object"""
        return self.backend.load_obj_data(obj)

    def set_backend(self, backend_name, **db_kwargs):
        """Change backend to the class given by backend_name initiated with db_kwargs"""
        if not isinstance(backend_name, str):
            # Then we assume that it is the backend itself, not the backend name
            self.backend = backend_name
        elif backend_name in BACKEND_CLASSES:
            BackendClass = BACKEND_CLASSES[backend_name]
            self.backend = BackendClass(**db_kwargs)
        else:
            raise NotImplementedError(
                f"ixdat doesn't recognize db_name = '{backend_name}'. If this is a new"
                "database backend, make sure it is added to the DATABASE_BACKENDS "
                "constant in ixdat.backends."
                "Or manually set it directly with DB.backend = <my_backend>"
            )
        return self.backend


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
    # TODO: restructure column_attrs, extra_column_attrs, extra_linkers so that they are
    #  sufficient to fully define the tables in an SQL backend
    column_attrs = None  # THIS SHOULD BE OVERWRITTEN IN INHERITING CLASSES
    extra_column_attrs = None  # THIS CAN BE OVERWRITTEN IN INHERITING CLASSES
    extra_linkers = None  # THIS CAN BE OVERWRITTEN IN INHERITING CLASSES
    # TODO: derive child_attrs somehow from the above class attributes, and have it in
    #   a way where it's easy to tell which id goes with which attribute, i.e. s_ids
    #   goes with series_list
    child_attrs = None  # THIS SHOULD BE OVERWRITTEN IN CLASSES WITH DATA REFERENCES

    def __init__(self, backend=None, **self_as_dict):
        """Initialize a Saveable object from its dictionary serialization

        This is the default behavior, and should be overwritten using an argument-free
        call to super().__init__() in inheriting classes.

        Args:
            self_as_dict: all key-word arguments are by default set to object attributes
        """
        for attr, value in self_as_dict:
            setattr(self, attr, value)
        if self_as_dict and not self.column_attrs:
            self.column_attrs = {attr: attr for attr in self_as_dict.keys()}
        self._backend = None
        self.backend = backend  # backend's setter will look up the backend name
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
                    f"{self!r} comes from {self.backend_name} "
                    "but did not get an id from its backend."
                )
        return self._id

    @property
    def short_identity(self):
        """short_identity is the backend if different from the active backend and the id

        FIXME: The overloaded return here is annoying and dangerous, but necessary for
          `Measurement.from_dict(m.as_dict())` to work as a copy, since the call to
          `fill_object_list` has to specify where the objects represented by
          PlaceHolderObjects live. Note that calling save() on a Saveable object will
          turn the backends into DB.backend, so this will only give id's when saving.
        This is (usually) sufficient to tell if two objects refer to the same thing,
        when used together with the class attribute table_name
        """
        if self.backend is DB.backend:
            return self.id
        return self.backend, self.id

    @property
    def full_identity(self):
        """The full immutable object short_identity as (str, str, str, int)

        Specifically: (backend_type, backend.address, table_name, id)
        """
        return self.backend_type, self.backend.address, self.table_name, self.id

    @property
    def backend(self):
        """The backend the Saveable object was loaded from or last saved to."""
        if not self._backend:
            self._backend = database_backends["none"]
        return self._backend

    @backend.setter
    def backend(self, new_backend):
        """"""
        new_backend = new_backend or DB.new_object_backend
        if isinstance(new_backend, str):
            if new_backend in database_backends:
                new_backend = database_backends[new_backend]
            elif new_backend in BACKEND_CLASSES:
                new_backend = BACKEND_CLASSES[new_backend]()
            else:
                print(f"WARNING! {self} given unrecognized backend = {new_backend}")
        self._backend = new_backend

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
        self.backend = backend

    def get_main_dict(self, exclude=None):
        """Return dict: serializition only of the row of the object's main table

        Args:
            exclude (list): List of attribute names to leave out of the dict
        """
        exclude = exclude or []
        if self.column_attrs is None:
            raise DataBaseError(
                f"{self!r} can't be serialized because the class "
                f"{self.__class__.__name__} hasn't defined column_attrs"
            )
        self_as_dict = {  # FIXME: probably better as loop, fix with table definitions.
            attr: getattr(self, attr)
            for attr in self.column_attrs
            if attr not in exclude
        }
        return self_as_dict

    def as_dict(self, exclude=None):
        """Return dict: serialization of the object main and auxiliary tables"""

        # So that a new object initiated using the dictionary returned here
        #   can find objects referenced by identities (i.e. s_ids for series_list),
        #   we need to make sure they are saved in a real backend. Before adding the
        #   id's to a list.
        # FIXME: There would be a more precise and elegant way to do this if the id's
        #   and corresponding attributes could be connected through class attributes.
        #   To do with table definitions.
        if self.child_attrs:
            for child_object_list_name in self.child_attrs:
                for child_obj in getattr(self, child_object_list_name):
                    if child_obj.backend is database_backends["none"]:
                        database_backends["memory"].save(child_obj)

        exclude = exclude or []
        self_as_dict = self.get_main_dict(exclude=exclude)
        if self.extra_column_attrs:
            # FIXME: comprehension best as loop. Will be redone with proper table defs.
            aux_tables_dict = {
                table_name: {
                    attr: getattr(self, attr) for attr in extras if attr not in exclude
                }
                for table_name, extras in self.extra_column_attrs.items()
            }
            for aux_dict in aux_tables_dict.values():
                self_as_dict.update(**aux_dict)
        if self.extra_linkers:
            # FIXME: comprehension best as loop. Will be redone with proper table defs.
            linker_tables_dict = {
                (table_name, linked_table_name): {attr: getattr(self, attr)}
                for table_name, (linked_table_name, attr) in self.extra_linkers.items()
                if attr not in exclude
            }
            for linked_attrs in linker_tables_dict.values():
                self_as_dict.update(**linked_attrs)

        return self_as_dict

    def __eq__(self, other):
        """Return whether self is functionally equivalent to other

        This means that everything in the dictionary representation is either equal
        or close enough, and that every owned Saveable object in self and other are
        equal by the same condition.

        FIXME: as_dict() should perhaps be an ordered dict. That way we could ensure
            that the order of the checks, in general, and in particular of the
            property names further down, is intentional to keep cheap result determining
            comparisons first and expensive ones last, for performance reasons
        """
        if self is other:
            # If they're actually the same object of course they're equal.
            return True
        if self.__class__ is not other.__class__:
            # If they're not the same class, they are not equal
            return False
        # Otherwise we compare their dictionary representations.
        self_as_dict = self.as_dict()
        other_as_dict = other.as_dict()
        if not len(self_as_dict) == len(other_as_dict):
            # If they don't have the same number of items, they are not equal:
            return False
        if self.extra_linkers:
            linker_id_names = [
                id_name
                for (
                    linker_table_name,
                    (linked_table_name, id_name),
                ) in self.extra_linkers.items()
            ]  # FIXME: This will be made much simpler with coming metaprogramming
        else:
            linker_id_names = []
        for key in self_as_dict:
            # Here we go through the values
            if key not in other_as_dict:
                # other.as_dict() must have all the keys of self.as_dict() to be equal
                return False
            if key in linker_id_names:
                # So, we don't want the linker names.
                # FIXME: Right now there's no way to figure out what the property (named
                #   in child_attrs) is called from the id_name :( ... so we have to
                #   pass here and do a new loop with the child_attrs.
                pass

            if not thing_is_close(self_as_dict[key], other_as_dict[key]):
                # Then the values aren't close (for floats and np arrays) or aren't
                # equal (for all else)
                return False

        # Now we have to go through the owned Saveable objects:
        if self.child_attrs:
            for object_list_name in self.child_attrs:
                object_list = getattr(self, object_list_name)
                other_object_list = getattr(other, object_list_name)
                # These two object lists need to have every corresponding element equal:
                for object, other_object in zip(object_list, other_object_list):
                    if object != other_object:
                        return False

        # If False hasn't been returned yet, then self and other are functionally equal.
        return True

    # This is necessary, because overriding __eq__ means that __hash__ is set to None
    # https://docs.python.org/3/reference/datamodel.html#object.__hash__
    # On the other hand, many Saveable objects are mutable, so maybe shouldn't have hash
    __hash__ = object.__hash__

    def save(self, db=None):
        """Save self and return the id. This sets self.backend_name and self.id"""
        db = db or self.db
        return db.save(self)

    @classmethod
    def get_all_column_attrs(cls):
        """List all attributes of objects of cls that correspond to table columns"""
        all_attrs = cls.column_attrs
        if cls.extra_column_attrs:
            for table, attrs in cls.extra_column_attrs.items():
                all_attrs = all_attrs.union(attrs)
        if cls.extra_linkers:
            for table, (ref_table, attr) in cls.extra_linkers.items():
                all_attrs.add(attr)
        return all_attrs

    @classmethod
    def from_dict(cls, obj_as_dict):
        """Return an object built from its serialization."""
        return cls(**obj_as_dict)

    @classmethod
    def get(cls, i, backend=None):
        """Open an object of cls given its id (the table is cls.table_name)"""
        old_backend = DB.backend
        DB.set_backend(backend or old_backend)
        obj = DB.get(cls, i)  # gets it from the requested backend.
        DB.set_backend(old_backend)
        return obj

    def load_data(self, db=None):
        """Load the data of the object, if ixdat in its laziness hasn't done so yet"""
        db = db or self.db
        return db.load_obj_data(self)


class PlaceHolderObject:
    """A tool for ixdat's laziness, instances sit in for Saveable objects."""

    def __init__(self, identity, cls, backend=None):
        """Initiate a PlaceHolderObject with info for loading the real obj when needed

        Args:
            identity (int or tuple): The id (principle key) of the object represented OR
                the short identity, i.e. a tuple of the id and the backend. In the later
                case, identity[1] overrides a backend if given
            cls (class): Class inheriting from Saveable and thus specifiying the table
            backend (Backend, optional): by default, placeholders objects must live in
                the active backend. This is the case if loaded with get().
        """
        if isinstance(identity, int):
            i = identity
        else:
            backend, i = identity
        self.id = i
        self.cls = cls
        if not backend:  #
            backend = DB.backend
        if not backend or backend == "none" or backend is database_backends["none"]:
            raise DataBaseError(f"Can't make a PlaceHolderObject with backend={backend}")
        self.backend = backend

    def get_object(self):
        """Return the loaded real object represented by the PlaceHolderObject"""
        return self.cls.get(self.id, backend=self.backend)

    @property
    def short_identity(self):
        """Placeholder also has a short_identity to check equivalence without loading"""
        if self.backend is DB.backend:
            return self.id
        return self.backend, self.id


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
            be loaded. Must be specified if object_list is empty.
    """
    cls = cls or object_list[0].__class__
    object_list = object_list or []
    provided_series_ids = [s.id for s in object_list]
    if not obj_ids:
        return object_list
    for identity in obj_ids:
        if identity not in provided_series_ids:
            object_list.append(PlaceHolderObject(identity=identity, cls=cls))
    return object_list


def with_memory(function):
    """Decorator for saving all new Saveable objects initiated in the memory backend"""

    def function_with_memory(*args, **kwargs):
        DB.new_object_backend = "memory"
        to_return = function(*args, **kwargs)
        DB.new_object_backend = "none"
        return to_return

    return function_with_memory
