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

from contextlib import contextmanager

from .exceptions import DataBaseError
from .backends import BACKEND_CLASSES, database_backends
from .tools import deprecate, thing_is_close


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

    @contextmanager
    def temporary_backend(self, backend=None):
        """Switch to a backend for the duration of the `with` block, then switch back.

        Loading an object needs its backend to be the active one, since building
        its lazily loaded related objects happens through that global setting. Using a
        `with`-block guarantees the old backend is restored even if loading fails
        partway through, so a failed load can't leave later saves/loads pointed at
        the wrong database.
        """
        old_backend = self.backend
        self.set_backend(backend or old_backend)
        try:
            yield self.backend
        finally:
            self.set_backend(old_backend)

    def get(self, cls, i, backend=None):
        """Select and return object of Saveable class cls with id=i from the backend"""
        with self.temporary_backend(backend) as selected_backend:
            return selected_backend.get(cls, i)

    def load(self, cls, name, backend=None):
        """Select and return object of Saveable class cls with name=name from backend"""
        with self.temporary_backend(backend) as selected_backend:
            return selected_backend.load(cls, name)

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
    return DB.set_backend(db_name, **db_kwargs)


def get_database_name():
    """Return the name of the class of which the database backend is an instance"""
    return DB.backend.__class__.__name__


class Relationship:
    """Describe how a Saveable object refers to other Saveable objects.

    The key in a class's ``relationships`` dictionary is the object attribute,
    while this description connects it to the stored id attribute and table(s).

    Args:
        linked_table (str): The table containing the related object.
        id_attr (str): The attribute which returns the related object's saved id,
            or its ordered ids when ``many=True``.
        many (bool): Whether the object attribute contains an ordered list.
        storage_table (str or None): The table which stores the relationship. A
            relationship containing many objects requires its own table of
            connections. A relationship containing one object uses the owner's main
            table when this is left out.
        save_related (bool): Whether saving the owner also sends the related object
            or objects through the backend's save process first. This gives new
            related objects saved ids. A save with ``force=True`` also updates related
            objects which are already saved.
    """

    def __init__(
        self,
        linked_table,
        id_attr,
        *,
        many=False,
        storage_table=None,
        save_related=True,
    ):
        if many and not storage_table:
            raise ValueError(
                "A relationship containing many objects needs a storage_table."
            )
        self.linked_table = linked_table
        self.id_attr = id_attr
        self.many = many
        self.storage_table = storage_table
        self.save_related = save_related


def same_short_identity(first, second):
    """Return whether two ``(backend, id)`` references reach the same row.

    A short identity contains a live backend object. Two backend objects can use
    separate connections to the same storage, so comparing the tuples directly
    can report a difference for references which reach the same row. The table is
    supplied by the surrounding relationship or object class.
    """
    first_backend, first_id = first
    second_backend, second_id = second
    return first_id == second_id and first_backend.shares_storage_with(second_backend)


class Saveable:
    """Base class for table-representing classes implementing database functionality.

    This enables seamless interoperability between database tables and ixdat classes.
    Classes inheriting from this need to provide just a bit of info to define the
    corresponding table, and then saving and loading should just work.

    At a minimum, the `table_name` and `column_attrs` class attributes need to be
    overwritten in inheriting classes to define the name and columns of the main
    corresponding table. Sub-sub classes can use `extra_column_attrs` to add extra
    columns via an auxiliary table without changing the main table name.

    A class can also say what Python type of value each of its columns holds, with
    `column_types`, and describe references to other Saveable objects with
    `relationships`. SQL backends use this to build their tables (see
    :module:`~ixdat.backends.relational`).

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
            containing subclass attributes. Definitions are merged over the class
            ancestry, so each class only declares the tables and columns it adds.
        extra_linkers (dict): Older form of relationship metadata, kept so existing
            external Saveable classes continue to work.
        column_types (dict): {attr: python_type} giving the type of the value stored
            in a column, for the columns where it matters. Supported types are
            ``int``, ``float``, ``str``, ``dict``, ``list``, ``tuple``, and
            ``numpy.ndarray``. A column left out here has no fixed type. Definitions
            are merged over the inheriting classes, so a class only declares the
            columns it adds itself.
        column_references (dict): Older form of single-reference metadata, kept so
            existing external Saveable classes continue to work.
        relationships (dict): {object_attr: Relationship} connecting a Python object
            attribute to its stored id attribute. This also says whether the
            relationship holds one object or many, which table stores its ids, and
            whether its related objects join the same save process.

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
    extra_linkers = None  # LEGACY RELATIONSHIP DESCRIPTION
    # every ixdat table has a name. Inheriting classes add the types of the columns
    # they introduce themselves; get_column_types() merges them back together:
    column_types = {"name": str}
    column_references = None  # LEGACY SINGLE-REFERENCE DESCRIPTION
    relationships = None  # THIS CAN BE OVERWRITTEN IN CLASSES WITH REFERENCES
    child_attrs = None  # LEGACY LIST OF RELATED OBJECT ATTRIBUTES

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
        """Return ``(backend, id)`` so a reference says where its object lives.

        Use :func:`same_short_identity` when comparing references from the same
        table. Separate backend objects may connect to the same storage.
        """
        return self.backend, self.id

    @property
    def full_identity(self):
        """Return the storage address, table, and id as an immutable tuple.

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
        main_column_attrs = self.get_main_column_attrs()
        if self.column_attrs is None:
            raise DataBaseError(
                f"{self!r} can't be serialized because the class "
                f"{self.__class__.__name__} hasn't defined column_attrs"
            )
        self_as_dict = {  # FIXME: probably better as loop, fix with table definitions.
            attr: getattr(self, attr)
            for attr in main_column_attrs
            if attr not in exclude
        }
        return self_as_dict

    def as_dict(self, exclude=None):
        """Return dict: serialization of the object main and auxiliary tables"""

        # Save unsaved related objects in memory before reading their id attributes.
        # This gives a dictionary copy enough information to find those objects again.
        for related_obj in self.iter_related_objects():
            if related_obj.backend is database_backends["none"]:
                database_backends["memory"].save(related_obj)

        exclude = exclude or []
        self_as_dict = self.get_main_dict(exclude=exclude)
        for attrs in self.get_extra_column_attrs().values():
            for attr in attrs:
                if attr not in exclude:
                    self_as_dict[attr] = getattr(self, attr)
        for _, attr in self.get_extra_linkers().values():
            if attr not in exclude:
                self_as_dict[attr] = getattr(self, attr)

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
        linker_id_names = {id_name for _, id_name in self.get_extra_linkers().values()}
        relationship_id_names = {
            relationship.id_attr for relationship in self.get_relationships().values()
        }
        for key in self_as_dict:
            # Here we go through the values
            if key not in other_as_dict:
                # other.as_dict() must have all the keys of self.as_dict() to be equal
                return False
            if key in linker_id_names or key in relationship_id_names:
                continue

            if not thing_is_close(self_as_dict[key], other_as_dict[key]):
                # Then the values aren't close (for floats and np arrays) or aren't
                # equal (for all else)
                return False

        compared_attrs = set()
        for object_attr, relationship in self.get_relationships().items():
            compared_attrs.add(object_attr)
            object_list = self._related_object_list(object_attr, relationship)
            other_object_list = other._related_object_list(object_attr, relationship)
            if len(object_list) != len(other_object_list):
                return False
            if any(
                obj != other_obj
                for obj, other_obj in zip(object_list, other_object_list)
            ):
                return False

        # Compare objects declared with the older child_attrs metadata as well.
        for object_attr in self.child_attrs or ():
            if object_attr in compared_attrs:
                continue
            object_list = self._legacy_child_list(object_attr)
            other_object_list = other._legacy_child_list(object_attr)
            if len(object_list) != len(other_object_list):
                return False
            if any(
                obj != other_obj
                for obj, other_obj in zip(object_list, other_object_list)
            ):
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

    def _related_object_list(self, object_attr, relationship):
        """Return one relationship's value as a list for traversal/comparison."""
        value = getattr(self, object_attr)
        if relationship.many:
            return list(value or ())
        return [] if value is None else [value]

    def _legacy_child_list(self, object_attr):
        """Return an older child_attrs value as a list."""
        value = getattr(self, object_attr)
        if value is None:
            return []
        if isinstance(value, (list, tuple, set)):
            return list(value)
        return [value]

    def iter_related_objects(self):
        """Yield objects which should be saved before this object."""
        relationship_attrs = set()
        for object_attr, relationship in self.get_relationships().items():
            relationship_attrs.add(object_attr)
            if relationship.save_related:
                yield from self._related_object_list(object_attr, relationship)
        for object_attr in self.child_attrs or ():
            if object_attr not in relationship_attrs:
                yield from self._legacy_child_list(object_attr)

    @classmethod
    def get_all_column_attrs(cls):
        """List all attributes of objects of cls that correspond to table columns"""
        all_attrs = set(cls.get_main_column_attrs())
        for attrs in cls.get_extra_column_attrs().values():
            all_attrs.update(attrs)
        for _, attr in cls.get_extra_linkers().values():
            all_attrs.add(attr)
        return all_attrs

    @classmethod
    def get_main_column_attrs(cls):
        """Return columns stored in the object's main table."""
        attrs = set(cls.column_attrs or ())
        for relationship in cls.get_relationships().values():
            if not relationship.many and relationship.storage_table is None:
                attrs.add(relationship.id_attr)
        return attrs

    @classmethod
    def get_extra_column_attrs(cls):
        """Return extension-table columns merged over ``cls``'s ancestry.

        Each class declares only the extension tables and columns it introduces.
        Columns declared for the same table are combined. Ancestors with a different
        main table describe a separate persistence model and are left out.
        """
        merged = {}
        for ancestor in reversed(cls.__mro__):
            if getattr(ancestor, "table_name", None) != cls.table_name:
                continue
            for table_name, attrs in (
                ancestor.__dict__.get("extra_column_attrs") or {}
            ).items():
                merged.setdefault(table_name, set()).update(attrs)
        for relationship in cls.get_relationships().values():
            if not relationship.many and relationship.storage_table:
                merged.setdefault(relationship.storage_table, set()).add(
                    relationship.id_attr
                )
        return merged

    @classmethod
    def get_extra_linkers(cls):
        """Return old and new connection-table definitions for ``cls``.

        The two-item tuples are the format used by the older ``extra_linkers`` API.
        New relationships are converted to that format for existing backend code.
        """
        merged = {}
        for ancestor in reversed(cls.__mro__):
            if getattr(ancestor, "table_name", None) == cls.table_name:
                merged.update(ancestor.__dict__.get("extra_linkers") or {})
        for relationship in cls.get_relationships().values():
            if relationship.many:
                merged[relationship.storage_table] = (
                    relationship.linked_table,
                    relationship.id_attr,
                )
        return merged

    @classmethod
    def get_relationships(cls):
        """Return Relationship definitions merged over ``cls``'s ancestry."""
        return cls._merged_class_dict("relationships")

    @classmethod
    def get_column_types(cls):
        """Return {attr: type_name} for all of cls's columns which have a fixed type

        Each class in cls's ancestry contributes the columns it declares itself in
        `column_types`, with the more specific class winning where they disagree.
        Merging lets a class name only the columns it adds. A class inheriting from
        two table-defining classes gets the column types of both.
        """
        return cls._merged_class_dict("column_types")

    @classmethod
    def get_column_references(cls):
        """Return {id_attr: table_name} for ids which point to another table.

        This combines older ``column_references`` declarations with relationships
        containing one object.
        """
        references = cls._merged_class_dict("column_references")
        for relationship in cls.get_relationships().values():
            if not relationship.many:
                references[relationship.id_attr] = relationship.linked_table
        return references

    @classmethod
    def _merged_class_dict(cls, attr):
        """Return the dict class attribute `attr` merged over cls's ancestry"""
        merged = {}
        for ancestor in reversed(cls.__mro__):
            # __dict__, rather than getattr, so that each ancestor contributes only
            # what it declares itself, and the reversed order decides the winner:
            merged.update(ancestor.__dict__.get(attr) or {})
        return merged

    @classmethod
    def from_dict(cls, obj_as_dict):
        """Return an object built from its serialization."""
        return cls(**obj_as_dict)

    @classmethod
    def get(cls, i, backend=None):
        """Open an object of cls given its id (the table is cls.table_name)"""
        return DB.get(cls, i, backend=backend)

    @classmethod
    def load(cls, name, backend=None):
        """Open the most recently saved object of cls with the given name"""
        return DB.load(cls, name, backend=backend)

    @deprecate(
        "0.3.0",
        "`load_data` takes the backend to load the data from. Pass it as `backend=`.",
        "0.4.0",
        kwarg_name="db",
    )
    def load_data(self, backend=None, db=None):
        """Load the data of the object, if ixdat in its laziness hasn't done so yet

        Args:
            backend (Backend): The backend to load the data from. By default, the
                backend this object came from, which is not necessarily the active
                one (e.g. if it was loaded with `get(..., backend=...)`).
            db (Backend): DEPRECATED alias for `backend`.
        """
        return (backend or db or self.backend).load_obj_data(self)


class PlaceHolderObject:
    """A tool for ixdat's laziness, instances sit in for Saveable objects."""

    def __init__(self, identity, cls, backend=None):
        """Initiate a PlaceHolderObject with info for loading the real obj when needed

        Args:
            identity (int or tuple): A local integer id or the ``(backend, id)``
                returned by ``short_identity``. A backend in the tuple takes priority
                over the separate ``backend`` argument.
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
        """Return an identity that can be checked without loading the real object."""
        return self.backend, self.id


def fill_object_list(object_list, obj_ids, cls=None):
    """Add PlaceHolderObjects to object_list for any unrepresented obj_ids.

    Args:
        object_list (list of objects or None): The objects already known,
            in a list. This is the list to be appended to. If None, an empty
            list will be appended to.
        obj_ids (list or None): Local ids or ``(backend, id)`` references to ensure
            are represented. Missing references become PlaceHolderObjects.
        cls (Saveable class): the class remembered by any PlaceHolderObjects
            added to the object_list, so that eventually the right object will
            be loaded. Must be specified if object_list is empty.
    """
    cls = cls or object_list[0].__class__
    object_list = object_list or []
    if not obj_ids:
        return object_list
    for identity in obj_ids:
        if isinstance(identity, int):
            backend, i = DB.backend, identity
        else:
            backend, i = identity
        if not any(
            same_short_identity(obj.short_identity, (backend, i)) for obj in object_list
        ):
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
