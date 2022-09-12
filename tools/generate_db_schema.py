from collections import defaultdict
from pprint import pprint, pformat

import numpy as np
import pyperclip

from ixdat import db
from ixdat.db import Saveable


def collect_db_classes(*, verbose=False):
    """Return all DB table defining classes

    Args:
        verbose (bool): The to print information out about the collected information

    Returns:
        list(cls), set((str, OwnedObjectList)): Returns data for the primary tables
        and linker tables. The primary tables is returned as a list of the classes
        that define them and the linker tables as a tuple of
        (class_that_owns_the_owned_object_list, owned_object_list).

    """
    # Collect classes that map to the same DB table
    table_to_class_mapping = defaultdict(list)
    for cls in db.ALL_SAVABLE_CLASSES:
        table_to_class_mapping[cls.table_name].append(cls)

    primary_table_classes = _extract_primary_classes(
        table_to_class_mapping, verbose=verbose
    )
    linker_tables = _extract_linker_tables(primary_table_classes, verbose=verbose)
    return primary_table_classes, linker_tables


def _extract_primary_classes(table_to_class_mapping, *, verbose=False):
    """Return the set of primary classes from all Saveable classes, formed from the
    `table_name -> [class]` mapping in `table_to_class_mapping`

    """
    # Collect the table that actually defines the table, which is the highest class in
    # terms of inheritance (super) of the list of class that all map to the same table
    primary_table_classes = set()
    for table_name, classes in table_to_class_mapping.items():
        top_class = None
        for cls in classes:
            if not top_class:
                top_class = cls
            else:
                if cls in top_class.__mro__:
                    top_class = cls
        primary_table_classes.add(top_class)

        if verbose:
            print("")
            print(f"Table --- {top_class.table_name} ---")
            print(
                f"Defined by {top_class.__name__}(Parent table class: "
                f""
                # top_class.parent_table_class may be None
                f"{getattr(top_class.parent_table_class, '__name__', None)})"
            )

            columns_lines = pformat(top_class.columns).split("\n")
            print(f"Columns: {columns_lines[0]}")
            for line in columns_lines[1:]:
                print("        ", line)
            if len(classes) > 1:
                print("Derived, non-table inducing, sub-classes:")
                for cls in classes:
                    if cls is top_class:
                        continue
                    print("  ", cls.__name__)

    # Sort the primary tables by their distance to Saveable (in terms of inheritance),
    # and then by name, to make sure that tables that are higher on the inheritance
    # chain are generated first
    sorted_primary_table_classes = sorted(
        primary_table_classes,
        key=lambda cls_: (cls_.__mro__.index(Saveable), cls_.__name__),
    )

    return sorted_primary_table_classes


def _extract_linker_tables(primary_table_classes, *, verbose=False):
    """Return the linker tables defined in the set of `primary_table_classes`
    as a set of (source_class, OwnedObjectList) pairs

    """
    linker_tables = set()  # Of (OwningClass, OwnedObjectList)
    for cls in primary_table_classes:
        for owned_object_list in cls.owned_object_lists:
            linker_tables.add((cls, owned_object_list))

    if verbose:
        print("\nLinker tables:")
        for source_class, owned_object_list in linker_tables:
            print(f"{source_class.__name__: <20} --->", owned_object_list)

    return linker_tables


# Type translation for BDML
_type_translation = {
    int: "INTEGER",
    str: "TEXT",
    float: "REAL",
    np.ndarray: "BLOB",
    dict: "JSON",
}
_id_type = _type_translation[int]


def generate_dbdiagramio_DBML(primary_table_classes, linker_tables):
    """Generate DB schema source code for dbdiagram.io"""
    schema, table_name_to_id_column_name = _generate_DBML_for_primary_tables(
        primary_table_classes
    )

    schema += _generate_DBML_for_linker_tables(
        linker_tables, primary_table_classes, table_name_to_id_column_name
    )

    return schema


def _generate_DBML_for_primary_tables(primary_table_classes):
    """Return DBML for primary tables and table_name -> id_column_name mapping"""
    schema = ""
    # Generate schema for primary classes
    table_name_to_id_column_name = {}
    for cls in primary_table_classes:
        schema += f"table {cls.table_name}{{\n"
        # If the class has a `parent_table_class`, then it is an expansion of a parent
        # class and the id from the parent class is re-used
        if cls.parent_table_class:
            if isinstance(cls.parent_table_class, tuple):  # Multiple inheritance case
                for parent_class in cls.parent_table_class:
                    id_column_name = table_name_to_id_column_name[
                        parent_class.table_name
                    ]
                    schema += (
                        f"  {parent_class.table_name.rstrip('s')}_id {_id_type} "
                        f"[pk, ref: - {parent_class.table_name}."
                        f"{id_column_name}]\n"
                    )
            else:  # Single super class
                id_column_name = f"{cls.parent_table_class.table_name.rstrip('s')}_id"
                table_name_to_id_column_name[cls.table_name] = id_column_name
                schema += (
                    f"  {id_column_name} {_id_type} "
                    f"[pk, ref: - {cls.parent_table_class.table_name}.id]\n"
                )
        else:
            table_name_to_id_column_name[cls.table_name] = "id"

        for column in cls.columns:
            schema += f"  {column.name} {_type_translation[column.ctype]}"
            if column.name == "id":
                schema += " [pk, increment]\n"
            else:
                schema += "\n"
        schema += "}\n\n"

    return schema, table_name_to_id_column_name


def _generate_DBML_for_linker_tables(
    linker_tables, primary_table_classes, table_name_to_id_column_name
):
    """Return DBML for `linker_tables`, using also the set of `primary_table_classes` and
    the table_name -> id_column_name mapping

    """
    schema = ""
    # Generate schema for linker tables
    table_name_to_primary_class = {cls.table_name: cls for cls in primary_table_classes}
    for parent_class, owned_object_list in linker_tables:
        # Broken, so skip for now
        if "component" in owned_object_list.joining_table_name:
            continue

        # Owner foreign key
        schema += f"table {owned_object_list.joining_table_name}{{\n"

        def link(cls_):
            """Return column DBML for link class"""
            # Determin the name of the column the foreign key points to
            foreign_key_column_name = table_name_to_id_column_name[cls_.table_name]
            if foreign_key_column_name == "id":
                column_name = cls_.table_name.rstrip("s") + "_id"
            else:
                column_name = foreign_key_column_name

            # Generate link text on the form:
            #   calibration_id INTEGER [pk, ref: > ms_calibrations.calibration_id]
            link_text = (
                f"  {column_name} {_id_type} [pk, ref: > "
                f"{cls_.table_name}.{foreign_key_column_name}]\n"
            )
            return link_text

        # Generate columns for the two foreign keys in a linker class e.g:
        #   calibration_id INTEGER [pk, ref: > ms_calibrations.calibration_id]
        #   ms_cal_result_id INTEGER [pk, ref: > ms_cal_results.id]
        owned_cls = table_name_to_primary_class[owned_object_list.owned_object_class]
        for cls in (parent_class, owned_cls):
            schema += link(cls)
        schema += "}\n\n"

    return schema


def main():
    """Generate all DB defining table information, convert to DBML and print and copy
    to clipboard

    """
    primary_table_classes, linker_tables = collect_db_classes(verbose=False)
    schema = generate_dbdiagramio_DBML(primary_table_classes, linker_tables)
    pyperclip.copy(schema)
    print(schema)


if __name__ == "__main__":
    main()
