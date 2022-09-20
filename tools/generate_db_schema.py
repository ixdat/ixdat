from collections import defaultdict
from pprint import pformat

import numpy as np
import pyperclip
from ixdat.db import TABLE_CLASSES, LINKERS

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
    schema = _generate_DBML_for_primary_tables(primary_table_classes)

    schema += _generate_DBML_for_linker_tables(linker_tables)

    return schema


def _generate_DBML_for_primary_tables(primary_table_classes):
    """Return DBML for primary tables and table_name -> id_column_name mapping"""
    schema = ""
    # Generate schema for primary classes
    for cls in primary_table_classes:
        schema += f"table {cls.table_name}{{\n"

        for column in cls.columns:
            schema_line = f"  {column.name} {_type_translation[column.ctype]}"
            dbml_column_specs = []
            if column.name == cls.primary_key:
                dbml_column_specs.append("pk")
            if column.foreign_key:
                fk_table, fk_column = column.foreign_key
                dbml_column_specs.append(f"ref: - {fk_table}.{fk_column}")

            if dbml_column_specs:
                schema_line += " [" + ", ".join(dbml_column_specs) + "]"
            schema += schema_line + "\n"
        schema += "}\n\n"

    return schema


def _generate_DBML_for_linker_tables(linker_tables):
    """Return DBML for `linker_tables`, using also the set of `primary_table_classes` and
    the table_name -> id_column_name mapping

    """
    schema = ""
    # Generate schema for linker tables
    for parent_class, owned_object_list in linker_tables:

        # linker table name
        schema += f"table {owned_object_list.joining_table_name}{{\n"

        # foreign key to the parent table
        top_cls = parent_class
        while top_cls.parent_table_class:
            top_cls = top_cls.parent_table_class
        parent_table_name = top_cls.table_name
        parent_id_column_name = (
            owned_object_list.parent_object_id_column_name
            or parent_table_name.rstrip("s") + "_id"
        )
        parent_link = f"  {parent_id_column_name} int [ref:> {parent_table_name}.id]"

        # foreign key to the owned object table
        owned_table_name = owned_object_list.owned_object_table_name
        owned_id_column_name = (
            owned_object_list.owned_object_id_column_name
            or owned_table_name.rstrip("s") + "_id"
        )
        owned_link = f"  {owned_id_column_name} int [ref:> {owned_table_name}.id]"

        # order column
        order_line = "  order int"

        # primary key column
        pk_lines = [
            "  indexes {",
            f"    ({parent_id_column_name}, {owned_id_column_name}) [pk]",
            "  }",
        ]

        # put it together:
        schema += (
            "\n".join([parent_link, owned_link, order_line] + pk_lines + ["}"]) + "\n\n"
        )

    return schema


def main():
    """Generate all DB defining table information, convert to DBML and print and copy
    to clipboard

    """
    schema = generate_dbdiagramio_DBML(
        primary_table_classes=TABLE_CLASSES, linker_tables=LINKERS
    )
    pyperclip.copy(schema)
    print(schema)


if __name__ == "__main__":

    from ixdat.techniques.ec_ms import ECMSCalibration

    cols = ECMSCalibration.get_full_columns()
    print(cols)
    owned = ECMSCalibration.get_full_owned_object_lists()
    print(owned)

    main()
