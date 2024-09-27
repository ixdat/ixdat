class CinfData_Plugin:
    """Class implement direct database read access using external module Cinfdata"""

    def __init__(
        self,
    ):
        self.cinfdata = None
        self._managed_cinfdata_object = None
        self._context_manager_kwargs = None

    def activate(self):
        from cinfdata import Cinfdata

        self.DB = Cinfdata

    def connect(self, setup_name=None, grouping_column=None):
        """setup_name (str): The name of the table inside the database
        grouping_column (str): Either the 'timestamp'/'comment' or 'Comment' column"""
        return self.DB(setup_name=setup_name, grouping_column=grouping_column)

    def __call__(self, **kwargs):
        """**kwargs: setup_name (str) and grouping_column (str) (see connect())"""
        self._context_manager_kwargs = kwargs
        return self

    def __enter__(self):
        if self._managed_cinfdata_object:
            self._context_manager_kwargs = None
            raise RuntimeError(
                "Using the cinfdata plugin as a context manager more "
                "than once at the same time is not supported"
            )

        self._managed_cinfdata_object = self.connect(**self._context_manager_kwargs)
        return self._managed_cinfdata_object

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._managed_cinfdata_object.connection.close()
        self._managed_cinfdata_object = None
        self._context_manager_kwargs = None
