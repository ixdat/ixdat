def pytest_addoption(parser):
    """Add a program argument option for running external tests with pytest.

    note: this is not related to invoke
    """
    parser.addoption(
        "--external",
        action="store_true",
        dest="external",
        default=False,
        help="enable external tests",
    )


def pytest_configure(config):
    """Set tests marked as external not to run by default."""
    if not config.option.external:
        setattr(config.option, "markexpr", "not external")
