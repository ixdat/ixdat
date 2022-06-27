def pytest_addoption(parser):
    parser.addoption(
        "--external",
        action="store_true",
        dest="external",
        default=False,
        help="enable external tests",
    )


def pytest_configure(config):
    if not config.option.external:
        setattr(config.option, "markexpr", "not external")
